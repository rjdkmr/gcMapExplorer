#!/usr/bin/env python3
#
# Author: Rajendra Kumar
#
# This file is part of gcMapExplorer
# Copyright (C) 2016-2017  Rajendra Kumar, Ludvig Lizana, Per Stenberg
#
# gcMapExplorer is a free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# gcMapExplorer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with gcMapExplorer.  If not, see <http://www.gnu.org/licenses/>.
#
#=============================================================================

import numpy as np
from multiprocessing import Process
import re, sys, os, shutil
import json
import logging

from . import ccmapHelpers as cmh
from . import ccmap as cmp


def remove_by_scores(peaks, ccmap):
	boxes = dict()
	c = 0
	for p in peaks:
		box = TADBOX(p)
		boxes[c] = box
		c += 1

	for i in range(len(boxes)):
		if boxes[i] is None:
			continue
		si0 = boxes[i].startingPoint[0]
		si1 = boxes[i].startingPoint[1]
		ei0 = boxes[i].endPoint[0]
		ei1 = boxes[i].endPoint[1]
		for j in range(i):
			if boxes[j] is None:
				continue
			sj0 = boxes[j].startingPoint[0]
			sj1 = boxes[j].startingPoint[1]
			ej0 = boxes[j].endPoint[0]
			ej1 = boxes[j].endPoint[1]

			if (abs(si0 - sj0) <=4 and abs(si1 - sj1) <=4) or (abs(ei0 - ej0) <=4 and abs(ei1 - ej1) <=4) :
				if abs(boxes[i].get_score(ccmap, 1.0) - boxes[j].get_score(ccmap, 1.0)) <= 0.05:
					if boxes[i].size > boxes[j].size:
						boxes[j] = None
					else:
						boxes[i] = None
						break

	new_peaks = []
	scores = []
	for key in boxes:
		if boxes[key] is not None:
			new_peaks.append(boxes[key].Vertex)
			scores.append(boxes[key].get_score(ccmap, 1.0))

	idx = np.argsort(scores)
	percent = int(len(idx)*50/100)

	return np.asarray(new_peaks)


def remove_peaks_from_diagonals(peaks, includeCells=4):
	new_peaks = []
	for p in peaks:
		if abs(p[0] - p[1]) <= includeCells:
			continue
		if p[0] > p[1]:
			continue
		new_peaks.append(p)
	return np.asarray(new_peaks)

class TADBOX(object):
	def __init__(self, vertex):
		self.Vertex = vertex
		self.size = None
		self.startingPoint = None
		self.endPoint = None
		self.set_start_end_and_size()
		self.edgeIndex = []
		self.get_edge_index()
		self.midPoint = self.get_mid_point()
		self.score = None

	def get_edge_index(self):
		iVertex = self.startingPoint[0] + self.size
		jVertex = self.startingPoint[1] + self.size
		for i in range(self.startingPoint[0], iVertex):
			if abs(i - jVertex) <= 4:
				continue
			self.edgeIndex.append([i, jVertex])

		for j in range(self.startingPoint[1], jVertex):
			if abs(j - iVertex) <= 4:
				continue
			self.edgeIndex.append([iVertex, j])

		self.edgeIndex = np.asarray(self.edgeIndex).T

	def set_start_end_and_size(self):
		minp = min(self.Vertex[0], self.Vertex[1])
		maxp = max(self.Vertex[0], self.Vertex[1])
		self.startingPoint = (minp, minp)
		self.endPoint = (maxp, maxp)
		self.size = maxp - minp + 1

	def get_mid_point(self):
		i = int ( (self.startingPoint[0] + self.size)/2 )
		self.midPoint = [ i, i ]

	def calc_score(self, ccmap, pfactor):
		self.score = -1 * pfactor * np.log(ccmap.matrix[self.edgeIndex]).mean()

	def get_score(self, ccmap, pfactor):
		if self.score is None:
			self.calc_score(ccmap, pfactor)
		return self.score

	def copy(self, del_index=False):
		box = TADBOX(self.startingPoint, self.size)
		if del_index:
			box.edgeIndex = None
		return box

'''
class TADFINDER(object):
	def __init__(self, ccmap, pfactor=1, mc_steps=100000, minsize=100, maxsize=1e5, expand=1000):
		self.hicimap = ccmap
		self.check_input_box_size(minsize, maxsize, expand)
		self.OutPutBoxes = None
		self.minBoxSize = minsize
		self.maxBoxSize = maxsize
		self.expansion = 1000
		self.boxsize = minsize
		self.moving_box = None
		self.pfactor = pfactor
		self.mc_steps = mc_steps

	def check_input_box_size(self, minsize, maxsize, expand):
		if self.hicimap.shape[0] < maxsize:
			raise ValueError ('Input maximum box size is larger than HiC-map size.')

		if self.hicimap.shape[0] < minsize:
			raise ValueError ('Input minimum box size is larger than HiC-map size.')


	def minimize_box_shape_local(self):
		""" Find a box with minimum energy at given starting point for a given box size
		"""
		old_score = self.moving_box.get_score(self.hicimap, self.pfactor)
		new_box = self.moving_box
		for i in range(1, self.expansion):
			if self.moving_box.size + i >= self.ccmap.shape[0]:
				break
			temp_box = TADBOX(self.moving_box.startingPoint, self.moving_box.size + i)
			new_score = temp_box.get_score(ccmap, self.pfactor)
			if new_score < old_score:
				new_box = temp_box
				oldscore = new_score
		self.OutPutBoxes.append(new_box.cop(del_index=True))

	def move_box(self):
		"""	Move box according to Metropolis Monte Carlo for a given box size
		"""
		if moving_box is not None:
			old_score = self.moving_box.get_score(self.hicimap, self.pfactor)
		else:
			old_score = 100000

		s = np.random.random_integers(int(self.boxsize/2), high=self.hicimap.shape[0]-int(self.boxsize/2))

		new_box = TADBOX((s, s), self.boxsize)
		new_score = new_box.get_score(self.hicimap, self.pfactor)
		moved = False
		if new_score < old_score:
			self.moving_box = new_box
			moved = True
		else:
			p = np.exp( -1 * self.pfactor * (new_score - old_score) )
			if p > random.random_sample():
				self.moving_box = new_box
				moved = True


	def do_MC_sampling(self):
		""" Do monte carlo sampling for the given box size
		"""
		for i in range(self.steps):
			if self.move_box():
				self.minimize_box_shape_local()

	def find(self, forSize=None):
		""" Do monte carlo sampling for all box sizes
		"""
		if forSize is None:
			for b in range(self.minBoxSize, self.maxBoxSize, self.expansion):
				self.do_MC_sampling()
				self.boxsize = b + self.expansion
		else:
			self.boxsize = forSize
			self.do_MC_sampling()

	def remove_duplicates(self):
		new_boxes = []
		for i in range(len(self.OutPutBoxes)):
			if i == 0:
				new_boxes.append(self.OutPutBoxes[i])
			else:
				for b in new_boxes:
					if b.startingPoint == self.OutPutBoxes[i].startingPoint and b.size == self.OutPutBoxes[i].size:
						continue
					else:
						new_boxes.append(self.OutPutBoxes[i])
						break

		self.OutPutBoxes = new_boxes

	def write_output_file(self, filename):
		fout = open(filename, 'w')
		for b in self.OutPutBoxes:
			fout.write('{0} {1} {2}' .format(b.startingPoint[0], b.size, np.exp( -1 * (b.score/self.pfactor) )))

	def run(self, outfile, forSize=None):
		self.find(forSize)
		self.remove_duplicates()
		self.write_output_file(outfile)
'''

from math import ceil
from multiprocessing import cpu_count

def _get_chunks(shape, ncpu):
	"""Split the array into equal sized chunks based on the number of
	available processors. The last chunk in each dimension absorbs the
	remainder array elements if the number of CPUs does not divide evenly into
	the number of array elements.

	Examples
	--------
	>>> _get_chunks((4, 4), 4)
	((2, 2), (2, 2))
	>>> _get_chunks((4, 4), 2)
	((2, 2), (4,))
	>>> _get_chunks((5, 5), 2)
	((2, 3), (5,))
	>>> _get_chunks((2, 4), 2)
	((1, 1), (4,))
	"""
	chunks = []
	nchunks_per_dim = int(ceil(ncpu ** (1./len(shape))))

	used_chunks = 1
	for i in shape:
		if used_chunks < ncpu:
			regular_chunk = i // nchunks_per_dim
			remainder_chunk = regular_chunk + (i % nchunks_per_dim)

			if regular_chunk == 0:
				chunk_lens = (remainder_chunk,)
			else:
				chunk_lens = ((regular_chunk,) * (nchunks_per_dim - 1) + (remainder_chunk,))
		else:
			chunk_lens = (i,)

		chunks.append(chunk_lens)
		used_chunks *= nchunks_per_dim
	return tuple(chunks)


def apply_parallel(function, array, chunks=None, depth=0, mode=None, extra_arguments=(), extra_keywords={}):
	"""Map a function in parallel across an array.

	Split an array into possibly overlapping chunks of a given depth and
	boundary type, call the given function in parallel on the chunks, combine
	the chunks and return the resulting array.

	Parameters
	----------
	function : function
	    Function to be mapped which takes an array as an argument.
	array : numpy array
	    Array which the function will be applied to.
	chunks : int, tuple, or tuple of tuples, optional
	    A single integer is interpreted as the length of one side of a square
	    chunk that should be tiled across the array.  One tuple of length
	    ``array.ndim`` represents the shape of a chunk, and it is tiled across
	    the array.  A list of tuples of length ``ndim``, where each sub-tuple
	    is a sequence of chunk sizes along the corresponding dimension. If
	    None, the array is broken up into chunks based on the number of
	    available cpus. More information about chunks is in the documentation
	    `here <https://dask.pydata.org/en/latest/array-design.html>`_.
	depth : int, optional
	    Integer equal to the depth of the added boundary cells. Defaults to
	    zero.
	mode : {'reflect', 'symmetric', 'periodic', 'wrap', 'nearest', 'edge'}, optional
	    type of external boundary padding.
	extra_arguments : tuple, optional
	    Tuple of arguments to be passed to the function.
	extra_keywords : dictionary, optional
	    Dictionary of keyword arguments to be passed to the function.

	Notes
	-----
	Numpy edge modes 'symmetric', 'wrap', and 'edge' are converted to the
	equivalent `dask` boundary modes 'reflect', 'periodic' and 'nearest',
	respectively.
	"""
	import dask.array as da

	if chunks is None:
		shape = array.shape
		try:
			ncpu = cpu_count()
		except NotImplementedError:
			ncpu = 4
		chunks = _get_chunks(shape, ncpu)

	if mode == 'wrap':
		mode = 'periodic'
	elif mode == 'symmetric':
		mode = 'reflect'
	elif mode == 'edge':
		mode = 'nearest'

	def wrapped_func(arr):
		return function(arr, *extra_arguments, **extra_keywords)

	darr = da.from_array(array, chunks=chunks)
	return darr.map_overlap(wrapped_func, depth, boundary=mode).compute()
