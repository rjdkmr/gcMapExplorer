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
import os, re
import string
import random

from gcMapExplorer.config import getConfig

# get configuration
config = getConfig()

dtype_npBINarray = 'float32'

def sorted_nicely( inputList ):
	""" Sorts the given given list in the way that is expected.

	Parameters
	----------
	inputList : list
		The input list to be sorted.

	Returns
	-------
	outputList : list
		The sorted list

	"""
	convert = lambda text: int(text) if text.isdigit() else text
	alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
	return sorted(inputList, key = alphanum_key)

def locate_significant_digit_after_decimal(value):
	""" Get location at which significant digit start after decimal

	Parameters
	----------
	value : float
		Input value

	Returns
	-------
	value : int
		Number of zeros after which digit start in a small decimal number.

	"""
	location = 0
	base = 0
	while( base <=0 ):
		base = int( abs(value) * 10 ** location  )
		location = location + 1

	return location-1

def kth_diag_indices(k, a):
	""" Get diagonal indices of 2D array 'a' offset by 'k'

	Parameters
	----------
	k : int
		Diagonal offset

	a : numpy.ndarray
		Input numpy 2D array

	Returns
	-------
	indices : tuple of two numpy.ndarray
		It contain indences of elements that are at the diagonal offset by 'k'.

	"""
	rows, cols = np.diag_indices_from(a)
	if k < 0:
		return rows[:k], cols[-k:]
	elif k > 0:
		return rows[k:], cols[:-k]
	else:
		return rows, cols

def detectOutliers1D(points, thresh=3.5):
	"""Returns a boolean array with ``True`` if points are outliers and False
	otherwise.

	Parameters
	----------
	points : numpy.ndarray
		An numobservations by numdimensions array of observations

	thresh : float
		The modified z-score to use as a threshold. Observations with
	    a modified z-score (based on the median absolute deviation) greater
	    than this value will be classified as outliers.


	Returns
	-------
	outBool : numpy.ndarray
		A numobservations-length boolean array.

	References
	----------
	    Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
	    Handle Outliers", The ASQC Basic References in Quality Control:
	    Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.

	"""

	if len(points.shape) == 1:
		points = points[:,None]
	#print(points)
	median = np.median(points, axis=0)
	diff = np.sum((points - median)**2, axis=-1)
	diff = np.sqrt(diff)
	med_abs_deviation = np.median(diff)

	modified_z_score = 0.6745 * diff / med_abs_deviation

	return modified_z_score > thresh

def detectOutliersMasked1D(points, thresh=3.5):
	"""Returns a masked array where outliers are masked with preserved input mask.

	Parameters
	----------
	points : numpy.ma.ndarray
		An numobservations by numdimensions array of observations

	thresh : float
		The modified z-score to use as a threshold. Observations with
	    a modified z-score (based on the median absolute deviation) greater
	    than this value will be classified as outliers.

	Returns
	-------
	maskArray : numpy.ma.ndarray
		A numobservations-length masked array.

	References
	----------
	    Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
	    Handle Outliers", The ASQC Basic References in Quality Control:
	    Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.

	"""

	inputPoint = points.compressed().copy()
	idx = detectOutliers1D(inputPoint, thresh=thresh)
	inputPoint[idx] = 0.0

	result = np.zeros(points.shape, dtype=points.dtype)
	result[~points.mask] = inputPoint[:]
	result = np.ma.masked_values(result, 0.0)

	return result

def resolutionToBinsize(resolution):
	"""Return the bin size from the resolution unit

	It is a convenient function to convert resolution unit to binsize.
	It has a support of base (b), kilobase (kb), megabase (mb) and
	gigabase (gb) unit. It also convert decimal resolution unit as shown below
	in examples.

	Parameters
	----------

	resolution:	str
		resolution in b, kb, mb or gb.

	Returns
	-------
	binsize : int
		bin size

	Examples
	--------
		>>> resolutionToBinsize('1b')
		1
		>>> resolutionToBinsize('10b')
		10
		>>> resolutionToBinsize('1kb')
		1000
		>>> resolutionToBinsize('16kb')
		16000
		>>> resolutionToBinsize('1.23kb')
		1230
		>>> resolutionToBinsize('1.6mb')
		1600000
		>>> resolutionToBinsize('1.457mb')
		1457000

	"""
	factor = None
	base = None
	if re.search('\d+b', resolution) is not None:
		factor = 1
	elif re.search('kb', resolution) is not None:
		factor = 1000
	elif re.search('mb', resolution) is not None:
		factor = 1000000
	elif re.search('gb', resolution) is not None:
		factor = 1000000000
	else:
		raise ValueError(' \'{0}\' keyword of resolution not implemented.' .format(resolution))

	if re.search('(\d+.\d+)', resolution) is not None:
		base = float( re.search('(\d+.\d+)', resolution).group() )
	else:
		base = int( re.search('(\d+)', resolution).group() )

	return int(base*factor)

def binsizeToResolution(binsize):
	"""Return the resolution unit from the bin size

	It is a convenient function to convert binsize into resolution unit.
	It has a support of base (b), kilobase (kb), megabase (mb) and
	gigabase (gb) unit. It also convert binsize to decimal resolution unit as
	shown below in examples.

	Parameters
	----------

	binsize:	int
		bin size

	Returns
	-------
	resolution : str
		resolution unit

	Examples
	--------
		>>> binsizeToResolution(1)
		'1b'
		>>> binsizeToResolution(10)
		'10b'
		>>> binsizeToResolution(10000)
		'10kb'
		>>> binsizeToResolution(100000)
		'100kb'
		>>> binsizeToResolution(125500)
		'125.5kb'
		>>> binsizeToResolution(1000000)
		'1mb'
		>>> binsizeToResolution(1634300)
		'1.6343mb'

	"""

	suffix = { 0:'b', 1: 'kb', 2: 'mb', 3: 'gb' }
	rest = None
	count = 0
	while (1):
		rest = binsize/(1000**count)
		if rest >= 1000:
			count += 1
		else:
			break

	if rest - int(rest) == 0:
		resolution = str(int(rest)) + suffix[count]
	else:
		resolution = str(rest) + suffix[count]

	return resolution


def getRandomName(size=10, chars=string.ascii_letters + string.digits ):
    """ Random name generator

    Parameters
    ----------
    size : int
        Number of alphabets in the name.

    Returns
    -------
    name : str
        Randomly generated name.

    """
    s = ''.join(random.choice(chars) for _ in range(size))
    return s

class MapNotFoundError(Exception):
	def __init__(self, value):
		super(MapNotFoundError, self).__init__()
		self.value = value

	def __str__(self):
		return repr(self.value)

class ResolutionNotFoundError(Exception):
	def __init__(self, value):
		super(ResolutionNotFoundError, self).__init__()
		self.value = value

	def __str__(self):
		return repr(self.value)
