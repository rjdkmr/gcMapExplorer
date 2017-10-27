
.. |markov chain link| raw:: html

	<a href="https://en.wikipedia.org/wiki/Markov_chain#Example" target="_blank">Markov-chain example</a>


statDist module
===============

This module contains functions and methods for calculating stationary
distribution by assuming markov-chain.

This module contains functions to calculate probability transition matrix and
subsequently to calculate stationary distribution.

Example
-------
Below is a simple example to calculate stationary distribution.
It consists of three major steps:

  * Calculate median-subtracted normalized matrix
  * Calculate probability transition matrix
  * Calculate stationary distribution


.. code-block:: python

    # median-subtracted normalized matrix, stype should be 'o-e', Other arguments can be changed.
    gmlib.normalizer.normalizeGCMapByMCFS('input_raw.gcmap', 'mcfs_O-E.gcmap', stats='median', stype='o-e')

    # Probability transition matrix, calculate matrix at '40kb' resolution.
    gmlib.statDist.transitionProbabilityMatrixForGCMap('mcfs_O-E.gcmap', 'prob_mat.gcmap', '40kb')

    # Stationary distribution at '40kb' resolution
    gmlib.statDist.statDistrByEigenDecompForGCMap('out_prob_mat.gcmap', 'stat_dist.h5', '40kb')



Summary
-------

.. currentmodule:: gcMapExplorer.lib

.. autosummary::
    statDist.calculateTransitionProbabilityMatrix
    statDist.transitionProbabilityMatrixForCCMap
    statDist.transitionProbabilityMatrixForGCMap
    statDist.statDistrByEigenDecompForCCMap
    statDist.statDistrByEigenDecompForGCMap
    statDist.stationaryDistributionByEigenDecomp


Documentation
-------------

.. automodule:: gcMapExplorer.lib.statDist
	:members:
