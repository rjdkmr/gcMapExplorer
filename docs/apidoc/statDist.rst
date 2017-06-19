
.. |markov chain link| raw:: html

	<a href="https://en.wikipedia.org/wiki/Markov_chain#Example" target="_blank">Markov-chain example</a>


statDist module
===============

This module contains functions and methods for calculating stationary
distribution by assuming markov-chain.

This module contains functions to calculate probablity transition matrix and
subsequently to calculate stationary distribution.

Example
-------
Below is a simple example to calculate stationary distribution.
It consists of three major steps:

  * Calcualte median-substracted normalized matrix
  * Calculate probablity transition matrix
  * Calculate stationary distribution


.. code-block:: python

    # median-substracted normalized matrix, stype should be 'o-e', Other arguments can be changed.
    gmlib.normalizer.normalizeGCMapByMCFS('input_raw.gcmap', 'mcfs_O-E.gcmap', stats='median', stype='o-e')

    # Probablity transition matrix, calculate matrix at '40kb' resolution.
    gmlib.statDist.transitionProbablityMatrixForGCMap('mcfs_O-E.gcmap', 'prob_mat.gcmap', '40kb')

    # Stationary distribution at '40kb' resolution
    gmlib.statDist.statDistrByEigenDecompForGCMap('out_prob_mat.gcmap', 'stat_dist.h5', '40kb')



Summary
-------

.. currentmodule:: gcMapExplorer.lib

.. autosummary::
    statDist.calculateTransitionProbablityMatrix
    statDist.transitionProbablityMatrixForCCMap
    statDist.transitionProbablityMatrixForGCMap
    statDist.statDistrByEigenDecompForCCMap
    statDist.statDistrByEigenDecompForGCMap
    statDist.stationaryDistributionByEigenDecomp


Documentation
-------------

.. automodule:: gcMapExplorer.lib.statDist
	:members:
