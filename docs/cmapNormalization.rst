Normalization of Hi-C maps
--------------------------

To normalize the Hi-C maps, several methods are implemented.

* **Iterative Correction** (IC) [1]_
  This method normalize the raw contact map by removing biases from
  experimental procedure. This is an method of matrix balancing, however, in the
  normalized, sum of rows and columns are **not** equal to one.
* **Knight-Ruiz Matrix Balancing** (KR) [2]_
  The Knight-Ruiz (KR) matrix balancing is a fast algorithm to normalize a
  symmetric matrix. A doubly stochastic matrix is obtained after this
  normalization. In this matrix, sum of rows and columns are equal to one.
* **Vanilla-Coverage** (VC) [3]_
  This method was first used for inter-chromosomal map. Later it was used for
  intra-chromosomal map by `Rao et al., 2014 <http://dx.doi.org/10.1016/j.cell.2014.11.021>`_.
  This is a simple method where at first each element is divided by sum of
  respective row and subsequently divided by sum of respective column.
* **Median Contact Frequency Scaling** (MCFS)
  This method can be used to normalize contact map using Median contact values
  for particular distance between two locations/coordinates. At first, Median
  distance contact frequency for each distance is calculated. Subsequently,
  the observed contact frequency is divided by median contact frequency obtained
  for distance between the two locations.

**To perform these normalizations, following tools are implemented:**

.. toctree::
    cmapNormalizer : A GUI application for normalization <commands/cmapNormalizer>
    normKR : Normalization by Knight-Ruiz matrix balancing <commands/normKR>
    normVC : Normalization by Vanilla-Coverage <commands/normVC>
    normIC : Normalization by Iterative Correction <commands/normIC>
    normMCFS : Scale maps using Median/Mean Contact Frequency <commands/normMCFS>


References
~~~~~~~~~~

.. [1] Imakaev *et al*. Iterative correction of Hi-C data reveals hallmarks of
    chromosome organization. `Nature Methods 9, 999–1003 (2012).  <https://doi.org/10.1038/nmeth.2148>`_

.. [2] Knight P and D. Ruiz. A fast algorithm for matrix balancing.
    `IMA J Numer Anal (2013) 33 (3): 1029-1047. <https://doi.org/10.1093/imanum/drs019>`_

.. [3] Lieberman-Aiden *et al*. Comprehensive mapping of long-range interactions reveals
    folding principles of the human genome.
    `Science (2009) 326 : 289-293. <https://doi.org/10.1126/science.1181369>`_
