################
Marching Squares
################

|License|_ |LastCommit|_

.. |License| image:: https://img.shields.io/github/license/NaokiHori/MarchingSquares
.. _License: https://opensource.org/licenses/MIT

.. |LastCommit| image:: https://img.shields.io/github/last-commit/NaokiHori/MarchingSquares/main
.. _LastCommit: https://github.com/NaokiHori/MarchingSquares/commits/main

.. image:: https://github.com/NaokiHori/MarchingSquares/blob/main/.github/thumbnail.png
   :width: 100%

********
Overview
********

`Marching squares <https://en.wikipedia.org/wiki/Marching_squares>`_ for wall-bounded and periodic boundaries.

**********
Motivation
**********

On one hand just for fun (99%).
On the other hand (1%) I wanted to extend `contour in Matplotlib <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.contour.html>`_ for periodic domains.

**********
Dependency
**********

   * C compiler

   * GNU make (not essential but recommended)

   * Gnuplot (not essential but recommended)

***********
Quick start
***********

Sample dataset is attached (in ``data`` directory) for stand-alone execution:

.. code-block::

   $ make
   $ ./a.out

Result is written to ``clusters.dat``, which is an ASCII file containing coordinates of each polygonal chain.
They can be visualised by ``Gnuplot``, e.g.,

.. code-block::

   $ plot 'clusters.dat' u 2:3 w lp

***
API
***

The main function of this library is ``make_clusters``, which is declared as

.. code-block::

   int make_clusters(
       const bool periods[2],
       const double lengths[2],
       const size_t sizes[2],
       const double threshold,
       const double * restrict xs,
       const double * restrict ys,
       const double * restrict values,
       size_t * restrict nclusters,
       cluster_t *** restrict clusters
   );

Parameters:

* ``periods``

   Boundary condition of two directions, periodic (``true``) or not (``false``)

* ``lengths``

   Physical sizes

* ``sizes``

   Number of grid points

* ``threshold``

   Threshold to draw contour lines

* ``xs``

   ``x`` coordinate (length of this vector should be ``sizes[0]``)

* ``ys``

   ``y`` coordinate (length of this vector should be ``sizes[1]``)

* ``values``

   Two-dimensional scalar field

* ``nclusters``

   Result, number of clusters

* ``clusters``

   Result, all polygonal chains

