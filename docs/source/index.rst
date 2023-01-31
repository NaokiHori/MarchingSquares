################
Marching Squares
################

UNDER CONSTRUCTION.

Although the algorithm in this project is fairly straightforward and general, I think it is worthwhile to note some minor but non-trivial stuffs.

**********************
Scalar field and cells
**********************

In this example, periodic boundary is assumed in :math:`x` direction, while :math:`y` direction is wall-bounded.
The input scalar field (two-dimensional array) is defined at each block dot.
The number of cells, which are formed by the dashed lines, is decreased by :math:`1` in both directions.
In :math:`x` direction, however, because of the periodicity, cells formed by the right and left bounds should be considered, giving one additional cell.

Thus, the number of cells to be considered leads

.. math::

   \begin{cases}
      \text{Wall-bounded} & N_i - 1, \\
      \text{Periodic    } & N_i,
   \end{cases}

where :math:`i = x, y`.

This can be seen in the code here:

.. myliteralinclude:: /../../src/cluster.c
   :language: c
   :tag: number of cells (NOT size of "values") in each direction

.. figure:: images/scalar_field_and_cells/result.png
   :width: 100%

   Relation between input scalar field and cells for marching squares.

***************
Cell and arrows
***************

.. figure:: images/cell_and_arrows/result.png
   :width: 100%

   Possible arrow types in a cell.

There are in total :math:`18` arrow types, which are implemented here:

.. myliteralinclude:: /../../src/cluster.c
   :language: c
   :tag: possible arrow types

