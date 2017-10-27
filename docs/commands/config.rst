config
======

**Description:**
To print configuration file and clean scratch directory

This can be used to print configuration file and its content.

**Usage:**

.. code-block:: bash

    gcMapExplorer config [-h] [-cs] [-pc]



``-cs,--clean-scratch``
-----------------------
Clean scratch directory
gcMapexplorer try to remove all temporary files present in scracth
directory. However, due to crash and errors, sometimes these files cannot be removed.
Therefore, by this simple command, all temporary files from the scracth directory
is removed.

.. note:: Temporary files could be huge and therefore it is advised to run this command periodically.
.. note:: Do not use this command when any of the gcMapExplorer tools or module are in execution process.

``-pc,--print-config``
----------------------
Print configuration file. It will print location configuration file and its content.
