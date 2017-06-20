
config module
=============

It contains functions that handle configuration of gcMapExplorer. gcMapExplorer
uses some default settings and options. This can be read and changed through
these modules.

Configuraion file structure
---------------------------

::

    Configuration
        ├─────────── Dirs
        │             └──────── WorkingDirectory
        │
        └─────────── Programs
                      ├──────── bigWigInfo
                      └──────── bigWigToWig


Examples
--------

.. code-block:: python3

    import gcMapExplorer

    gcMapExplorer.config.cleanScratch()   # Clean default scratch directory

    # Change scratch directory
    gcMapExplorer.config.updateConfig('Dirs', 'WorkingDirectory', 'Path/to/new/scratch/directory')

    # Set path to bigWigInfo program
    gcMapExplorer.config.updateConfig('Programs', 'bigWigInfo', 'Path/to/bigWigInfo')

    # Set path to bigWigToWig program
    gcMapExplorer.config.updateConfig('Programs', 'bigWigToWig', 'Path/to/bigWigToWig')

    # Print current configuration file content
    gcMapExplorer.config.printConfig()

    # Get configuration
    config = gcMapExplorer.config.getConfig()

    # Get scratch directory
    print(config['Dirs']['WorkingDirectory'])


Summary
-------
.. currentmodule:: gcMapExplorer.config

.. autosummary::
    updateConfig
    getConfig
    printConfig
		cleanScratch

.. automodule:: gcMapExplorer.config
	:members: updateConfig, getConfig, printConfig, cleanScratch
