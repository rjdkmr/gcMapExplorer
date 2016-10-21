.. |appdirs| raw:: html

  <a href="https://pypi.python.org/pypi/appdirs" target="_blank"> appdirs </a>

.. |numpy| raw:: html

  <a href="http://www.numpy.org/" target="_blank"> numpy </a>

.. |scipy| raw:: html

  <a href="http://www.scipy.org/" target="_blank"> scipy </a>

.. |matplotlib| raw:: html

  <a href="http://matplotlib.org/" target="_blank"> matplotlib </a>

.. |h5py| raw:: html

  <a href="http://www.h5py.org/" target="_blank"> h5py </a>

.. |Cython| raw:: html

  <a href="http://cython.org/" target="_blank"> Cython </a>

.. |PyQt5| raw:: html

  <a href="https://www.riverbankcomputing.com/software/pyqt/download5" target="_blank"> PyQt5 </a>

.. |Homebrew| raw:: html

  <a href="http://brew.sh/" target="_blank"> Homebrew </a>

.. |WinPython| raw:: html

  <a href="https://winpython.github.io/" target="_blank"> WinPython3-Qt5 </a>

Requirements and Installation
=============================


Requirements
------------

gcMapExplorer is written in Python3, therefore, it **requires Python3** for
installation. It also requires several external Python packages.

**Package Required during installation:** It has to be installed before gcMapExplorer installation.

* |Cython|

**Package required after installation:** These packages are installed automatically during gcMapExplorer installation.

* |appdirs|
* |numpy|
* |scipy|
* |matplotlib|
* |h5py|

**Package required to install manually:**

* |PyQt5| - It needs to be installed manually. In case of **Python-3.5**, it can be installed automatically from PyPI.

****

Installation Steps on Linux
---------------------------

1. Python3 is available through package managers such as **yum** (Fedora, CentOs), **YaST** (OpenSuse) and **apt-get**
   (Ubuntu, Linux Mint). For example on ubuntu: run ``sudo apt-get install python3`` command to install Python3.

2. Install Cython by ``pip3 install Cython`` command.

3. Similar to Python3, PyQt5 is available through package managers. For example on ubuntu: run ``sudo apt-get install python3-pyqt5`` command
   to install Python3.

4. Install **gcMapExplorer** by ``pip3 install gcMapExplorer`` command.


****


Installation Steps on MacOS
---------------------------

1. Python3 is available through |Homebrew| package manager. After installing Homebrew, run ``brew install python3`` command to install Python3.

2. Install Cython by ``pip3 install Cython`` command.

3. Similar to Python3, PyQt5 is available through |Homebrew|. Run ``brew install pyqt5 --with-python3`` command to install pyqt5.

4. Install **gcMapExplorer** by ``pip3 install gcMapExplorer`` command.


****

Installation Steps on Windows OS
--------------------------------

1. Download and install |WinPython|. Note that WinPython should include PyQt5.

2. Open WinPython directory (Default is in C:/ drive) and click on **"WinPython Command Prompt"**. It will open a command prompt terminal.

3. Run ``pip3 install gcMapExplorer`` in command prompt terminal to install **gcMapExplorer**

.. note::
  To execute gcMapExplorer command, simple command prompt terminal (from Start Menu) might not work.
  Use **"WinPython Command Prompt"** present in WinPython directory to launch or execute gcMapExplorer.
