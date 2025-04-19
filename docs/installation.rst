.. include:: links.rst

.. _Installation:

Installation
============

PhyloTracer requires Python v3.0+, ete3, HyDe, pandas, numpy, tqdm, time, pypdf4, matplotlib, pyqt5 Python modules (listed below).

Miniconda
---------

We recommend using a Python distribution such as Miniconda to make it easier
to manage modules and environments.

.. code:: bash

  # Get Miniconda for your operating system (Mac or Linux)
  # Answer yes to the questions the Installer asks
  # These commands will download Python 3.6 for Mac OSX
  curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
  bash Miniconda3-latest-MacOSX-x86_64.sh

Required Python Packages
------------------------

With Miniconda installed, we can use ``pip`` (or ``conda``) to install all of the Python
modules that PhyloTracer requires (ete3, HyDe, pandas, numpy, tqdm, time, pypdf4, matplotlib, pyqt5).

.. code:: bash

  # Install packages with pip
  conda install ete3, HyDe, pandas, numpy, tqdm, time, pypdf4, matplotlib, pyqt5

Installing PhyloTracer
----------------------

There are two ways that PhyloTracer can be installed once you have all of the required Python modules:
(1) install from PyPI using pip or (2) clone from GitHub and install manually.

To install from PyPI, all we need to do is type the following command:

.. code:: bash

  #installing from PyPI
  pip install pPhyloTracer

Next, we'll take a look at how to install PhyloTracer by cloning it from GitHub.
The commands below take you through every step to accomplish this:

.. code:: bash

    #A convenient one-click installation by using conda (https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html) with the following commands:
    git clone https://github.com/YiyongZhao/PhyloTracer.git
    cd PhyloTracer
    conda env create -f environment.yml
    conda activate phylotracer

    #Alternatively, a convenient one-click installation by using pip (the package installer for Python) with the following commands:
    chmod +x install_packages.sh
    bash install_package.sh

    #Reminder for potential visualization issues: qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found and this application failed to start because no Qt platform plugin could be initialized. Reinstalling the application may fix this problem.
    #Alternative available platform plugins include: eglfs, linuxfb, minimal, minimalegl, offscreen, vnc, wayland-egl, wayland, wayland-xcomposite-#egl, wayland-xcomposite-glx, webgl, xcb. before running PhyloTracer, please execute the following bash command:
    export QT_QPA_PLATFORM=linuxfb
