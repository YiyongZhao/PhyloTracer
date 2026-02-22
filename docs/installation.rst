.. include:: links.rst

.. _Installation:

Installation
============

PhyloTracer requires **Python 3.8–3.12**, ete3, HyDe, pandas, numpy, tqdm, pypdf, matplotlib, and pyqt5.

.. warning::

   **Python 3.13 is NOT supported.** PhyloTracer depends on ``ete3``, which uses
   the ``cgi`` module that was removed in Python 3.13. Please use Python 3.8–3.12.

Quick Install (Recommended)
---------------------------

.. code:: bash

  # Create a conda environment with Python 3.12
  conda create -n phylotracer python=3.12 -y
  conda activate phylotracer

  # Install PyQt5 (required by ete3 for tree visualization)
  conda install -c conda-forge pyqt=5 -y

  # Install from PyPI
  pip install PhyloTracer

  # For headless/server environments
  export QT_QPA_PLATFORM=offscreen

  # Verify installation
  PhyloTracer --help

Miniconda
---------

We recommend using a Python distribution such as Miniconda to make it easier
to manage modules and environments.

.. code:: bash

  # Get Miniconda for your operating system (Mac or Linux)
  # Answer yes to the questions the Installer asks
  # These commands will download Python 3 for Mac OSX
  curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
  bash Miniconda3-latest-MacOSX-x86_64.sh

Required Python Packages
------------------------

With Miniconda installed, we can use ``pip`` (or ``conda``) to install all of the Python
modules that PhyloTracer requires (ete3, HyDe, pandas, numpy, tqdm, pypdf, matplotlib, pyqt5).

.. code:: bash

  # Install packages with conda
  conda install ete3 pandas numpy tqdm matplotlib pyqt5

Installing PhyloTracer
----------------------

There are two ways that PhyloTracer can be installed once you have all of the required Python modules:
(1) install from PyPI using pip or (2) clone from GitHub and install manually.

To install from PyPI, all we need to do is type the following command:

.. code:: bash

  # Installing from PyPI
  pip install PhyloTracer

Next, we'll take a look at how to install PhyloTracer by cloning it from GitHub.
The commands below take you through every step to accomplish this:

.. code:: bash

    # A convenient one-click installation by using conda
    # (https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html)
    git clone https://github.com/YiyongZhao/PhyloTracer.git
    cd PhyloTracer
    conda env create -f environment.yml
    conda activate PhyloTracer

    # Alternatively, a convenient one-click installation by using pip
    chmod +x install_packages.sh
    bash install_packages.sh

Docker
------

PhyloTracer is also available as a Docker image for fully reproducible environments:

.. code:: bash

    # Build the image
    docker build -t phylotracer .

    # Run PhyloTracer
    docker run --rm -v $(pwd):/data phylotracer --help

Troubleshooting
---------------

**ImportError: cannot import name 'NodeStyle' from 'ete3':**

This happens when PyQt5 is not installed. ``ete3`` requires PyQt5 for tree
visualization components (``NodeStyle``, ``TreeStyle``, ``TextFace``, etc.).
Install it via conda:

.. code:: bash

    conda install -c conda-forge pyqt=5 -y

**Qt platform plugin error:**

If you encounter ``qt.qpa.plugin: Could not load the Qt platform plugin "xcb"``
when running on a headless server, set the following environment variable:

.. code:: bash

    export QT_QPA_PLATFORM=offscreen

**Python 3.13 error (``ModuleNotFoundError: No module named 'cgi'``):**

This is caused by ``ete3`` using the deprecated ``cgi`` module which was removed
in Python 3.13. The solution is to downgrade to Python 3.12:

.. code:: bash

    conda create -n phylotracer python=3.12 -y
    conda activate phylotracer
    pip install PhyloTracer
