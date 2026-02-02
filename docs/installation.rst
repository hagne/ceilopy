Installation
=============

Dependencies
------------
Dependencies are not installed automatically. Make sure to install the
runtime requirements (conda example below). If you plan to run the production
CLI, also install ``productomator``.

.. code-block:: bash

    conda install -c conda-forge numpy pandas xarray matplotlib netcdf4 python-magic

``python-magic`` needs the system ``libmagic`` library (conda installs it
automatically; with pip you may need to install it via your OS package manager,
e.g. ``brew install libmagic`` on macOS).

Optional (production CLI):

.. code-block:: bash

    pip install productomator

Install from github
-------------------

Clone from Github, enter directory and install with git: 

.. code-block:: bash

    $ git clone https://github.com/hagne/ceilopy
    $ cd ceilopy
    $ pip install .

Build the documentation
-----------------------

Docs are built with Sphinx + nbsphinx.

.. code-block:: bash

    $ cd docs
    $ pip install sphinx nbsphinx sphinx-rtd-theme
    $ make html

Open the documentation
----------------------

.. code-block:: bash

    $ open _build/html/index.html
    
