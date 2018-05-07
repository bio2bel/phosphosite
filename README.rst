Bio2BEL Phosphosite |build| |coverage| |documentation|
======================================================
Converts post-translational modifications and associated mutations to BEL

Installation |pypi_version| |python_versions| |pypi_license|
------------------------------------------------------------
``bio2bel_phosphosite`` can be installed easily from `PyPI <https://pypi.python.org/pypi/bio2bel_phosphosite>`_ with the
following code in your favorite terminal:

.. code-block:: sh

    $ python3 -m pip install bio2bel_phosphosite

or from the latest code on `GitHub <https://github.com/bio2bel/phosphosite>`_ with:

.. code-block:: sh

    $ python3 -m pip install git+https://github.com/bio2bel/phosphosite.git@master

Setup
-----
PhosphoSitePlus can be downloaded and populated from either the Python REPL or the automatically installed command line
utility.

Python REPL
~~~~~~~~~~~
.. code-block:: python

    >>> import bio2bel_phosphosite
    >>> phosphosite_manager = bio2bel_phosphosite.Manager()
    >>> phosphosite_manager.populate()

Command Line Utility
~~~~~~~~~~~~~~~~~~~~
.. code-block:: bash

    bio2bel_phosphosite populate

.. |build| image:: https://travis-ci.org/bio2bel/phosphosite.svg?branch=master
    :target: https://travis-ci.org/bio2bel/phosphosite
    :alt: Build Status

.. |coverage| image:: https://codecov.io/gh/bio2bel/phosphosite/coverage.svg?branch=master
    :target: https://codecov.io/gh/bio2bel/phosphosite?branch=master
    :alt: Coverage Status

.. |documentation| image:: https://readthedocs.org/projects/phosphosite/badge/?version=latest
    :target: http://phosphosite.readthedocs.io
    :alt: Documentation Status

.. |climate| image:: https://codeclimate.com/github/bio2bel/phosphosite/badges/gpa.svg
    :target: https://codeclimate.com/github/bio2bel/phosphosite
    :alt: Code Climate

.. |python_versions| image:: https://img.shields.io/pypi/pyversions/bio2bel_phosphosite.svg
    :alt: Stable Supported Python Versions

.. |pypi_version| image:: https://img.shields.io/pypi/v/bio2bel_phosphosite.svg
    :alt: Current version on PyPI

.. |pypi_license| image:: https://img.shields.io/pypi/l/bio2bel_phosphosite.svg
    :alt: MIT License
