README
==========

|license| |DL| |release| |PYPI_version| |Python_version| |workflows| |fork| |stars|

.. |license| image:: https://img.shields.io/github/license/kousuke-nakano/cifcomb
.. |release| image:: https://img.shields.io/github/release/kousuke-nakano/cifcomb/all.svg
.. |DL| image:: https://img.shields.io/pypi/dm/cifcomb
.. |Python_version| image:: https://img.shields.io/pypi/pyversions/cifcomb
.. |fork| image:: https://img.shields.io/github/forks/kousuke-nakano/cifcomb?style=social
.. |stars| image:: https://img.shields.io/github/stars/kousuke-nakano/cifcomb?style=social
.. |workflows| image:: https://github.com/kousuke-nakano/cifcomb/actions/workflows/cifcomb-pytest.yml/badge.svg
.. |PyPI_version| image:: https://badge.fury.io/py/cifcomb.svg

Synthesized compounds search from cif database (cifcomb).

Installation
------------

cifcomb can be obtained from PyPI

.. code-block:: console

    pip install cifcomb


Development
^^^^^^^^^^^

If you prefer to install from source,
instead follow the procedures below.

.. code-block:: console

    git clone https://github.com/kousuke-nakano/cifcomb.git
    cd cifcomb
    pip install -e .

Quick use
---------

Preparation of a CIF data file (from CIF files)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You should prepare a CIF data file from your cif files

.. code-block::

    cifcomb -m cif_dir1 cif_dir2 ....

Search a CIF file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You should search a CIF in the generated database.

.. code-block::

    cifcomb -c MgO TiO2 ....

Get a CIF file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You should copy a CIF with an ID in the generated database.

.. code-block::

    cifcomb -g 102345


Additional information
^^^^^^^^^^^^^^^^^^^^^^

For additional information, you can use the help command:

.. code-block:: console

    cifcomb -h

or you can refer to the documentation.


How to release
--------------

Work on the devel or on a new branch

.. code-block:: console
    
    git merge <new branch> devel # if you work on a new branch.
    git push origin devel # A GitHub Action triggers pytests.

Check the next-version version

.. code-block:: console

    # Confirm the version number via `setuptools-scm`
    python -m setuptools_scm
    e.g., 1.1.4.dev28+gceef293.d20221123 -> <next-version> = v1.1.4 or v1.1.4-alpha(for pre-release)

Add and push with the new tag

.. code-block:: console

    # Push with tag
    git tag <next-version>  # e.g., git tag v1.1.4  # Do not forget "v" before the version number!
    git push origin devel --tags  # or to the new branch

Send a pull request to the master branch on GitHub. After the pull request is approved and the devel branch is merged to the master, a GitHub Action checks if the automatic deploy works using test-pyPI (if the commit is tagged correctly, e.g., v1.1.0).

Finally, do a new release with a release note on GitHub. The new release trigggers an implemented GitHub Action that automatically uploads the package to PyPI (if the commit is tagged correctly, e.g., v1.1.0).

Contributing to the project
---------------------------

If you want to contribute to the project, report a bug, or ask for
a new feature, please `raise an issue <https://github.com/kousuke-nakano/cifcomb/issues>`_.
