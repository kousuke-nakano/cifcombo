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

Synthesized compositions search from cif database.

Installation via PyPI
----------------------------------------------------------

[In progress] cifcomb can be obtained from PyPI

.. code-block:: console

    pip install cifcomb


Installation from source
----------------------------------------------------------

If you prefer to install from source,
instead follow the procedures below.

.. code-block:: console

    git clone https://github.com/kousuke-nakano/cifcomb.git
    cd cifcomb
    pip install -e .

Quick use
----------------------------------------------------------

Preparation of a CIF data file (from CIF files)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You should prepare a CIF data file from your cif files

.. code-block::

    cifcomb -m cif_dir1 cif_dir2 ....

The generation process is parallelized using `joblib` library. Nevertheless, it takes a very long time if the specified CIF database is huge (e.g., It takes ~ 2 days for the entire COD/Crystallography Open Database [http://www.crystallography.net/cod/] CIF files with 256 cores of Intel(R) Xeon Phi(TM) CPU 7210 @ 1.30GHz).

If one wants to use a prepared CIF database, plz. contact to the developer Kosuke Nakano [kousuke_1123@icloud.com].

Search a CIF file from the generated database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You should search a CIF that can be sum of the input compositions from the generated database.

.. code-block::

    %cifcomb -s CaO TiO2
    1*Ca1 Ti1 O3 = 1*Ca1O1 + 1*Ti1O2, cifid=1000022, SG=Pnma
    ...

Get a CIF file from the generated database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You should copy a CIF with an ID in the generated database.

.. code-block::

    cifcomb -g 1000022
    # copying 1000022.cif to the current directory.

Decompose an input composition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can decompose a target composition with the input compositions

.. code-block::

    %cifcomb -s CaO TiO2 -t CaTiO3
    1*CaTiO3 = 1*CaO + 1*TiO2


Additional information
^^^^^^^^^^^^^^^^^^^^^^

For additional information, you can use the help command:

.. code-block:: console

    cifcomb -h

or you can refer to the documentation.


How to release
----------------------------------------------------------

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
