.. index:: improper_style none

improper_style none command
===========================

Syntax
""""""

.. code-block:: LAMMPS

   improper_style none

Examples
""""""""

.. code-block:: LAMMPS

   improper_style none

Used in input scripts:

  .. parsed-literal::

       examples/rdf-adf/in.spce
       examples/VISCOSITY/in.cos.1000SPCE
       examples/PACKAGES/agni/in.adatom
       examples/PACKAGES/agni/in.vacancy

Description
"""""""""""

Using an improper style of none means improper forces and energies are
not computed, even if quadruplets of improper atoms were listed in the
data file read by the :doc:`read_data <read_data>` command.

See the :doc:`improper_style zero <improper_zero>` command for a way to
calculate improper statistics, but compute no improper interactions.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`improper_style zero <improper_zero>`

Default
"""""""

none
