.. index:: fix broken/bonds

fix broken/bonds command
===================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID broken/bonds N keyword attribute1 attribute2 ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* broken/bonds = style name of this fix command
* N = prepare data for output every this many timesteps
* one or more attributes may be appended

  .. parsed-literal::

       possible attributes = id1 id2 time x y z xstore ystore zstore

  .. parsed-literal::

          id1, id2 = IDs of 2 atoms in the bond
          time = the time the bond broke
          x, y, z = the center of mass position of the 2 atoms when the bond broke
          x, y, z = the inintial center of mass position of the 2 atoms

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all broken/bonds 1000 id1 id2 time
   dump 1 all local 100 broken_bonds.dat f_1[1] f_1[2] f_2[3]

Description
"""""""""""

This compute tracks bonds that are broken bonds for bond styles in
the :doc:`BPM modle <HowTo_bpm>` and records data. Must be used in 
conjunction with either the :doc:`bpm/rotational <bond_bpm_rotational>` 
or :doc:`bpm/simple <bond_bpm_simple>` bond styles. 
Data is accumulated over a span of *N* timesteps after which it is cleared
The number of datums generated, aggregated across all processors, equals 
the number of broken bonds. bonds are only included if both
atoms are in the included in the specified fix group. 

----------

Restart, fix_modify, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.  
None of the :doc:`fix_modify <fix_modify>` options are
relevant to this fix.  
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

Output info
"""""""""""

This compute calculates a local vector or local array depending on the
number of input values.  The length of the vector or number of rows in
the array is the number of broken bonds.  If a single input is
specified, a local vector is produced.  If two or more inputs are
specified, a local array is produced where the number of columns = the
number of inputs.  The vector or array can be accessed by any command
that uses local values from a compute as input.  See the :doc:`Howto output <Howto_output>` doc page for an overview of LAMMPS output
options.

The vector or array values will be doubles that correspond to the
specified attribute.

Restrictions
""""""""""""

This fix can only be used if LAMMPS was built with the BPM
package.  See the :doc:`Build package <Build_package>` doc page for more
info.

Must be used in conjuction with either the :doc:`bpm/rotational <bond_bpm_rotational>` 
or :doc:`bpm/simple <bond_bpm_simple>` bond styles. 

Related commands
""""""""""""""""

:doc:`bond_style bpm/rotational <bond_bpm_rotational>` and 
:doc:`bond_style bpm/simple <bond_bpm_simple>`

Default
"""""""

none
