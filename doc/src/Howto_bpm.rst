Bonded particle models
===============

Bonded particle models are used to simulate mesoscale solids.
Solids are constructed as a collection of particles which each
represent a coarse-grained region of space much larger than the
atomistic scale. Particles within a solid region are then connected 
by a network of bonds to provide solid elasticity. 

Unlike traditional bonds in molecular dynamics, the equilibrium 
bond length can vary between bonds. Bonds store the reference state. 
This includes setting the equilibrium length equal to the initial 
distance between the two particles but can also include data on the 
bond orientation for rotational models. This produces a stress free 
initial state. Furthermore, bonds are allowed to break under large 
strains producing fracture.

Bonds can be created using a :doc:`read data <read_data>` 
or :doc:`create bond <create_bond>` command. Alternatively, a 
:doc:`molecule <molecule>` template with bonds can be used with 
:doc:`fix deposit <fix_deposit>` or :doc:`fix pour <fix_pour>` to
create solid grains.
In this implementation, bonds store their reference state when they 
are first computed in the setup of a simulation run. Data is then
preserved across run commands and is written to :doc:`binary restart files <restart>` 
such that restarting the system will not reset the reference state of a bond.

As bonds can be broken between neighbor list builds, :doc:`special_bonds <special_bonds>` 
are incompatible with bonded particle models. Therefore, special bonds
must be turned off using an :doc:`atom modify <atom_modify>` command.
This has the extra benefit of saving memory as these systems often have
very dense networks of bonds.

Currently there are two types of bonds included in this package. The first
bond style, :doc:`bond bpm/simple <bond_bpm_simple>`, only applies pairwise, 
central body forces. Point particles must have :doc:`bond atom style <atom_style>`
and may be thought of as nodes in a spring network. Alternatively,
the second bond style, :doc:`bond bpm/rotational <bond_bpm_rotational>`, 
resolves tangential forces and torques arising with the shearing, bending, 
and twisting of the bond due to rotation or displacement of particles.
Particles are similar to those used in the :doc:`granular package <Howto_granular>`,
:doc:`atom style sphere <atom_style>`. However, they must also track the
current orientatio. of particles and therefore use a derived :doc:`atom style sphere/bpm <atom_style>`.
This also requires a unique integrator :doc:`fix nve/sphere/bpm <fix_nve_sphere_bpm>`
which numerically integrates orientation similar to :doc:`fix nve/asphere <fix_nve_asphere>`.

Unlike :doc:`bond quartic <bond_quartic>`, the pairwise interaction is 
not subtracted during the bond computation. Certain pair styles included 
in this package are designed to skip force computation for bonded particles.
This includes :doc:`pair gran/hertz/history/bpm <pair_gran>` for models with 
rotational degrees of freedom and :doc:`pair bpm/simple <pair_bpm_simple>` 
for models with point particles. Alternatively, any other pair style can be
used but forces will just be overlaid.

In these simulations, bonds are often the main focus. Therefore, several
extra computes and fixes are provided to monitor bonds in the system.
First is :doc:`compute nbond/atom <compute_nbond_atom>`, this compute
simply tallies the current number of bonds per atom. Second is 
:doc:`fix bond/broken <fix_bond_broken>` which records instances of bond 
breakage for output.

