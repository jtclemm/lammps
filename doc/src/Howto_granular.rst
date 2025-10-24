
Granular models
===============

Granular systems are composed of finite size particles with rotational
degrees of freedom, unlike atomic/point particles. This means an applied
torque can induce rotation at an angular velocity. This is commonly known as
the discrete element method or DEM. Typically, these particles are represented
as a single atom in LAMMPS and are spherical (or circular disks in 2D) with a
defined diameter but aspherical particles can also be modeled. This page
describes various methods available in the GRANULAR package to simulate
granular materials as well as other relevant computes and fixes.

----------

.. contents:: Outline of topics

----------


.. _basics

Package basics - atom and pair styles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In most granular simulations, atom style :doc:`sphere <atom_style>` is
required as it includes both a per-atom diameter and density. To
integrate both the translational and rotational degrees of freedom,
`fix nve/sphere <fix_nve_sphere>` is commonly required. Exceptions include
simulations of rigid overlapping spheres, bonded particles, or some barostats
as described below.

.. note::

   By default, for 2d systems, granular particles are still modeled
   as 3d spheres, not 2d discs (circles), meaning their moment of inertia
   will be the same as in 3d.  If you wish to model granular particles in
   2d as 2d discs, see the note on this topic on the :doc:`Howto 2d <Howto_2d>`
   doc page, where 2d simulations are discussed.

To calculate contact forces and torques, four pair styles are available

* :doc:`pair_style gran/history <pair_gran>`
* :doc:`pair_style gran/no_history <pair_gran>`
* :doc:`pair_style gran/hertzian <pair_gran>`
* :doc:`pair_style granular <pair_granular>`

The *gran/X* pair styles were the first DEM pair styles implemented and
are largely superseded by the *granular* pair style which allows one to
pick a wider ranger of contact models. More detail is found in the doc
page of pair *granular*, but essentially a contact model is broken into
several submodels that describe some component of the interaction: normal
forces, tangential/rolling/twisting friction, damping, and heat conduction.
Each submodel style can be mixed and matched. New submodel styles can be
modularly added as described in the :doc:`modifying granular sub-models
<Modify_gran_sub_mod>` page.

.. _fixes

Fixes to control and drive systems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are many ways to drive and control a granular system through
various fixes. Here is a brief, although certainly not comprehensive,
list of some fixes that are more relevant to granular systems.

First, some fixes can apply various types of forces directly to
particles. This includes general fixes such as :doc:`fix addforce
<fix_addforce>` and :doc:`fix setforce <fix_setforce>` which can
apply or constrain forces (but not torques) on individual particles.
Two more granular-specific alternatives are :doc:`fix gravity <fix_gravity>`
and :doc:`fix freeze <fix_freeze>`. Akin to addforce, fix gravity
will apply additional forces to each particle but it accounts for
the mass of a particle. Like fix setforce, fix freeze can be used
to remove forces on particles and freeze them in space. However,
fix freeze not only removes translational but also rotational
velocities. This can be useful in freezing particles that make up
walls or boundary conditions. Note that the :doc:`create_atoms
<create_atoms>` command has the ability to read an STL file to
place atoms along an arbitrary surface mesh. For computational
efficiency, you can eliminate needless pairwise computations
between frozen atoms using :doc:`neigh_modify <neigh_modify>` exclude

Another method to create walls/boundaries is :doc:`fix wall/gran
<fix_wall_gran>` which generates flat planar walls. In a similar
fashion to pair granular, one can specify a wide range of contact
models for this wall. Another option is :doc:`fix wall/gran/region
<fix_wall_gran_region>` which creates walls using the :doc:`region
<region>` command which supports many geometric primitive shapes.
Note that these commands have options to translate or rotate walls
in time.

For simulations with periodic boundary conditions, there are multiple
fixes that apply global deformation, either by specifying strain and/or
a target stress tensor. :doc:`Fix deform <fix_deform>` is the standard
fix for applying a defined strain or strain rate on the simulation cell.
This is generalized in :doc:`fix deform/pressure <fix_deform_pressure>`
which allows one to specify both strain- or stress-based deformation
conditions or mixtures of the two. When specifying stress-driven
deformation, a linear controller is used to move the system towards the
target stress state. Other barostats common in Molecular Dynamics are
also available including a granular-specific :doc:`Nose-Hoover NPH
barostat <fix_nph_sphere>` and a :doc:`Berendsen barostat <fix_berendsen>`.

To apply drag to particles, one might consider using
* :doc:`fix damping/cundall <fix_damping_cundall>`
* :doc:`fix viscous <fix_viscous>`
* :doc:`fix langevin <fix_langevin>`
Fix damping/Cundall applies velocity-independent forces and torques
to remove kinetic energy from the system. Fix viscous and fix
langevin can be used to apply equivalent velocity-dependent damping
forces, although fix langevin also includes an optional thermostat
and optional damping torques.

:doc:`Fix pour <fix_pour>` is useful fix for pouring non-overlapping
atoms into a system. It is designed to account for the effects of fix
gravity and a gradual fill rate such that atoms are deposited in a
gradually rising zone.

.. _thermal

Thermal evolution
^^^^^^^^^^^^^^^^^

By default, atoms do not store thermal quantities and thermal
evolution is not integrated. To enable this feature, one must
first add the temperature and heatflow atom variables to atoms
using :doc:`fix property/atom <fix_property_atom>`. Then,
:doc:`fix heat/flow <fix_heat_flow>` must be added to integrate
temperature. Optional conduction between particles is enabled
using the heat/flow options in the granular pair style.

To control the heat flow into/out of the system, the granular
wall fixes have an option to set a boundary temperature and
:doc:`fix add/heat <fix_add_heat>` directly adds/removes heat
to particles in the defined group.

.. _computes

Granular specific computes
^^^^^^^^^^^^^^^^^^^^^^^^^^

LAMMPS includes many computes which are helpful in a wide
range of scenarios and systems. In this section, a few of
the more granular-specific computes are listed and described.
However, it is recommended that new users still browse the
list of all available computes in LAMMPS as many of these are
still very useful in granular simultations.

Typically, when calculating the temperature or kinetic
energy of particles, most computes only account for
translational degrees of freedom. To measure rotational
contributions, there are the computes:
* :doc:`compute erotate/sphere <compute_erotate_sphere>`
* :doc:`compute temp/sphere <compute_temp_sphere>`

To measure the statistics of contacts, several computes provide
granular-specific support and account for the finite size of
particles. These include
* :doc:`compute contact/atom <compute_contact_atom>` which calculates particle coordination number
* :doc:`compute rattlers <compute_nonaffine_rattlers>` which identifies undercoordinated rattlers in the system
* :doc:`compute fabric <compute_fabric>` which calculates the fabric tensor from the network of contacts

Lastly, :doc:`compute nonaffine/displacement
<compute_nonaffine_displacement>` calculates a local strain tensor
and nonaffine deviations from the tensor from the change in contact
networks over time. This compute can calculate D2min as well.

.. _advanced

Advanced features and options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To efficiently simulate high levels of polydispersity, LAMMPS includes
improved algorithms for building the neighborlist. This can be activated
using the *multi* option in the :doc:`neighbor command <neighbor>` and
the associated information on its documentation page.

To simulate aspherical granular particles, there are two primary options.
The first is the :doc:`rigid package <fix_rigid>` which combines an
assortment of spherical particles into a single rigid body. This method
is also known as the multi-sphere approach in the literature. The second
is the :doc:`Bonded Particle Model (BPM) package <Howto_BPM>` which can
connect an assortment of spherical particles with cohesive bonds to
create an elastic body with (optional) fracture or plasticity. This package
is capable of reproducing the cohesive beam model or bonded discrete
element method. To determine which approach is better suited for an
application, it is recommended that one consider the importance of
absolute rigidity versus elasticity and benchmark the two options. The
relative costs can depend on many factors including the size of a solid
body relative to a processor's domain (fix rigid performance degrades
with very large grains), whether bodies are hollow or solid, and the
complexity of pair force calculations relative to bond force calculations.

In simulations with significant spatial variations in particle density,
:doc:`fix balance <fix_balance>` can be used to adjust the domain decomposition
in an attempt to balance the number of particles, and therefore the
computational load, on each processor. For instance in a simulation of a
conical hopper, load balancing might increase the domain size for processors
assigned to the bottom of the cone to avoid idle time.
