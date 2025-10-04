Granular models
===============

Granular system are composed of finite size particles with rotational
degrees of freedom as opposed to point particles. This means they have
an angular velocity and torque can be imparted to them to cause them
to rotate as is typical in the discrete element method (DEM). Typically,
these are spheres (or disks in 2D) with a defined diameter but aspherical
bodies are also a possibility. This page describes various methods available
to simulate granular materials and relevant computes and fixes.

----------

For most granular simulations, atom style :doc:`sphere <atom_style>` is
required as it includes both a per-atom diameter and density. To
integrate both the translational and rotational degrees of freedom,
`fix nve/sphere <fix_nve_sphere>` is typically required (unless using
a rigid set of overlapping spheres as described below).

To calculate contact forces and torques, a four pair styles are available

* :doc:`pair_style gran/history <pair_gran>`
* :doc:`pair_style gran/no_history <pair_gran>`
* :doc:`pair_style gran/hertzian <pair_gran>`
* :doc:`pair_style granular <pair_granular>`

The *gran/X* pair styles were the first DEM pair styles implemented and
are largely superseded by the *granular* pair style which includes a
wider range of contact models (e.g. rolling/twisting friction, cohesive
contacts, and heat transfer).

To model heat conduction, one must add the temperature and heatflow
atom variables with:

* :doc:`fix property/atom <fix_property_atom>`

a temperature integration fix

* :doc:`fix heat/flow <fix_heat_flow>`.

----------

There are

* :doc:`fix gravity <fix_gravity>`
* :doc:`fix damping/cundall <fix_damping_cundall>`
* :doc:`fix freeze <fix_freeze>`
* :doc:`fix pour <fix_pour>`

This compute

* :doc:`compute erotate/sphere <compute_erotate_sphere>`

calculates rotational kinetic energy which can be :doc:`output with thermodynamic info <Howto_output>`.
The compute

* :doc:`compute fabric <compute_fabric>`

calculates various versions of the fabric tensor for granular and non-granular pair styles.


These commands implement fix options specific to granular systems:

* :doc:`fix wall/gran <fix_wall_gran>`
* :doc:`fix wall/gran/region <fix_wall_gran_region>`

The fix style *freeze* zeroes both the force and torque of frozen
atoms, and should be used for granular system instead of the fix style
*setforce*\ .

and a heat conduction option defined in both

* :doc:`pair_style granular <pair_granular>`
* :doc:`fix wall/gran <fix_wall_gran>`

For computational efficiency, you can eliminate needless pairwise
computations between frozen atoms by using this command:

* :doc:`neigh_modify <neigh_modify>` exclude

.. note::

   By default, for 2d systems, granular particles are still modeled
   as 3d spheres, not 2d discs (circles), meaning their moment of inertia
   will be the same as in 3d.  If you wish to model granular particles in
   2d as 2d discs, see the note on this topic on the :doc:`Howto 2d <Howto_2d>`
   doc page, where 2d simulations are discussed.

To add custom granular contact models, see the
:doc:`modifying granular sub-models page <Modify_gran_sub_mod>`.
