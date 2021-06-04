.. index:: bond_style bpm/rotational

bond_style bpm/rotational command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   bond_style bpm/rotational

Examples
""""""""

.. code-block:: LAMMPS

   bond_style bpm/rotational
   bond_coeff 1 

Description
"""""""""""

The *bpm/rotational* bond style computes forces and torques based 
on deviations from the initial reference state of the two atoms.
The reference state is stored by each bond when it is first computed
in the setup of a run. Data is then preserved across run commands and
is written to :doc:`binary restart files <restart>` such that restarting
the system will not reset the reference state of a bond.

Forces include a normal and tangential component. The base normal force
has a magnitude of

.. math::

   F_r = K_r (r - r_0)
   
where :math:`K_r` is a stiffness and :math:`r` is the current distance and 
:math:`r_0` is the initial distance between the two particles. In 
compression, a multiplicative stiffining term can be added such that

.. math::

   F_r = K_r (r - r_0) e^{c (r0/r - 1.0)}
   
where :math:`c` controls the strength of the stiffening.

A tangential force is applied perpendicular to the normal direction
which is proportional to the tangential shear displacement with a stiffness
of :math:`K_s`. This tangential force also induces a torque.
In addition, bending and twisting torques are also applied to particles
which are proportional to angular bending and twisting displacements with 
stiffnesses of :math`K_b` and :math:`K_t', respectively.
Details on the calculations of shear displacements and angular displacements 
can be found in :ref:`(Wang) <Wang>` and :ref:`(Wang and Mora) <WangMora>`.

Bonds will break under sufficient stress. A breaking criteria is calculated

.. math::

   B = \alpha(r, r_0) \frac{F_r}{F_{r,c}} + \frac{F_s}{F_{s,c}} + 
       \frac{\tau_b}{\tau_{b,c}} + \frac{\tau_t}{\tau_{t,c}}
   
where :math:`F_s` is the magnitude of the shear force and 
:math:`\tau_b` and :math:`\tau_t` are the magnitudes of the bending and
twisting forces, respectively. The corresponding variables :math:`F_{r,c}`
:math:`F_{s,c}`, :math:`\tau_{b,c}`, and :math:`\tau_{t,c}` are critical
limits to each force or torque. The term :math:`\alpha` is simply one in
extension and zero in compression such that the normal force component 
does not contribute to the breaking criteria in compression.
If :math:`B` is ever equal to or exceeds one, the bond will break. 
This is done by setting by setting its type to 0 such that forces and 
torques are no longer computed.

After computing the base magnitudes of the forces and torques, they are 
all multiplied by an extra factor :math:`w` to smoothly interpolate 
forces and torques to zero as the bond breaks. This term is calculated 
as :math:`w = (1.0 - B)^2`.

Finally, additional damping forces and torques are applied to the two
particles. A force is applied proportional to the difference in the
normal velocity of particles using a similar construction as 
dissipative particle dynamics (:ref:`(Groot) <Groot1>`):

.. math::

   F_D = - \gamma_v w (\hat{r} \bullet \vec{v})
   
where :math:`\gamma_v` is the damping strength, :math:`\hat{r}` is the 
radial normal vector, and :math:`\vec{v}` is the velocity difference
between the two particles. Torques are also applied to damp relative
angular velocties between the two particlces. This includes damping 
the relative twisting velocity and the relative tangential angular velocity.
Both torques are the product of :math:`\gamma_\omega w` and the respective
components of the angular velocity.

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`K_r`           (force/distance units)
* :math:`K_s`           (force/distance units)
* :math:`K_t`           (force units)
* :math:`K_b`           (force units)
* :math:`F_{r,c}`       (force units)
* :math:`F_{s,c}`       (force units)
* :math:`\tau_{b,c}`    (force*distance units)
* :math:`\tau_{t,c}`    (force*distance units)
* :math:`\gamma_v`      (force/velocity units)
* :math:`\gamma_\omega` (distance*force/seconds/radians units)
* :math:`c`             (dimensionless)

As bonds can be broken between neighbor list builds, :doc:`special_bonds <special_bonds>` 
will not work with this bond style. Unlike :doc:`bond quartic <bond_quartic>`,
the pairwise interaction is not subtracted during the bond computation.
Certain pair styles, namely :doc:`pair gran/hertz/history/bpm <pair_gran>`,
are designed to skip force computation for bonded particles. Alternatively,
pair forces can be overlaid.

Note that when bonds are dumped to a file via the :doc:`dump local <dump>` command, bonds with type 0 are not included.  The
:doc:`delete_bonds <delete_bonds>` command can also be used to query the
status of broken bonds or permanently delete them, e.g.:

.. code-block:: LAMMPS

   delete_bonds all stats
   delete_bonds all bond 0 remove


----------

Restart
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This bond style writes the reference state of each bond to 
:doc:`binary restart files <restart>`. Loading a restart
file will properly resume bonds.

Restrictions
""""""""""""

This bond style can only be used if LAMMPS was built with the BPM
package.  See the :doc:`Build package <Build_package>` doc page for more
info.

The *bpm/rotational* style requires :doc:`special_bonds <special_bonds>`
be turned off using the :doc:`special_bonds <special_bonds>` command. 

The *bpm/rotational* style requires :doc:`atom style sphere/bpm <atom_style>`.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`fix broken/bonds <fix_broken_bonds>`,
:doc:`fix nve/sphere/bpm <fix_nve_sphere_bpm>`

Default
"""""""

none


.. _Wang:

**(Wang)** Wang, Acta Geotechnica, 4,
p 117-127 (2009).

.. _WangMora:

**(Wang and Mora)** Wang, Mora, Advances in Geocomputing, 
119, p 183-228 (2009).

.. _Groot1:

**(Groot)** Groot and Warren, J Chem Phys, 107, 4423-35 (1997).
