.. _sec-hawt:

Model bodies and constraints for a HAWT
---------------------------------------

In this section we describe our approach to modeling a horiztonal axis wind turbine (HAWT), including
the specific "bodies" that include massless nodes with six DOF (position
and orientation), 6-DOF rigid bodies, rigid connections, and flexible
beams. In the model descirbed here, the tower and blades are flexible
beams, with all other components being rigid.
Table `1 <#table:hawt-cs>`__ defines the subscripts associated with
coordinate systems used to define the model, which are illustrated in
Figure `1 <#fig:geom>`__. For example, the tower top is defined by
:math:`\underline{x}^\mathrm{r}_\mathrm{tt}`, :math:`\underline{u}_\mathrm{tt}`,
:math:`\underline{\underline{R}}^\mathrm{r}_\mathrm{tt}`, and
:math:`\underline{\underline{R}}_\mathrm{tt}`, which is the reference
position, displacement, reference orientation, and relative rotation,
respectively. The current orientation, denoted by a c-superscript, is
:math:`\underline{\underline{R}}^\mathrm{c}_\mathrm{tt} = \underline{\underline{R}}_\mathrm{tt} \underline{\underline{R}}^\mathrm{r}_\mathrm{tt}`.

         +--------------+-------------------------------------------------------+
         | subscript    | description                                           |
         +==============+=======================================================+
         | g            | global/inertial                                       |
         +--------------+-------------------------------------------------------+
         | tb           | tower base; first node in tower finite element        |
         +--------------+-------------------------------------------------------+
         | tt           | tower top; last node in tower finite element          |
         +--------------+-------------------------------------------------------+
         | yb           | yaw-base; bearing center of mass                      |
         +--------------+-------------------------------------------------------+
         | sb           | shaft base (does not rotate with rotor)               |
         +--------------+-------------------------------------------------------+
         | a            | azimuth (rotates with rotor)                          |
         +--------------+-------------------------------------------------------+
         | h            | hub center of mass (rotates with rotor)               |
         +--------------+-------------------------------------------------------+
         | c\ :math:`i` | coned CS for blade :math:`i` (rotates with rotor;     |
         |              | does not pitch)                                       |
         +--------------+-------------------------------------------------------+
         | b\ :math:`i` | blade CS for blade :math:`i` (rotates with rotor;     |
         |              | pitches)                                              |
         +--------------+-------------------------------------------------------+
         | nm           | Nacelle center of mass                                |
         +--------------+-------------------------------------------------------+

.. figure:: images/hawt.png
   :alt: A plot of data
   :width: 50%
   :align: center

   Schematic showing the coordinate systems and elements for a horizontal axis wind turbine.  Only one blade is shown.


Constraints
~~~~~~~~~~~

In this section we define the constraints for our HAWT model. We start
at the base and move up through the connections. We assume the tower
base is fixed, for which the constraints can be written

.. math:: 
   
   \underline{\Phi}_\mathrm{tb} = 
   \underline{\Phi}^\mathrm{BC6}\left(\underline{u}_\mathrm{tb}, \underline{\underline{R}}_\mathrm{tb}, \underline{0}, 
   \underline{\underline{I}} \right) \in \mathbb{R}^6

The tower is represented by a flexible beam.
The constraint between the tower top and the yaw base is complicated by
the yaw controller and the motion of the beam model for the tower top.
We apply the yaw controller to the yaw base plate in the reference
configuration through the Rigid66 constraint:

Keeping the yaw displacement tied to tower top displacement, the
constraints between the tower-top and yaw baseplate are

.. math::

   \underline{\Phi}_\mathrm{yb-tt} = \underline{\phi}^\mathrm{Rigid66}\left(
   \underline{u}_\mathrm{yb}, \underline{\underline{R}}_\mathrm{yb}, 
   \underline{u}_\mathrm{tt}, \underline{\underline{R}}_\mathrm{tt}, \underline{\underline{R}}_\mathrm{yc}, 
   \underline{r}^\mathrm{r}_{\mathrm{yb-tt}} \right) \in \mathbb{R}^6 


The shaft base and yaw CS’s are assumed to be rigidly connected. Here,
shaft-base relative rotation is the same as the yaw relative rotation:

.. math::

   \underline{\Phi}_\mathrm{sb-yb} = 
   \underline{\Phi}^\mathrm{Rigid66}\left(
   \underline{u}_\mathrm{sb}, \underline{\underline{R}}_\mathrm{sb}, 
   \underline{u}_\mathrm{yb}, \underline{\underline{R}}_\mathrm{yb}, \underline{\underline{I}},
   \underline{r}^\mathrm{r}_{\mathrm{sb-yb}} \right) \in \mathbb{R}^6

The azimuth CS is free to rotate about the
:math:`\widehat{x}_\mathrm{sb}` axis and is tied to the location of the
shaft-base CS.

.. math::

   \underline{\Phi}_\mathrm{a-sb} = 
   \underline{\Phi}^\mathrm{Revolute}\left(\underline{u}_\mathrm{a}, \underline{\underline{R}}_\mathrm{a},
   \underline{\underline{R}}_\mathrm{a}^\mathrm{r},
   \underline{u}_\mathrm{sb}, \underline{\underline{R}}_\mathrm{sb},  
   \underline{\underline{R}}_\mathrm{sb}^\mathrm{r}
   \right) \in \mathbb{R}^5 


The hub coordinate system (for a rigid shaft) has the same relative
rotation as the azimuth CS

.. math::

   \underline{\Phi}_\mathrm{h-a} = \underline{\Phi}^\mathrm{Rigid66}\left(
   \underline{u}_\mathrm{h}, \underline{\underline{R}}_\mathrm{h},
   \underline{u}_\mathrm{a}, \underline{\underline{R}}_\mathrm{a}, \underline{\underline{I}},
   \underline{r}^\mathrm{r}_{\mathrm{h-a}} \right) \in \mathbb{R}^6


The cone CS for blade :math:`i` rotates with the hub and the cone CS
displacement is equal to that of the hub, but does not pitch.

.. math::

   \underline{\Phi}_{\mathrm{c}i\mathrm{-h}} = \underline{\Phi}^\mathrm{Rigid66}\left(
   \underline{u}_{\mathrm{c}i}, \underline{\underline{R}}_{\mathrm{c}i},
   \underline{u}_\mathrm{h}, \underline{\underline{R}}_\mathrm{h}, \underline{\underline{I}},
   \underline{r}^\mathrm{r}_{\mathrm{c}i-\mathrm{h}} \right) \in \mathbb{R}^6

The blade-root CS for blade :math:`i` is offset from the cone CS and it
pitches.

.. math::

   \underline{\Phi}_{\mathrm{b}i-\mathrm{c}i} = \underline{\Phi}^\mathrm{Rigid66}\left(
   \underline{u}_{\mathrm{b}i}, \underline{\underline{R}}_{\mathrm{b}i},
   \underline{u}_{\mathrm{c}i}, \underline{\underline{R}}_{\mathrm{c}i}, \underline{\underline{R}}_{\mathrm{pc}i},
   \underline{r}^\mathrm{r}_{\mathrm{b}i-\mathrm{c}i} \right) \in \mathbb{R}^6

The nacelle-mass and yaw-base CS’s are assumed to be rigidly connected.

.. math::

   \underline{\Phi}_\mathrm{nm-yb} = \underline{\Phi}^\mathrm{Rigid66}\left(
   \underline{u}_\mathrm{nm}, \underline{\underline{R}}_\mathrm{nm},
   \underline{u}_\mathrm{yb}, \underline{\underline{R}}_\mathrm{yb}, \underline{\underline{I}},
   \underline{r}^\mathrm{r}_{\mathrm{nm}-\mathrm{yb}} \right) \in \mathbb{R}^6

Constraint Gradient
~~~~~~~~~~~~~~~~~~~

Constraint gradient matrices are described in the Constraint Library section of the Theory documentation.

Data layout
~~~~~~~~~~~

Tangent stiffness matrix

.. figure:: images/tangent.png
   :alt: A plot of data layout
   :width: 50%
   :align: center

   Schematic illustrating the nonzero entries of the tangent stiffness matriix.

.. figure:: images/constraint-grad.png
   :alt: A plot of data layout
   :width: 50%
   :align: center

   Schematic illustrating the nonzero entries of the constraint-gradient matriix.


