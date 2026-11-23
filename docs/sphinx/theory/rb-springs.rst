.. _sec-rb-springs:

Rigid body with three springs
-----------------------------

This section provides model details for a rigid-body solver with three
massless taut-line "mooring lines", represented by geometrically
nonlinear springs. The model configuration is appropriate for modeling
the DeepCWind floating platform.

The generalized model degrees of freedom are

.. math::

   \underline{q} = \begin{bmatrix} 
     \underline{u}_\mathrm{RB} \\
     \underline{\underline{R}}_\mathrm{RB} \\
     \underline{u}_\mathrm{sp1} \\
     \underline{u}_\mathrm{sp2} \\
     \underline{u}_\mathrm{sp3}
    \end{bmatrix} \quad
    \underline{v} = \begin{bmatrix} 
     \dot{\underline{u}}_\mathrm{RB} \\
     \underline{\omega}_\mathrm{RB} \\
     \dot{\underline{u}}_\mathrm{sp1} \\
     \dot{\underline{u}}_\mathrm{sp2} \\
     \dot{\underline{u}}_\mathrm{sp3}
    \end{bmatrix}

where :math:`\underline{q} \in \mathbb{R}^3 \times \mathrm{SO}(3) \times \mathbb{R}^{18}` and :math:`\underline{v} \in \mathbb{R}^{24}`,
:math:`\underline{u}_\mathrm{RB} \in \mathbb{R}^3` and
:math:`\underline{\underline{R}}_\mathrm{RB} \in \mathrm{SO}(3)` are rigid
body displacement and orientation, respectively,
:math:`\underline{u}_{\mathrm{sp}i} \in \mathbb{R}^6` is spring i end
displacement (2 nodes with 3 DOF each; see :ref:`sec-spring`) and

.. math::

   \underline{u}_{\mathrm{sp}i} = 
   \begin{bmatrix}
   \underline{u}_{\mathrm{sp}i1} \\ \underline{u}_{\mathrm{sp}i2} 
   \end{bmatrix}

The system input parameters include the following:

- :math:`\underline{\underline{M}}^*_\mathrm{RB} \in \mathbb{R}^{6\times 6}`: Rigid
  body mass matrix in material coordinate system

- :math:`\underline{x}^\mathrm{r}_\mathrm{RB} \in \mathbb{R}^3`: Rigid body reference
  position

- :math:`\underline{\underline{R}}^\mathrm{r}_\mathrm{RB} \in \mathrm{SO}(3)`: Rigid
  body reference orientation

- :math:`\underline{u}^\mathrm{init}_\mathrm{RB} \in \mathbb{R}^3`: Rigid body
  initial displacement

- :math:`\underline{\underline{R}}^\mathrm{init}_\mathrm{RB} \in \mathrm{SO}(3)`:
  Rigid body initial orientation

- :math:`\dot{\underline{u}}^\mathrm{init}_\mathrm{RB} \in \mathbb{R}^3`: Rigid body
  initial velocity

- :math:`\underline{\omega}^\mathrm{init}_\mathrm{RB} \in \mathbb{R}^3`: Rigid body
  initial angular velocity

- :math:`\underline{x}_{\mathrm{sp,r}i} \in \mathbb{R}^6`: Spring
  :math:`i` reference location

- :math:`k_{\mathrm{sp}i} \in \mathbb{R}`: Spring :math:`i` spring
  constant

- :math:`L_{\mathrm{sp,r}i} \in \mathbb{R}`: Spring :math:`i` unstretched
  length

and where

.. math::

   \underline{x}^{\mathrm{sp,r}i} =
   \begin{bmatrix}
   \underline{x}^{\mathrm{sp,r}i}_1 \\ \underline{x}^{\mathrm{sp,r}i}_2
   \end{bmatrix}

The constraints for each spring can be written

.. math::

   \begin{aligned}
    \underline{u}_1^{\mathrm{sp}i} &=  \underline{u} + \underline{\underline{R}} \underline{r}^{\mathrm{sp,r}i}_1
   - \underline{r}^{\mathrm{sp,r}i}  = \underline{u} + \left(\underline{\underline{R}}-\underline{\underline{I}}\right) 
   \underline{r}^{\mathrm{sp,r}i}_1 \\
   \underline{u}_2^{\mathrm{sp}i} &= 0
   \end{aligned}

where :math:`\underline{r}^{\mathrm{sp,r}i}_1 
= \underline{x}^{\mathrm{sp,r}i}_1 - \underline{x}^\mathrm{r}`.

The full-system residual required for time integration can be written as follows (see Eq. :eq:`residual`):

.. math::

   \underline{R} = \begin{bmatrix}
   \underline{\underline{M}}\, \dot{\underline{v}} +\underline{g} - \underline{f} + \underline{\underline{B}}^T\underline{\lambda} \\
   \underline{\Phi}
   \end{bmatrix} 
   \in \mathbb{R}^{42}

where :math:`\underline{\lambda} \in \mathbb{R}^{18}`, :math:`\underline{\underline{B}} \in \mathbb{R}^{18 \times 24}`

.. math::

   \underline{\underline{M}} 
   = \begin{bmatrix}
   \underline{\underline{M}}_\mathrm{RB} & \underline{\underline{0}}_{6\times 18}\\
   \underline{\underline{0}}_{18\times 6} & \underline{\underline{0}}_{18\times 18}
   \end{bmatrix} \in \mathbb{R}^{24\times 24}
   \,,\quad \underline{g} = \begin{bmatrix} 
   \underline{g}_\mathrm{RB} \\ 
   \underline{g}_\mathrm{sp1} \\
   \underline{g}_\mathrm{sp2} \\
   \underline{g}_\mathrm{sp3} 
   \end{bmatrix} \in \mathbb{R}^{24}
   \,,\quad
   \underline{f}
   = 
   \begin{bmatrix}
   \underline{f}_\mathrm{RB}\\
   \underline{0}_{18}\\
   \end{bmatrix} 
   \in \mathbb{R}^{24}

In the above, :math:`\underline{\underline{M}}_\mathrm{RB} \in \mathbb{R}^{6 \times 6}` and :math:`\underline{g}_\mathrm{RB}\in \mathbb{R}^6` 
are defined in :ref:`sec-rigid`, :math:`\underline{f}_\mathrm{RB}\in \mathbb{R}^6` is the force applied to the rigid-body center of mass, 
and :math:`\underline{g}_{\mathrm{sp}i}\in\mathbb{R}^6` is defined in :ref:`sec-spring`. The constraints are given by

.. math::

   \underline{\Phi} =
   \begin{bmatrix}
   \underline{\Phi}^\mathrm{Rigid63}_\mathrm{RB-sp1,1}\left(\underline{u}_\mathrm{RB}, 
   \underline{\underline{R}}_\mathrm{RB}, \underline{u}_\mathrm{sp1,1},\underline{r}^r_\mathrm{RB-sp1,1} \right) \\
   \underline{\Phi}^\mathrm{BC3}_\mathrm{sp1,2}\left(\underline{u}_\mathrm{sp1,2},\underline{0}\right) \\
   \underline{\Phi}^\mathrm{Rigid63}_\mathrm{RB-sp2,1}\left(\underline{u}_\mathrm{RB}, 
   \underline{\underline{R}}_\mathrm{RB}, \underline{u}_\mathrm{sp2,1},\underline{r}^r_\mathrm{RB-sp2,1} \right) \\
   \underline{\Phi}^\mathrm{BC3}_\mathrm{sp2,2}\left(\underline{u}_\mathrm{sp2,2},\underline{0}\right) \\
   \underline{\Phi}^\mathrm{Rigid63}_\mathrm{RB-sp3,1}\left(\underline{u}_\mathrm{RB}, 
   \underline{\underline{R}}_\mathrm{RB}, \underline{u}_\mathrm{sp3,1},\underline{r}^r_\mathrm{RB-sp3,1} \right) \\
   \underline{\Phi}^\mathrm{BC3}_\mathrm{sp3,2}\left(\underline{u}_\mathrm{sp3,2},\underline{0}\right) 
   \end{bmatrix} \quad
   \in \mathbb{R}^{18}

The specific entries of :math:`\underline{\underline{B}}` are given by

.. math::

   \begin{aligned}
   \underline{\underline{B}} =
   \begin{bmatrix}
   & \left[ \underline{\underline{B}}_\mathrm{RB-sp1,1}^\mathrm{Rigid63}  \right]_{(1:3,1:6)}
   & \left[\underline{\underline{B}}_\mathrm{RB-sp1,1}^\mathrm{Rigid63}\right]_{(1:3,7:9)}  & \underline{\underline{0}} & \underline{\underline{0}} 
   & \underline{\underline{0}} & \underline{\underline{0}} & \underline{\underline{0}} \\
   %
   \underline{\underline{0}} &  \underline{\underline{0}}
   & \underline{\underline{0}} & \underline{\underline{B}}^\mathrm{BC3}_\mathrm{sp1,2} & \underline{\underline{0}} 
   & \underline{\underline{0}} & \underline{\underline{0}} & \underline{\underline{0}} \\
   %
   & \left[ \underline{\underline{B}}_\mathrm{RB-sp2,1}^\mathrm{Rigid63}  \right]_{(1:3,1:6)}
   & \underline{\underline{0}}  & \underline{\underline{0}} & \left[\underline{\underline{B}}_\mathrm{RB-sp2,1}^\mathrm{Rigid63}\right]_{(1:3,7:9)} 
   & \underline{\underline{0}} & \underline{\underline{0}} & \underline{\underline{0}} \\
   %
   \underline{\underline{0}} &  \underline{\underline{0}}
   & \underline{\underline{0}} & \underline{\underline{0}} & \underline{\underline{0}} 
   &  \underline{\underline{B}}^\mathrm{BC3}_\mathrm{sp2,2} & \underline{\underline{0}} & \underline{\underline{0}} \\
   %
   & \left[ \underline{\underline{B}}_\mathrm{RB-sp3,1}^\mathrm{Rigid63}  \right]_{(1:3,1:6)}
   & \underline{\underline{0}} & \underline{\underline{0}} & \underline{\underline{0}}
   & \underline{\underline{0}} & \left[\underline{\underline{B}}_\mathrm{RB-sp3,1}^\mathrm{Rigid63}\right]_{(1:3,7:9)}  & \underline{\underline{0}}\\
   %
   \underline{\underline{0}} &  \underline{\underline{0}}
   & \underline{\underline{0}} & \underline{\underline{0}} & \underline{\underline{0}} 
   & \underline{\underline{0}} & \underline{\underline{0}} & \underline{\underline{B}}^\mathrm{BC3}_\mathrm{sp3,2} \\
   \end{bmatrix}
   \end{aligned}

The full matrices :math:`\underline{\underline{C}}`, :math:`\underline{\underline{K}}`, and :math:`\underline{\underline{K}}^\Phi` for the iteration matrix Eq. :eq:`iteration` in :ref:`sec-gen-alpha` are given by

.. math::

   \underline{\underline{C}} = \begin{bmatrix}
   \underline{\underline{C}}^\mathrm{rb} & \underline{\underline{0}}_{6 \times 18}\\
   \underline{\underline{0}}_{18 \times 6} & \underline{\underline{0}}_{18 \times 18}
   \end{bmatrix} 
   \in \mathbb{R}^{24 \times 24}

.. math::

   \underline{\underline{K}} = \begin{bmatrix}
   \underline{\underline{K}}^\mathrm{rb} & \underline{\underline{0}}_{6 \times 6} & \underline{\underline{0}}_{6 \times 6} & \underline{\underline{0}}_{6 \times 6}\\
   \underline{\underline{0}}_{6 \times 6} & \underline{\underline{K}}^\mathrm{sp1} & \underline{\underline{0}}_{6 \times 6} & \underline{\underline{0}}_{6 \times 6}\\
   \underline{\underline{0}}_{6 \times 6} & \underline{\underline{0}}_{6 \times 6} & \underline{\underline{K}}^\mathrm{sp2} & \underline{\underline{0}}_{6 \times 12}\\
   \underline{\underline{0}}_{6 \times 6} & \underline{\underline{0}}_{6 \times 6} &  \underline{\underline{0}}_{6 \times 6} & \underline{\underline{K}}^\mathrm{sp3}
   \end{bmatrix}
   \in \mathbb{R}^{24 \times 24}

.. math::

   \underline{\underline{K}}^\Phi = \begin{bmatrix}
   \underline{\underline{0}}_{3 \times 3} & 
   \underline{\underline{0}}_{3 \times 3} & 
   \underline{\underline{0}}_{3 \times 18} \\
   \underline{\underline{0}}_{3 \times 3} & 
   \left( 
    \widetilde{\lambda}_2  \widetilde{ \underline{\underline{R}} \underline{r}_1^{\mathrm{sp,r1}}} +  
    \widetilde{\lambda}_4  \widetilde{ \underline{\underline{R}} \underline{r}_1^{\mathrm{sp,r2}}} +  
    \widetilde{\lambda}_6  \widetilde{ \underline{\underline{R}} \underline{r}_1^{\mathrm{sp,r3}}}\right)  &
   \underline{\underline{0}}_{3 \times 18} \\
   \underline{\underline{0}}_{18 \times 3} & 
   \underline{\underline{0}}_{18 \times 3} & 
   \underline{\underline{0}}_{18 \times 18} 
   \end{bmatrix} \in \mathbb{R}^{24 \times 24}

.. note::

   The term :math:`\underline{\underline{K}}^\Phi` is include here for
   completeness, but is not in the current Kynema-FMB implementation, and
   nonlinear-system convergence has been satisfactory.


