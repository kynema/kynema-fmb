.. _sec-constraint-types:

Constraint Library
==================

This section catalogs the constraint types available in Kynema-FMB. For each constraint we
give the constraint equation :math:`\underline{\Phi}`, its gradient
:math:`\underline{\underline{B}} = \partial \underline{\Phi}/\partial \underline{q}` used to
assemble the residual and iteration matrix (see :ref:`sec-gen-alpha`), and, where
:math:`\underline{\underline{B}}` is not constant with respect to the system state, the
additional iteration-matrix term :math:`\underline{\underline{K}}^\Phi` that arises from the
constraint force :math:`\underline{\underline{B}}^T \underline{\lambda}`. Unless noted
otherwise, :math:`\alpha` denotes the base/parent node of a constraint and, for two-node
constraints, :math:`\beta` denotes the dependent/child node.

BC3
---

Translational boundary condition for a 3-DOF node :math:`\alpha`:

.. math::

   \underline{\Phi}^\mathrm{BC3}_\alpha \left(\underline{u}_\alpha,\underline{u}^\mathrm{p}_\alpha\right) = \begin{bmatrix}
     \underline{u}_\alpha - \underline{u}^\mathrm{p}_\alpha 
   \end{bmatrix}
   \quad \in \mathbb{R}^3 

where :math:`\underline{u}^\mathrm{p}_\alpha` is the prescribed position of node :math:`\alpha`. 
Taking the variation of the translational constraint, the constraint gradient can be expressed as

.. math::

   \delta \underline{\Phi}^\mathrm{BC3}_\alpha = \underline{\underline{B}}^\mathrm{BC3}_\alpha \delta \underline{u}_\alpha 


where

.. math::
   \underline{\underline{B}}^\mathrm{BC3}_\alpha = \underline{\underline{I}}
   \quad \in \mathbb{R}^{3\times3}

Since :math:`\underline{\Phi}^\mathrm{BC3}_\alpha` is linear in :math:`\underline{u}_\alpha`,
:math:`\underline{\underline{B}}^\mathrm{BC3}_\alpha` is constant with respect to the system
state, and the associated additional iteration-matrix term in Eq. :eq:`iteration` vanishes:

.. math::

   \underline{\underline{K}}_\alpha^\mathrm{BC3,\Phi} = \underline{\underline{0}}
   \quad \in \mathbb{R}^{3\times3}


BC6
---

Translation and rotation boundary condition with prescribed position :math:`\underline{u}^\mathrm{p}_\alpha \in \mathbb{R}^3` and orientation :math:`\underline{\underline{R}}^\mathrm{p}_\alpha \in SO(3)`
for a 6-DOF node :math:`\alpha` with configuration in :math:`\mathbb{R}^3 \times SO(3)`:

.. math::

   \underline{\Phi}^\mathrm{BC6}_\alpha \left(\underline{u}_\alpha, \underline{\underline{R}}_\alpha,\underline{u}^\mathrm{p}_\alpha, \underline{\underline{R}}_\alpha^\mathrm{p}\right)= \begin{bmatrix}
     \underline{u}_\alpha - \underline{u}^\mathrm{p}_\alpha
     \\ \mathrm{axial}\left( \underline{\underline{R}}_\alpha \left(\underline{\underline{R}}^\mathrm{p}_\alpha\right)^T \right)
   \end{bmatrix}
   \quad \in \mathbb{R}^6

where we have used the axial vector to constrain rotation following
[@Sonneville-Bruls:2012].

Taking the variation of the translational and rotational constraint, the constraint gradient can be expressed as

.. math::

   \delta \underline{\Phi}^\mathrm{BC6}_\alpha = \underline{\underline{B}}^\mathrm{BC6}_\alpha
   \begin{bmatrix}
     \delta \underline{u}_\alpha \\
     \underline{\delta \theta}_\alpha
   \end{bmatrix}

where

.. math::

   \underline{\underline{B}}^\mathrm{BC6}_\alpha = \begin{bmatrix}
     \underline{\underline{I}} & \underline{\underline{0}} \\
     \underline{\underline{0}} & \mathrm{AX}\left( \underline{\underline{R}}_\alpha \left(\underline{\underline{R}}^\mathrm{p}_\alpha\right)^T \right)
   \end{bmatrix}
   \quad \in \mathbb{R}^{6\times6}

where :math:`\underline{\underline{R}}^\mathrm{p}_\alpha` is the prescribed rotation of rotaiton of node :math:`\alpha`, and we leveraged the following identity

.. math:: 

   \mathrm{axial}\left({ \underline{\underline{A}} \widetilde{a}}\right) =  
   \frac{1}{2} \left(\mathrm{tr} \left( \underline{\underline{A^T}} \right) \underline{\underline{I}} -  \underline{\underline{A}}^T \right) \underline{a}, 
   \quad \forall \underline{a} \in \mathbb{R}^3 \ \text{and} \ \underline{\underline{A}} \in \mathbb{R}^{3\times3}

and we introduced the operator

.. math:: \mathrm{AX}\left( \underline{\underline{A}} \right) =  \frac{1}{2} \left(\mathrm{tr} \left( \underline{\underline{A}} \right) \underline{\underline{I}} -  \underline{\underline{A}} \right)

Since :math:`\underline{\underline{B}}^\mathrm{BC6}_\alpha` is not constant with respect to
:math:`\underline{\theta}_\alpha`, there is an additional term that arises in
:math:`\underline{\underline{K}}^\Phi` in Eq. :eq:`iteration`. The associated constraint force is

.. math::

   \underline{F}^\mathrm{BC6}_\alpha = \left(\underline{\underline{B}}^\mathrm{BC6}_\alpha\right)^T \underline{\lambda}_\alpha,
   \quad \underline{\lambda}_\alpha = \begin{bmatrix} \underline{\lambda}_{\alpha,1} \\ \underline{\lambda}_{\alpha,2} \end{bmatrix} \in \mathbb{R}^6, 
   \quad \underline{\lambda}_{\alpha,1}, \underline{\lambda}_{\alpha,2} \in \mathbb{R}^3

Taking the variation, we obtain

.. math:: \delta \underline{F}^\mathrm{BC6}_\alpha =
   \underline{\underline{K}}_\alpha^\mathrm{BC6,\Phi}
   \begin{bmatrix}
   \delta \underline{u}_\alpha \\
   \underline{\delta \theta}_\alpha
   \end{bmatrix}
   + \left(\underline{\underline{B}}^\mathrm{BC6}_\alpha\right)^T \delta \underline{\lambda}_\alpha

where

.. math::

   \underline{\underline{K}}_\alpha^\mathrm{BC6,\Phi}
   =
   \begin{bmatrix}
    \underline{\underline{0}} & \underline{\underline{0}} \\
    \underline{\underline{0}} & \mathrm{AX} \left(  \underline{\underline{R}}_\alpha \left(\underline{\underline{R}}^\mathrm{p}_\alpha\right)^T 
    \widetilde{\lambda}_{\alpha,2} \right)
   \end{bmatrix}
   \quad \in \mathbb{R}^{6\times6}


Rigid63
--------

Rigid connection between a 6-DOF node :math:`\alpha` and a 3-DOF node :math:`\beta`,
i.e., :math:`\beta` is rigidly embedded in the frame of :math:`\alpha`.

.. math::

   \underline{\Phi}^\mathrm{Rigid63}_{\alpha-\beta} \left( \underline{u}_\alpha, \underline{\underline{R}}_\alpha, 
   \underline{u}_\beta, \underline{r}^\mathrm{r}_{\alpha-\beta}\right) = 
   \underline{u}_\beta - \underline{u}_\alpha
   - \left( \underline{\underline{R}}_\alpha - \underline{\underline{I}} \right) \underline{r}^\mathrm{r}_{\alpha-\beta}
   \quad \in \mathbb{R}^3

where :math:`\underline{r}^\mathrm{r}_{\alpha-\beta} = \underline{x}^\mathrm{r}_\beta - \underline{x}^\mathrm{r}_\alpha`
is the reference relative position vector from node :math:`\alpha` to node :math:`\beta`.

Taking the variation of the rigid connection constraint, the constraint gradient can be expressed as

.. math::

   \delta \underline{\Phi}^\mathrm{Rigid63}_{\alpha-\beta} = \underline{\underline{B}}^\mathrm{Rigid63}_{\alpha-\beta}
   \begin{bmatrix}
     \delta \underline{u}_\alpha \\
     \underline{\delta \theta}_\alpha \\
     \delta \underline{u}_\beta
   \end{bmatrix}

where

.. math::

   \underline{\underline{B}}^\mathrm{Rigid63}_{\alpha-\beta} = \begin{bmatrix}
     -\underline{\underline{I}} & \widetilde{\underline{\underline{R}}_\alpha \underline{r}^\mathrm{r}_{\alpha-\beta}} & \underline{\underline{I}}
   \end{bmatrix}
   \quad \in \mathbb{R}^{3\times9}

The associated constraint force is

.. math::

   \underline{F}^\mathrm{Rigid63}_{\alpha-\beta} = \left(\underline{\underline{B}}^\mathrm{Rigid63}_{\alpha-\beta}\right)^T \underline{\lambda}_{\alpha-\beta},
   \quad \underline{\lambda}_{\alpha-\beta} \in \mathbb{R}^3

Taking the variation, we obtain

.. math:: \delta \underline{F}^\mathrm{Rigid63}_{\alpha-\beta} =
   \underline{\underline{K}}_{\alpha-\beta}^\mathrm{Rigid63,\Phi}
   \begin{bmatrix}
   \delta \underline{u}_\alpha \\
   \underline{\delta \theta}_\alpha \\
   \delta \underline{u}_\beta
   \end{bmatrix}
   + \left(\underline{\underline{B}}^\mathrm{Rigid63}_{\alpha-\beta}\right)^T \delta \underline{\lambda}_{\alpha-\beta}

where, since :math:`\underline{\underline{B}}^\mathrm{Rigid63}_{\alpha-\beta}` depends on the
state only through :math:`\underline{\underline{R}}_\alpha`, only the
:math:`\underline{\delta \theta}_\alpha` column is non-zero:

.. math::

   \underline{\underline{K}}_{\alpha-\beta}^\mathrm{Rigid63,\Phi}
   =
   \begin{bmatrix}
    \underline{\underline{0}} & \underline{\underline{0}} & \underline{\underline{0}} \\
    \underline{\underline{0}} & -\widetilde{\lambda}_{\alpha-\beta} \widetilde{\underline{\underline{R}}_\alpha \underline{r}^\mathrm{r}_{\alpha-\beta}} & \underline{\underline{0}} \\
    \underline{\underline{0}} & \underline{\underline{0}} & \underline{\underline{0}}
   \end{bmatrix}


Rigid66
--------

Rigid connection between a 6-DOF node :math:`\alpha` and a 6-DOF node :math:`\beta`
with a prescribed relative rotation applied to node :math:`\alpha` in its referenced configuration, :math:`\underline{\underline{R}}^\mathrm{p}_{\alpha-\beta}\in SO(3)`
between the two frames:

.. math::

   \underline{\Phi}^\mathrm{Rigid66}_{\alpha-\beta}\left( \underline{u}_\alpha, \underline{\underline{R}}_\alpha, \underline{u}_\beta, \underline{\underline{R}}_\beta, 
   \underline{\underline{R}}^\mathrm{p}_{\alpha-\beta}, \underline{r}^\mathrm{r}_{\alpha-\beta} \right) = 
   \begin{bmatrix}
     \underline{u}_\beta - \underline{u}_\alpha - \left( \underline{\underline{R}}_\alpha - \underline{\underline{I}} \right) \underline{r}^\mathrm{r}_{\alpha-\beta}
     \\ \mathrm{axial}\left( \underline{\underline{R}}_\alpha \left(\underline{\underline{R}}^\mathrm{p}_{\alpha-\beta}\right)^T \underline{\underline{R}}_\beta^T \right)
   \end{bmatrix}
   \quad \in \mathbb{R}^6

where :math:`\underline{r}^\mathrm{r}_{\alpha-\beta} = \underline{x}^\mathrm{r}_\beta - \underline{x}^\mathrm{r}_\alpha`
is the reference relative position vector from node :math:`\alpha` to node :math:`\beta`. 
The rotational constraint is derived from the relation:

.. math::

   \underline{\underline{R}}_\alpha^\mathrm{c} =  \underline{\underline{R}}_\beta \underline{\underline{R}}^\mathrm{p}_{\alpha-\beta} \underline{\underline{R}}_\alpha^\mathrm{r}
  \quad \Rightarrow \quad \underline{\underline{R}}_\alpha = \underline{\underline{R}}_\beta \underline{\underline{R}}^\mathrm{p}_{\alpha-\beta}

Taking the variation of the rigid connection constraint, the constraint gradient can be expressed as

.. math::

   \delta \underline{\Phi}^\mathrm{Rigid66}_{\alpha-\beta} = \underline{\underline{B}}^\mathrm{Rigid66}_{\alpha-\beta}
   \begin{bmatrix}
     \delta \underline{u}_\alpha \\
     \underline{\delta \theta}_\alpha \\
     \delta \underline{u}_\beta \\
     \underline{\delta \theta}_\beta
   \end{bmatrix}

where

.. math::

   \underline{\underline{B}}^\mathrm{Rigid66}_{\alpha-\beta} = \begin{bmatrix}
     -\underline{\underline{I}} & \widetilde{\underline{\underline{R}}_\alpha \underline{r}^\mathrm{r}_{\alpha-\beta}} & \underline{\underline{I}} & \underline{\underline{0}} \\
     \underline{\underline{0}} & -\mathrm{AX}\left( \underline{\underline{R}}_\beta \underline{\underline{R}}^\mathrm{p}_{\alpha-\beta} \underline{\underline{R}}_\alpha^T \right) &
     \underline{\underline{0}} & \mathrm{AX}\left( \underline{\underline{R}}_\alpha \left(\underline{\underline{R}}^\mathrm{p}_{\alpha-\beta}\right)^T \underline{\underline{R}}_\beta^T \right)
   \end{bmatrix}
   \quad \in \mathbb{R}^{6\times12}

Proceeding as for :math:`\underline{\underline{B}}^\mathrm{BC6}_\alpha` above, the constraint
force is

.. math::

   \underline{F}^\mathrm{Rigid66}_{\alpha-\beta} = \left(\underline{\underline{B}}^\mathrm{Rigid66}_{\alpha-\beta}\right)^T \underline{\lambda}_{\alpha-\beta},
   \quad \underline{\lambda}_{\alpha-\beta} = \begin{bmatrix} \underline{\lambda}_{\alpha-\beta,1} \\ \underline{\lambda}_{\alpha-\beta,2} \end{bmatrix} \in \mathbb{R}^6

Taking the variation, we obtain

.. math:: \delta \underline{F}^\mathrm{Rigid66}_{\alpha-\beta} =
   \underline{\underline{K}}_{\alpha-\beta}^\mathrm{Rigid66,\Phi}
   \begin{bmatrix}
   \delta \underline{u}_\alpha \\
   \underline{\delta \theta}_\alpha \\
   \delta \underline{u}_\beta \\
   \underline{\delta \theta}_\beta
   \end{bmatrix}
   + \left(\underline{\underline{B}}^\mathrm{Rigid66}_{\alpha-\beta}\right)^T \delta \underline{\lambda}_{\alpha-\beta}

where

.. math::

   \begin{aligned}
   \underline{\underline{K}}_{\alpha-\beta}^\mathrm{Rigid66,\Phi}
   =
   \begin{bmatrix}
    \underline{\underline{0}} & \underline{\underline{0}} & \underline{\underline{0}} & \underline{\underline{0}}\\
    \underline{\underline{0}} & -\mathrm{AX} \left(\underline{\underline{R}}_\beta \underline{\underline{R}}^\mathrm{p}_{\alpha-\beta} \underline{\underline{R}}_\alpha^T \widetilde{\lambda}_{\alpha-\beta,2}  \right)
   -\widetilde{\lambda}_{\alpha-\beta,1} \widetilde{ \underline{\underline{R}}_\beta \underline{r}^\mathrm{r}_{\alpha-\beta} }
   & \underline{\underline{0}} & -\mathrm{AX2}\left( \underline{\lambda}_{\alpha-\beta,2}, \underline{\underline{R}}_\alpha \left(\underline{\underline{R}}^\mathrm{p}_{\alpha-\beta}\right)^T \underline{\underline{R}}_\beta^T \right) \\
    \underline{\underline{0}} & \underline{\underline{0}} & \underline{\underline{0}} & \underline{\underline{0}}\\
    \underline{\underline{0}} & \mathrm{AX2} \left( \underline{\lambda}_{\alpha-\beta,2}, \underline{\underline{R}}_\beta \underline{\underline{R}}^\mathrm{p}_{\alpha-\beta} \underline{\underline{R}}_\alpha^T\right)
   & \underline{\underline{0}} & \mathrm{AX} \left( \underline{\underline{R}}_\alpha \left(\underline{\underline{R}}^\mathrm{p}_{\alpha-\beta}\right)^T \underline{\underline{R}}_\beta^T  \widetilde{\lambda}_{\alpha-\beta,2} \right)
   \end{bmatrix}
   \quad \in \mathbb{R}^{12\times12}
   \end{aligned}

and where we introduced the operator

.. math::

   \begin{aligned}
   \mathrm{AX2}\left(\underline{a},\underline{\underline{A}}\right) = \frac{1}{2}
   \begin{bmatrix}
   a_1 (A_{23}- A_{32}) & -a_1 A_{13}-a_2 A_{32}-a_3 A_{33} & a_1 A_{12}+a_2 A_{22} + a_3 A_{23} \\
   a_2 A_{23} + a_1 A_{31} + a_3 A_{33} & a_2 (A_{31}-A_{13}) & -a_1 A_{11}-a_2 A_{21} - a_3 A_{13} \\
   -a_1 A_{21}-a_2 A_{22}-a_3 A_{32} & a_1 A_{11}+a_2 A_{12}+a_3 A_{31} & a_3 (A_{12}-A_{21})
   \end{bmatrix}
   \end{aligned}

such that

.. math::

  \mathrm{axial}\left({ \widetilde{a} \widetilde{b} \underline{\underline{A}}}\right) = \mathrm{AX2}\left(\underline{a},\underline{\underline{A}}\right) \underline{b}, 
   \quad \forall \underline{a}, \underline{b} \in \mathbb{R}^3, \ \underline{\underline{A}} \in \mathbb{R}^{3\times3}


Revolute
--------

Revolute joint between two collocated 6-DOF nodes :math:`\alpha` and :math:`\beta`, free to
rotate relative to one another about the :math:`\widehat{x}_\alpha` axis. Options exist for
prescribed torque, though that functionality is external to the constraint definition.

.. math::

   \underline{\Phi}^\mathrm{Revolute}_{\alpha-\beta} \left( \underline{u}_\alpha, \underline{\underline{R}}_\alpha, 
   \underline{\underline{R}}_\alpha^\mathrm{r}, \underline{u}_\beta, \underline{\underline{R}}_\beta, 
   \underline{\underline{R}}_\beta^\mathrm{r} \right) = 
   \begin{bmatrix}
     \underline{u}_\beta - \underline{u}_\alpha
     \\ \widehat{z}_\beta^T \widehat{x}_\alpha \\
     \widehat{y}_\beta^T \widehat{x}_\alpha
   \end{bmatrix} =
   \begin{bmatrix}
     \underline{u}_\beta - \underline{u}_\alpha
     \\ \widehat{z}_\beta^{\mathrm{r}T}  \underline{\underline{R}}_\beta^T  \underline{\underline{R}}_\alpha \widehat{x}_\alpha^\mathrm{r} \\
     \widehat{y}_\beta^{\mathrm{r}T}  \underline{\underline{R}}_\beta^T  \underline{\underline{R}}_\alpha \widehat{x}_\alpha^\mathrm{r}
   \end{bmatrix}
   \quad \in \mathbb{R}^5

Taking the variation of the revolute-joint constraint, the constraint gradient can be expressed as

.. math::

   \delta \underline{\Phi}^\mathrm{Revolute}_{\alpha-\beta} = \underline{\underline{B}}^\mathrm{Revolute}_{\alpha-\beta}
   \begin{bmatrix}
     \delta \underline{u}_\alpha \\
     \underline{\delta \theta}_\alpha \\
     \delta \underline{u}_\beta \\
     \underline{\delta \theta}_\beta
   \end{bmatrix}

where

.. math::

   \underline{\underline{B}}^\mathrm{Revolute}_{\alpha-\beta} = \begin{bmatrix}
     -\underline{\underline{I}} & \underline{\underline{0}} & \underline{\underline{I}} & \underline{\underline{0}} \\
     \underline{0}^T & -\widehat{z}_\beta^T \widetilde{x}_\alpha & \underline{0}^T & \widehat{z}_\beta^T \widetilde{x}_\alpha \\
     \underline{0}^T & -\widehat{y}_\beta^T \widetilde{x}_\alpha & \underline{0}^T & \widehat{y}_\beta^T \widetilde{x}_\alpha
   \end{bmatrix}
   \quad \in \mathbb{R}^{5\times12}

The associated constraint force is

.. math::

   \underline{F}^\mathrm{Revolute}_{\alpha-\beta} = \left(\underline{\underline{B}}^\mathrm{Revolute}_{\alpha-\beta}\right)^T \underline{\lambda}_{\alpha-\beta},
   \quad \underline{\lambda}_{\alpha-\beta} = \begin{bmatrix} \underline{\lambda}_{\alpha-\beta,1} \\ \lambda_{\alpha-\beta,2} \\ \lambda_{\alpha-\beta,3} \end{bmatrix} \in \mathbb{R}^5
   \quad \underline{\lambda}_{\alpha-\beta,1} \in \mathbb{R}^3

Taking the variation, we obtain

.. math:: \delta \underline{F}^\mathrm{Revolute}_{\alpha-\beta} =
   \underline{\underline{K}}_{\alpha-\beta}^\mathrm{Revolute,\Phi}
   \begin{bmatrix}
   \delta \underline{u}_\alpha \\
   \underline{\delta \theta}_\alpha \\
   \delta \underline{u}_\beta \\
   \underline{\delta \theta}_\beta
   \end{bmatrix}
   + \left(\underline{\underline{B}}^\mathrm{Revolute}_{\alpha-\beta}\right)^T \delta \underline{\lambda}_{\alpha-\beta}

where

.. math::

   \begin{aligned}
   \underline{\underline{K}}_{\alpha-\beta}^\mathrm{Revolute,\Phi}
   =
   \begin{bmatrix}
    \underline{\underline{0}} & \underline{\underline{0}} & \underline{\underline{0}} & \underline{\underline{0}}\\
    \underline{\underline{0}} & \lambda_{\alpha-\beta,2} \widetilde{z}_\beta \widetilde{x}_\alpha
           + \lambda_{\alpha-\beta,3} \widetilde{y}_\beta \widetilde{x}_\alpha
          & \underline{\underline{0}} &
            -\lambda_{\alpha-\beta,2} \widetilde{x}_\alpha \widetilde{z}_\beta
           - \lambda_{\alpha-\beta,3} \widetilde{x}_\alpha \widetilde{y}_\beta
          \\
    \underline{\underline{0}} & \underline{\underline{0}} & \underline{\underline{0}} & \underline{\underline{0}}\\
    \underline{\underline{0}} & -\lambda_{\alpha-\beta,2} \widetilde{z}_\beta \widetilde{x}_\alpha
           - \lambda_{\alpha-\beta,3} \widetilde{y}_\beta \widetilde{x}_\alpha
         &    \underline{\underline{0}} &
            \lambda_{\alpha-\beta,2} \widetilde{x}_\alpha \widetilde{z}_\beta
           + \lambda_{\alpha-\beta,3} \widetilde{x}_\alpha \widetilde{y}_\beta
   \end{bmatrix}
   \end{aligned}

.. note::

   The :math:`\underline{\underline{K}}^\Phi` terms are included here for completeness, but are
   not in the current Kynema-FMB implementation, and nonlinear-system convergence has been
   satisfactory.

.. container:: references csl-bib-body hanging-indent
   :name: refs

   .. container:: csl-entry
      :name: ref-Sonneville-Bruls:2012

      Sonneville, V., and O. Brüls. 2012. “Formulation of Kinematic
      Joints and Rigidity Constraints in Multibody Dynamics Using a Lie
      Group Approach.” In *Proceedings of the 2nd Joint International
      Conference on Multibody System Dynamics*. Stuttgart, Germany.