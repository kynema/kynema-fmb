---
title:  'Kynema-FMB: A Performance-Portable Flexible-Multibody-Dynamics Solver'
tags:
  - flexible multibody dynamics
  - differential-algebraic equations
  - geometrically exact beam theory
  - fluid-structure interaction
authors:
  - name: Michael A. Sprague
    orcid: 0000-0003-0519-0379
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 1
  - name: Derek Slaughter
#   equal-contrib: true
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: David Dement
    affiliation: 2
  - name: Faisal H. Bhuiyan
    orcid: 0009-0001-1029-0908
    affiliation: 1
  - name: Michael Kuhn
    orcid: 0000-0002-5506-7777
    affiliation: 1
  - name: Paul S. Crozier
    corresponding: false # (This is how to denote the corresponding author)
    affiliation: 2
affiliations:
 - name: National Laboratory of the Rockies
   index: 1
   ror: '036266993'
 - name: Sandia National Laboratories
   index: 2
   ror: '01apwpt12'
date: 12 December 2025
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal. ->
---

# Summary

Kynema-FMB [@kynema:2022] is an open-source performance-portable flexible
multibody (FMB) dynamics solver designed for time-domain simulations.
Kynema-FMB was designed with a narrow focus, namely to provide a lightweight,
fast, accurate FMB dynamics solver for coupling to computational-fluid-dynamics
(CFD) codes, especially the broader suite of Kyenma CFD codes (formerly called
ExaWind) [@Sprague-etal:2020;@Sharma2023;@Kuhn-etal:2025], for
fluid-structure-interaction (FSI) simulations.  Kynema-FMB is equipped to model
systems that can be represented as a collection of beams and rigid bodies that
are connected through constraints.  Degrees of freedom are defined in the
inertial/global frame of reference and include displacements and rotations
(formally as rotation matrices, but stored as quaternions).  The underlying
formulation is built on a Lie-group time integrator designed for index-3
differential-algebraic equations that is second-order accurate in time
[@Bruls-etal:2012].  Beam models are based on geometrically exact beam theory
[@Reissner:1973;@Bauchau:2011] and are discretized as high-order spectral
finite elements similar to those in BeamDyn [@Wang-etal:2017].  The governing
equations for many FMB systems constitute a highly nonlinear system of
constrained partial-differential equations.  Kynema-FMB uses analytical
Jacobians in the nonlinear-system solves in each time step.  Linear systems use
sparse storage and several third-party sparse-linear-system solvers are
enabled. Ill conditioning of linear systems is mitigated with preconditioning
described in @Bottasso-etal:2008.  Kynema-FMB is integrated with a simple
open-source controller [@abbas2022reference]. There is an application
programming interface (API) for coupling to geometry-resolved CFD (like that in
@Sharma2023) and actuator-force CFD (like that in @Kuhn-etal:2025).  In the
latter, for actuator-line models, Kynema-FMB includes an internal blade-element
aerodynamics solver that depends on user-provided lookup tables for
coefficients of lift and drag, i.e., aerodynamic polars.  Kynema-FMB is written
in C\texttt{++} and leverages Kokkos and Kokkos-Kernels [@KokkosEcosystem2021]
as its performance portability layer enabling simulations on both CPU and GPU
systems. The repository is equipped with extensive automated testing at the
unit and regression/system levels.

# Statement of need

The motivating use case for Kynema-FMB was to provide a modern
structural-dynamics code for FSI coupling to the Kynema suite of open-source
high-fidelity CFD solvers.  Kynema-FMB fills the need for an open-source FMB
structural dynamics code that is high fidelity, lightweight, robust, fast, and
capable of running on different computer architectures, including GPUs from
multiple manufacturers.  Kynema-FMB is well suited for simulating structures
that can be modeled as rigid bodies and slender deformable structures (i.e.,
beams), including offshore oil and gas infrastructure (floating platforms and
mooring lines), power transmission lines, and rotor craft.

Kynema is written in C\texttt{++} to facilitate integration with the
C\texttt{++}-based Kynema CFD codes.  While Kynema-FMB's primary target is
coupling to the Kynema CFD suite, it is designed to be readily coupled with
other fluid-dynamics codes as well.  A driving need was also to have a code
base that was developed using modern-software-development best practices to
ensure the code is future proof and has minimal maintenance requirements.  Best
practices include test-driven-development, a multi-platform build system,
nightly testing, thorough documentation, and extensive automated hierarchical
testing.

We note that the Kynema CFD codes were first established with coupling to
OpenFAST [@openfast;@Jonkman:2013], which
includes models for structural dynamics, aerodynamics, hydrodynamics, control
systems, drivetrain, etc.  While OpenFAST was critical to the establishment and
success of Kynema CFD codes, it faced challenges as a component in the Kynema
suite due to its limited unit and regression testing and that it is Fortran
code, which complicates creation and maintenance of APIs for coupling to
C\texttt{++} codes.

# State of the field

There are several established open-source, C\texttt{++}-based FMB dynamics
codes, including MBDyn [@Masarati-etal:2014], Exudyn [@Gerstmayr:2023] (which
is a follow on to HOTINT [@Gerstmayr:2009]), and Chrono (part of ProjectChrono)
[@Chrono2016;@projectChronoGithub].  Kynema shares features with several of
these codes, e.g., similar time-integration (ExuDyn and Chrono), use of the
open-source Eign library for linear algebra (Chrono and HOTINT), and
three-dimensional geometrically nonlinear beams (Chrono and MBDyn).  Kynema-FMB
stands apart in that it is:

  - equipped with an extensive hierarchical automated test suite (currently 448
    tests) that is run nightly and reported on CDash,

  - built with Kokkos/Kokkos-Kernels enabling shared-memory CPU-based computing
  as well as GPU-based computing (on GPU hardware from multiple vendors),

  - equipped with high-order spectral finite elements for beam discretization,
    and

  - equipped for Spack-based [@Gamblin-etal:2015] installation in addition to
    CMake.

# Software design

## Programming language and models

Kynema-FMB is written in C\texttt{++} with tight integration of the Kokkos
performance portability library. This approach allows a single code base to
achieve near-optimal performance when run on CPU or any GPU platform, rather
than requiring separate code paths. Optimized math routines, such as those
covered by BLAS and LAPACK packages or sparse linear solvers, are obtained from
specialized third-party libraries to ensure state-of-the-art performance.

## Accessibility

Kynema-FMB provides a user-friendly, high level API for developers to use to
define and run their structural dynamics problems and to couple Kynema-FMB with
other codes. This approach decouples users of Kynema-FMB from the low level
details of its implementation and improves the speed at which developers can
define and execute their problem. Advanced users are able to use the lower
level Kynema-FMB APIs directly or to define their own interfaces, should the
high level APIs not address their needs directly.

## Robustness

Kynema-FMB prioritizes software robustness through a comprehensive unit and
regression test suite, which are run through a Continuous Integration process.
Beyond that, Kynema-FMB continuously runs a variety of static and dynamic
analysis tools to identify potential bugs. Linters and manual code review on
every change help to ensure consistent, well designed, and sustainable
software.

## Performance

Kynema-FMB is performance-focused software. Core data structures are designed
to provide optimal cache usage and all algorithms are written to best take
advantage of on-chip resources. Kynema-FMB performs optimally on both CPU and
GPU, using hierarchical parallelism and other techniques to ensure performance
portability for problems of all sizes.

# Research impact statement

As described in @Sprague-etal:2025, Kynema-FMB has been 

  - successfully coupled to a blade-element aerodynamics solver and controller;
  simulations were compared against the OpenFAST simulations,

  - successfully coupled to a two-phase Kynema FSI CFD simulation of a geometry
  resolved floating platform where Kynema-FMB solved the rigid body motion.

As described in the Kynema-FMB documentation [@kynema-doc], it has also
been coupled to the Kynema-SGF structured-grid CFD code for actuator-line FSI
simulations of the NREL 5-MW turbine [@jonkman_definition_2009] in a turbulent
atmospheric boundary layer, with results compared against those simulated with
OpenFAST.  It was demonstrated that Kynema-FMB could provide stable, accurate
solutions wherein the structure and CFD models used the same time step ($\Delta
t_\mathrm{CFD}$), whereas OpenFAST required a much smaller time step ($\Delta
t_\mathrm{CFD}/20$) for stability. 

# AI usage disclosure

Generative AI was not used in software creation, online documentaiton, or paper authoring. 

# Acknowledgements

This work was authored in part by the National Laboratory of the Rockies for
the U.S. Department of Energy (DOE), operated under Contract No.
DE-AC36-08GO28308.  This work was also authored in part by Sandia National
Laboratories, a multimission laboratory managed and operated by National
Technology \& Engineering Solutions of Sandia, LLC, a wholly owned subsidiary
of Honeywell International Inc., for the DOE National Nuclear Security
Administration (NNSA) under contract DE-NA0003525.  Funding provided by the DOE
Office of Critical Minerals and Energy Innovation Integrated Energy Systems
Office.  Additional funding provided by the DOE Office of Science, Offices of
Advanced Scientific Computing Research (ASCR) and Biological and Environmental
Research through the FLOWMAS Energy Earthshot Research Center, and through the
ASCR Exascale Computing Project.  The views expressed in the article do not
necessarily represent the views of the DOE or the U.S. Government. The U.S.
Government retains and the publisher, by accepting the article for publication,
acknowledges that the U.S. Government retains a nonexclusive, paid-up,
irrevocable, worldwide license to publish or reproduce the published form of
this work, or allow others to do so, for U.S.  Government purposes.

The authors thank Prof.\ Olivier Brüls for his assistance in understanding the
Lie-group integrator, and Alan Williams and Rafael Mudafort for their early
assistance with the development plan and infrastructure.

# References


