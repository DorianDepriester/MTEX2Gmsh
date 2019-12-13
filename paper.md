---
title: 'Polycrystalline FEM representation: a tool
for generating 2D meshes from EBSD data'
tags:
  - Polycrystals
  - FEM
  - Grain Boundary
  - EBSD
  - mesh
author:
  - name: Depriester, Dorian
  orcid: 0000-0002-2881-8942
  affiliation: 1
affiliation:
  -name: MSMP laboratory (EA 7350), Ecole Nationale Supérieure d'Arts et Métiers, 2 cours des Arts et Métiers - 13617
Aix-en-Provence, France
date: 13 december 2019
bibliography: mtex2Gmsh.bib

# Summary
Micromechanics of polycrystalline aggregates consists in evaluating the thermo-mechanical behaviour of the aggregates at their grain scale. 
It is of great interest for investigating the influence of grain morphology (size and shape) and crystallographic texture on the macroscopic behaviour. 
If the investigated material is subjected to macroscopic deformation, the local strain can be obtained either experimentally, 
thanks to full-field measurement methods such as microgrid technique [1] or Digital Image Correlation (DIC) [2],
or thanks to numerical simulation of the microstructure. 
The latter needs to take into account the mechanical heterogeneities (due to the different constituents) and the anisotropy of each phase,
depending on its crystalline orientation.

In order to perform Finite Element Analysis (FEA) on a polycrystal, one needs to first generate a mesh based on either EBSD or reconstructed grains. 
In this mesh, the Graind Boundaries (GBs) must be accurately described since they play an important role in the overall behaviour of aggregates. 
Indeed, it is known that GBs increase the energy of the materials. 
The interfacial energy between two adjacent grains due to their boundary depends, among other parameters, on their misorientation and on the surface
normal of the boundary [4]. 
In addition, Zhong et al. [5] mentioned that the GB curvature is one of the most important properties of a microstructure. For instance, the driving
force for grain growth depends on the local curvature of the GBs.

This paper aims to propose a customisable and robust method for generating a 2D mesh based on EBSD with smooth and accurate definition of the GBs.
