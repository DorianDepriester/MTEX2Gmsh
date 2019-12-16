---
title: 'MTEX2Gmsh: a tool for generating 2D meshes from EBSD data'
tags:
  - Polycrystals
  - FEM
  - Grain Boundary
  - EBSD
  - mesh
  - Matlab
authors:
  - name: Dorian Depriester
    orcid: 0000-0002-2881-8942
    affiliation: 1
  - name: Régis Kubler
    orcid: 0000-0001-7781-5855
    affiliation: 1
affiliations:
  - name: MSMP laboratory (EA 7350), Ecole Nationale Supérieure d'Arts et Métiers, 2 cours des Arts et Métiers - 13617 Aix-en-Provence, France
    index: 1
date: 13 december 2019
bibliography: paper.bib

---

# Summary

Micromechanics of polycrystalline aggregates consists in evaluating the thermo-mechanical behaviour of the aggregates at their grain scale. It is of great interest for investigating the influence of grain morphology (size and shape) and crystallographic texture on the macroscopic behaviour. If the investigated material is subjected to macroscopic deformation, the local strain can be obtained either experimentally, thanks to full-field measurement methods such as microgrid technique [@Allais:1994] or Digital Image Correlation (DIC) [@Hild:2002], or thanks to numerical simulation of the microstructure. The latter needs to take into account the mechanical heterogeneities (due to the different constituents) and the anisotropy of each phase, depending on its crystalline orientation.

Orientation Imaging Microscopy (OIM), usually made from Electron Backscatter Diffraction (EBSD) data, is now widely used in the literature for investigating the grain morphology and local crystal orientations in crystalline materials. Raw EBSD data can be considered as matrices of measurements of crystallographic data: each dot contains information about the phase and its orientation at the corresponding position.

In order to perform Finite Element Analysis (FEA) on a polycrystal, one needs to first generate a mesh based on either EBSD or reconstructed grains. In this mesh, the Graind Boundaries (GBs) must be accurately described since they play an important role in the overall behaviour of aggregates. Indeed, it is known that GBs increase the energy of the materials. The interfacial energy between two adjacent grains due to their boundary depends, among other parameters, on their misorientation and on the surface normal of the boundary [@Priester:2012]. In addition, @Zhong:2017 mentioned that the GB curvature is one of the most important properties of a microstructure. For instance, the driving force for grain growth depends on the local curvature of the GBs.

@Latypov:2016 have proposed a program to generate regular pseudo-3D mesh, consisting in brick elements with only one element in thickness. Nevertheless, because of the regular structure, the GBs can not be smoothly described. In addition, the element size must be constant, possibly resulting in a huge number of elements, depending on the size and the spatial resolution of the orientation map. @Dancette:2016 have proposed the following method to generate a conforming mesh with smooth GBs:

* computation of the GBs based on a proper criterion;
* grain reconstruction using a graph theory -based method;
* spline interpolation of the GBs;
* meshing.

The criterion used by the previous authors for defining the GBs, called weight in the context of graph theory, was speciffcally designed for cubic phases. The geometry was meshed using the Gmsh software [@Geuzaine:2009].

The proposed software, named [`MTEX2Gmsh`](https://github.com/DorianDepriester/mtex2Gmsh), provides a customisable and robust method for generating a 2D mesh based on EBSD with smooth and accurate definition of the GBs. It is based on the MTEX toolbox for Matlab [@Bachmann:2011] and the Gmsh software. Figure 1 schematically illustrates the proposed algorithm. The program allows to mesh the volume with a couple of options, such as:

* increasing element size with increasing distance from the grains boundaries;
* element type (tetrahedron, wedge or brick elements);
* nesting the Region of Interest (ROI) into a larger medium.

This sofware comes with an Abaqus plugin for importing the mesh and allocating the phase and Euler Angles of each grain.


![Schematic representation of the algorithm used in `MTEX2Gmsh`: 1) once the grains are reconstructed thank to MTEX [@Bachmann:2011], the algorithm fetches all triple junctions (TJ) in the whole map; 2) each grain boundary is divided into TJ-to-TJ segments; 3) all those segments are smoothed using Bspline approximation; 3) this decription of the grains can be converted into Gmsh-readable files, allowing to mesh the whole region efficiently. The Bspline approximation results in very accurate difinitions of the grain boundaries, with limited serration (usually introduced by the EBSD resolution) and limited number of elements.](GraphicalAbstract.png)


# References
