# MTEX2Gmsh
This toolbox for Matlab allows to generate meshes from EBSD data. It is intended to perform Finite Element Analysis (FEA) at grain scale on polycrystal imaged by EBSD. It is based on [MTEX](http://mtex-toolbox.github.io/) and [Gmsh](http://gmsh.info/).

## Usage
This toolbox defines the class named `gmshGeo`. Once the grains are computed using MTEX, an instance of `gmshGeo` can be constructed. This object can be used to generate a Gmsh-readable file, in order to mesh it and perform FEA. For further details, please look at the [Wiki](https://github.com/DorianDepriester/mtex2Gmsh/wiki).

## Example
The file named [`Example.m`](https://github.com/DorianDepriester/mtex2Gmsh/blob/master/Example.m) provides a minimum working example of the usage of the `gmshGeo` class. The picture below illustrates a mesh obtained from an example dataset, provided along with MTEX.

![Mesh with size gradient](https://user-images.githubusercontent.com/7643875/57290884-a1abd180-70be-11e9-8383-8675d1221f93.png)


## Reference
The corresponding paper is still under review. In the meantime, please cite this code using the related DOI:
[![DOI](https://zenodo.org/badge/137471547.svg)](https://zenodo.org/badge/latestdoi/137471547)

or use the following entry in BibTeX :
```
@misc{mtex2Gmsh,
    author       = {Dorian Depriester},
    title        = {{MTEX2Gmsh}},
    month        = march,
    year         = 2019,
    doi          = {10.5281/zenodo.2595052},
    version      = {V1.0},
    publisher    = {Zenodo},
    url          = {https://github.com/DorianDepriester/mtex2Gmsh}
}
```
