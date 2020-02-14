# MTEX2Gmsh
This toolbox for Matlab allows to generate meshes from EBSD data. It is intended to perform Finite Element Analysis (FEA) at grain scale on polycrystal imaged by EBSD. It is based on [MTEX](http://mtex-toolbox.github.io/) and [Gmsh](http://gmsh.info/).

## How it works
This toolbox defines the class named `gmshGeo`. Once the grains are computed using MTEX, an instance of `gmshGeo` can be constructed. This object can be used to generate a Gmsh-readable file, in order to mesh it and perform FEA.

## Documentation and examples
Visit the corresponding [site](https://doriandepriester.github.io/MTEX2Gmsh/) to see examples and full documentation. Alternatively, you can check out the [``docs/Examples``](https://github.com/DorianDepriester/MTEX2Gmsh/tree/master/docs/Examples) folder.

## Reference
The corresponding paper is still under review. In the meantime, please cite this code using the related DOI:
[![DOI](https://zenodo.org/badge/137471547.svg)](https://zenodo.org/badge/latestdoi/137471547)

or use the following entry in BibTeX :
```
@misc{mtex2Gmsh,
  author       = {Dorian Depriester},
  title        = {{MTEX2Gmsh}},
  month        = jan,
  year         = 2020,
  publisher    = {Zenodo},
  version      = {2.0.0},
  doi          = {10.5281/zenodo.3628143},
  url          = {https://doi.org/10.5281/zenodo.3628143}
}
```

## Bug report
Please, use the [Issue](https://github.com/DorianDepriester/MTEX2Gmsh/issues) tab to report any bug or whish for new feature.
