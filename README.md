# MTEX2Gmsh
This toolbox for Matlab allows to generate meshes from EBSD data. It is intended to perform Finite Element Analysis (FEA) at grain scale on polycrystal imaged by EBSD. It is based on [MTEX](http://mtex-toolbox.github.io/) and [Gmsh](http://gmsh.info/).

## How it works
This toolbox defines the class named `gmshGeo`. Once the grains are computed using MTEX, an instance of `gmshGeo` can be constructed. This object can be used to generate a Gmsh-readable file, in order to mesh it and perform FEA.

## Requirements
This toolbox has been designed for MATLAB R2013b, but it may work or newer versions. In addition, the following are required:
- The [MTEX toolbox](https://mtex-toolbox.github.io/) (v 5.1.1 or newer) should be installed in your MATLAB session;
- The [Gmsh software](http://gmsh.info/) (v 4.5.2 or newer) should be installed on your computer (at least its binary should accessible).

It works on both Windows and Unix-like plateform (Linux and Mac OS).

### :warning: Linux user
When running the ``mesh`` command, you may stumble on the error below:

    /MATLAB/sys/os/glnxa64/libstdc++.so.6: version `GLIBCXX_3.4.21' not found (required by gmsh)
    
If so, instead of running 

    matlab
    
run

    LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libstdc++.so.6" matlab

## Documentation and examples
Visit the corresponding [site](https://doriandepriester.github.io/MTEX2Gmsh/) to see examples and [full documentation](https://doriandepriester.github.io/MTEX2Gmsh/html/index.html). Alternatively, you can check out the [``docs/Examples``](https://github.com/DorianDepriester/MTEX2Gmsh/tree/master/docs/Examples) folder.

## Unit test
The aforementioned examples can be easily reproduced. In addition, the reader can check out the reproductibility of minimal example on [Code Ocean](https://codeocean.com/capsule/8758800/tree/v1).

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

## Contribute
You can easily edit the present code so that it fits your needs (as long as this edit complies with the MIT licence). You are also welcome to contribute. In this case, please read [``CONTRIBUTING.md``](CONTRIBUTING.md).
