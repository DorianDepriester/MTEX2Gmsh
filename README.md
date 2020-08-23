# MTEX2Gmsh [![DOI](https://joss.theoj.org/papers/10.21105/joss.02094/status.svg)](https://doi.org/10.21105/joss.02094)

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
The aforementioned examples can be easily reproduced. In addition, the reader can check out the reproductibility of minimal example on [Code Ocean](https://codeocean.com/capsule/8758800/tree/v2).

## Reference
If you use this work, please cite the following paper: Depriester et al., (2020). MTEX2Gmsh: a tool for generating 2D meshes from EBSD data. *Journal of Open Source Software*, 5(52), 2094, https://doi.org/10.21105/joss.02094

## Bug report
Please, use the [Issue](https://github.com/DorianDepriester/MTEX2Gmsh/issues) tab to report any bug or whish for new feature.

## Contribute
You can easily edit the present code so that it fits your needs (as long as this edit complies with the MIT licence). You are also welcome to contribute. In this case, please read [``CONTRIBUTING.md``](CONTRIBUTING.md).
