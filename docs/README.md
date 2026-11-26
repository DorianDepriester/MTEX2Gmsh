# MTEX2Gmsh
Generates conforming meshes from grains reconstructed by MTEX.

## Purpose of this toolbox
In order to evaluate the thermo-mechanical behaviour of crystalline materials (such as metals or ceramics) at microscopic scale, one usually perform numerical simulation at grain scale using the Finite Element Method. In order to proceed, one must first create a mesh which is representative of the real material.

The microstructure of crystalline materials is usually made from Electron Backscattered Diffraction (EBSD) technique. Thus, this toolbox is designed to generate meshes from EBSD in a robust and accurate way.

## Examples
[![Example: aachen.m](./Examples/aachen.png)](./Examples/aachen.png)
[``aachen.m``](Examples/aachen.m)


[![Example: titanium_medium.m](./Examples/titanium_medium.png)](./Examples/titanium_medium.png)
[``titanium_medium.m``](Examples/titanium_medium.m)

[![Example: twins.m](./Examples/twins.png)](./Examples/twins.png)
[``twins.m``](Examples/twins.m)

## Documentation
### Online
You can navigate the documentation [here](html/index.html). Its sidebar menu is generated from [``html/helptoc.xml``](html/helptoc.xml) (the same file used by the MATLAB help browser): after editing it, run ``helptoc2yml`` to regenerate ``_data/toc.yml``.

Note: each page in ``html/`` must start with the Jekyll front matter (``layout: doc`` + ``title``) to be displayed with the sidebar. Since MATLAB ``publish`` overwrites the HTML files and strips this front matter, re-add it after each republication (or re-apply it to every file in ``html/``) before committing.

### From Matlab
The ``gmshGeo`` class is the core of this toolbox. For comprehensive details about it, just type

    help gmshGeo
    
The following command will print all the ``gmshGeo`` methods:

    methods gmshGeo
    
For details about a given method (let say ``plot``):

    help gmshGeo/plot

## Reference
For further details, check out the corresponding paper [[1]](#1). If you use this project, please cite it as follows:

````
@article{MTEX2Gmsh,
  doi = {10.21105/joss.02094},
  url = {https://doi.org/10.21105/joss.02094},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {52},
  pages = {2094},
  author = {Dorian Depriester and R\'egis Kubler},
  title = {MTEX2Gmsh: a tool for generating 2D meshes from EBSD data},
  journal = {Journal of Open Source Software}
}
````

<a id="1">[1]</a> Depriester et al., (2020). MTEX2Gmsh: a tool for generating 2D meshes from EBSD data. *Journal of Open Source Software*, 5(52), 2094, https://doi.org/10.21105/joss.02094
