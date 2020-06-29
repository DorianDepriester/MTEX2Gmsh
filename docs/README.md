## Purpose of this toolbox
In order to evaluate the thermo-mechanical behaviour of crystalline materials (such as metals or ceramics) at microscopic scale, one usually perform numerical simulation at grain scale using the Finite Element Method. In order to proceed, one must first create a mesh which is representative of the real material.

The microstructure of cytalline materials is usually made from Electron Backscattered Diffraction (EBSD) technique. Thus, this toolbox is designed to generate meshes from EBSD in a robust and accurate way.

## Examples
[![Example: aachen.m](./Examples/aachen.png)](./Examples/aachen.png)
[``aachen.m``](Examples/aachen.m)


[![Example: titanium_medium.m](./Examples/titanium_medium.png)](./Examples/titanium_medium.png)
[``titanium_medium.m``](Examples/titanium_medium.m)

[![Example: twins.m](./Examples/twins.png)](./Examples/twins.png)
[``twins.m``](Examples/twins.m)

## Documentation
### Online
You can navigate the documentation [here](html/index.html).

### From MATLAB
#### Full documentation
Once the toolbox is installed on your Matlab session, open the documentation of the present toolbox by typing:

    doc
    
Then, click on "MTEX2Gmsh toolbox", under the _Supplemental Software_ section (bottom right).

#### Help functions
The ``gmshGeo`` class is the core of this toolbox. For comprehensive details about it, just type

    help gmshGeo
    
The following command will print all the ``GmshGeo`` methods:

    methods gmshGeo
    
For details about a given method (let say ``plot``):

    help gmshGeo/plot
