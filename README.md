# Geometry-Generation-of-Runners-in-Salome-for-Casting
This script can give a very quick sketch or first version of a runner system for use in die-casting. Currently, and perhaps for a long time (because it is VERY difficult), the generation of such a runner system is only possible in the XY-plane. The runner has the typical trapezoidal shape. As cross-sectional areas of such runners are a very important factor, this is dealt with setting the area of the first segment and the area progression per segment (That is: if the area of the inlet trapezoid is 50mm² and the area_factor_per_seg = 0.9, the outlet of this segment will be 45mm²). A continous channel can be created via supplying a list of coordinates for the kinks, the area of the first inlet, area_factor_per_seg, scale_factor_XY_connector and a counter for naming the geometry objects. The kinks will then be connected with trapezoidal channels and so-called 'connectors' that fill the wedge-shaped voids in the channel with frustrums. These connectors can also be scaled via scale_factor_XY_connector to prevent any undeformable surfaces that may arise when, for example, arise when area_factor_per_seg is very small and the channel is very short (somehow Salome cannot deal with very short extruded channels). Some important numbers are printed to the Python shell in Salome. Run this script via laoding it in Salome.
Real casting trees can be created with adding up several continuous channels. Such an example is also in the script. As Salome still uses PYthon 3.6 we need to keep that in mind (no f-strings for example).

The following image shows an example for a continuous channel:

![Bildschirmfoto vom 2024-01-11 11-47-19](https://github.com/emefff/Geometry-Generation-of-Runners-in-Salome-for-Casting/assets/89903493/5f39e7bd-4482-46ef-98c6-a40bae138d3c)

A casting tree example might look like the follwing image. As alyways, the branching from one channel into two is problematic to get right (cross section areas must add up!):

![Bildschirmfoto vom 2024-01-11 11-49-04](https://github.com/emefff/Geometry-Generation-of-Runners-in-Salome-for-Casting/assets/89903493/5a017af1-96c7-4183-9825-0531dc8c583b)

emefff@gmx.at
