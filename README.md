# Geometry-Generation-of-Complete-Runner-Systems-in-Salome-for-Die-Casting-and-Thixomolding
This script can give a very quick sketch or first version of a runner system for use in die-casting. Currently, and perhaps for a long time (because it is VERY difficult), the generation of such a runner system is only possible in the XY-plane. The runner has the typical trapezoidal shape with rounded top edges. As cross-sectional areas of such runners are a very important factor, this is dealt with setting the area of the first segment and the area progression per segment (That is: if the area of the inlet trapezoid is 50mm² and the area_factor_per_seg = 0.9, the outlet of this segment will be 45mm²). A continous channel can be created via supplying a list of coordinates for the kinks, the area of the first inlet, area_factor_per_seg, scale_factor_XY_connector and a counter for naming the geometry objects. The kinks will then be connected with trapezoidal channels and so-called 'connectors' that fill the wedge-shaped voids in the channel with frustrums. These connectors can also be scaled via scale_factor_XY_connector to prevent any undeformable surfaces that may arise when, for example, arise when area_factor_per_seg is very small and the channel is very short (very short channels and steep area_progression lead to weard extruded bodies). Some important numbers are printed to the Python shell in Salome. Due to the areas being calculated before applying 3D-fillet, the areas in the console and those of the model differ a bit. Run this script via loading it in Salome.
Real casting trees can be created with adding up several continuous channels. Such an example is also in the script. As Salome still uses Python 3.6 we need to keep that in mind (no f-strings etc.).

The following image shows an example for a continuous channel generated from a list of points:

![Bildschirmfoto vom 2024-01-11 11-47-19](https://github.com/emefff/Geometry-Generation-of-Runners-in-Salome-for-Casting/assets/89903493/5f39e7bd-4482-46ef-98c6-a40bae138d3c)

The latest version of this script can also create inlet fans according to NADCA (we prefer those with constant angle). It is very important that such inlet fans have constant area or a defined area progression from first to last face. Just connecting the more quadratic first face with the last elongated and thin face is a very common mistake in runner design. We have seen it MANY TIMES from other designers, but it creates an area maximum somewhere near the center. The excess air present in this volume is blown into the part! A complete runner system with such ingates may look like this, the dummy part is just a 2.5mm plate with a step:

![Bildschirmfoto vom 2024-01-17 19-19-22](https://github.com/emefff/Geometry-Generation-of-Complete-Runner-Systems-in-Salome-for-Die-Casting-and-Thixomolding/assets/89903493/e2f816dc-0feb-4966-84db-db1e2009bcb0)

The following image shows the typical shape of such an ingate fan with constant area progression and constant opening angle. The shape is computed from 25 cross sections (is of course variable), however Salome does some interpolations on the curves for viewing. The nonlinear height is a typical feature of these types an ingate fans:

![Bildschirmfoto vom 2024-01-15 16-14-26](https://github.com/emefff/Geometry-Generation-of-Runners-in-Salome-for-Die-Casting/assets/89903493/8efb3a9b-bd0e-45fc-bb20-67f38fe204ca)

The included video shows how easy corrections can be made (however, it shows an older version of the model), both inner runners and ingates are moved +-2.5mm to the center. As almost everything is parametrized, corrections only take seconds instead of hours (in case the design engineer did not use any parametrization). Here is a screenshot after the corrections:

![Bildschirmfoto vom 2024-01-17 19-22-48](https://github.com/emefff/Geometry-Generation-of-Complete-Runner-Systems-in-Salome-for-Die-Casting-and-Thixomolding/assets/89903493/d405889a-6dc3-4ecc-a378-d1443e173644)


Your design engineer won't be happy, but that's the way it is. This script may also be suitable for injection molding, but we are absolutely no experts for that. Use at your own risk.

emefff@gmx.at
