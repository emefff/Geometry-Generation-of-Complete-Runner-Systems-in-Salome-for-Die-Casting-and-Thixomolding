#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.8.0 with dump python functionality
###

import sys
import salome
import numpy as np

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/mario')

###
### GEOM component
###

# import GEOM
from salome.geom import geomBuilder
# from salome.geom import geomtools as gt
from salome.geom.geomtools import GeomStudyTools as gst
import math
# import SALOMEDS

geompy = geomBuilder.New()

print(50*"*")
print("dir(geompy) = ", dir(geompy), "\n")
print(50*"*")
print("dir(gst) = ", dir(gst), "\n")
print(50*"*")

# origin and axes
O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )


def create_trapezoidal_channel(base_length, height, top_length, vector_direction,\
                        channel_length, length_increase, new_base_center, \
                        rot_deg_around_z, scale_factor_XY_connector, name_object):
    """
    Creates a channel with trapezoid cross section with rounded top and a 
    connector with same area (scale_factor_XY_connector==1) or changed area 
    (scale_factor_XY_connector!=1). A name_object can be passed, it will be the
    name of the object in the Salome tree. The channel is ectruded along an axis 
    (most of the time it will be the X-axis). It is rotated
    by an angle rot_deg_around_z and translated together with its connector 
    to new_base_center (the initial base_center is always O, the origin at 
    [0, 0, 0]). The channel length is channel_length, the length_increase is 
    a factor needed for scaling of the end of the channel (usually it is 
    calculated from the area_factor, which is much more useful).     

    Parameters
    ----------
    base_length : TYPE float
        DESCRIPTION. length of the base of the trapezoid = a
    height : TYPE float
        DESCRIPTION. height of trapezoid = h
    top_length : TYPE float
        DESCRIPTION. Length of the top edge of the trapezoid = c
    vector_direction : TYPE list of floats
        DESCRIPTION. direction of ectrusion, usually we will use an axis like
                     the X-axis
    channel_length : TYPE float
        DESCRIPTION. length of the extruded channel
    length_increase : TYPE float
        DESCRIPTION. used for scaling the end of the channel
    new_base_center : TYPE list of 3 floats
        DESCRIPTION. new_base_center in the beginning of the channel trapezoid #
                     base, we translate the channel to this location
    rot_deg_around_z : TYPE float 
        DESCRIPTION. angle in degrees the channel is rotated around the z-axis.
    scale_factor_XY_connector : TYPE float
        DESCRIPTION. Factor that scales the connector in x and y to prevent
                     undeformable surfaces that sometimes occur. 1 will leave it
                     the same size, <1 will decrease the size (mostly not useful)
                     and >1 will increase the connector.
    name_object : TYPE string
        DESCRIPTION. name for the object in the Salome tree.

    Returns
    -------
    None. The functions itself does not return anything but of course, geometries 
          are generated in the Salome tree.

    """
    # Vertex_5 = geompy.MakeVertex(0, base_length/2, 0) # base center
    # Vertex_6 = geompy.MakeVertex(new_base_center[0], new_base_center[1], new_base_center[2]) # new base center
    area_increase = length_increase**2
        
    center_of_mass_y = height/3*(base_length+2*top_length)/(base_length+top_length) # when extruding the rhomboid with length_increase, this is the center in y
    # but this results in the other end also being extruded down into the XY-plane --> we need
    # an additional rotation around X-axis
    deflection_of_end = center_of_mass_y * (length_increase - 1)
    rot_deg_around_y = math.asin(deflection_of_end / channel_length) * 180 / math.pi
    
    # we need a box to cut the face, we want to do a revolution of half the face
    # to make the connections between the channels, we will call these 'Connectors'
    Box_1 = geompy.MakeBoxDXDYDZ(200, 200, 200)
    Box_2 = geompy.MakeTranslation(Box_1, 0, -200, 0) # we us this to cut the Face_1
    
    vector = geompy.MakeVectorDXDYDZ(vector_direction[0], vector_direction[1], vector_direction[2]) # we will always us the X-axis but more is possible
    diff_base_top = base_length - top_length
    Vertex_1 = geompy.MakeVertex(0, -base_length/2, 0)
    Vertex_2 = geompy.MakeVertex(0, base_length/2, 0)
    Vertex_3 = geompy.MakeVertex(0, base_length/2 - diff_base_top/2, height)
    Vertex_4 = geompy.MakeVertex(0, diff_base_top/2 - base_length/2, height)
    Line_1 = geompy.MakeLineTwoPnt(Vertex_1, Vertex_2)
    Line_2 = geompy.MakeLineTwoPnt(Vertex_2, Vertex_3)
    Line_3 = geompy.MakeLineTwoPnt(Vertex_3, Vertex_4)
    Line_4 = geompy.MakeLineTwoPnt(Vertex_4, Vertex_1)
    Face_1 = geompy.MakeFaceWires([Line_1, Line_2, Line_3, Line_4], 1)
    Extrusion_1 = geompy.MakePrismVecH(Face_1, vector, channel_length, length_increase)
    Extrusion_2 = geompy.Rotate(Extrusion_1, OY, -rot_deg_around_y*math.pi/180.0)
    Extrusion_3 = geompy.Rotate(Extrusion_2, OZ, rot_deg_around_z*math.pi/180.0)
    Extrusion_4 = geompy.TranslateDXDYDZ(Extrusion_3, new_base_center[0], new_base_center[1], new_base_center[2])
    Fillet_1 = geompy.MakeFillet(Extrusion_4, 0.5, geompy.ShapeType["EDGE"], [17, 24])
    # geompy.addToStudy( Vertex_1, 'Vertex_1' )
    # geompy.addToStudy( Vertex_2, 'Vertex_2' )
    # geompy.addToStudy( Vertex_3, 'Vertex_3' )
    # geompy.addToStudy( Vertex_4, 'Vertex_4' )
    # geompy.addToStudy( Line_1, 'Line_1' )
    # geompy.addToStudy( Line_2, 'Line_2' )
    # geompy.addToStudy( Line_3, 'Line_3' )
    # geompy.addToStudy( Line_4, 'Line_4' )
    # geompy.addToStudy( Face_1, name_object+'___Face_1' )
    geompy.addToStudy(Fillet_1, name_object)
    
    # generate the connector, a revolution of a half trapezoid
    # to avoid undeformable surfaces, we scale it a little bit
    Face_1_Cut = geompy.MakeCutList(Face_1, [Box_2], True)
    Revolution_1 = geompy.MakeRevolution(Face_1_Cut, OZ, 360*math.pi/180.0)
    # scale_factor_XY_connector = 1.01
    Scale_1 = geompy.MakeScaleAlongAxes(Revolution_1, O, scale_factor_XY_connector, scale_factor_XY_connector, 1)
    Revolution_2 = geompy.TranslateDXDYDZ(Scale_1, new_base_center[0], new_base_center[1], new_base_center[2])
    Fillet_1 = geompy.MakeFillet(Revolution_2, 0.5, geompy.ShapeType["EDGE"], [5])
    # geompy.addToStudy( Revolution_1, name_object+"___Revolution_1")
    geompy.addToStudy( Fillet_1, name_object+"___Connector")
    
    
    area_inlet = (base_length + top_length) / 2 * height
    area_outlet = area_inlet * area_increase
    print("Inlet area of ", name_object, " is", area_inlet)
    print("Outlet area of ", name_object, " is", area_outlet)
    print("Area factor in ", name_object, " is", area_increase)
    
    
def dot_prod(a, b):
    """
    Calculates the dot product of two vectors a and b

    Parameters
    ----------
    a : TYPE list of coords
        DESCRIPTION. first vector
    b : TYPE list of coords
        DESCRIPTION. second vector

    Returns
    -------
    TYPE float
        DESCRIPTION. dot product of a and b.

    """
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


def create_continuous_channel(base_centers, area_first, area_factor_per_seg, \
                              scale_factor_XY_connector, counter):
    """
    Creates a continuous channel from the input parameters. Currently only working
    in the XY-plane, that means Z must be 0 in all base_centers. area_first is the 
    cross sectional area of the first trapezoid, area_factor_per_seg changes the outlet
    area to inlet_area*area_factor_per_seg. For the first segment, that is 
    area_first*area_factor_per_seg. Every channel segment gets a unique name with 
    consecutive numbers in Salome via the number in counter. scale_factor_XY_connector
    lets you scale the connector a bit to prevent surfaces that are not deformable 
    in casting (every surface of a casting must be deformable during die opening,
    that's why every surface must have a draft angle). Currently the shape of the
    trapezoid is hard-coded to c=0.67*a and h=0.67*a. Connectors are frustums with
    rounded tops, channels are rounded at their tops too (currently hard-coded).
       
    Parameters
    ----------
    base_centers : TYPE list of lists of float-triples
        DESCRIPTION. coordinates of the kinks in the cont. channel
    area_first : TYPE float
        DESCRIPTION. area of your first trapezoid
    area_factor_per_seg : TYPE float
        DESCRIPTION. number that is used to calculate the cross section area
                     progression in every segment. For example: if the inlet area
                     to a segment is 50mm² and area_factor_per_seg is 0.9, the 
                     outlet area of this section will be 0.9*50 = 45mm².
    scale_factor_XY_connector : TYPE float
        DESCRIPTION. To avoid undeformable surfaces in the kinks, we can increase
                     the connector (a frustum with rounded top)
    counter : TYPE int
        DESCRIPTION. a number used for naming the created channels and connectors.

    Returns
    -------
    None. The functions itself does not return anything but of course, geometries 
          are generated in the Salome tree.

    """
    for i,point in enumerate(base_centers):
        # all faces are simple trapezoids with area (a+c)/2*h
        # with a being the base_length, c the top_length and h the height 
        # area_first = area of the first face
        # area_factor_per_seg = the area increase per segment
        
        top_length_factor = 0.67 # the ratio that determines the length of c = top_length_factor * base_length
        base_length_first = math.sqrt(3 * area_first / (1 + top_length_factor) )
        top_length_first = top_length_factor * base_length_first
        height_first = base_length_first * 2 / 3
           
        length_increase_per_seg = math.sqrt(area_factor_per_seg) # the length increase
       
        if i == 0: # for the first we have to do some things differently
            vector0 = [1, 0, 0] # x-axis
            vector1 = [base_centers[i+1][0]-base_centers[i][0],\
                       base_centers[i+1][1]-base_centers[i][1],\
                       base_centers[i+1][2]-base_centers[i][2] ]
            channel0_length = math.sqrt(dot_prod(vector0, vector0))
            channel1_length = math.sqrt(dot_prod(vector1, vector1))
            
            # we need to check in which direction we will rotate
            if base_centers[i+1][1] <= base_centers[i][1]:
                angle_rot_XY = - math.acos( dot_prod(vector1, vector0)/( channel0_length*channel1_length) )
                # print(i, " ...", point, "..", vector0," --> ", vector1, " ########1  ", angle_rot_XY*180/math.pi)
            else:
                angle_rot_XY = math.acos( dot_prod(vector1, vector0)/( channel0_length*channel1_length) )
                # print(i, " ...", point, "..", vector0," --> ", vector1, " ########2  ", angle_rot_XY*180/math.pi) 
            
            # we generate the first channel and its connector
            create_trapezoidal_channel(base_length_first, height_first, top_length_first,\
                                [1,0,0], channel1_length, length_increase_per_seg,\
                                base_centers[i], angle_rot_XY*180/math.pi, \
                                scale_factor_XY_connector, 'Channel_'+str(counter)+"___"+str(i) )
            # print(i, " ...", vector0," --> ", vector1, " ########0  ", angle_rot_XY*180/math.pi)
            
        print("")
        if 0 < i < len(base_centers) - 1:
            base_length1 = base_length_first * length_increase_per_seg**i
            top_length1 = top_length_factor * base_length1
            height1 = base_length1 * 2 / 3

            vector0 = [1, 0, 0]
            vector1 = [ base_centers[i+1][0]-base_centers[i][0],\
                        base_centers[i+1][1]-base_centers[i][1],\
                        base_centers[i+1][2]-base_centers[i][2] ]
            channel0_length = math.sqrt(dot_prod(vector0, vector0))
            channel1_length = math.sqrt(dot_prod(vector1, vector1))
            
            # we need to check in which direction we will rotate
            if base_centers[i+1][1] <= base_centers[i][1]:
                angle_rot_XY = - math.acos( dot_prod(vector1, vector0)/( channel0_length*channel1_length) )
                # print(i, " ...", point, "..", vector0," --> ", vector1, " ########1  ", angle_rot_XY*180/math.pi)
            else:
                angle_rot_XY = math.acos( dot_prod(vector1, vector0)/( channel0_length*channel1_length) )
                # print(i, " ...", point, "..", vector0," --> ", vector1, " ########2  ", angle_rot_XY*180/math.pi) 
                           
            # new base center
            base_center_new = base_centers[i]
            # print("base_center_new = ", base_center_new)
            
            # we generate a channel and its connector that connects to its predecessor
            create_trapezoidal_channel(base_length1, height1, top_length1,\
                                [1,0,0], channel1_length, length_increase_per_seg,\
                                base_center_new, angle_rot_XY*180/math.pi, \
                                scale_factor_XY_connector, 'Channel_'+str(counter)+"___"+str(i) )
    

######################## CREATE A TYPICAL RUNNER ##############################
# let's generate a typical casting tree, remember when we seprate into two channels
# we need to halve the inlets, or accordingly, choose the first runner with double
# the outlet section area
# we supply a list of coords where our kinks are.
# level 1
point_list1 = [ [0, 0, 0], [50, 70, 0], [90, 90, 0] ]
point_list8 = [ [89.5, 89.6, 0], [100, 90, 0] ]
# level 2
point_list2 = [ [98, 90, 0], [110, 50, 0] ]
point_list3 = [ [98, 90, 0], [110, 130, 0] ]
point_list9 = [ [110, 50, 0], [115.5, 50, 0] ]
point_list10 = [ [110, 130, 0], [115.5, 130, 0] ]
# level 3
point_list4 = [ [112.5, 50, 0], [140, 30, 0], [145, 30, 0] ]
point_list5 = [ [112.5, 50, 0], [140, 70, 0], [145, 70, 0] ]
point_list6 = [ [115, 130, 0], [140, 110, 0], [145, 110, 0] ]
point_list7 = [ [115, 130, 0], [140, 150, 0], [145, 150, 0] ]

# create the segments of the runner
# level1
create_continuous_channel(point_list1, 100, 0.9, 1.03, 1)
create_continuous_channel(point_list8, 87, 0.51, 1.03, 8)
# level 2
create_continuous_channel(point_list2, 45, 0.8, 1.03, 2)
create_continuous_channel(point_list3, 45, 0.8, 1.03, 3)
create_continuous_channel(point_list9, 42, 0.65, 1.03, 9)
create_continuous_channel(point_list10, 42, 0.65, 1.03, 10)
# level 3
create_continuous_channel(point_list4, 28, 0.9, 1.03, 4)
create_continuous_channel(point_list5, 28, 0.9, 1.03, 5)
create_continuous_channel(point_list6, 28, 0.9, 1.03, 6)
create_continuous_channel(point_list7, 28, 0.9, 1.03, 7)

####################### CREATE A SINGLE CHANNEL ###############################
# with a continuous channel, we do not need a scale_factor_XY_connector other than 1.0
point_list15 = [[0, 0, 0], [50, 25, 0], [100, 30, 0], [150, -50, 0] ]
create_continuous_channel(point_list15, 100, 0.9, 1.0, 15)

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
