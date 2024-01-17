"""
Created on Mon Jan 15 13:49:22 2024

@author: emefff
"""

#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.8.0 with dump python functionality
###

import sys
import salome
# import numpy as np

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


def create_fan(base_center, angle_rot_XY, area_first, length_fan, height_ingate,\
               area_progression_fan_total, counter):
    """
    Creates a fan with constant opening angle according to NADCA. These fans can
    have constant area progression (area_progression_fan_total = 1) or a slightly
    decreasing area (area_progression_fan_total < 1). The location of the base_center
    of the fan (or better its inlet face) is determined by base_center, with angle_rot_XY
    it can be turned by this angle (unit: degrees), area_first sets the approximate
    area of the first face (the inlet face of the fan itself, attention with the naming
    the inlet at the cast part is the last face of the fan :-) ) length_fan sets
    the total length of the fan, height_ingate sets the ingate height at the part
    (this is the last face of the fan), counter is just for naming the fan in 
    Salome.
    
    This function can also create a 'dumb' fan with straight segment edges, but 
    these are commented out at the moment. They are not very useful, because 
    straight segment edges cannot be used with fillets.

    Parameters
    ----------
    base_center : TYPE list of float
        DESCRIPTION. coords of the base_center of the first face of the fan.
    angle_rot_XY : TYPE float
        DESCRIPTION. rotation angle in degrees
    area_first : TYPE float
        DESCRIPTION. area of the first fan face
    length_fan : TYPE float
        DESCRIPTION. total length of the fan
    height_ingate : TYPE float
        DESCRIPTION. height of the ingate, the last face of the fan.
    area_progression_fan_total : TYPE float
        DESCRIPTION. total area progressoin of the fan. if = 1 then area_first
                     is also the area of the last segment. if != 1 then the area
                     of the last segment is area_progression_fan_total * area_first
    counter : TYPE int
        DESCRIPTION. number for naming the fans in Salome

    Returns
    -------
    None. Technically this funtions returns None but a geom object is generated
          in the Salome tree.

    """
   
    # we set a number of segments, the more the better but more expensive
    number_segs_fan = 25
    number_segs_fan += 1 # we need to add 1 o get the actual number entered
    length_per_seg = length_fan / number_segs_fan
                      
    # let's hardcode some lists we need later
    height_per_seg_list = []
    base_centers_per_seg_near_origin_list = []
    base_lengths_per_seg_list = []
    areas_per_seg_list = []
    top_lengths_per_seg_list = []
    
    # determine the area range and steps
    area_ingate = area_progression_fan_total * area_first
    area_progr_last = area_progression_fan_total
    area_progr_per_seg = math.exp(math.log(area_progr_last) / (number_segs_fan-1))
    length_factor_per_seg = math.sqrt(area_progr_per_seg)
    # length_factor_per_seg = 1
    # print("area_progression_total = ", area_progression_fan_total)
    # print("area_progression_per_seg = ", area_progr_per_seg)
    # print("length_factor_per_seg = ", length_factor_per_seg)
        
    # determine the base length ranges, height ranges and top_length ranges
    c_factor_first = 2/3
    base_length_first = math.sqrt( 2 * area_first / (2/3 + 4/9) )      # base_length_first = a_first follows from c_first = 2/3 * a_first and h_first = 2/3*a_first
    c_factor_ingate = 0.98 # toplength of ingate is c_factor_ingate * base_length_ingate
    base_length_ingate = 2 * area_ingate / height_ingate / (1 + c_factor_ingate)
    height_first = base_length_first * 2 / 3
    base_length_step_per_seg = (base_length_first - base_length_ingate) / (number_segs_fan-1)
    
    top_length_first = base_length_first * c_factor_first
    top_length_ingate = base_length_ingate * c_factor_ingate

    print("")
    print("Creating fan nr ", counter)
    print("base_length_first = ", base_length_first,"height_first = ", height_first, ".... area first = ", (base_length_first + top_length_first)/2 * height_first)
    print("base_length_part_ingate = ", base_length_ingate, ".... area part ingate = ", (base_length_ingate + top_length_ingate)/2 * height_ingate)
    print(40*"*")

    top_length_step_per_seg = (base_length_ingate*c_factor_ingate - base_length_first*c_factor_first) / (number_segs_fan-1)
   
    faces_vars_list = []
    
    for i in range(number_segs_fan):
        channel_seg_length = i * length_per_seg
        base_centers_per_seg_near_origin_list.append([channel_seg_length, 0, 0])
        
        area = area_first * area_progr_per_seg** i
        areas_per_seg_list.append(area)
        
        base_length = base_length_first - i * base_length_step_per_seg
        base_lengths_per_seg_list.append(base_length)
        
        top_length = top_length_first + i * top_length_step_per_seg
        top_lengths_per_seg_list.append(top_length)
        
        # we have to calculate h here
        height = 2 * area / (base_length + top_length)
        height_per_seg_list.append(height)
        
        # we test if every value is correct by calculcating the area per segment
        # area_test = (base_length + top_length) / 2 * height
        # print(area_test) ... ok!
        #print("i = ", i,"channel_seg_length", channel_seg_length, "area = ", area, "base_length = ", base_length, \
        #      "top_length = ", top_length, "height = ", height)
        object_id = create_trapezoidal_fan_seg(base_lengths_per_seg_list[i], \
                height_per_seg_list[i], top_lengths_per_seg_list[i], \
                base_centers_per_seg_near_origin_list[i], \
                'Fan_trapezoid_'+str(counter)+"___"+str(i))
        faces_vars_list.append([object_id, 'Fan_trapezoid_'+str(counter)+"___"+str(i)])
    
    
    radius = 0.5
 
    # FOR A SMOOTH SOLID BODY, WE NEED TO BUILD CURVES THROUGH ALL CORNERS AND MAKE A SOLID FROM THESE
    # THEN WE CAN APPLY FILLET WITHOUT PROBLEMS
    # CREATE 4 CURVES TROUGH SOME CORNERS OF ALL THE Fan_segments
    # first curve
    vertex9_vars_list = []
    for i in range(len(faces_vars_list)):
        vertex9_var_name_str = faces_vars_list[i][1]+'_vertex_'+str(i)
        vertex9 = geompy.GetSubShape(faces_vars_list[i][0], [9])
        vertex9_vars_list.append([vertex9, vertex9_var_name_str])
    # we need to make a list of all vars in vertex9_vars[i][0]
    vertex9_curve_list = []
    for i in range(len(vertex9_vars_list)):
        vertex9_curve_list.append(vertex9_vars_list[i][0])
    Curve_1 = geompy.MakeInterpol(vertex9_curve_list, False, False)
    ##geompy.addToStudy( Curve_1, 'Curve_1' )
    # second curve
    vertex7_vars_list = []
    for i in range(len(faces_vars_list)):
        vertex7_var_name_str = faces_vars_list[i][1]+'_vertex_'+str(i)
        vertex7 = geompy.GetSubShape(faces_vars_list[i][0], [7])
        vertex7_vars_list.append([vertex7, vertex7_var_name_str])
    # we need to make a list of all vars in vertex7_vars[i][0]
    vertex7_curve_list = []
    for i in range(len(vertex7_vars_list)):
        vertex7_curve_list.append(vertex7_vars_list[i][0])
    Curve_2 = geompy.MakeInterpol(vertex7_curve_list, False, False)
    ##geompy.addToStudy( Curve_2, 'Curve_2' )
    # third curve
    vertex4_vars_list = []
    for i in range(len(faces_vars_list)):
        vertex4_var_name_str = faces_vars_list[i][1]+'_vertex_'+str(i)
        vertex4 = geompy.GetSubShape(faces_vars_list[i][0], [4])
        vertex4_vars_list.append([vertex4, vertex4_var_name_str])
    # we need to make a list of all vars in vertex4_vars[i][0]
    vertex4_curve_list = []
    for i in range(len(vertex4_vars_list)):
        vertex4_curve_list.append(vertex4_vars_list[i][0])
    Curve_3 = geompy.MakeInterpol(vertex4_curve_list, False, False)
    ##geompy.addToStudy( Curve_3, 'Curve_3' )
    # fourth curve
    vertex5_vars_list = []
    for i in range(len(faces_vars_list)):
        vertex5_var_name_str = faces_vars_list[i][1]+'_vertex_'+str(i)
        vertex5 = geompy.GetSubShape(faces_vars_list[i][0], [5])
        vertex5_vars_list.append([vertex5, vertex5_var_name_str])
    # we need to make a list of all vars in vertex5_vars[i][0]
    vertex5_curve_list = []
    for i in range(len(vertex5_vars_list)):
        vertex5_curve_list.append(vertex5_vars_list[i][0])
    Curve_4 = geompy.MakeInterpol(vertex5_curve_list, False, False)
    ##geompy.addToStudy( Curve_4, 'Curve_4' )
    
    # NOW WE HAVE TO BUILD WIRES FOR EVERY FACE WE NEED FOR THE SHELL AND SUBSEQUENTLY, THE SOLID
    # OBVIOUSLY, WE NEED 4 WIRES FOR 6 FACES OF THE SOLID, BECAUSE WE ALREADY HAVE 2 (FIRST AND LAST TRAPEZOID)
    # wire 1
    # for this wire 1 we need the edge_8 of the first and the last trapezoid
    First_edge_8 = geompy.GetSubShape(faces_vars_list[0][0], [8])
    Last_edge_8 = geompy.GetSubShape(faces_vars_list[-1][0], [8])
    Wire_1 = geompy.MakeWire([First_edge_8, Last_edge_8, Curve_1, Curve_2], 0.001)
    ##geompy.addToStudy( Wire_1, 'Wire_1' )
    # wire 2
    # for this wire 2 we need the edge_6 of the first and the last trapezoid
    First_edge_6 = geompy.GetSubShape(faces_vars_list[0][0], [6])
    Last_edge_6 = geompy.GetSubShape(faces_vars_list[-1][0], [6])
    Wire_1_edge_9 = geompy.GetSubShape(Wire_1, [9])
    Wire_2 = geompy.MakeWire([First_edge_6, Last_edge_6, Curve_4, Wire_1_edge_9], 0.001)
    ##geompy.addToStudy( Wire_2, 'Wire_2' )
    # wire 3
    # for this wire 3 we need the edge_10 of the first and the last trapezoid
    First_edge_10 = geompy.GetSubShape(faces_vars_list[0][0], [10])
    Last_edge_10 = geompy.GetSubShape(faces_vars_list[-1][0], [10])
    Wire_1_edge_5 = geompy.GetSubShape(Wire_1, [5])
    Wire_3 = geompy.MakeWire([First_edge_10, Last_edge_10, Curve_3, Wire_1_edge_5], 0.001)
    ##geompy.addToStudy( Wire_3, 'Wire_3' )
    # wire 4
    # for this wire 4 we need the edge_3 of the first and the last trapezoid
    First_edge_3 = geompy.GetSubShape(faces_vars_list[0][0], [3])
    Last_edge_3 = geompy.GetSubShape(faces_vars_list[-1][0], [3])
    Wire_2_edge_9 = geompy.GetSubShape(Wire_2, [9])
    Wire_3_edge_5 = geompy.GetSubShape(Wire_3, [5])
    Wire_4 = geompy.MakeWire([First_edge_3, Last_edge_3, Wire_2_edge_9, Wire_3_edge_5], 0.001)
    ##geompy.addToStudy( Wire_4, 'Wire_4' )
    
    # WE BUILD ALL $ REMAINING FACES
    Face_1 = geompy.MakeFaceWires([Wire_1], 0)
    Face_2 = geompy.MakeFaceWires([Wire_2], 0)
    Face_3 = geompy.MakeFaceWires([Wire_3], 0)
    Face_4 = geompy.MakeFaceWires([Wire_4], 0)
    
    # WE BUILD A SHELL WITH THE 4 NEW FACES AND THE FIRST TRAPEZOID AND THE LAST
    First_Trapezoid = faces_vars_list[0][0]
    Last_Trapezoid = faces_vars_list[-1][0]
    
    Shell_1 = geompy.MakeShell([First_Trapezoid, Last_Trapezoid, Face_1, Face_2, Face_3, Face_4])
    Solid_1 = geompy.MakeSolid([Shell_1])
    Fillet_1 = geompy.MakeFillet(Solid_1, radius, geompy.ShapeType["EDGE"], [15, 19])
    
    Rotation_1 = geompy.MakeRotation(Fillet_1, OZ, angle_rot_XY*math.pi/180.0)
    Translation_1 = geompy.TranslateDXDYDZ(Rotation_1, base_center[0], base_center[1], base_center[2])
    geompy.addToStudy( Translation_1, 'Fan_'+str(counter) )

    # PRINT SOME LISTS FOR CHECKING
    # print("")
    # print("base_lengths_list = ", base_lengths_per_seg_list)
    # print("")
    # print("top_lengths_list = ", top_lengths_per_seg_list)
    # print("")
    # print("heights_list = ", height_per_seg_list)
    # print("")
    # print("areas_list = ", areas_per_seg_list)
    # print("")
        

def create_trapezoidal_fan_seg(base_length, height, top_length, \
                        base_center, \
                        name_object):
    """
    Creates a trapezoidal fan segment for the constant cross-section fan.     

    Parameters
    ----------
    base_length : TYPE float
        DESCRIPTION. base length of the trapezoid
    height : TYPE float
        DESCRIPTION. height of the trapezoid
    top_length : TYPE float
        DESCRIPTION. top length of the trapezoid
    base_center : TYPE list of floats
        DESCRIPTION. coords of the base center of the fan segment or fan face
    name_object : TYPE str
        DESCRIPTION. name of the generated object.

    Returns
    -------
    Face_2 : TYPE Salome object
        DESCRIPTION. the returned Salome object is needed to keep track of it. 

    """
    # we only create a trapezoid and return it
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
    Face_2 = geompy.TranslateDXDYDZ(Face_1, base_center[0], base_center[1], base_center[2])
    
    area_inlet = (base_length + top_length) / 2 * height
    # print("")
    # print("Inlet area of ", name_object, " is", area_inlet)
    # print("*****************************************")
    return Face_2


def create_trapezoidal_channel(base_length, height, top_length, vector_direction,\
                        channel_length, length_factor, base_center, \
                        rot_deg_around_z, scale_factor_XY_connector, connector_at_end, \
                        name_object):
    """
    Creates a channel with trapezoid cross section with rounded top and a 
    connector with same area (scale_factor_XY_connector==1) or changed area 
    (scale_factor_XY_connector!=1). A name_object can be passed, it will be the
    name of the object in the Salome tree. The channel is ectruded along an axis 
    (most of the time it will be the X-axis). It is rotated
    by an angle rot_deg_around_z and translated together with its connector 
    to base_center (the initial base_center is always O, the origin at 
    [0, 0, 0]). The channel length is channel_length, the length_factor is 
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
    length_factor : TYPE float
        DESCRIPTION. used for scaling the end of the channel
    base_center : TYPE list of 3 floats
        DESCRIPTION. base_center in the beginning of the channel trapezoid #
                     base, we translate the channel to this location
    rot_deg_around_z : TYPE float 
        DESCRIPTION. angle in degrees the channel is rotated around the z-axis.
    scale_factor_XY_connector : TYPE float
        DESCRIPTION. Factor that scales the connector in x and y to prevent
                     undeformable surfaces that sometimes occur. 1 will leave it
                     the same size, <1 will decrease the size (mostly not useful)
                     and >1 will increase the connector.
    connector_at_end: boolean
        DESCRIPTION. if True, additional connector (frustum) will be placed at 
                     end of channel.
    name_object : TYPE string
        DESCRIPTION. name for the object in the Salome tree.

    Returns
    -------
    None. The functions itself does not return anything but of course, geometries 
          are generated in the Salome tree.

    """
    # Salome cannot deal with very short channels, so we strech the channel by a
    # factor and shrink it back later. This is useful especially if channels are
    # short and the area progression is harsh (for example length 10mm and length_factor=0.5)
    
    stretch_factor_X = 100
    radius = 0.5
    
    # Vertex_5 = geompy.MakeVertex(0, base_length/2, 0) # base center
    # Vertex_6 = geompy.MakeVertex(base_center[0], base_center[1], base_center[2]) # new base center
    area_increase = length_factor**2
        
    center_of_mass_y = height/3*(base_length+2*top_length)/(base_length+top_length) # when extruding the rhomboid with length_factor, this is the center in y
    # but this results in the other end also being extruded down into the XY-plane --> we need
    # an additional rotation around X-axis
    deflection_of_end = center_of_mass_y * (length_factor - 1)
    rot_deg_around_y = math.asin(deflection_of_end / channel_length) * 180 / math.pi / stretch_factor_X
    
    # we need a box to cut the face, we want to do a revolution of half the face
    # to make the connections between the channels, we will call these 'Connectors'
    Box_1 = geompy.MakeBoxDXDYDZ(200, 200, 200)
    
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
    # extrude much longer than necessary due to Salome not being able to handle very short channels
    Extrusion_1 = geompy.MakePrismVecH(Face_1, vector, channel_length*stretch_factor_X, length_factor)
    Extrusion_2 = geompy.Rotate(Extrusion_1, OY, -rot_deg_around_y*math.pi/180.0)

    # scaling back with 1/stretch_factor_X
    Scale_1 = geompy.MakeScaleAlongAxes(Extrusion_2, O, 1/stretch_factor_X, 1, 1)

    Extrusion_3 = geompy.Rotate(Scale_1, OZ, rot_deg_around_z*math.pi/180.0)
    Extrusion_4 = geompy.TranslateDXDYDZ(Extrusion_3, base_center[0], base_center[1], base_center[2])
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
    
    # generate the connector(s), a revolution of a half trapezoid
    # to avoid undeformable surfaces, we scale it a little bit
    Vertex_1 = geompy.MakeVertex(channel_length, 0, 0) # this is the endpoint in the non-rotated and translated channel segment
    Face_1_Cut = geompy.MakeCutList(Face_1, [Box_1], True)
    Revolution_1 = geompy.MakeRevolution(Face_1_Cut, OZ, 360*math.pi/180.0)
    Scale_1 = geompy.MakeScaleAlongAxes(Revolution_1, O, scale_factor_XY_connector, scale_factor_XY_connector, 1)

    Revolution_2 = geompy.TranslateDXDYDZ(Scale_1, base_center[0], base_center[1], base_center[2])
    Fillet_1 = geompy.MakeFillet(Revolution_2, radius, geompy.ShapeType["EDGE"], [5])
    geompy.addToStudy( Fillet_1, name_object+"___Connector_1")
    if connector_at_end == True: # we also put a connector at the end of the channel if connector_at_end = True
        Scale_2 = geompy.MakeScaleAlongAxes(Revolution_1, O, scale_factor_XY_connector, scale_factor_XY_connector, 1) # we take the same connector for the end
        Translation_2 = geompy.TranslateDXDYDZ(Scale_2, channel_length, 0, 0) # shift it to the end
        Rotation_2 = geompy.Rotate(Translation_2, OZ, rot_deg_around_z*math.pi/180.0) # rotate it 
        Translation_3 = geompy.TranslateDXDYDZ(Rotation_2, base_center[0], base_center[1], base_center[2]) # translate to end point
        Fillet_2 = geompy.MakeFillet(Translation_3, radius, geompy.ShapeType["EDGE"], [5]) # for whatever reason, we cannot just take Scale_1 ?? WTF
        geompy.addToStudy( Fillet_2, name_object+"___Connector_2") # is larger than Fille_1 because one is generated with the inlet area, the other with
        # the possibliy smaller outlet area (if length_factor < 1), not pretty

    area_inlet = (base_length + top_length) / 2 * height
    area_outlet = area_inlet * area_increase
    print("")
    print("Creating ", name_object, " and" , name_object+"___Connector")
    print("Inlet area of ", name_object, " is", area_inlet)
    print("Outlet area of ", name_object, " is", area_outlet)
    print("Area factor in ", name_object, " is", area_increase)
    print(40*"*")
    
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
                              scale_factor_XY_connector, connector_at_end, counter):
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
    connector_at_end : boolean
        DESCRIPTION. if True, additional connector will be placed at end of channel.
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
           
        length_factor_per_seg = math.sqrt(area_factor_per_seg) # the length increase
        
        if 0 <= i < len(base_centers) - 1:
            base_length1 = base_length_first * length_factor_per_seg**i
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
                                [1,0,0], channel1_length, length_factor_per_seg,\
                                base_center_new, angle_rot_XY*180/math.pi, \
                                scale_factor_XY_connector, connector_at_end, 'Channel_'+str(counter)+"___"+str(i) )
    

###############################################################################
######################## CREATE A TYPICAL RUNNER ##############################
# let's generate a typical casting tree, remember when we seprate into two channels
# we need to halve the inlets, or accordingly, choose the first runner with double
# the outlet section area
# we supply a list of coords where our kinks are.
# level 1
origin = [0, 0, 0]
kink01  = [25, 50, 0]
kink02  = [55, 80, 0]
kink03 =  [80, 90, 0]
kink04 =  [80, 90, 0]
point_list01 = [ origin, kink01, kink02, kink03 ]
# point_list08 = [ kink03, kink04 ]
# level 2
kink05 = [80, 90, 0]
kink06 = [112.5, 50, 0]
kink07 = [112.5, 130, 0] 
kink08 = [112.5, 50, 0]
kink09 = [112.5, 130, 0]
point_list02 = [ kink05, kink06 ]
point_list03 = [ kink05, kink07 ]

# level 3
kink10 = [112.5, 50, 0]
kink11 = [140, 30, 0]
kink12 = [112.5, 50, 0]
kink13 = [140, 70, 0]
kink14 = [115, 130, 0]
kink15 = [140, 110, 0]
kink16 = [140, 150, 0]
kink17 = [145, 30, 0]
kink18 = [145, 70, 0]
kink19 = [145, 110, 0]
kink20 = [145, 150, 0]

point_list04 = [ kink10, kink11, kink17 ]
point_list05 = [ kink12, kink13, kink18 ]
point_list06 = [ kink07, kink15, kink19 ]
point_list07 = [ kink07, kink16, kink20 ]

# create the segments of the runner
# level1
create_continuous_channel(point_list01, 100, 0.9, 1.03, True, 1)

# level 2
create_continuous_channel(point_list02, 45, 1, 1.03, True, 2)
create_continuous_channel(point_list03, 45, 1, 1.03, True, 3)

# level 3
create_continuous_channel(point_list04, 28, 0.9, 1.03, False, 4)
create_continuous_channel(point_list05, 28, 0.9, 1.03, False, 5)
create_continuous_channel(point_list06, 28, 0.9, 1.03, False, 6)
create_continuous_channel(point_list07, 28, 0.9, 1.03, False, 7)

####################### CREATE A SINGLE CHANNEL ###############################
# with a continuous channel, we do not need a scale_factor_XY_connector other than 1.0
# point_list015 = [[10, 10, 0], [50, 25, 0], [100, 30, 0], [150, -50, 0] ]
# create_continuous_channel(point_list015, 100, 0.9, 1.0, True, 15)

################### CREATE THE INLET FANS FOR THE RUNNER ######################
create_fan(kink17, 0, 22.68, 25, 1, 0.9, 1) # for area_first just chose last area of channel!
create_fan(kink18, 0, 22.68, 20, 1.2, 0.9, 2) # values are printed in Salome Python shell
create_fan(kink19, 0, 22.68, 20, 1.2, 0.9, 3)
create_fan(kink20, 0, 22.68, 25, 1, 0.9, 4)

# ###############################################################################
# ######################### MAKING SOME CORRECTIONS #############################
# # MAKING CORRECTIONS IS A PIECE OF CAKE
# the inner runner+fans are shifted by +-2.5mm to the center
# level 1
origin = [0, 0, 0]
kink01  = [25, 50, 0]
kink02  = [55, 80, 0]
kink03 =  [80, 90, 0]
kink04 =  [80, 90, 0]
point_list01 = [ origin, kink01, kink02, kink03 ]
# point_list08 = [ kink03, kink04 ]
# level 2
kink05 = [80, 90, 0]
kink06 = [112.5, 50, 0]
kink07 = [112.5, 130, 0] 
kink08 = [112.5, 50, 0]
kink09 = [112.5, 130, 0]
point_list02 = [ kink05, kink06 ]
point_list03 = [ kink05, kink07 ]

# level 3
kink10 = [112.5, 50, 0]
kink11 = [140, 30, 0]
kink12 = [112.5, 50, 0]
kink13 = [140, 72.5, 0]
kink14 = [115, 130, 0]
kink15 = [140, 107.5, 0]
kink16 = [140, 150, 0]
kink17 = [145, 30, 0]
kink18 = [145, 72.5, 0]
kink19 = [145, 107.5, 0]
kink20 = [145, 150, 0]

point_list04 = [ kink10, kink11, kink17 ]
point_list05 = [ kink12, kink13, kink18 ]
point_list06 = [ kink07, kink15, kink19 ]
point_list07 = [ kink07, kink16, kink20 ]

# create the segments of the runner
# level1
create_continuous_channel(point_list01, 100, 0.9, 1.03, True, 1)

# level 2
create_continuous_channel(point_list02, 45, 1, 1.03, True, 2)
create_continuous_channel(point_list03, 45, 1, 1.03, True, 3)

# level 3
create_continuous_channel(point_list04, 28, 0.9, 1.03, False, 4)
create_continuous_channel(point_list05, 28, 0.9, 1.03, False, 5)
create_continuous_channel(point_list06, 28, 0.9, 1.03, False, 6)
create_continuous_channel(point_list07, 28, 0.9, 1.03, False, 7)

####################### CREATE A SINGLE CHANNEL ###############################
# with a continuous channel, we do not need a scale_factor_XY_connector other than 1.0
# point_list015 = [[10, 10, 0], [50, 25, 0], [100, 30, 0], [150, -50, 0] ]
# create_continuous_channel(point_list015, 100, 0.9, 1.0, True, 15)

################### CREATE THE INLET FANS FOR THE RUNNER ######################
create_fan(kink17, 0, 22.68, 25, 1, 0.9, 1) # for area_first just chose last area of channel!
create_fan(kink18, 0, 22.68, 20, 1.2, 0.9, 2) # values are printed in Salome Python shell
create_fan(kink19, 0, 22.68, 20, 1.2, 0.9, 3)
create_fan(kink20, 0, 22.68, 25, 1, 0.9, 4)


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()