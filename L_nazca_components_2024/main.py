# -*- coding: utf-8 -*-
"""
Created on Thu May  2 2024
Last edit  Thu Jun 27 2024

Loads Nazca and makes layers, components and other useful functions for
photonic circuit creation. Includes markers, straight and bent waveguides, beam
splitters, racetrack and ring resonators of straight, symmetric and pulley design,
recctangle and circular spiral delay lines and thermo-optic and electro-optic
phase shifters. Most components have two associated functions. The one that ends with
'_' (B) might be less intuitive, so the one without '_' (A) is more user-friendly.
'A' always calls 'B'. 'B' calls Nazca and makes the component's cell. Call function
'A' to get the cell, and then put the cell with '.put(x,y,a)'. Most points are
'pins', that is, they have an x coordinate, y coordinate and an angle a in
degrees. Angle 0 points to the right in exported files. To export a file,
we recommend 'nd.export_gds()', where nd is the Nazca module.

@author: LeoKrisK
"""
# All lengths and coordinates are in microns (μm).
# Cells are Nazca cells.

import numpy as np
import nazca as nd
import nazca.demofab as demo
from nazca.interconnects import Interconnect
import random

###############################################################################
############################# Abbreviations ###################################
###############################################################################
# bs    - beam splitter
# ps    - phase shifter
# i     - input
# o     - output

# xs    - cross-section
# ls    - layers
# w     - width
# l     - length (usually used for straight, horizontal lines)
# hi    - high / upper / top / above
# lo    - low / lower / bottom / below
# m     - mid / middle
# co    - coupled (region)
# tk    - (race)track
# s     - srt / strt / straight
# b     - bend / bent
# p     - pt(s) / point(s)  or  pos / position
# r     - radius  or  real
# a     - ang / angle
# os    - offset
# block - component of circuit
# con   - conductor / conductive
# d     - delta

# NxN   - related to N times N unitary matrix circuit
# NxNm  - NxN unitary matrix circuit with Mach-Zehnder interferometer design
# oa    - one_all
# aoa   - all_one_all

### For pins:
# 'a0' is the first input pin
# 'a1' is the second input pin, etc.
# 'b0' is the first output pin
# 'b1' is the second output pin, etc.

### btype can be:
# 'sbend' - an s-bend
# 'sin'   - a sine bend
# 'cobra' - a cobra ben
# 'bsb'   - bend-straight-bend
# More can be added
###############################################################################
############################# Define Layers ###################################
###############################################################################
def load_all_layers(all_layers):
    """
    Adds all layers in 'all_layers' dict.
    """
    for l_str, l_code in all_layers.items():
        nd.add_layer(name=l_str, layer=l_code)

def make_xsection(xs_name, xs_ls):
    """
    Creates a cross-section named 'xs_name', which has the layers found in
    xs_ls[xs_name] (a list of layer names).

       Parameters
    xs_name (str): The xsection's name.
    layers (dict of list of strings): Dictionary. The key is the xs name, the
        value is a list of the layers in the xs.
    """
    layers = xs_ls[ xs_name ] # Gets the specific layer list
    nd.add_xsection(name=xs_name)
    for l_str in layers:
        nd.add_layer2xsection(xsection=xs_name, layer=l_str)

###############################################################################
########################### Generic Functions #################################
###############################################################################
def get_n_str():
    """
    Returns a random number that increases with each call. Used to avoid
    generating components with the same name. By default, each cell gets a
    number added at the end.
    """
    return ' '+str(random.random())

def pythagoras(p1, p2):
    """
    Takes two points p1, p2 of the form (x,y) or (x,y,a) and returns the
    Euclidean distance between them.
    """
    return ((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)**0.5

def atan(p1, p2):
    """
    Takes two points p1, p2 of the form (x,y) or (x,y,a) and returns
    arctan(delta_y / delta_x).
    """
    ang = None
    if p2[0] != p1[0]:
        ang = np.arctan((p2[1] - p1[1])/(p2[0] - p1[0]))
    else:
        ang = np.pi/2
    return ang

def polar_to_cart(r, ang):
    """
    Converts polar coordinates (r,angle) to Cartesian coordinates (x,y).
    """
    x, y = r * np.cos(ang), r * np.sin(ang)
    return x, y

def add_pos(p1, p2, minus=False, add_a=True):
    """
    Adds together points p1 and p2 of the form (x,y) or (x,y,a) component-wise,
    or, if minus=True, takes p1-p2. If add_a=False, returns p1[2] as the angle
    instead of p1[2] +- p2[2].
    """
    p_sum = None
    if len(p1) == 3 and len(p2) == 3:
        if not minus:
            if add_a:
                p_sum = (p1[0]+p2[0], p1[1]+p2[1], p1[2]+p2[2])
            else:
                p_sum = (p1[0]+p2[0], p1[1]+p2[1], p1[2])
        else:
            if add_a:
                p_sum = (p1[0]-p2[0], p1[1]-p2[1], p1[2]-p2[2])
            else:
                p_sum = (p1[0]-p2[0], p1[1]-p2[1], p1[2])
    else: # Assumes p1 and p2 have (at least) two components
        if not minus:
            p_sum = (p1[0]+p2[0], p1[1]+p2[1])
        else:
            p_sum = (p1[0]-p2[0], p1[1]-p2[1])
    return p_sum

def rotate_pos(pos, rot_a, pos0=(0,0,0), rotate_a=False):
    """
    Rotates point pos=(x,y,a) by some angle rot_a (in degrees) around pos0.
    Mostly used to change between coordinate systems. If rotate_a=True, adds
    rot_a to angle a. If rotate_a=False, a stays unchanged.
    """
    # Code inspired by LyleScott on
    # https://gist.github.com/LyleScott/e36e08bfb23b1f87af68c9051f985302
    rot_a_rad = rot_a/180*np.pi # Converts degrees to radians
    x, y, a    = pos
    x0, y0, a0 = pos0
    
    xx = x0 + (x-x0) * np.cos(rot_a_rad) + (y-y0) * np.sin(rot_a_rad)
    yy = y0 - (x-x0) * np.sin(rot_a_rad) + (y-y0) * np.cos(rot_a_rad)
    
    aa = a
    if rotate_a:
        aa += rot_a
    return (xx,yy,aa)

def check_unitary(A, tolerance=1e-14, st=False):
    """
    Checks if A is a unitary matrix of dimension NxN. Returns True if each
    element of A*A_dagger and A_dagger*A is identity to within the error given
    by 'tolerance'. If st=True, may also print further information.
    """
    N = int(np.sqrt(A.size))
    A_dagg = np.conjugate(np.transpose(A))
    Ident_1 = np.matmul(A, A_dagg)
    Ident_2 = np.matmul(A_dagg, A)    
    almost_unit = True
    for i in range(N):
        for j in range(N):
            if i != j: # Tests off-diagonal elements
                if abs(Ident_1[i,j]) > tolerance:
                    almost_unit = False
                if abs(Ident_2[i,j]) > tolerance:
                    almost_unit = False

            if almost_unit == True and abs(Ident_1[i,j] - Ident_2[i,j]) > tolerance:
                # Compares Ident_1 and Ident_2
                almost_unit = False
    
    is_unit = almost_unit
    
    scale_fact = Ident_1[0,0]
    Ident_1_nor, Ident_2_nor = Ident_1, Ident_2
    if scale_fact != 0:
        # 'Normalizing' the matrices
        Ident_1_nor = Ident_1 / scale_fact
        Ident_2_nor = Ident_2 / scale_fact
    
    for i in range(N):
        # Tests if the normalized matrices are identity
        if abs(Ident_1_nor[i,i]-1) > tolerance or abs(Ident_2_nor[i,i]-1) > tolerance:
            almost_unit = False
    for i in range(N):
        # Tests if the original matrices are identity
        if abs(Ident_1[i,i]-1) > tolerance or abs(Ident_2[i,i]-1) > tolerance:
            is_unit = False
    if is_unit:
        almost_unit = False
    
    if st:
        if not is_unit:
            print('A * A_dagger = ')
            print(Ident_1)
            print('A_dagger * A = ')
            print(Ident_2)
        if abs(scale_fact-1) > tolerance and almost_unit and scale_fact != 0:
            print('A is unitary if the entire matrix is divided by',
                  np.sqrt(scale_fact), '.')
    return is_unit

###############################################################################
########################## Component Functions ################################
###############################################################################
def make_plus_(layers, x_wide, y_wide, x_tall=None, y_tall=None, st=False):
    """
    Creates and returns a cell of a marker that looks like a plus sign.

       Parameters
    layers (list of strings) : The layers.
    x_wide (float) : Length in x of the wide rectangle.
    y_wide (float) : Length in y of the wide rectangle.
    x_tall (float, optional) : Length in x of the tall rectangle.
    y_tall (float, optional) : Length in y of the tall rectangle.
        If x_tall and y_tall are not defined, makes the plus symmetric.
    st (bool, optional) : If True, prints notes and warnings.
    """
    if x_tall == None or y_tall == None:
        x_tall = y_wide
        y_tall = x_wide
    with nd.Cell(name='Plus Marker'+get_n_str()) as PM_cell:
        pts_x = nd.geometries.rectangle(length=x_wide, height=y_wide)
        pts_y = nd.geometries.rectangle(length=x_tall, height=y_tall)
        
        # So the origin is in the middle of the plus
        pos1 = ( (x_tall-x_wide)/2,      0, 0 )
        pos2 = ( 0,      (y_wide-y_tall)/2, 0 )
        
        for l_str in layers:
            elm1 = nd.Polygon(points=pts_x, layer=l_str).put(pos1)
            elm2 = nd.Polygon(points=pts_y, layer=l_str).put(pos2)
    return PM_cell

def make_tp_(layers, w_in, w_out, l, st=False):
    """
    Creates and returns a taper cell.

       Parameters
    layers (list of strings) : The layers.
    w_in  (float) : Width at the input i (usually left).
    w_out (float) : Width at the output o (usually right).
    l (float) : Length from i to o (usually in x).
    st (bool, optional) : If True, prints notes and warnings.
    """
    with nd.Cell('Taper'+get_n_str()) as T_cell:
        pts = nd.geometries.taper(length=l, width1=w_in, width2=w_out)
        for l_str in layers:
            elm1 = nd.Polygon(points=pts, layer=l_str).put(0,0,180)
        
        # Adds pins (manually)
        nd.Pin("a0").put(0,0,0)
        nd.Pin("b0").put(-l,0,180)
    return T_cell

def make_strt_object_(xs, w, l, is_con=False, st=False):
    """
    Creates and returns a straight waveguide (or conductor) cell.

       Parameters
    xs (str)  : Cross-section, specifying the layers.
    w (float) : Width (usually in y).
    l (float) : Length (usually in x).
    is_con (bool, optional) : If True, changes name to 'Straight Conductor'.
    st (bool, optional) : If True, prints notes and warnings.
    """
    ic = Interconnect(xs=xs)
    cname = 'Straight Waveguide'
    if is_con:
        cname = 'Straight Conductor'
    with nd.Cell(name=cname+get_n_str()) as SO_cell:
        elm1 = ic.strt(length=l, width=w).put(0)
        # Adds pins
        nd.Pin("a0", pin=elm1.pin['a0']).put()
        nd.Pin("b0", pin=elm1.pin['b0']).put()
    return SO_cell

def mbo_update_name_(cname, btype, st=False):
    """
    Adds the bend type to the name of a bent waveguide (or conductor).

       Parameters
    cname (str) : Generic name of the object.
    btype (str) : Bend type.
    st (bool, optional) : If True, prints notes and warnings.
    """
    if btype.lower() == 'sbend':
        cname = 'S-Bend ' + cname
    elif btype.lower() == 'sin':
        cname = 'Sine ' + cname
    elif btype.lower() == 'cobra':
        cname = 'Cobra ' + cname
    elif btype.lower() == 'bsb':
        cname = 'Bend-Straight-Bend ' + cname
    return cname

def mbo_get_xsb_xgl_(xs, w, r, pos2, st=False):
    """
    Defines x_sb, the length in x of the sbend waveguide, and the extra strt gap
    length x_gl added later.
    """
    ic = Interconnect(xs=xs)
    x2, y2, a2 = pos2
    x_sb = 0
    if abs(y2) >= abs(r*2):
        x_sb = 2*r
    else:
        # Calls function in Nazca to get the precise sbend length
        sb_solve = ic._sbend_solve(radius=r, width=w, pin=(0,0,0), offset=y2)
        x_sb = sb_solve[1]['xya'][0]
    
    x_gl = abs(x2) - x_sb
    return x_sb, x_gl

def mbo_get_pos2_(p1, p2, st=False):
    """
    Get pos2 of a bent object, given input point p1 and output point p2.
    """
    pos2 = add_pos(p2, p1, minus=True)
    pos2 = rotate_pos(pos2, p1[2])
    return pos2

def mbo_axy_arc(xy, r, ang, da, accr, st=False):
    """
    Calculates the points along an arc to accuracy accr, adds them to xy, a
    list of points on the curve.
    """
    ang_new = ang
    if da != 0:
        da_step = da*accr
        i_max = int(da/da_step)
        for i in range(i_max):
            ang_new += da_step # Add small angle
            dx, dy = polar_to_cart(abs(r*da_step/180*np.pi), ang_new/180*np.pi)
            # Add point
            xy.append(( xy[-1][0]+dx, xy[-1][1]+dy ))
            
        ang += da # The angle changes
    return xy, ang

def mbo_get_all_xy_(xs, w, r, pos2, btype, ang=0, st=False):
    """
    Calculates all points in a bent object. Code inspired by Nazca's
    'interconnect' and 'mask_elements' modules.
    """
    ic = Interconnect(xs=xs)
    x2, y2, a2 = pos2
    if y2 == 0 and a2 == 0:
        return [(pos2[0], pos2[1])], [0] # If strt line
    
    l = 0
    xy = [(0,0)]
    x_sb, x_gl = 0, 0
    accr = 0.1/r # Define accr
    if accr > 0.01:
        accr = 0.01
    
    if btype.lower() == 'sbend':
        parse, result = ic._sbend_solve(radius=r, width=w, pin=(0,0,0), offset=y2)
        geo = result['geo']

        el_list=[]
        new_geo=[]
        for i, tube in enumerate(geo):
            if tube[0] == 's':
                if abs(tube[1][0])<1.e-5:
                    continue
                dr, temp_dump = tube[1]
                dx, dy = polar_to_cart(dr, ang/180*np.pi)
                xy.append(( xy[-1][0]+dx, xy[-1][1]+dy ))
            elif tube[0] == 'b':
                if abs(tube[1][0])<1.e-5:
                    continue
                da, r, temp_dump = tube[1]
                xy, ang = mbo_axy_arc(xy, r, ang, da, accr, st)

        for i, tube in enumerate(geo):
            if tube[0] == 's':
                l += tube[1][0]
            elif tube[0] == 'b':
                angle = tube[1][0]
                radius = tube[1][1]
                l += abs(radius * angle/180*np.pi)
        x_sb, x_gl = mbo_get_xsb_xgl_(xs, w, r, pos2, st)
        if x_sb != x2: # The potential extra gap length
            l += x_gl
            xy_new = [(0,0)]
            for x, y in xy:
                x = x + x_gl/2
                xy_new.append((x,y))
            xy = xy_new
            xy.append((xy[-1][0]+x_gl/2, xy[-1][1]))

    elif btype.lower() == 'sin':
        grow = ((0,0),(0,0),0,0)
        xya = (x2, y2, 0)
        acc = 0.001
        xy = nd.mask_elements.curve2polyline(nd.generic_bend.sinebend_point,
                                             xya, acc, (x2, y2))
        xy = nd.util.polyline2edge(xy, 0, grow=grow, line=True)
    elif btype.lower() == 'cobra':
        grow = ((0,0),(0,0),0,0)
        xya = pos2
        offset1 = nd.mask_elements.__get_offset(xs, w, 0)
        offset2 = nd.mask_elements.__get_offset(xs, w, 0)
        xya = (xya[0], xya[1]-offset1, xya[2])
        xya = (xya[0] - offset2*np.sin(np.radians(xya[2])), xya[1] +
               offset2*np.cos(np.radians(xya[2])), xya[2])
        acc = 0.001
        
        A, B, L, Rmin = nd.generic_bend.gb_coefficients(xya, radius1=0, radius2=0)
        xy = nd.mask_elements.curve2polyline(nd.generic_bend.gb_point, xya, acc,
                                             (A, B, L))
        xy = nd.util.polyline2edge(xy, 0, grow=grow, line=True)
    elif btype.lower() == 'bsb':
        parse, curves = ic._bend_strt_bend_p2p_solve((0,0,0), pos2, radius1=r,
                                                 radius2=r, width=w)
        radius1 = r
        radius2 = r
        width = w
        variations = [curves]

        for curve in variations:
            param   = curve['solution']
            Ltot  = param['Ltot']
            L     = param['L']
            b     = param['angle1']
            e     = param['angle2']
            
            da = b/np.pi*180 # Bend
            xy, ang = mbo_axy_arc(xy, r, ang, da, accr, st)
            
            dr = L # Straight
            dx, dy = polar_to_cart(dr, ang/180*np.pi)
            xy.append(( xy[-1][0]+dx, xy[-1][1]+dy ))
            
            da = e/np.pi*180 # Bend
            xy, ang = mbo_axy_arc(xy, r, ang, da, accr, st)
    return xy

def mbo_get_l_(xs, w, r, pos2, btype, st=False):
    """
    Calculates the total length through a bent object. Code inspired by Nazca's
    'interconnect' and 'mask_elements' modules.
    """
    ic = Interconnect(xs=xs)
    x2, y2, a2 = pos2
    if y2 == 0 and a2 == 0:
        return x2 # Strt line
    
    l = 0
    xy = []
    x_sb, x_gl = 0, 0
    if btype.lower() == 'sbend':
        parse, result = ic._sbend_solve(radius=r, width=w, pin=(0,0,0), offset=y2)
        geo = result['geo']
        
        for i, tube in enumerate(geo):
            if tube[0] == 's':
                l += tube[1][0]
            elif tube[0] == 'b':
                angle = tube[1][0]
                radius = tube[1][1]
                l += abs(radius * angle/180*np.pi)
        x_sb, x_gl = mbo_get_xsb_xgl_(xs, w, r, pos2, st)
        if x_sb != x2: # The potential extra gap length
            l += x_gl
    elif btype.lower() == 'sin':
        xy, temp_dump = mbo_get_all_xy_(xs, w, r, pos2, btype, st)
    elif btype.lower() == 'cobra':
        xy, temp_dump = mbo_get_all_xy_(xs, w, r, pos2, btype, st)
    elif btype.lower() == 'bsb':
        bsb_solve = ic._bend_strt_bend_p2p_solve((0,0,0), pos2, radius1=r,
                                                 radius2=r, width=w)
        l = bsb_solve[1]['solution']['Ltot']
    if xy != []:
        for i in range(len(xy)-1): # Finds length between all xy points
            p1 = xy[i]
            p2 = xy[i+1]
            l += pythagoras(p1, p2)
    return l

def make_b_object_(xs, w, r, pos2, btype, is_con=False, st=False):
    """
    Creates and returns a bent waveguide (or conductor) cell. The cell usually
    constists of two bends: a first bend away from the horizontal, then
    a second bend back onto the horizontal.

       Parameters
    xs (str)  : Cross-section, specifying the layers.
    w (float) : Width (usually in y).
    r (float) : Bending radius of the two bends. Not used for 'cobra' and 'sin'.
    pos2 (tuple of floats) : Position (x,y,a) of the output, when the input
        is the origin (0,0,0).
    btype (str) : Bend type.
    is_con (bool, optional) : If True, changes name to 'Conductor'.
    st (bool, optional) : If True, prints notes and warnings.
    """
    ic = Interconnect(xs=xs)
    cname = 'Waveguide'
    if is_con:
        cname = 'Conductor'
    cname = mbo_update_name_(cname, btype, st)

    x2, y2, a2 = pos2
    if y2 == 0 and a2 == 0:
        # Ends function if straight line
        return make_strt_object_(xs, w, x2, is_con=is_con, st=st)
    
    x_sb, x_gl = 0, 0
    add_pins = True
    with nd.Cell(name=cname+get_n_str()) as BW_cell:
        
        # More bend types can always be added
        
        if btype.lower() == 'sbend':
            x_sb, x_gl = mbo_get_xsb_xgl_(xs, w, r, pos2, st)
            if x_sb == x2: # If nothing needs to be added in x
                elm_bend = ic.sbend(radius=r, width=w, pin=(0,0,0), offset=y2).put(0)
            else: # Adds straight lines before and after
                p_ang = 0
                if x2 < 0:
                    y2   = -y2 # If negative gap, rotates pin angles by 180 deg
                    x_gl = -x_gl
                    p_ang = 180
                
                elm_si = ic.strt(length=(x_gl/2), width=w).put(0)
                pos_b = ( x_gl/2, 0, p_ang )                # For elm_bend
                elm_bend = ic.sbend(radius=r, width=w, pin=(0,0,0), offset=y2).put(pos_b)
                pos_o = elm_bend.pin['b0'].move(0,0,p_ang)  # For elm_so
                elm_so = ic.strt(length=(x_gl/2), width=w).put(pos_o)
                if st:
                    print('NOTE: Straight lines have been added to make sbend '+
                          'curve(s) connect.')
                # Adds sbend pins
                add_pins = False
                nd.Pin("a0", pin=elm_si.pin['a0']).put()
                nd.Pin("b0", pin=elm_so.pin['b0']).put()
            
        elif btype.lower() == 'sin':
            elm_bend = ic.sinebend(width=w, pin=(0,0,0), distance=x2, offset=y2).put(0)
        elif btype.lower() == 'cobra':
            elm_bend = nd.mask_elements.cobra(xya=pos2, width1=w, width2=w, xs=xs).put(0)
        elif btype.lower() == 'bsb':
            elm_bend = ic.bend_strt_bend_p2p((0,0,0), pos2, radius=r, width=w).put(0)
        else:
            raise ValueError("Invalid value of btype. Should be sbend, sin, "+
                             "cobra or bsb.")
        
        if add_pins: # Adds pins if not added already
            nd.Pin("a0", pin=elm_bend.pin['a0']).put()
            nd.Pin("b0", pin=elm_bend.pin['b0']).put()
    return BW_cell

def make_bo_pair_(xs, w1, w2, r1, r2, pos2, pos3, pos4, btype1, btype2,
        is_con=False, st=False):
    """
    Creates and returns a cell with a pair of bent waveguides (or conductors).
    The pair consists of object 1 (usually the upper one) and 2 (the lower one).

       Parameters
    xs (str)   : Cross-section, specifying the layers.
    w1 (float) : Width of object 1 (usually in y).
    w2 (float) : Width of object 2 (usually in y).
    r1 (float) : Bending radius of object 1.
    r2 (float) : Bending radius of object 2.
    pos2 (tuple of floats) : Position (x,y,a) of object 1's output, when
        object 1's input is the origin (0,0,0).
    pos3 (tuple of floats) : Position (x,y,a) of object 2's input, when
        object 1's input is (0,0,0).
    pos4 (tuple of floats) : Position (x,y,a) of object 2's output, when
        object 1's input is (0,0,0).
    btype1 (str) : Bend type of object 1.
    btype2 (str) : Bend type of object 2.
    is_con (bool, optional) : If True, changes name to 'Conductor'.
    st (bool, optional) : If True, prints notes and warnings.
    """
    cname = 'Waveguide'
    if is_con:
        cname = 'Conductor'
    cname = mbo_update_name_(cname, btype1, st) + ' Pair' # Only btype1 in name
    with nd.Cell(name=cname+get_n_str()) as BWP_cell:
        elm_bend_1 = make_b_object_(xs, w1, r1, pos2, btype1, is_con, st).put(0)
        elm_bend_2 = make_b_object_(xs, w2, r2, pos4, btype2, is_con, st).put(pos3)
            
        # Adds pins
        nd.Pin("a0", pin=elm_bend_1.pin['a0']).put()
        nd.Pin("a1", pin=elm_bend_2.pin['a0']).put()
        nd.Pin("b0", pin=elm_bend_1.pin['b0']).put()
        nd.Pin("b1", pin=elm_bend_2.pin['b0']).put()
    return BWP_cell

def make_sqrc_(xs, w, x_os, y_os, r, pos, st=False):
    """
    Creates and returns a squircle cell (rectangle with circular arcs as corners).
    If rectangle sides are zero, we have a circle.
    
       Parameters
    xs (str)  : Cross-section, specifying the layers.
    w (float) : Width.
    x_os (float) : The rectangle side lengths in x.
    y_os (float) : The rectangle side lengths in y.
        If x_os and y_os are zero, we have a circle. If only one is zero, we
        have a racetrack.
    r (float) : Bending radius of all corners. Also defines size of corners
        and radius of circle.
    pos (tuple of floats) : Position (x,y,a) of the squircle, when
        the center of the bottom left circular arc is the origin (0,0,0).
    st (bool, optional)   : If True, prints notes and warnings.
    """
    ic = Interconnect(xs=xs)
    cname = 'Squircle'
    if x_os == 0 and y_os == 0:
        cname = 'Circle'
    with nd.Cell(name=cname+get_n_str()) as S_cell:
        x0, y0, a0 = pos
        r = r-w/2 # To behave like a Nazca circle (i.e., radius doesn't include the width)
        pos_arc1 = ( x0-r,        y0,          a0-90  )
        pos_arc2 = ( x0 + x_os,   y0-r,        a0     )
        pos_arc3 = ( x0,          y0+r + y_os, a0+180 )
        pos_arc4 = ( x0+r + x_os, y0   + y_os, a0+90  )
        
        elm_arc1 = nd.mask_elements.bend(radius=r, width=w, angle=90, xs=xs).put(pos_arc1)
        elm_arc2 = nd.mask_elements.bend(radius=r, width=w, angle=90, xs=xs).put(pos_arc2)
        elm_arc3 = nd.mask_elements.bend(radius=r, width=w, angle=90, xs=xs).put(pos_arc3)
        elm_arc4 = nd.mask_elements.bend(radius=r, width=w, angle=90, xs=xs).put(pos_arc4)
        if x_os != 0: # Horizontal lines
            elm_s1 = ic.strt_p2p(pin1=elm_arc1.pin['b0'], pin2=elm_arc2.pin['a0'],
                                 width=w).put()
            elm_s3 = ic.strt_p2p(pin1=elm_arc3.pin['a0'], pin2=elm_arc4.pin['b0'],
                                 width=w).put()
        if y_os != 0: # Vertical lines
            elm_s2 = ic.strt_p2p(pin1=elm_arc1.pin['a0'], pin2=elm_arc3.pin['b0'],
                                 width=w).put()
            elm_s4 = ic.strt_p2p(pin1=elm_arc2.pin['b0'], pin2=elm_arc4.pin['a0'],
                                 width=w).put()
    return S_cell

def make_bs_(xs, w0, w1, w2, l0, l1, l2, r1, r2, pos2, pos3, btype1, btype2,
             gap=None, flip=False, st=False):
    """
    Creates and returns a beam splitter cell. Consists of an input wg and two
    output wgs, waveguide 1 and 2. If flipped, instead has two inputs and one
    output.

       Parameters
    xs (str) : Cross-section, specifying the layers.
    w0 (float) : Width of the strt input waveguide (usually the left).
    w1 (float) : Width of output waveguide 1 (usually the upper right one).
    w2 (float) : Width of output waveguide 2 (usually the lower right one).
    l0 (float) : Length of input waveguide (usually in x).
    l1 (float) : The extra length of the strt part of output waveguide 1.
    l2 (float) : The extra length of the strt part of output waveguide 2.
    r1 (float) : Bending radius of output waveguide 1.
    r2 (float) : Bending radius of output waveguide 2.
    pos2 (tuple of floats) : Position (x,y,a) of the output of the bent part of
        output waveguide 1.
    pos3 (tuple of floats) : Position (x,y,a) of the output of the bent part of
        output waveguide 2.
    btype1 (str) : Bend type of output waveguide 1.
    btype2 (str) : Bend type of output waveguide 2.
    gap (float, optional) : Distance between the coupling regions of the
        output waveguides. NoneType means no gap. Positive means output wg 1
        moves up, negative means output wg 2 moves down.
    flip (bool, optional) : If True, reflects beam splitter across y axis.
    st (bool, optional)   : If True, prints notes and warnings.
    """
    ic = Interconnect(xs=xs)
    with nd.Cell(name='Beam Splitter'+get_n_str()) as B_cell:
        if gap != None:
            if gap > 0:
                x2, y2, a2 = pos2
                pos2 = (x2, y2+gap+(w1+w2)/2, a2) # Accounts for w
            if gap < 0:
                x3, y3, a3 = pos3
                pos3 = (x3, y3+gap-(w1+w2)/2, a3) # Accounts for w

        if flip:
            x2, y2, a2 = pos2
            pos2 = (x2, -y2, a2) # Switches output wg 1 and 2
            x3, y3, a3 = pos3
            pos3 = (x3, -y3, a3)
        
        elm_s_i  = ic.strt(length=l0, width=w0).put(0)
        elm_s_o1 = ic.strt(length=l1, width=w1).put(pos2)
        elm_s_o2 = ic.strt(length=l2, width=w2).put(pos3)
        
        xi, yi, ai = elm_s_i.pin['b0'].xya()
        xo1, yo1, ao1 = pos2
        xo2, yo2, ao2 = pos3
        pos_elm_so1 = ( xo1-xi, yo1, 0 )
        pos_elm_so2 = ( xo2-xi, yo2, 0 )
        
        pos_elm_si1 = (xi, yi, ai)
        pos_elm_si2 = (xi, yi, ai)
        if gap != None:
            k = 1
            if flip:
                k = -1
            if gap > 0:
                x, y, a = pos_elm_si1
                pos_elm_si1 = (x, y+k*(gap+(w1+w2)/2), a) # Accounts for w
                x, y, a = pos_elm_so1
                pos_elm_so1 = (x, y-k*(gap+(w1+w2)/2), a) # Accounts for w
            if gap < 0:
                x, y, a = pos_elm_si2
                pos_elm_si2 = (x, y+k*(gap-(w1+w2)/2), a) # Accounts for w
                x, y, a = pos_elm_so2
                pos_elm_so2 = (x, y-k*(gap-(w1+w2)/2), a) # Accounts for w
        is_con = False
        elm_bend_1 = make_b_object_(xs, w1, r1, pos_elm_so1, btype1, is_con,
                                st).put(pos_elm_si1)
        elm_bend_2 = make_b_object_(xs, w2, r2, pos_elm_so2, btype2, is_con,
                                st).put(pos_elm_si2)
            
        # Adds pins
        if not flip:
            nd.Pin("a0", pin=elm_s_i.pin['a0']).put()
            nd.Pin("b0", pin=elm_s_o1.pin['b0']).put()
            nd.Pin("b1", pin=elm_s_o2.pin['b0']).put()
        else:
            nd.Pin("a0", pin=elm_s_o1.pin['b0']).put()
            nd.Pin("a1", pin=elm_s_o2.pin['b0']).put()
            nd.Pin("b0", pin=elm_s_i.pin['a0']).put()
    return B_cell

def make_sym_bs_(xs, w, l_in, l_out, r, pos2, btype, gap=None, flip=False, st=False):
    """
    Creates and returns a beam splitter cell, where the output waveguides are
    identical but go in opposite directions (a symmetric beam splitter). All
    widths are the same. If flipped, instead has two inputs and one output.

       Parameters
    xs (str)  : Cross-section, specifying the layers.
    w (float) : Width of all waveguides.
    l_in  (float) : Length of input waveguide (usually in x).
    l_out (float) : The extra length of the strt parts of the output waveguides.
    r (float) : Bending radii of the output waveguides.
    pos2 (tuple of floats) : Position (x,y,a) of the output of the bent part of
        output waveguide 1.
    btype (str) : Bend type of the output waveguides.
    gap (float, optional) : Distance between the coupling regions of the
        output waveguides. NoneType means no gap. Positive means output wg 1
        moves up, negative means output wg 2 moves down.
    flip (bool, optional) : If True, reflects beam splitter across y axis.
    st (bool, optional)   : If True, prints notes and warnings.
    """
    with nd.Cell(name='Symmetric Beam Splitter'+get_n_str()) as SB_cell:
        x2, y2, a2 = pos2
        pos3 = ( x2, -y2, -a2 )
        # Calls the general BS function
        sym_bs = make_bs_(xs, w, w, w, l_in, l_out, l_out, r, r, pos2, pos3,
                          btype, btype, gap, flip, st).put(0)
                
        # Adds pins
        nd.Pin("a0", pin=sym_bs.pin['a0']).put()
        nd.Pin("b0", pin=sym_bs.pin['b0']).put()
        if not flip:
            nd.Pin("b1", pin=sym_bs.pin['b1']).put() # Has a second output
        else:
            nd.Pin("a1", pin=sym_bs.pin['a1']).put() # Has a second input
    return SB_cell

def make_bump_wg_(xs, w, l, r_hi, r_lo, ang, pos_arc, pos_strt2, btype1, btype2,
                  flip=False, st=False):
    """
    Creates and returns a cell of a waveguide with a symmetric bump. From input
    to output, the bump consists of a bent btype wg, a circular arc, an optional
    flat plateau (usually in x), a circular arc and a bent btype wg.

       Parameters
    xs (str)  : Cross-section, specifying the layers.
    w (float) : Width.
    l (float) : Length of flat plateau, the middle of the bump. Usually zero.
    r_hi (float) : Bending radii of the circular arcs.
    r_lo (float) : Bending radii of the bent i/o waveguides.
    ang (float)  : Maximum angle of each circular arc. The arcs go through
        [0,ang]. By default, the zero angle is straight up, since the arcs
        start from there and go outwards and down.
    pos_arc (tuple of floats) : The position of the first circular arc (usually
        the left one).
    pos_strt2 (tuple of floats) : Position (x,y,a) of the output of the output
        waveguide (usually to the right).
    btype1 (str) : Bend type of the input waveguide.    
    btype2 (str) : Bend type of the output waveguide.
    flip (bool, optional) : If True, reflects the wg across the x axis.
    st (bool, optional)   : If True, prints notes and warnings.
    """
    ic = Interconnect(xs=xs)
    with nd.Cell(name='Bump Waveguide'+get_n_str()) as BW_cell:
        xa, ya, aa = pos_arc
        pos_arc_flip = ( xa,   ya, aa+180 )
        pos_arc      = ( xa+l, ya, aa     ) # Adds the potential plateau length
        
        # Makes arcs and plateau
        elm_arc1 = nd.mask_elements.bend(radius=r_hi, width=w, angle= ang, xs=xs).put(pos_arc_flip)
        elm_arc2 = nd.mask_elements.bend(radius=r_hi, width=w, angle=-ang, xs=xs).put(pos_arc)
        if l != 0: # If plateau
            elm_s_mid = ic.strt_p2p(elm_arc1.pin['a0'], elm_arc2.pin['a0'], width=w).put()
        
        # Makes btype bends
        xa1, ya1, aa1 = elm_arc1.pin['b0'].xya()
        xa2, ya2, aa2 = elm_arc2.pin['b0'].xya()
        x2, y2, a2    = pos_strt2
        
        pos_b2_i = ( x2-xa2, y2-ya2, aa2-a2  )
        pos_b2_o = ( x2,     y2,     a2+180  )
        pos_b1_o = ( xa1,    ya1,    aa1+180 )
        
        is_con = False
        elm_bend_1 = make_b_object_(xs, w, r_lo, pos_b1_o, btype1, is_con, st).put(0)
        elm_bend_2 = make_b_object_(xs, w, r_lo, pos_b2_i, btype2, is_con, st).put(pos_b2_o)
        # Adds pins
        if not flip:
            nd.Pin("a0", pin=elm_bend_1.pin['a0']).put()
            nd.Pin("b0", pin=elm_bend_2.pin['a0']).put()
        else:
            nd.Pin("a0", pin=elm_bend_2.pin['a0']).put() # Flips if negative bump
            nd.Pin("b0", pin=elm_bend_1.pin['a0']).put()            
    return BW_cell

def mrtr_update_name_(cname, rtr_type, x_os, y_os, is_ring=False, st=False):
    """
    Updates name of racetrack resonator (straight, symmetric or pulley).

       Parameters
    cname (str)    : Generic name of the resonator.
    rtr_type (str) : Type of racetrack (straight, symmetric, pulley).
    x_os (float) : The rectangle side lengths in x.
    y_os (float) : The rectangle side lengths in y.
        If x_os and y_os are zero, we have a circle. If only one is zero, we
        have a racetrack.
    is_ring (bool, optional) : If True, changes name to 'Ring Resonator'.
    st (bool, optional) : If True, prints notes and warnings.
    """
    if   x_os == 0 and y_os != 0:
        cname = 'Vertical ' + cname   # Updates name if vertical racetrack
    elif x_os != 0 and y_os == 0:
        cname = 'Horizontal ' + cname # Updates name if horizontal racetrack
    
    if x_os == 0 and y_os == 0 and not is_ring and st:
        print('NOTE: Racetrack resonator cell created with x_os=0 and y_os=0. '+
              'This is identical to a ring resonator. If intentional, '+
              'set is_ring=True.')
    if (x_os != 0 or y_os != 0) and is_ring and st:
        print('NOTE: Ring resonator cell with non-zero x_os and/or y_os detected. '+
              'x_os and y_os do nothing in a ring resonator.')
    
    if is_ring: # Updates cname if ring
        cname = 'Ring Resonator'
    
    if   rtr_type == 'srt': # Name depending on type of coupled region
        cname = 'Straight ' + cname
    elif rtr_type == 'sym':
        cname = 'Symmetric ' + cname
    elif rtr_type == 'pul':
        cname = 'Pulley ' + cname
    else:
        raise ValueError("Invalid value of rtr_type. Should be srt, sym or pul.")
    return cname
 
def make_rtr_(xs, rtr_type, w_co, w_tk, x_os_tk, y_os_tk, r_tk, l_co, r_co, ang,
              pos_arc, pos_strt2, pos_tk, btype1, btype2, flip=False,
              is_ring=False, st=False):
    """
    Creates and returns a racetrack resonator cell (straight, symmetric or
    pulley).

       Parameters
    xs (str) : Cross-section, specifying the layers.
    rtr_type (str) : Type of track/ring (straight, symmetric, pulley).
    w_co (float) : Width of the coupled region.
    w_tk (float) : Width of the track/ring.
    x_os_tk (float) : The rectangle side lengths of the track/ring in x.
    y_os_tk (float) : The rectangle side lengths of the track/ring in y.
    r_tk (float) : Bending radius of the track/ring. Also defines size of
        corners or radius of circle.
    l_co (float) : Length of the wg coupled region (the part closest to the
        track/ring).
    r_co (float) : Bending radius of the bent i/o waveguides of the coupled
        region. Used for make_bump_wg_(), when rtr_type != 'srt'.
    ang (float)  : Maximum angle of each circular arc of the coupled region.
        The arcs go through [0,ang]. By default, the zero angle is straight up,
        since the arcs start from there and go outwards and down. Used for
        make_bump_wg_(), when rtr_type != 'srt'.
    pos_arc (tuple of floats) : The position of the left circular arc of
        the coupled region. Used for make_bump_wg_(), when rtr_type != 'srt'.
    pos_strt2 (tuple of floats) : Position (x,y,a) of the output of the output
        waveguide (usually the right) of the coupled region. Used for
        make_bump_wg_(), when rtr_type != 'srt'.
    pos_tk (tuple of floats) : Position (x,y,a) of the track/ring, when the
        center of the bottom left circular arc is the origin (0,0,0).
    btype1 (str) : Bend type of the input waveguide of the coupled region. Used
        for make_bump_wg_(), when rtr_type != 'srt'.
    btype2 (str) : Bend type of the output waveguide of the coupled region. Used
        for make_bump_wg_(), when rtr_type != 'srt'.
    flip (bool, optional)    : If True, reflects the resonator across the x axis.
    is_ring (bool, optional) : If True, changes name to 'Ring Resonator'.
    st (bool, optional)      : If True, prints notes and warnings.
    """
    ic = Interconnect(xs=xs)
    cname = 'Racetrack Resonator'
    cname = mrtr_update_name_(cname, rtr_type, x_os_tk, y_os_tk, is_ring, st)
    if is_ring:
        x_os_tk = y_os_tk = 0

    with nd.Cell(name=cname+get_n_str()) as RR_cell:
        r_hi = 0
        # Makes the bump/strt waveguide (the coupling region)
        elm_co_wg = None
        if rtr_type == 'sym':
            r_hi = r_tk
        elif rtr_type == 'pul':
            ang = -ang
            l_co += x_os_tk # So that l=0 gives the expected result for pulley resonators
            r_hi = abs(pos_arc[1] - pos_tk[1]) # To make the track and arc concentric
        
        if rtr_type == 'srt':
            elm_co_wg = ic.strt(length=l_co, width=w_co).put(0)
        else:
            elm_co_wg = make_bump_wg_(xs, w_co, l_co, r_hi, r_co, ang,
                    pos_arc, pos_strt2, btype1, btype2, flip=False, st=st).put(0)
                    # No flip, otherwise elm_co_wg is flipped twice
        
        # Makes the track/ring
        elm_tk = make_sqrc_(xs, w_tk, x_os_tk, y_os_tk, r_tk, pos_tk, st).put(0)
        
        # Adds pins
        if not flip:
            nd.Pin("a0", pin=elm_co_wg.pin['a0']).put()
            nd.Pin("b0", pin=elm_co_wg.pin['b0']).put()
        else:
            nd.Pin("a0", pin=elm_co_wg.pin['b0']).put()
            nd.Pin("b0", pin=elm_co_wg.pin['a0']).put()
    return RR_cell

###############################################################################
####################### Spiral Delay Line Functions ###########################
###############################################################################

def mdl_get_l_(xs, w, dr, x_os, y_os, r_lo, n_spirals, btype, r_b, st=False):
    """
    Calculates the total length of a spiral delay line (rectangle or circular).
    """
    l = 0
    l += 2*n_spirals*(x_os + y_os) # Strt lengths (rectangle parts)

    r = r_lo - w/2
    for i in range(n_spirals): # Circular arc lengths (circular parts)
        l += 2*np.pi*r
        r += dr
    
    ## Gets pos2_rot, used below in _mbo_get_l()
    y_os2 = y_os + dr
    pos_arc4 = ( r_lo+x_os/2, y_os2/2, +90 )
    
    pos_m2 = (pos_arc4[0]-r_lo, pos_arc4[1]-dr+r_lo, pos_arc4[2]+90)
    pos_m1 = (-pos_m2[0], -pos_m2[1], 0)

    p2 = ( pos_m2[0]-pos_m1[0], pos_m2[1]-pos_m1[1], 0 )
    p2_rot = rotate_pos(p2, pos_m1[2])
    pos2_rot = (p2_rot[0], p2_rot[1], pos_m2[2]-pos_m1[2]+180)
    
    l_mid = mbo_get_l_(xs, w, r_b, pos2_rot, btype, st)
    l += l_mid
    return l    

def mdl_update_xflips_(x, dx, flip, l, l_t1, l_t2, st=False):
    """
    Finds new values of x and dx for spiral delay line parameter optimization.
    x can be x_os, y_os or r_lo.
    """
    if abs(l_t1 - l) < abs(l_t2 - l):
        flip = 0
    else:
        dx = -dx # Flips (jumped over the optimal x value)
        flip += 1
    x += dx
    if flip > 0:
        dx = dx/2 # More precise step dx, if flipped
    return x, dx, flip

def mdl_optimize_size_(xs, w, dr, l, btype, x_os=None, y_os=None, r_lo=None,
                       r_b=None, n_spirals=2, dx=10, dy=10, drl=10, max_steps=100,
                       l_error=0.01, xy_same=False, st=False):
    """
    For spiral delay line: finds locally optimal values of x_os, y_os and r_lo
    so the total length is within 'l_error' of l. If one or more of x_os, y_os
    or r_lo are given, no optimization occurs for those parameters. Reaches the
    local critical point by using a simple form of gradient descent, in steps
    of dx, dy and drl.
    """
    rb_change = False
    if r_b == None:
        rb_change = True # Can change r_b
    
    # Cases when no optimization
    if x_os != None:
        dx = 0
    if y_os != None:
        dy = 0
    if r_lo != None:
        drl = 0

    if r_lo == None:
        r_lo = l/20
    
    if rb_change:
        r_b = r_lo/2 - w
        if r_b < w:
            r_b = w
    
    if   x_os == None and y_os == None: # Initial guesses
        x_os = r_lo
        y_os = r_lo
    elif x_os == None and y_os != None:
        x_os = y_os
    elif x_os != None and y_os == None:
        y_os = x_os
    
    if xy_same:
        if dy != 0:
            dx = dy
            x_os = y_os
        else:
            dy = dx
            y_os = x_os
        if st:
            print('NOTE: May override initial x_os or y_os value due to '+
                  'xy_same=True.')

    flip_x  = 0
    flip_y  = 0
    flip_rl = 0
    i = 0
    optimal_dl = False
    while not optimal_dl:
        if x_os < 0:
            x_os = 0
        if y_os < 0:
            y_os = 0
        if r_lo < 0:
            r_lo = 0
        
        i += 1
        if xy_same: # Treats x_os and y_os identically
            if dx != 0:
                y_os = x_os
                l_t1 = mdl_get_l_(xs, w, dr, x_os+dx, y_os+dx, r_lo, n_spirals,
                                  btype, r_b, st)
                l_t2 = mdl_get_l_(xs, w, dr, x_os-dx, y_os-dx, r_lo, n_spirals,
                                  btype, r_b, st)
                x_os, dx, flip_x = mdl_update_xflips_(x_os, dx,
                                flip_x, l, l_t1, l_t2, st)
                y_os = x_os
                dy = dx
        else: # Optimize x_os and y_os independently
            if dx != 0:
                l_t1 = mdl_get_l_(xs, w, dr, x_os+dx, y_os, r_lo, n_spirals,
                                  btype, r_b, st)
                l_t2 = mdl_get_l_(xs, w, dr, x_os-dx, y_os, r_lo, n_spirals,
                                  btype, r_b, st)
                x_os, dx, flip_x = mdl_update_xflips_(x_os, dx,
                                    flip_x, l, l_t1, l_t2, st)
            if dy != 0:
                l_t1 = mdl_get_l_(xs, w, dr, x_os, y_os+dy, r_lo, n_spirals,
                                  btype, r_b, st)
                l_t2 = mdl_get_l_(xs, w, dr, x_os, y_os-dy, r_lo, n_spirals,
                                  btype, r_b, st)
                y_os, dy, flip_y, = mdl_update_xflips_(y_os, dy,
                                    flip_y, l, l_t1, l_t2, st)

        if drl != 0:
            l_t1 = mdl_get_l_(xs, w, dr, x_os, y_os, r_lo+drl, n_spirals,
                              btype, r_b, st)
            l_t2 = mdl_get_l_(xs, w, dr, x_os, y_os, r_lo-drl, n_spirals,
                              btype, r_b, st)
            r_lo, drl, flip_rl = mdl_update_xflips_(r_lo,
                        drl, flip_rl, l, l_t1, l_t2, st)
            if rb_change:
                r_b = r_lo/2 - w
                if r_b < w:
                    r_b = w
        
        if x_os < 0 or y_os < 0 or r_lo < 0:
            raise ValueError('No solution found for delay line. '+
                             'gap, n_spirals and r_b might be too large.')
        if i >= max_steps:
            if st:
                print('NOTE: Optimization limit max_steps reached for '+
                      'delay line. Increase max_steps for better precision.')
            optimal_dl = True
        
        # Gets the length at each step
        l_r = mdl_get_l_(xs, w, dr, x_os, y_os, r_lo, n_spirals, btype, r_b, st)
        if abs(l_r - l) < l_error:
            optimal_dl = True # Reached goal
    return x_os, y_os, r_lo, l_r, r_b

def make_delay_line_(xs, w, gap, btype, l=None, x_os=None, y_os=None, r_lo=None,
                     r_b=None, n_spirals=None, dx=10, dy=10, drl=10, max_steps=100,
                     l_error=0.01, xy_same=False, flip=False, is_circ=False,
                     st=False):
    """
    Creates and returns a spiral delay line cell (rectangle or circular),
    consisting of two spiral objects and a center waveguide. Each spiral is a
    collection of circular arcs and straight waveguides. If l is not given,
    x_os, y_os and r_lo all need to be given. If l is given, x_os, y_os and
    r_lo will be calculated automatically. Therefore, out of x_os, y_os and
    r_lo, define none of them, or at most two.

       Parameters
    xs (str)  : Cross-section, specifying the layers.
    w (float) : Width.
    gap (float) : Gap size between neighboring spirals.
    btype (str) : Bend type of the center waveguide connecting the two spiral objects.
    l (float, optional) : Length of entire delay line.
    x_os (float, optional) : The rectangle side lengths of the spiral in x.
    y_os (float, optional) : The rectangle side lengths of the spiral in y.
    r_lo (float, optional) : The smallest bending radius of both spiral parts.
    r_b (float, optional)  : Bending radius of the center waveguide. If not
        defined, becomes r_lo/2 - w.
    n_spirals (int, optional) : Number of spirals (number of iterations).
    dx (float, optional)  : In optimizing x_os, y_os and r_lo, starting
        correction for x_os.
    dy (float, optional)  : In optimizing x_os, y_os and r_lo, starting
        correction for y_os.
    drl (float, optional) : In optimizing x_os, y_os and r_lo, starting
        correction for r_lo.
    max_steps (int, optional) : Max number of steps in optimizing x_os, y_os and r_lo.
    l_error (float, optional) : The largest allowed error between the total
        desired length l and actual length of the delay line.
    xy_same (bool, optional) : If True, x_os and y_os are equal (in optimization).
    flip (bool, optional) : If True, reflects the delay line across the x axis.
    is_circ (bool, optional) : If True, creates a circular delay line.
    st (bool, optional)   : If True, prints notes and warnings.
    """
    ic = Interconnect(xs=xs)
    dr = gap + w

    if x_os == 0 and y_os == 0 and not is_circ and st:
        print('NOTE: Rectangle Delay Line cell created with x_os=0 and y_os=0. '+
          'This is identical to a Circular Spiral Delay Line. If intentional, '+
          'set is_circ=True.')
    if ((x_os != 0 and x_os !=None) or (y_os != 0 and y_os != None)) and is_circ and st:
        print('NOTE: Circular Spiral Delay Line cell with non-zero x_os and/or '+
          'y_os. x_os and y_os do nothing when is_circ=True.')
    if is_circ:
        x_os = y_os = 0
    
    l_r = None
    if l != None:
        if x_os == None or y_os == None or r_lo == None:
            if n_spirals == None:
                n_spirals = int((l/(4*gap))**0.5)

            x_os, y_os, r_lo, l_r, r_b = mdl_optimize_size_(xs, w, dr, l, btype, x_os,
                y_os, r_lo, r_b, n_spirals, dx, dy, drl, max_steps, l_error, xy_same, st)
        else:
            if is_circ:
                raise ValueError('Circular delay line is overdetermined. When l '+
                                 'is defined, r_lo should not be defined.')
            else:
                raise ValueError('Rectangle delay line is overdetermined. When l '+
                                 'is defined, at least one of x_os, y_os and r_lo '+
                                 'should not be defined.')
    else:
        if x_os == None or y_os == None or r_lo == None:
            if is_circ:
                raise ValueError('Circular delay line is underdetermined. When l '+
                                 'is undefined, r_lo needs to be defined.')
            else:
                raise ValueError('Rectangle delay line is underdetermined. When l '+
                                 'is undefined, x_os, y_os and r_lo need to be '+
                                 'defined.')

    if n_spirals == None:
        n_spirals = 2

    sdname = 'Spiral Delay Line'+get_n_str()
    if is_circ:
        sdname = 'Circular ' + sdname
    else:
        sdname = 'Rectangle ' + sdname
    with nd.Cell(name=sdname) as RSDL_cell:
        if r_b == None:
            r_b = r_lo
        r = r_lo
        y_os2 = y_os + dr # y_os of the next spiral
        pos_m1 = None
        pos_m2 = None
        ang1, ang2, ang3, ang4 = 90,90,90,90
        for i in range(n_spirals): # Code similar to make_sqrc_()
            pos_arc1 = ( -r-x_os/2, -y_os2/2,   -90  )
            pos_arc2 = ( +x_os/2,   -r-y_os2/2, 0    )
            pos_arc3 = ( -x_os/2,   +r+y_os2/2, +180 )
            pos_arc4 = ( +r+x_os/2, +y_os2/2,   +90  )
            
            pos_arc1 = (pos_arc1[0], pos_arc1[1]+dr, pos_arc1[2])
            pos_arc4 = (pos_arc4[0], pos_arc4[1]-dr, pos_arc4[2])
            
            if flip:
                # Flip horizontally. Because of pins, this later becomes a vertical flip
                pos_arc1 = (-pos_arc1[0], pos_arc1[1], pos_arc1[2])
                pos_arc2 = (-pos_arc2[0], pos_arc2[1], pos_arc2[2]+180)
                pos_arc3 = (-pos_arc3[0], pos_arc3[1], pos_arc3[2]+180)
                pos_arc4 = (-pos_arc4[0], pos_arc4[1], pos_arc4[2])
                ang1 = -90
                ang2 = -90
                ang3 = -90
                ang4 = -90
            
            if r != 0:
                elm_arc1 = nd.mask_elements.bend(radius=r, width=w, angle=ang1, xs=xs).put(pos_arc1)
                elm_arc2 = nd.mask_elements.bend(radius=r, width=w, angle=ang2, xs=xs).put(pos_arc2)
                elm_arc3 = nd.mask_elements.bend(radius=r, width=w, angle=ang3, xs=xs).put(pos_arc3)
                elm_arc4 = nd.mask_elements.bend(radius=r, width=w, angle=ang4, xs=xs).put(pos_arc4)

            if x_os != 0: # Horizontal lines
                pt = elm_arc1.pin['b0'].xya()
                pos_s1_new = (pt[0], pt[1]-dr, pt[2]) # Position of the next spiral
                pt = elm_arc4.pin['b0'].xya()
                pos_s3_new = (pt[0], pt[1]+dr, pt[2]+180)
                elm_s1 = ic.strt_p2p(pin1=pos_s1_new, pin2=elm_arc2.pin['a0'],
                                width=w).put()
                elm_s3 = ic.strt_p2p(pin1=elm_arc3.pin['a0'], pin2=pos_s3_new,
                                width=w).put()
            if y_os != 0: # Vertical lines
                elm_s2 = ic.strt_p2p(pin1=elm_arc1.pin['a0'], pin2=elm_arc3.pin['b0'],
                                width=w).put()
                elm_s4 = ic.strt_p2p(pin1=elm_arc2.pin['b0'], pin2=elm_arc4.pin['a0'],
                                width=w).put()
            
            if i==0:
                # Saves positions in the middle of the spiral, used for the center wg
                if r != 0:
                    pos_m1 = elm_arc1.pin['b0'].xya()
                    pos_m2 = elm_arc4.pin['b0'].xya()
                else:
                    pos_arc4 = ( r_lo+x_os/2, y_os2/2, +90 )
                    pos_m2 = (pos_arc4[0]-r_lo, pos_arc4[1]-dr+r_lo, pos_arc4[2]+90)
                    pos_m1 = (-pos_m2[0], -pos_m2[1], 0)
            r += dr
        
        # Makes the center wg
        p2 = ( pos_m2[0]-pos_m1[0], pos_m2[1]-pos_m1[1], 0 )
        p2_rot = rotate_pos(p2, pos_m1[2])
        pos2_rot = (p2_rot[0], p2_rot[1], pos_m2[2]-pos_m1[2]+180)
        is_con = False
        elm_mid = make_b_object_(xs, w, r_b, pos2_rot, btype, is_con, st).put(pos_m1)

        # Adds pins
        pos_a0, pos_b0 = None, None
        if x_os != 0:
            pos_a0 = elm_s1.pin['a0']
            pos_b0 = elm_s3.pin['b0']
        else:
            pos_a0 = elm_arc2.pin['a0']
            pos_b0 = elm_arc3.pin['a0']
            
        nd.Pin("a0").put(pos_a0)
        nd.Pin("b0").put(pos_b0)

    if st:
        # Prints delay line information
        if l_r == None:
            l_r = mdl_get_l_(xs, w, dr, x_os, y_os, r_lo, n_spirals, btype, r_b, st)
        extra_str = " strt x and y offsets x_os = "+str(x_os)+", y_os = "+str(y_os)
        if is_circ:
            extra_str = ""
        print("NOTE: Created '"+sdname+"' cell with exact length l_r = "+str(l_r)+
              ", minimum radius r_lo = "+str(r_lo)+ ", number of spirals n_spirals = "+
              str(n_spirals)+","+extra_str+", bend radius r_b = "+str(r_b)+" and "+
              " maximum radius r_max = "+str(r)+". "+
              "For better precision, decrease l_error = "+str(l_error)+".")
    return RSDL_cell

###############################################################################
###################### Electric Conductor Functions ###########################
###############################################################################

def make_pad_(layers, x, y, st=False):
    """
    Creates and returns a rectangular conductor pad cell. Used at inputs and
    outputs of conductors.

       Parameters
    layers (list of strings) : The layers.
    x (float) : Length in x.
    y (float) : Length in y.
    st (bool, optional) : If True, prints notes and warnings.
    """
    with nd.Cell(name='Conducting Pad'+get_n_str()) as CP_cell:
        pts_pad = nd.geometries.rectangle(length=x, height=y)
        for l_str in layers:
            pos1 = (-x/2, -y/2, 0) # So the center is at the origin
            elm1 = nd.Polygon(points=pts_pad, layer=l_str).put(pos1)
    return CP_cell

def mps_all_cons_(all_xs, all_w, all_r, all_p, all_btype, st=False):
    """
    Creates all of the (bent/strt) parts of a single phase shifter, given
    'all_' lists. Lists go from conductor input to output.
    
       Parameters
    all_xs (list of strings) : All cross-sections, specifying the layers.
    all_w (list of floats) : All widths.
    all_r (list of floats) : All bending radii. Not used for 'cobra' and 'sin'.
    all_p (list of tuples) : All points (x,y,a) through the phase shifer.
        We connect bent btype objects between all neighbouring points.
    all_btype (list of strings) : All bend types.
    st (bool, optional) : If True, prints notes and warnings.
    """
    con_num = len(all_p)
    elms = []
    for i in range(con_num - 1):
        p1 = all_p[i]   # Current position
        p2 = all_p[i+1] # Next position
        xs = all_xs[i]
        w = all_w[i]
        r = all_r[i]
        btype = all_btype[i]
        pos2 = mbo_get_pos2_(p1, p2, st)
        elms.append( make_b_object_(xs, w, r, pos2,
                     btype, is_con=True, st=st).put(p1) )
    return elms

def make_to_ps_(xs, all_w, all_r, p0, all_btype, all_pts, layers_pad=['L1'],
        xy_pad=None, pads=False, flip=False, st=False):
    """
    Creates and returns a thermo-optic phase shifter cell (a conductor). 'all_'
    lists go from conductor input to output.

       Parameters
    xs (str) : Cross-section of conductor, specifying the layers.
    all_w (list of floats) : Widths between all points in all_pts.
    all_r (list of floats) : All bending radii. Not used for 'cobra' and 'sin'.
    p0 (tuple of floats)   : The starting position (x,y,a), the input.
    all_btype (list of strings) : All bend types.
    all_pts (list of floats) : All points (x,y,a) through the phase shifter,
        relative to p0, from input to output.
    layers_pad (str, optional) : Layers of pad.
    xy_pad (tuple, optional) : Lengths in x,y of input pad, then x,y of output
        pad. Has the form (x1, y1, x2, y2).
    pads (bool, optional) : If True, adds pads.
    flip (bool, optional) : If True, reflects the phase shifter across the x axis.
    st (bool, optional)   : If True, prints notes and warnings.
    """
    with nd.Cell(name='Thermo-Optic Phaser Shifter'+get_n_str()) as PS_cell:
        con_num = len(all_pts)
        all_xs = [xs] * con_num # Assumes the entire ps has the same xsection
        all_pts_r = []
        # Makes phase shifter
        for i in range(con_num):
            # Gets the correct coordinates, adding p0 to all_pts
            if flip:
                all_pts[i] = (all_pts[i][0], -all_pts[i][1], -all_pts[i][2])
            all_pts_r.append(  add_pos(p0, all_pts[i])  )
        
        elms = mps_all_cons_(all_xs, all_w, all_r,
            all_pts_r, all_btype, st)
        
        # Makes pads
        pos_i = elms[0].pin['a0']
        pos_o = elms[-1].pin['b0']
        if pads:
            x_pad1, y_pad1, x_pad2, y_pad2 = xy_pad
            pad1 = make_pad_(layers_pad, x_pad1, y_pad1, st).put(pos_i)
            pad2 = make_pad_(layers_pad, x_pad2, y_pad2, st).put(pos_o)        
        
        # Adds pins
        if flip:
            nd.Pin("a0", pin=pos_o).put()
            nd.Pin("b0", pin=pos_i).put()
        else:
            nd.Pin("a0", pin=pos_i).put()
            nd.Pin("b0", pin=pos_o).put()
    return PS_cell

def make_eo_ps_(xs, all_w1, all_w2, all_r1, all_r2, p0, all_btype1, all_btype2,
        all_pts1, all_pts2, layers_pad=['L1'], xy_pad=None, pads=False,
        flip=False, st=False):
    """
    Creates and returns an electro-optic phase shifter cell (a conductor). 'all_'
    lists go from conductor input to output. Has two parts: part 1 (usually
    the top part) and part 2 (usually the bottom part). Part 1 goes from
    input to middle, while part 2 goes from middle to output. Near the middle,
    the two parts are parallel, creating the desired electric field. There, one
    part is on the top and the other is on the bottom.
    
       Parameters
    xs (str) : Cross-section of conductor, specifying the layers.
    all_w1 (list of floats) : Widths of part 1 of the phase shifter.
    all_w2 (list of floats) : Widths of part 2 of the phase shifter.
    all_r1 (list of floats) : All bending radii of part 1. Not used for 'cobra'
        and 'sin'.
    all_r2 (list of floats) : All bending radii of part 2. Not used for 'cobra'
        and 'sin'.
    p0 (tuple of floats) : The starting position (x,y,a), the input.
    all_btype1 (list of strings) : All bend types of part 1.
    all_btype2 (list of strings) : All bend types of part 2.
    all_pts1 (list of floats) : All points (x,y,a) through part 1 of the
        phase shifter, relative to p0, from input to middle.
    all_pts2 (list of floats) : All points (x,y,a) through part 2 of the
        phase shifter, relative to p0, from middle to output.
    layers_pad (str, optional) : Layers of pad.
    xy_pad (tuple, optional) : Lengths in x,y of input pad, then x,y of output
        pad. Has the form (x1, y1, x2, y2).
    pads (bool, optional) : If True, adds pads.
    flip (bool, optional) : If True, reflects the phase shifter across the x axis.
    st (bool, optional)   : If True, prints notes and warnings.
    """
    with nd.Cell(name='Electro-Optic Phaser Shifter'+get_n_str()) as PS_cell:
        con_num1 = len(all_pts1)
        con_num2 = len(all_pts2)
        all_xs = [xs] * (con_num1 + con_num2)
        # Assumes the entire ps has the same xsection. Note: On average, the
        # addition makes all_xs double the needed size.
        
        # Makes phase shifter part 1
        all_pts_r = []
        for i in range(con_num1):
            # Gets the correct coordinates, adding p0 to all_pts1
            if flip:
                all_pts1[i] = (all_pts1[i][0], -all_pts1[i][1], -all_pts1[i][2])
            all_pts_r.append(  add_pos(p0, all_pts1[i])  )
        
        elms1 = mps_all_cons_(all_xs, all_w1, all_r1,
            all_pts_r, all_btype1, st)
        
        # Makes phase shifter part 2
        all_pts_r = []
        for i in range(con_num2):
            # Gets the correct coordinates, adding p0 to all_pts2
            if flip:
                all_pts2[i] = (all_pts2[i][0], -all_pts2[i][1], -all_pts2[i][2])
            all_pts_r.append(  add_pos(p0, all_pts2[i])  )
        
        elms2 = mps_all_cons_(all_xs, all_w2, all_r2,
            all_pts_r, all_btype2, st)
        
        # Makes pads
        pos_i = elms1[0].pin['a0']
        pos_o = elms2[-1].pin['b0']
        if pads:
            x_pad1, y_pad1, x_pad2, y_pad2 = xy_pad
            pad1 = make_pad_(layers_pad, x_pad1, y_pad1, st).put(pos_i)
            pad2 = make_pad_(layers_pad, x_pad2, y_pad2, st).put(pos_o)        
        
        # Adds pins
        if flip:
            nd.Pin("a0", pin=pos_o).put()
            nd.Pin("b0", pin=pos_i).put()
        else:
            nd.Pin("a0", pin=pos_i).put()
            nd.Pin("b0", pin=pos_o).put()
    return PS_cell

###############################################################################
######################## User-Friendly Functions ##############################
###############################################################################

def plus(x_wide, y_wide, x_tall=None, y_tall=None, layers=['L1'], st=False, put=False):
    """
    A marker in the shape of a plus sign.

       Parameters
    x_wide (float) : Length in x of the wide rectangle.
    y_wide (float) : Length in y of the wide rectangle.
    x_tall (float, optional) : Length in x of the tall rectangle.
    y_tall (float, optional) : Length in y of the tall rectangle.
        If x_tall and y_tall are not defined, makes the plus symmetric.
    layers (list of strings) : The layers.
    st (bool, optional)  : If True, prints notes and warnings.
    put (bool, optional) : If True, puts the cell.
    """
    plus = make_plus_(layers, x_wide, y_wide, x_tall, y_tall, st)
    if put:
        plus.put(0)
    return plus

def taper(w_in, w_out, l, layers=['L1'], st=False, put=False):
    """
    A taper of length l, defaults to being horizontal.

       Parameters
    w_in  (float) : Width at the input i (left).
    w_out (float) : Width at the output o (right).
    l (float) : Length from i to o (in x).
    layers (list of strings) : The layers.
    st (bool, optional)  : If True, prints notes and warnings.
    put (bool, optional) : If True, puts the cell.
    """
    taper = make_tp_(layers, w_in, w_out, l, st)
    if put:
        taper.put(0)
    return taper

def strt_wg(w, l, xs='XS1', st=False, put=False):
    """
    A straight waveguide of length l.

       Parameters
    w (float) : Width (in y).
    l (float) : Length (in x).
    xs (str, optional)   : Cross-section, specifying the layers.
    st (bool, optional)  : If True, prints notes and warnings.
    put (bool, optional) : If True, puts the cell.
    """
    s_WG = make_strt_object_(xs, w, l, False, st)
    if put:
        s_WG.put(0)
    return s_WG

def bent_wg(w, x_tot, y_tot, r, btype, xs='XS1', st=False, put=False):
    """
    A bent waveguide, bending away from the horizontal and then back towards
    the horizontal.

       Parameters
    w (float) : Width.
    x_tot (float) : Total displacement in x.
    y_tot (float) : Total displacement in y.
    r (float)   : Bending radius. Not used for 'cobra' and 'sin'.
    btype (str) : Bend type.
    xs (str, optional)   : Cross-section, specifying the layers.
    st (bool, optional)  : If True, prints notes and warnings.
    put (bool, optional) : If True, puts the cell.
    """
    # x_tot = 0 means completely vertical:
    pos2 = (x_tot, y_tot, 0)
    is_con = False
    bent_WG = make_b_object_(xs, w, r, pos2, btype, is_con, st)
    if put:
        bent_WG.put(0)
    return bent_WG

def bwgp_get_py2_3_4_(gap_in, gap_out, y_tot_b_hi, y_tot_b_lo, st=False):
    """
    Calculates p2_y, p3_y, p4_y, the y values for the bent wg pair. p2_y is
    the output of the upper wg. p3_y is the input of the lower wg. p4_y is the
    output of the lower wg. The input of the upper wg is the origin.
    """
    p2_y = None
    p3_y = None
    p4_y = None
    
    none_count = 0
    none_list = []
    if gap_in == None: # Checks that only one or zero of the four are defined
        none_count += 1
        none_list.append('gap_in')
    if gap_out == None: 
        none_count += 1
        none_list.append('gap_out')
    if y_tot_b_hi == None: 
        none_count += 1
        none_list.append('y_tot_b_hi')
    if y_tot_b_lo == None: 
        none_count += 1
        none_list.append('y_tot_b_lo')
    
    if none_count > 1:
        # Raises ValueError
        none_str = ''
        i_max = len(none_list) - 1
        for i in range(i_max + 1): # Formats the error nicely
            none_str += none_list[i]
            if i == i_max-1:
                none_str += ' and '
            elif i < i_max:
                none_str += ', '
        str_error = 'Invalid values for bent waveguide pair. ' + none_str +' are not defined. This is more than one. Only one variable should not be defined.'
        raise ValueError(str_error)
    elif none_count == 0:
        if st:
            print('NOTE: The four variables gap_in, gap_out, y_tot_b_hi and '+
                  'y_tot_b_lo are all defined. Only three are needed. Currently, '+
                  'gap_in is not being used.')
    
    ###### The values for the four valid cases
    if gap_out == None:
        p2_y = y_tot_b_hi
        p3_y = -gap_in
        p4_y = y_tot_b_lo
    elif y_tot_b_lo == None:
        p2_y = y_tot_b_hi
        p3_y = -gap_in
        p4_y = p2_y + (gap_in-gap_out)
    elif y_tot_b_hi == None:
        p3_y = -gap_in
        p4_y = y_tot_b_lo
        p2_y = p4_y + (gap_out-gap_in)
    else:
        p2_y = y_tot_b_hi
        p4_y = y_tot_b_lo
        p3_y = p2_y - gap_out - p4_y
    return p2_y, p3_y, p4_y

def generic_bent_wg_pair(w_hi, r_hi, x_tot_b_hi,
                         y_tot_b_hi=None, btype_hi='cobra', w_lo=1, r_lo=10,
                         x_tot_b_lo=20, y_tot_b_lo=None, btype_lo='cobra',
                         gap_in=None, gap_out=None, xs='XS1',
                         st=False, put=False): # Default values set so that y_tot_b_hi is not necessary
    """
    A pair of independent, bent waveguides. Of the four variables gap_in,
    gap_out, y_tot_b_hi and y_tot_b_lo, we only need to define three.

       Parameters
    w_hi (float) : Width of the upper waveguide.
    r_hi (float) : Bending radius of the upper wg. Not used for 'cobra' and 'sin'.
    x_tot_b_hi (float) : Total bend displacement in x of the upper wg.
    y_tot_b_hi (float) : Total bend displacement in y of the upper wg.
    btype_hi (str) : Bend type of the upper waveguide.
    
    w_lo (float) : Width of the lower waveguide.
    r_lo (float) : Bending radius of the lower wg. Not used for 'cobra' and 'sin'.
    x_tot_b_lo (float) : Total bend displacement in x of the lower wg.
    y_tot_b_lo (float) : Total bend displacement in y of the lower wg.
    btype_lo (str) : Bend type of the lower waveguide.
    
    gap_in (float)  : Gap at the input.
    gap_out (float) : Gap at the output.
    
    xs (str, optional)   : Cross-section, specifying the layers.
    st (bool, optional)  : If True, prints notes and warnings.
    put (bool, optional) : If True, puts the cell.
    """
    w1 = w_hi
    w2 = w_lo
    r1 = r_hi
    r2 = r_lo
    btype1 = btype_hi
    btype2 = btype_lo
    # x_tot = 0 means completely vertical:
    p2_x = x_tot_b_hi
    p4_x = x_tot_b_lo
    
    p2_y, p3_y, p4_y = bwgp_get_py2_3_4_(gap_in, gap_out,
                                         y_tot_b_hi, y_tot_b_lo, st)

    pos2 = (p2_x,   p2_y,             0)
    pos3 = (0,      p3_y-(w1 + w2)/2, 0) # Accounts for w
    pos4 = (p4_x,   p4_y,             0)
    is_con = False
    GDBWG = make_bo_pair_(xs, w1, w2, r1, r2, pos2, pos3, pos4, btype1, btype2,
        is_con, st)
    if put:
        GDBWG.put(0)
    return GDBWG

def squircle(w, x_os, y_os, r, xs='XS1', st=False, put=False):
    """
    A rectangle with circular corners, or a racetrack.

       Parameters
    w (float) : Width.
    x_os (float) : The rectangle side lengths in x.
    y_os (float) : The rectangle side lengths in y.
        If one of them is zero, we have a racetrack.
    r (float) : Bending radius of all corners. Also defines size of corners.
    
    xs (str, optional)   : Cross-section, specifying the layers.
    st (bool, optional)  : If True, prints notes and warnings.
    put (bool, optional) : If True, puts the cell.
    """
    # So the center is at the origin:
    pos = ( -x_os/2, -y_os/2, 0 )
    sqr = make_sqrc_(xs, w, x_os, y_os, r, pos, st)
    if put:
        sqr.put(0)
    return sqr

def circle(w, r, xs='XS1', st=False, put=False):
    """
    A circle.

       Parameters
    w (float) : Width.
    r (float) : Radius.
    xs (str, optional) : Cross-section, specifying the layers.
    st (bool, optional) : If True, prints notes and warnings.
    put (bool, optional) : If True, puts the cell.
    """
    circ = make_sqrc_(xs, w, 0, 0, r, (0,0,0), st)
    if put:
        circ.put(0)
    return circ

def generic_beamsplit(w_in, l_in, w_hi, l_hi, r_hi, x_tot_b_hi, y_tot_b_hi,
                      btype_hi, w_lo, l_lo, r_lo, x_tot_b_lo, y_tot_b_lo, btype_lo,
                      gap=None, flip=False, xs='XS1', st=False, put=False):
    """
    A beam splitter with independent, bent waveguides. Consists of an input wg
    and two output wgs, the uppwer and lower waveguides. If flipped, instead
    has two inputs and one output wg.

       Parameters
    w_in (float) : Width of the input waveguide.
    l_in (float) : Length of input waveguide (in x).

    w_hi (float) : Width of the upper waveguide.
    l_hi (float) : Extra, strt length (in x) at the upper output.
    r_hi (float) : Bending radius of the upper wg. Not used for 'cobra' and 'sin'.
    x_tot_b_hi (float) : Total bend displacement in x of the upper wg.
    y_tot_b_hi (float) : Total bend displacement in y of the upper wg.
    btype_hi (str) : Bend type of the upper waveguide.

    w_lo (float) : Width of the lower waveguide.
    l_lo (float) : Extra, strt length (in x) at the lower output.
    r_lo (float) : Bending radius of the lower wg. Not used for 'cobra' and 'sin'.
    x_tot_b_lo (float) : Total bend displacement in x of the lower wg.
    y_tot_b_lo (float) : Total bend displacement in y of the lower wg.
    btype_lo (str) : Bend type of the lower waveguide.
    gap (float, optional) : Distance between the coupling regions of the
        output waveguides. NoneType means no gap. Positive means upper output
        moves up, negative means lower output moves down.

    flip (bool, optional) : If True, reflects the beam splitter across the y axis.
    xs (str, optional)   : Cross-section, specifying the layers.
    st (bool, optional)  : If True, prints notes and warnings.
    put (bool, optional) : If True, puts the cell.
    """
    w0 = w_in
    w1 = w_hi
    w2 = w_lo
    l0 = l_in
    l1 = l_hi
    l2 = l_lo
    r1 = r_hi
    r2 = r_lo
    btype1 = btype_hi
    btype2 = btype_lo
    # x_tot = 0 means completely vertical:
    pos2_x = l0 + x_tot_b_hi
    pos3_x = l0 + x_tot_b_lo
    pos2 = (pos2_x, y_tot_b_hi, 0)
    pos3 = (pos3_x, y_tot_b_lo, 0)
    GBS = make_bs_(xs, w0, w1, w2, l0, l1, l2, r1, r2, pos2, pos3,
        btype1, btype2, gap, flip, st)
    if put:
        GBS.put(0)
    return GBS

def sym_beamsplit(w, l_in, l_out, x_tot_b, y_tot_b, r, btype, gap=None, flip=False,
                  xs='XS1', st=False, put=False):
    """
    A beam splitter whose top and bottom outputs are symmetrical.

       Parameters
    w_in (float)  : Width of all waveguides.
    l_in (float)  : Length of input waveguide (in x).
    l_out (float) : Extra, strt length (in x) at the outputs.
    x_tot_b (float) : Total bend displacement in x of the upper output wg.
    y_tot_b (float) : Total bend displacement in y of the upper output wg.
    r (float)   : Bending radii of the output wgs. Not used for 'cobra' and 'sin'.
    btype (str) : Bend type of the output wgs.
    gap (float, optional) : Distance between the coupling regions of the
        output waveguides. NoneType means no gap. Positive means upper output
        moves up, negative means lower output moves down.

    flip (bool, optional) : If True, reflects the beam splitter across the y axis.
    xs (str, optional)   : Cross-section, specifying the layers.
    st (bool, optional)  : If True, prints notes and warnings.
    put (bool, optional) : If True, puts the cell.
    """
    # x_tot = 0 means completely vertical:
    pos2_x = l_in + 2*r + x_tot_b
    pos2 = ( pos2_x, y_tot_b, 0 )

    symBS = make_sym_bs_(xs, w, l_in, l_out, r, pos2, btype, gap, flip, st)
    if put:
        symBS.put(0)
    return symBS

def bump_wg(w, l_bump, x_tot, x_bump_middle, bump_height, r_hi, r_lo, ang,
        btype_in, btype_out, flip=False, xs='XS1', st=False, put=False):
    """
    A waveguide with a symmetric bump. From input to output, the bump consists
    of a bent btype wg, a circular arc, an optional flat plateau (in x), a
    circular arc and a bent btype wg.
    
       Parameters
    w (float) : Width.
    l_bump (float) : Length of flat plateau, the middle of the bump. Usually zero.
    x_tot (float)  : Length in x from the input to the output.
    x_bump_middle (float) : Length in x from the input to the middle of the bump.
        Also works for negative values.
    bump_height (float) : Height of the bump's top, above the input position x.
    r_hi (float) : Bending radii of the circular arcs.
    r_lo (float) : Bending radii of the bent i/o waveguides.
    ang (float)  : Maximum angle of each circular arc. The arcs go through
        [0,ang]. By default, the zero angle is straight up, since the arcs
        start from there and go outwards and down.
    btype_in (str)  : Bend type of the input waveguide.
    btype_out (str) : Bend type of the output waveguide.
    flip (bool, optional) : If True, reflects the wg across the x axis.
    xs (str, optional)   : Cross-section, specifying the layers.
    st (bool, optional)  : If True, prints notes and warnings.
    put (bool, optional) : If True, puts the cell.
    """
    btype1 = btype_in
    btype2 = btype_out
    
    pos_arc_x = x_bump_middle - l_bump/2
    pos_arc   = ( pos_arc_x, bump_height, 0 )
    pos_strt2 = ( x_tot,     0,           0 )
    BW = make_bump_wg_(xs, w, l_bump, r_hi, r_lo, ang, pos_arc, pos_strt2,
        btype1, btype2, flip, st)
    if put:
        BW.put(0)
    return BW

def strt_racetrack_reson(gap, w_co, l_co, w_tk, x_os_tk, y_os_tk, r, flip=False,
        xs='XS1', is_ring=False, st=False, put=False):
    """
    A straight racetrack or ring resonator.

       Parameters
    gap (float)  : Gap size between coupling region and resonator.
    w_co (float) : Width of the (strt) coupling region.
    l_co (float) : Length of the (strt) coupling region.
    w_tk (float) : Width of the resonator.
    x_os_tk (float) : The rectangle side lengths of the track/ring in x.
    y_os_tk (float) : The rectangle side lengths of the track/ring in y.
        If one of them is zero, we have a racetrack.
    r (float) : Bending radius of the resonator. Also defines size of corners.
    flip (bool, optional) : If True, reflects the resonator across the x axis.
    xs (str, optional) : Cross-section, specifying the layers.
    is_ring (bool, optional) : If True, creates a ring resonator.
    st (bool, optional)  : If True, prints notes and warnings.
    put (bool, optional) : If True, puts the cell.
    """
    x_os = x_os_tk
    y_os = y_os_tk
    
    pos_tk_x = (l_co - x_os)/2
    # Above line makes coupled region symmetric about the middle of the track
    pos_tk_y = gap + r + w_co/2
    pos_tk   = ( pos_tk_x, pos_tk_y, 0 )
    RaR = make_rtr_(xs, 'srt', w_co, w_tk, x_os, y_os, r, l_co, 0, 0, 0, 0,
        pos_tk, 0, 0, flip, is_ring, st)
    if put:
        RaR.put(0)
    return RaR

def strt_ring_reson(gap, w_co, l_co, w_ring, r, flip=False, xs='XS1',
                    st=False, put=False):
    """
    A straight ring resonator.

       Parameters
    gap (float)  : Gap size between coupling region and resonator.
    w_co (float) : Width of the (strt) coupling region.
    l_co (float) : Length of the (strt) coupling region.
    w_ring (float) : Width of the ring resonator.
    r (float) : Bending radius of the ring resonator.
    flip (bool, optional) : If True, reflects the resonator across the x axis.
    xs (str, optional)   : Cross-section, specifying the layers.
    st (bool, optional)  : If True, prints notes and warnings.
    put (bool, optional) : If True, puts the cell.
    """
    RiR = strt_racetrack_reson(gap, w_co, l_co, w_ring, 0, 0,
            r, flip, xs, is_ring=True, st=st, put=put)
    if put:
        RiR.put(0)
    return RiR

def sym_racetrack_reson(gap, w_co, l_co, r_co, w_tk, x_os_tk,
        y_os_tk, r_tk, x_tot, bump_height, ang, btype_in, btype_out, flip=False,
        xs='XS1', is_ring=False, st=False, put=False):
    """
    A symmetric racetrack or ring resonator.

       Parameters
    gap (float)  : Gap size between coupling region and resonator.
    w_co (float) : Width of the (bent) coupling region.
    l_co (float) : Length of the (bent) coupling region, the flat plateau of
        the bump.
    r_co (float) : Bending radius of the bent i/o waveguides of the coupled
        region, at the bottom.
    w_tk (float) : Width of the resonator.
    x_os_tk (float) : The rectangle side lengths of the track/ring in x.
    y_os_tk (float) : The rectangle side lengths of the track/ring in y.
        If one of them is zero, we have a racetrack.
    r_tk (float) : Bending radius of the resonator. Also defines size of corners.
    x_tot (float)  : Length in x from the input to the output.
    bump_height (float) : Height of the plateau of the coupled region's bump,
        where the input/output is the zero level.
    ang (float) : Maximum angle of each circular (bump) arc. The arcs go through
        [0,ang]. By default, the zero angle is straight up, since the arcs
        start from there and go outwards and down.
    btype_in (str)  : Bend type of the input waveguide.    
    btype_out (str) : Bend type of the output waveguide.

    flip (bool, optional) : If True, reflects the resonator across the x axis.
    xs (str, optional)   : Cross-section, specifying the layers.
    is_ring (bool, optional) : If True, creates a ring resonator.
    st (bool, optional)  : If True, prints notes and warnings.
    put (bool, optional) : If True, puts the cell.
    """
    btype1 = btype_in
    btype2 = btype_out
    
    pos_arc_x = (x_tot - l_co)/2
    pos_tk_x  = (x_tot - x_os_tk)/2
    # Above code makes coupled region symmetric about the middle of the track
    pos_tk_y  = r_tk + bump_height + gap + w_co/2

    pos_arc   = ( pos_arc_x, bump_height, 0 )
    pos_strt2 = ( x_tot,     0,           0 )
    pos_tk    = ( pos_tk_x,  pos_tk_y,    0 )
    
    RaR = make_rtr_(xs, 'sym', w_co, w_tk, x_os_tk, y_os_tk, r_tk, l_co, r_co,
        ang, pos_arc, pos_strt2, pos_tk, btype1, btype2, flip, is_ring, st)
    if put:
        RaR.put(0)
    return RaR

def sym_ring_reson(gap, w_co, l_co, r_co, w_ring, r_ring, x_tot, bump_height, ang,
        btype_in, btype_out, flip=False, xs='XS1', st=False, put=False):
    """
    A symmetric ring resonator.

       Parameters
    gap (float)  : Gap size between coupling region and resonator.
    w_co (float) : Width of the (bent) coupling region.
    l_co (float) : Length of the (bent) coupling region, the flat plateau of
        the bump.
    r_co (float) : Bending radius of the bent i/o waveguides of the coupled
        region, at the bottom.
    w_ring (float) : Width of the ring resonator.
    r_ring (float) : Bending radius of the ring resonator.
    x_tot (float)  : Length in x from the input to the output.
    bump_height (float) : Height of the plateau of the coupled region's bump,
        where the input/output is the zero level.
    ang (float) : Maximum angle of each circular (bump) arc. The arcs go through
        [0,ang]. By default, the zero angle is straight up, since the arcs
        start from there and go outwards and down.
    btype_in (str)  : Bend type of the input waveguide.    
    btype_out (str) : Bend type of the output waveguide.

    flip (bool, optional) : If True, reflects the resonator across the x axis.
    xs (str, optional)   : Cross-section, specifying the layers.
    st (bool, optional)  : If True, prints notes and warnings.
    put (bool, optional) : If True, puts the cell.
    """
    RiR = sym_racetrack_reson(gap, w_co, l_co, r_co, w_ring, 0, 0, r_ring,
            x_tot, bump_height, ang, btype_in, btype_out, flip, xs,
            is_ring=True, st=st, put=put)
    if put:
        RiR.put(0)
    return RiR

def pulley_racetrack_reson(gap, w_co, l_co_extra, r_co, w_tk,
        x_os_tk, y_os_tk, r_tk, x_tot, bump_height, ang, btype_in,
        btype_out, flip=False, xs='XS1', is_ring=False, st=False, put=False):
    """
    A pulley racetrack or ring resonator.

       Parameters
    gap (float)  : Gap size between coupling region and resonator.
    w_co (float) : Width of the (bent) coupling region.
    l_co_extra (float) : Extra length of the (strt) coupling region, the bottom
        of the bump.
    r_co (float) : Bending radius of the bent i/o waveguides of the coupled
        region.
    w_tk (float) : Width of the resonator.
    x_os_tk (float) : The rectangle side lengths of the track/ring in x.
    y_os_tk (float) : The rectangle side lengths of the track/ring in y.
        If one of them is zero, we have a racetrack.
    r_tk (float)  : Bending radius of the resonator. Also defines the of corners.
    x_tot (float)  : Length in x from the input to the output.
    bump_height (float) : Height of the plateau of the coupled region's bump,
        where the input/output is the zero level.
    ang (float) : Maximum angle of each circular (bump) arc. The arcs go through
        [0,ang]. By default, the zero angle is straight up, since the arcs
        start from there and go outwards and down.
    btype_in (str)  : Bend type of the input waveguide.    
    btype_out (str) : Bend type of the output waveguide.

    flip (bool, optional) : If True, reflects the resonator across the x axis.
    xs (str, optional)   : Cross-section, specifying the layers.
    is_ring (bool, optional) : If True, creates a ring resonator.
    st (bool, optional)  : If True, prints notes and warnings.
    put (bool, optional) : If True, puts the cell.
    """
    l_co = l_co_extra
    btype1 = btype_in
    btype2 = btype_out
    
    pos_arc_x = (x_tot - x_os_tk - l_co)/2
    pos_tk_x  = (x_tot - x_os_tk)/2
    # Above lines makes coupled region symmetric about the middle of the track
    pos_tk_y  = r_tk + bump_height + gap + w_co/2
    
    pos_arc   = ( pos_arc_x, bump_height, 0 )
    pos_strt2 = ( x_tot,     0,           0 )
    pos_tk    = ( pos_tk_x,  pos_tk_y,    0 )
    
    RaR = make_rtr_(xs, 'pul', w_co, w_tk, x_os_tk, y_os_tk, r_tk, l_co, r_co,
        ang, pos_arc, pos_strt2, pos_tk, btype1, btype2, flip, is_ring, st)
    if put:
        RaR.put(0)
    return RaR

def pulley_ring_reson(gap, w_co, l_co_extra, r_co, w_ring, r_ring, x_tot,
        bump_height, ang, btype_in, btype_out, flip=False,
        xs='XS1', st=False, put=False):
    """
    A pulley ring resonator.

       Parameters
    gap (float)  : Gap size between coupling region and resonator.
    w_co (float) : Width of the (bent) coupling region.
    l_co_extra (float) : Extra length of the (strt) coupling region, the bottom
        of the bump.
    r_co (float) : Bending radius of the bent i/o waveguides of the coupled
        region.
    w_ring (float) : Width of the ring resonator.
    r_ring (float) : Bending radius of the ring resonator.
    x_tot (float)  : Length in x from the input to the output.
    bump_height (float) : Height of the plateau of the coupled region's bump,
        where the input/output is the zero level.
    ang (float) : Maximum angle of each circular (bump) arc. The arcs go through
        [0,ang]. By default, the zero angle is straight up, since the arcs
        start from there and go outwards and down.
    btype_in (str)  : Bend type of the input waveguide.    
    btype_out (str) : Bend type of the output waveguide.

    flip (bool, optional) : If True, reflects the resonator across the x axis.
    xs (str, optional)   : Cross-section, specifying the layers.
    st (bool, optional)  : If True, prints notes and warnings.
    put (bool, optional) : If True, puts the cell.
    """
    RiR = pulley_racetrack_reson(gap, w_co, l_co_extra, r_co, w_ring, 0, 0,
            r_ring, x_tot, bump_height, ang, btype_in, btype_out, flip, xs,
            is_ring=True, st=st, put=put)
    if put:
        RiR.put(0)
    return RiR

def rect_spiral_delay(w, gap, btype, l_tot=None, x_os=None, y_os=None, r_lo=None,
        r_b=None, n_spirals=2, dx_start=10, dy_start=10, drl_start=10, max_steps=100,
        l_error=0.01, xy_same=False, flip=False, xs='XS1', st=False, put=False):
    """
    A rectangle spiral delay line, consisting of two spiral objects and a center
    waveguide. Each spiral is a collection of circular arcs and straight
    waveguides, a squircle-like shape. If l_tot is not given, x_os, y_os and
    r_lo all need to be given. If l_tot is given, x_os, y_os and r_lo will
    be calculated automatically. Therefore, out of x_os, y_os and r_lo, define
    none of them, or at most two.

       Parameters
    w (float)   : Width.
    gap (float) : Gap size between neighboring spirals.
    btype (str) : Bend type of the center waveguide connecting the two spiral objects.
    l_tot (float, optional) : Length of the entire delay line.
    x_os (float, optional) : The rectangle side lengths of the spiral in x.
    y_os (float, optional) : The rectangle side lengths of the spiral in y.
    r_lo (float, optional) : The smallest bending radius of both spiral parts.
    r_b (float, optional)  : Bending radius of the center waveguide. If not
        defined, becomes r_lo/2 - w.
    n_spirals (int, optional) : Number of spirals (number of iterations).
    dx_start (float, optional)  : In optimizing x_os, y_os and r_lo, starting
        correction for x_os.
    dy_start (float, optional)  : In optimizing x_os, y_os and r_lo, starting
        correction for y_os.
    drl_start (float, optional) : In optimizing x_os, y_os and r_lo, starting
        correction for r_lo.
    max_steps (int, optional) : Max number of steps in optimizing x_os, y_os and r_lo.
    l_error (float, optional) : The largest allowed error between the total
        desired length l_tot and actual length of the delay line.
    xy_same (bool, optional) : If True, x_os and y_os are equal (in optimization).
    flip (bool, optional) : If True, reflects the delay line across the x axis.
    xs (str, optional)   : Cross-section, specifying the layers.
    st (bool, optional)  : If True, prints notes and warnings.
    put (bool, optional) : If True, puts the cell.
    """
    l = l_tot
    dx  = dx_start
    dy  = dy_start
    drl = drl_start
    
    is_circ = False
    RSD = make_delay_line_(xs, w, gap, btype, l, x_os, y_os, r_lo, r_b, n_spirals,
        dx, dy, drl, max_steps, l_error, xy_same, flip, is_circ, st)
    
    if put:
        RSD.put(0)
    return RSD

def circ_spiral_delay(w, gap, btype, l_tot=None, r_lo=None, r_b=None, n_spirals=2,
        drl_start=10, max_steps=100, l_error=0.01, flip=False, xs='XS1', st=False,
        put=False):
    """
    A circular spiral delay line, consisting of two spiral objects and a center
    waveguide. Each spiral is a collection of circular arcs, a squircle-like shape.
    Define either l_tot or r_lo.
    
       Parameters
    w (float)   : Width.
    gap (float) : Gap size between neighboring spirals.
    btype (str) : Bend type of the center waveguide connecting the two spiral objects.
    l_tot (float, optional) : Length of the entire delay line.
    r_lo (float, optional)  : The smallest bending radius of both spiral parts.
    r_b (float, optional)   : Bending radius of the center waveguide. If not
        defined, becomes r_lo/2 - w.
    n_spirals (int, optional) : Number of spirals (number of iterations).
    drl_start (float, optional) : In optimizing r_lo, starting correction for r_lo.
    max_steps (int, optional) : Max number of steps in optimizing r_lo.
    l_error (float, optional) : The largest allowed error between the total
        desired length l_tot and actual length of the delay line.
    flip (bool, optional) : If True, reflects the delay line across the x axis.
    xs (str, optional)   : Cross-section, specifying the layers.
    st (bool, optional)  : If True, prints notes and warnings.
    put (bool, optional) : If True, puts the cell.
    """
    l = l_tot
    drl = drl_start
    
    xy_same = False
    is_circ = True
    CSD = make_delay_line_(xs, w, gap, btype, l, 0, 0, r_lo, r_b, n_spirals,
        0, 0, drl, max_steps, l_error, xy_same, flip, is_circ, st)
    
    if put:
        CSD.put(0)
    return CSD

def pad(x, y, xs='XS1', st=False, put=False):
    """
    A rectangular conductor pad cell. Used at inputs and outputs of conductors
    (phase shifters).

       Parameters
    x (float) : Length in x.
    y (float) : Length in y.
    xs (str, optional)  : Cross-section, specifying the layers.
    st (bool, optional) : If True, prints notes and warnings.
    """
    Pad = make_pad_(xs, x, y, st)
    if put:
        Pad.put(0)
    return Pad

def thermo_ps(all_w, all_r, p0, all_btype, all_pts, xy_pad=None, pads=False,
        flip=False, xs='XS2', layers_pad=['L1'], st=False, put=False):
    """
    A thermo-optic phase shifter (a conductor). 'all_' lists go from input
    to output.

       Parameters
    all_w (list of floats) : Widths.
    all_r (list of floats) : All bending radii. Not used for 'cobra' and 'sin'.
    p0 (tuple of floats)   : The starting position (x,y,a), the input.
    all_btype (list of strings) : All bend types.
    all_pts (list of floats) : All points (x,y,a) through the phase shifter,
        relative to p0, from input to output.
    xy_pad (tuple, optional) : Lengths in x,y of first pad, then x,y of second
        pad. Has the form (x1, y1, x2, y2).
    pads (bool, optional) : If True, adds pads.
    flip (bool, optional) : If True, reflects the phase shifter across the x axis.
    xs (str, optional)   : Cross-section of conductor, specifying the layers.
    layers_pad (str, optional) : Layers of pad.
    st (bool, optional)  : If True, prints notes and warnings.
    put (bool, optional) : If True, puts the cell.
    """
    to_ps = make_to_ps_(xs, all_w, all_r, p0, all_btype, all_pts, layers_pad,
            xy_pad, pads, flip, st)
    if put:
        to_ps.put(0)
    return to_ps

def electro_ps(all_w1, all_w2, all_r1, all_r2, p0, all_btype1, all_btype2, all_pts1,
        all_pts2, xy_pad=None, pads=False, flip=False, xs='XS2', layers_pad=['L1'],
        st=False, put=False):
    """
    An electro-optic phase shifter (a conductor). 'all_' lists go from input
    to output. Has two parts: part 1 (usually the top part) and part 2 (usually
    the bottom part). Part 1 goes from input to middle, while part 2 goes from
    middle to output. Near the middle, the two parts are parallel, creating the
    desired electric field. There, one part is on the top and the other is on
    the bottom.

       Parameters
    all_w1 (list of floats) : Widths of part 1 of the phase shifter.
    all_w2 (list of floats) : Widths of part 2 of the phase shifter.
    all_r1 (list of floats) : All bending radii of part 1. Not used for 'cobra'
        and 'sin'.
    all_r2 (list of floats) : All bending radii of part 2. Not used for 'cobra'
        and 'sin'.
    p0 (tuple of floats) : The starting position (x,y,a), the input.
    all_btype1 (list of strings) : All bend types of part 1.
    all_btype2 (list of strings) : All bend types of part 2.
    all_pts1 (list of floats) : All points (x,y,a) through part 1 of the phase
        shifter, relative to p0, from input to output.
    all_pts2 (list of floats) : All points (x,y,a) through part 2 of the phase
        shifter, relative to p0, from input to output.
    xy_pad (tuple, optional) : Lengths in x,y of first pad, then x,y of second
        pad. Has the form (x1, y1, x2, y2).
    pads (bool, optional) : If True, adds pads.
    flip (bool, optional) : If True, reflects the phase shifter across the x axis.
    xs (str, optional)   : Cross-section of conductor, specifying the layers.
    layers_pad (str, optional) : Layers of pad.
    st (bool, optional)  : If True, prints notes and warnings.
    put (bool, optional) : If True, puts the cell.
    """
    eo_ps =  make_eo_ps_(xs, all_w1, all_w2, all_r1, all_r2, p0, all_btype1,
            all_btype2, all_pts1, all_pts2, layers_pad, xy_pad, pads, flip, st)
    if put:
        eo_ps.put(0)
    return eo_ps

def pap_update_values_(i_lo, i_lo2, i_hi, i_max, all_w, all_r, all_btype, all_pts,
            p0=None, w=None, r=None, btype=None, point=None, ori_type='zero',
            ori_ang=False, st=False):
    """
    Updates w, r, btype and point of the new phase shifter point that is added
    in ps_add_point(). Takes into account ori_type, which can be 'zero', 'p0',
    'prev', 'next', 'avg' or an integer beween 0 and index_max. w, r, btype
    and point can start as NoneTypes. The point then takes all of the properties
    from around it in the ps.
    """
    # Takes the value of the conductor that is 'split in half'
    if w == None:
        w = all_w[i_lo]
    if r == None:
        r = all_r[i_lo]
    if btype == None:
        btype = all_btype[i_lo]
    if point == None:
        # Average position
        pl = all_pts[i_lo2]
        ph = all_pts[i_hi]
        point = ( (pl[0]+ph[0])/2, (pl[1]+ph[1])/2, (pl[2]+ph[2])/2 )
    
    # Integrate ori_type
    ori_type_is_str = None
    try: # Tests if ori_type is a string
        str_test = ori_type.lower()
        ori_type_is_str = True
    except ValueError:
        ori_type_is_str = False
    
    if ori_type_is_str:
        if   ori_type.lower() == 'zero':
            point = point # No change
        elif ori_type.lower() == 'p0':
            # p0 is origin
            point = add_pos(point, p0, add_a=ori_ang)
        elif ori_type.lower() == 'prev':
            # Previous point is origin
            point = add_pos(point, all_pts[i_lo2], add_a=ori_ang)
        elif ori_type.lower() == 'next':
            # Next point is origin
            point = add_pos(point, all_pts[i_hi], add_a=ori_ang)
        elif ori_type.lower() == 'avg':
            # Origin is the average of previous and last point
            pl = all_pts[i_lo2]
            ph = all_pts[i_hi]
            p_avg = ( (pl[0]+ph[0])/2, (pl[1]+ph[1])/2, (pl[2]+ph[2])/2 )
            point = add_pos(point, p_avg, add_a=ori_ang)
        else:
            raise ValueError("ori_type should be 'zero', 'p0', 'prev', 'next', "+
                             "'avg' or a number.")
    else:
        if int(ori_type) >= 0 and int(ori_type) <= i_max:
            # Origin is the point at index int(ori_type)
            i = int(ori_type)
            point = add_pos(point, all_pts[i], add_a=ori_ang)
        else:
            raise ValueError("Invalid integer value of ori_type.")
    return w, r, btype, point

def ps_add_point(all_w, all_r, all_btype, all_pts, p0=None, index=0, w=None,
            r=None, btype=None, point=None, ori_type='zero', 
            ori_ang=False, ps_index=None, hi_or_lo=None, st=False):
    """
    Adds a point to a phase shifter. The position depends on point and
    ori_type. Returns updated lists all_w, all_r, all_btype, all_pts. Can
    also be used for aoa lists when designing an NxN unitary circuit.
    
       Parameters
    all_w (list of floats) : Widths between all points in all_pts.
    all_r (list of floats) : All bending radii. Not used for 'cobra' and 'sin'.
    all_btype (list of string) : All bend types.
    all_pts (list of tuples) : All points (x,y,a) through the phase shifter,
        relative to p0, from input to output.
    p0  (tuple, optional) : The starting position (x,y,a) of the ps, the input.
    index (int, optional) : Index in all_pts for the added point.
    w (float, optional)  : Width right after the added point.
    r (float, optional)  : Bending radius right after the added point.
    btype (str, optional) : Bend type right after the added point.
    point (tuple, optional) : Position (x,y,a) of added point.
    ori_type (str or int, optional) : Chooses the origin of 'point'. 'zero'
        is (0,0,0). 'p0' is p0. 'prev' is the previous point. 'next' is the 
        next point. 'avg' is the average between the previous and next. If an
        integer n, the position at n in all_pts is the origin.
    ori_ang (bool, optional) : If False, only the position (x,y) of the
        origin gets added to 'point'. If True, the angle a also gets added.
    ps_index (int, optional) : If NoneType, does nothing. If int, treats
        all lists with prefix 'all_' as 'aoa' lists of lists (see below), and
        picks out the list at index ps_index.
    hi_or_lo (str, optional) : If NoneType, does nothing. If 'hi', treats
        all lists with prefix 'all_' as 'aoa' lists of lists (see below), and
        picks out the aoa list at 0, corresponding to ps's with inputs/outputs
        above the circuit. If 'lo', picks the aoa list at 1, for ps's with i/o
        below the circuit.
        For example, all_w[0][ps_index] gives the real all_w of upper ps number
        'ps_index'. all_w[1][ps_index] gives the real all_w for the lower ps.
    st (bool, optional) : If True, prints notes and warnings.
    """
    if ps_index != None and hi_or_lo != None:
        # If aoa lists
        n = None
        if hi_or_lo.lower() == 'hi':
            n = 0
        elif hi_or_lo.lower() == 'lo':
            n = 1
        else:
            raise ValueError("hi_or_lo should be 'hi' or 'lo' or NoneType.")

        # Here, the starting lists are actually 'aoa_', not 'all_'
        # So we copy the lists to the correct variable names:
        aoa_w = all_w
        aoa_r = all_r
        aoa_btype = all_btype
        aoa_pts = all_pts
        
        # Get the correct 'all_' lists
        all_w = aoa_w[n][ps_index]
        all_r = aoa_r[n][ps_index]
        all_btype = aoa_btype[n][ps_index]
        all_pts = aoa_pts[n][ps_index]
        p0 = p0[n][ps_index]
    i_max = len(all_pts) - 1
    
    # Starts creating new 'all_' lists
    # Add all points before index 'index'
    all_w2 = all_w[:index]
    all_r2 = all_r[:index]
    all_btype2 = all_btype[:index]
    all_pts2 = all_pts[:index]
    
    # Define i_lo, i_lo2 and i_hi
    i_lo = index - 1 # Used for w, r, btype in _pap_update_values()
    if i_lo < 0:
        i_lo = 0
    i_lo2 = i_lo # Used for 'points' in _pap_update_values()
    if i_lo > i_max-1:
        i_lo = i_max-1
    if i_lo2 > i_max:
        i_lo2 = i_max
    i_hi = index # Used for 'points' in _pap_update_values()
    if i_hi > i_max:
        i_hi = i_max

    # Update values at point 'index'
    w, r, btype, point = pap_update_values_(i_lo, i_lo2, i_hi, i_max, all_w,
        all_r, all_btype, all_pts, p0, w, r, btype, point, ori_type,
        ori_ang, st)
    
    # Add new point at index 'index'
    all_w2.append(w)
    all_r2.append(r)
    all_btype2.append(btype)
    all_pts2.append(point)
    
    # Add all later points
    all_w2 += all_w[index:]
    all_r2 += all_r[index:]
    all_btype2 += all_btype[index:]
    all_pts2 += all_pts[index:]
    
    if ps_index != None:
        # Updates the aoa lists to have the changed ps at index ps_index
        aoa_w[n][ps_index] = all_w2
        aoa_r[n][ps_index] = all_r2
        aoa_btype[n][ps_index] = all_btype2
        aoa_pts[n][ps_index] = all_pts2
        return aoa_w, aoa_r, aoa_btype, aoa_pts
    else:
        return all_w2, all_r2, all_btype2, all_pts2

