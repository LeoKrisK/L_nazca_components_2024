# -*- coding: utf-8 -*-
"""
Created on Mon May 13 2024
Last edit  Thu Jun 27 2024

A template for creating a universal multiport interferometer design on a photonic
circuit platform, using our module 'L_nazca_components_2024', containing
'L_nazca_components_2024.NxN_circuit', which creates the circuit.

When having many modes, the phase shifter connections (conductors) overlap,
which is not desired. The current solution is to manually use
circ.NxN_swap_ps_io() to swap the input/output of two phase shifters, and
ncm.ps_add_point() to add more points in one conductor, such that the
conductors move around each other without touching.

@author: LeoKrisK
"""

import numpy as np
import nazca as nd
import nazca.demofab as demo
from nazca.interconnects import Interconnect
import random

import L_nazca_components_2024 as ncm
import L_nazca_components_2024.NxN_circuit as circ

## Dictionary where each layer name points to the corresponding gds-code
all_layers = {
    'L1': (1,0),  # 'Shallow_Guide' layer
    'L2': (5,0),  # 'Active' layer
    'L3': (10,0), # 'MetalDC' layer
    'L4': (3,10)  # 'Deep_stub' layer
    }

ncm.load_all_layers(all_layers)

## Dictionary of xsections (cross-sections), where the value is the
## corresponding layers of the xsection.
xs_ls = {}

### Makes 'XS1' xsection, which uses layer L1 #################################
xs1        = 'XS1'  # Name
xs_ls[xs1] = ['L1'] # List of the layers in XS1
ncm.make_xsection(xs1, xs_ls)

### Makes 'XS2' xsection, which uses layer L2 #################################
xs2        = 'XS2'  # Name
xs_ls[xs2] = ['L2'] # List of the layers in XS2
ncm.make_xsection(xs2, xs_ls)

### Makes 'XS3' xsection, which uses L1 *and* L2 ##############################
xs3        = 'XS3'        # Name
xs_ls[xs3] = ['L1', 'L2'] # List of the layers in XS2
ncm.make_xsection(xs3, xs_ls)

### Feel free to add more xsections below #####################################
#xs4        = ''   # Name
#xs_ls[xs4] = [''] # List of layers
#ncm.make_xsection(xs4, xs_ls)

#xs5        = ''   # Name
#xs_ls[xs5] = [''] # List of layers
#ncm.make_xsection(xs4, xs_ls)

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
############################# Default Values ##################################
st = False      # 'Show Text' - Print out errors, warnings etc.
xs = xs1        # Default cross-section, which is 'XS1'
w = 1.2         # Default width
btype = 'cobra' # Default bend type for bent objects (a cobra bend)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###################### NxN Unitary Circuit Variables ##########################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
N = 2                   # Number of dimension in matrix / number of modes
c_type = 'bell21clem'   # Type of circuit
# c_type can currently be 'clem16' or 'bell21clem', the design from
# "Optimal design for universal multiport interferometers" and the Clements
# design in "Further Compactifying Linear Optical Unitaries".
M = circ.NxNm_get_M(c_type, N, st) # Number of columns in circuit / number of rows in all matrices below

# In the below description, we assume c_type = 'clem16', the design from 2016.
# The NxN circuit consists of N rows and M=8*N+6 columns. All matrices with
# position information, of the form A[i][j], have dimensions MxN. In other
# words, for A[i][j], i specifies the column of the CIRCUIT, j the row of the
# CIRCUIT. Simultaneously, of course, i specifies the row of the MATRIX A,
# j the column of the MATRIX A, as defined by the notation A[i][j]. Wherever
# we say "columns" or "rows", we mean the physical circuit, not the matrix.

# The circuit has 8*N+4 columns since, for each row, the input has one taper,
# then one straight wg, then one bent wg. Similar for the output. These are
# independent of N, they add 4. 8*N comes from the blocks.

# The circuit might be said to be divided in N blocks. Each block has 8 columns.
# The block consists of MZI's stacked in a layer/column. The 'length' is thus
# one MZI, while the 'height' is however many MZI are used in parallel (so,
# however many are 'stacked'). If N=5, we have 2 stacked MZI in each block.
#
# There are always two types of blocks. We call them type 0 and type 1. They are
# slightly different for even and odd N. For odd N, type 0 has an MZI between
# modes (1,2), and between (3,4),..., until (N-2,N-1). Then the last mode is a
# straight wg, no MZI. For odd N, type 1, we instead start with the straight wg
# at mode 1, then have MZI's between (2,3), (4,5), until (N-1,N). For even N, type
# 0 blocks have N/2 MZI's, and no room for straight wgs. For even N type 1, the
# first and last modes are straight wgs, and we have N/2-1 MZI's in between.
#
# All circuits follow the block order: type 0, type 1, type 0,..., until we have
# N blocks. The first block is always type 0, at least in our scheme. Within
# each block, we first have (in column 1) a strt wg, then a bwg, then (in column
# 2), a strt wg that couples with a neighboring mode to create a 50:50 beam splitter,
# then a bwg. Then these four columns are then repeated to get 8 columns.
#
# NxNm_apy_block_type_0() and NxNm_apy_block_type_1() define the y values for blocks
# of type 0 and type 1. They call NxNm_apy_make_column(), which creates each column.
# The blocks are connected in series with NxNm_ap_y(). We define the x values
# with NxN_ap_x(). Everything is then combined with NxN_ap_tups().

###############################################################################
###############################################################################
################## WAVEGUIDE (wg) and TAPER (t) variables #####################
###############################################################################
###############################################################################
w0 = 1.2         # Default width of all wgs
wit0 = 0.5         # Default width at the input tapers (it)
wot0 = 0.5         # Default width at the ouput tapers (ot)
r0 = 100         # Default bending radius of all wg bends
btype0 = 'sbend' # Default bend type of wgs

## x distances (lengths):
lt0  = 200       # Default length of tapers
lio0 = 2000      # Default length of (straight) input/output wgs
lco  = 500       # Default length of coupling regions
l0   = 1000      # Default length of all other wgs
xb0  = 200       # How much each (small) bent wg (b) moves in x

## y distances (heights):
gap_large = 40           # The bigge rof the two y gaps, at the start and end (input and output)
gap_small = 5            # The smaller of the two y gaps in the middle (the beam splitter), for type 0 blocks
gap_large_t1 = gap_large # The bigger of the two y gaps, for type 1 blocks
gap_small_t1 = gap_small # The smaller of the two y gaps (the beam splitter), for type 1 blocks

def NxNm_get_wg_data(c_type, N, M, w0, wit0, wot0, btype0, r0, lt0, lio0, lco, l0,
        xb0, gap_large, gap_small, gap_large_t1, gap_small_t1, st=False):
    """
    Used in NxN unitary matrix circuit with Mach-Zehnder interferometer design.
    Returns lists and matrices with wg information.
    """
    # delta_x_list contains a list of the lengths in x of each element in the
    # wg circuit. In the MZI case, it looks like this:
    delta_x_list = None
    if   c_type.lower() == 'clem16':
        delta_x_list = [lio0] + (2*N)*[l0, xb0, lco, xb0] + [lio0]
        if N == 2:
            delta_x_list = [lio0] + 2*[l0, xb0, lco, xb0] + [lio0]
    elif c_type.lower() == 'bell21clem':
        delta_x_list = [lio0] + (N)*[0, xb0, lco, xb0, l0, xb0, lco, xb0] + [lio0]
        if N == 2:
            delta_x_list = [lio0] + [0, xb0, lco, xb0, l0, xb0, lco, xb0] + [lio0]
    
    # Gets all_pos, a matrix with all wg start and end points in the circuit
    # Assumes all wg have the same width w0
    all_pos = circ.NxNm_get_all_pos(c_type, N, w0, gap_large,
            gap_small, gap_large_t1, gap_small_t1, delta_x_list, st)
    all_lt = M*[N*[lt0]]
    all_pos = circ.NxN_ap_add_taper(all_pos, all_lt, st)
    
    # Between all points in all_pos, there is either a wg or a taper
    # This is captured by cell_type, a matrix that says if we have a wg or taper
    # Below, '' is a wg, while the input/output are tapers
    cell_type = [N*['taper']] + (M-3)*[N*['']] + [N*['taper']]

    all_w = [N*[wit0]] + (M-2)*[N*[w0]] + [N*[wot0]] # All wg/taper widths
    all_r = M*[N*[r0]]         # All wg bending radii
    all_btype = M*[N*[btype0]] # All wg bend types

    return all_pos, cell_type, all_w, all_r, all_btype

# Call NxN_get_wg_data()
all_pos, cell_type, all_w, all_r, all_btype = NxNm_get_wg_data(c_type, N, M, w0, wit0, wot0, btype0, r0, lt0, lio0, lco, l0, xb0, gap_large, gap_small, gap_large_t1, gap_small_t1, st)

###############################################################################
###############################################################################
####################### PHASE SHIFTER (ps) variables ##########################
###############################################################################
###############################################################################
###### Another name for the phase shifter is the CONDUCTOR
ps_type = 'to'      # Type of phase shifter (thermo-optic 'to' or electro-optic 'eo')
init_layer  = True  # If there is an initial later of phase shifters. Used with 'bell21clem'
final_layer = True  # If there is a final layer of phase shifters. Used with 'clem16' and 'bell21clem'
ps_flip_connections = False # If True, flips the order of ps i/o connections
pads = True         # If the conductors have pads

def NxNm_get_misc_ps_data(c_type, N, M, ps_flip_connections, init_layer=False,
                          final_layer=False, st=False):
    """
    Used in NxN unitary matrix circuit with Mach-Zehnder interferometer design.
    Returns various, extra floats, lists and matrices with ps information.
    """
    # Gets number of phase shifters:
    ps_num = circ.NxNm_ps_num(c_type, N, init_layer, final_layer, st)
    
    # If j is less than j_edge or greater than N-j_edge, the phase shifter auto-
    # matically connects to the closest side of the electrical i/o
    j_edge = int((N+1)/2) - 1.5
    
    ps_mid_ntype  = [0,1] * int(ps_num+1)
    # When phase shifers are 'near the middle', their ntype is taken
    # from this lists, which here is of the form [0,1,0,1,...].
    
    # ps_all_ij contains all (i,j) coordinates (in all_pos) where ps appear.
    # It is a list of two lists. The first list represents the ps's that
    # connect to the top of the circuit, the second list connect to the bottom.
    # Note: the ntype may change inside NxNm_ps_all_ij() if there are
    # too many ps on one side.
    ps_all_ij = circ.NxNm_ps_all_ij(c_type, N, M, init_layer, final_layer,
                        ps_mid_ntype, j_edge, ps_flip_connections, st)
    ps_num_hi = len(ps_all_ij[0]) # Number of ps connecting to above circuit
    ps_num_lo = len(ps_all_ij[1]) # Number of ps connecting to below circuit
    return ps_num, j_edge, ps_mid_ntype, ps_all_ij, ps_num_hi, ps_num_lo

# Call NxNm_get_misc_ps_data()
ps_num, j_edge, ps_mid_ntype, ps_all_ij, ps_num_hi, ps_num_lo = NxNm_get_misc_ps_data(c_type, N, M, ps_flip_connections, init_layer, final_layer, st)

###############################################################################
### Variables for ONE phase shifter (which are, by default, copied to all other
### ps's, but can be changed manually later down in the code) #################
ps_w0 = w/2       # Default width of conductors that do not interact to the wgs
ps_w1 = w         # Default width of conductors that interact to the wgs
ps_w2 = w         # Only for electro-optic: default width below the wg
ps_r = 30         # Default bending radius of conductors
ps_btype = 'bsb'  # Default bend type of conductors

ps_gap = 2*w      # Default gap size of ps active/interaction region (>0 is above wg)
ps_l = 500        # Default length of ps active region
ps_os = 0         # Default x offset of ps from starting x coordinate
ps_ang = 90       # Default angle at ps i/o. 90 means straight up/down

ps_xs = 'XS2'     # Default xsection of conductors (ps)
ps_layers_pad = ['L2'] # Default xsection of pads (for ps)
# Size of left pad
x_pad1 = 200
y_pad1 = 200
# Size of right pad
x_pad2 = 200
y_pad2 = 200
# Combine pad information
ps_xy_pad = (x_pad1, y_pad1, x_pad2, y_pad2)

ps_pi_ex = [[(0,-200,0),(45,-380,45)]]  # Extra points between input and active region, for one ps
ps_po_ex = [[(0,-200,0),(-45,-380,45)]] # Extra points between active region and output, for one ps

def NxN_circuit_xy(all_pos, init_layer=False, final_layer=False, st=False):
    """
    Used in NxN unitary matrix circuit with arbitrary design. Returns total
    length and height of the circuit.
    """
    # Gets total y of circuit
    circuit_height = all_pos[0][0][1] - all_pos[0][-1][1]

    # Gets total x of circuit
    x_start = all_pos[0][0][0]
    x_end = all_pos[-1][0][0]
    circuit_length = x_end - x_start

    return circuit_length, circuit_height

# Gets length and height of circuit
circuit_length, circuit_height = NxN_circuit_xy(all_pos, init_layer, final_layer, st)

##### Phase shift input/output variables ######################################
## x distances (lengths):
ps_x_hi = 221   # Starting x coordinate of ps i/o on the top left
ps_x_lo = 221   # Starting x coordinate of ps i/o on the bottom left
ps_dx_hi = 2800 # Change in x between two pads, top
ps_dx_lo = 2800 # Change in x between two pads, bottom

## y distances (heights):
ps_y_hi = 1300                        # Height (in y) of ps i/o on the top
ps_y_lo = -circuit_height - ps_y_hi  # Height (in y) of ps i/o on the bottom
# ps_y_lo reflects ps_y_hi so the distance above is the same as below
# There is no ps i/o change in y through the circuit, so ps_dy = 0

###############################################################################
##### Phase shifter variables, for ALL ps #####################################
# Copies from above for all ps

def NxNm_get_all_ps_data(ps_all_ij, ps_y_hi, ps_y_lo, ps_x_hi, ps_x_lo,
                ps_dx_hi, ps_dx_lo, ps_num_hi, ps_num_lo, ps_xs, ps_w1, ps_w2,
                ps_gap, ps_l, ps_os, ps_ang, ps_layers_pad, ps_xy_pad, ps_pi_ex,
                ps_po_ex, st=False):
    """
    Used in NxN unitary matrix circuit of arbitrary design. Returns lists and
    matrices with ps information. 'all_' lists are lists with one element per
    phase shifter.
    """
    ps_all_w1 = [ps_num_hi*[ps_w1]] + [ps_num_lo*[ps_w1]] # Width of the upper active region (for eo), or
    # width of the total active region (for to)
    ps_all_w2 = [ps_num_hi*[ps_w2]] + [ps_num_lo*[ps_w2]] # Width of the secondary (lower) active region (only for eo)
    ps_all_gap = [ps_num_hi*[ps_gap]] + [ps_num_lo*[ps_gap]] # All ps gaps
    ps_all_l   = [ps_num_hi*[ps_l]] + [ps_num_lo*[ps_l]]     # All ps lengths in x of the active regions
    ps_all_os  = [ps_num_hi*[ps_os]] + [ps_num_lo*[ps_os]]   # All ps offsets in x from the active region starting point
    ps_all_ang = [ps_num_hi*[ps_ang]] + [ps_num_lo*[ps_ang]] # All ps angles at the i/o
    
    # Gets ps_all_io: all ps i/o positions
    ps_all_io = circ.NxN_ps_all_io(ps_all_ij, ps_y_hi, ps_y_lo, ps_x_hi, ps_x_lo,
                                   ps_dx_hi, ps_dx_lo, ps_all_ang, st)
    
    ps_all_xs  = [ps_num_hi*[ps_xs]] + [ps_num_lo*[ps_xs]] # All ps cross-sections
    all_layers_pad = [ps_num_hi*[ps_layers_pad]] + [ps_num_lo*[ps_layers_pad]] # All cross-sections for all ps pads
    all_xy_pad = [ps_num_hi*[ps_xy_pad]] + [ps_num_lo*[ps_xy_pad]] # All size variables for all pads

    ps_aoa_pi_ex = [ps_num_hi*ps_pi_ex]+[ps_num_lo*ps_pi_ex] # Extra points between input and active region, for all ps
    ps_aoa_po_ex = [ps_num_hi*ps_po_ex]+[ps_num_lo*ps_po_ex] # Extra points between active region and output, for all ps

    return ps_all_io, ps_all_xs, ps_all_w1, ps_all_w2, ps_all_ang, ps_all_gap, ps_all_l, ps_all_os, all_layers_pad, all_xy_pad, ps_aoa_pi_ex, ps_aoa_po_ex

# Call NxNm_get_ps_data()
ps_all_io, ps_all_xs, ps_all_w1, ps_all_w2, ps_all_ang, ps_all_gap, ps_all_l, ps_all_os, all_layers_pad, all_xy_pad, ps_aoa_pi_ex, ps_aoa_po_ex = NxNm_get_all_ps_data(ps_all_ij, ps_y_hi, ps_y_lo, ps_x_hi, ps_x_lo, ps_dx_hi, ps_dx_lo, ps_num_hi, ps_num_lo, ps_xs, ps_w1, ps_w2, ps_gap, ps_l, ps_os, ps_ang, ps_layers_pad, ps_xy_pad, ps_pi_ex, ps_po_ex, st)

# If we do not want extra points between input and active region:
#ps_aoa_pi_ex = [[None],[None]]
# If we do not want extra points between active region and output:
#ps_aoa_po_ex = [[None],[None]]

# We can swap the ps i/o with the following function:
#ps_all_io = circ.NxN_swap_ps_io(
#    ps_all_io,
#    hi_or_lo = 'both',
#    index1 = 1,
#    index2 = 2,
#    io1 = 0,
#    io2 = 1,
#    st = st
#    )

# To swap many i/o, we can iterate the function in a loop:
#for i in range(3,8,2): # The 2 means i has values 3,5,7
#    ps_all_io = circ.NxN_swap_ps_io(
#        ps_all_io,
#        hi_or_lo = 'both',
#        index1 = i,
#        index2 = i+1,
#        io1 = 0,
#        io2 = 0,
#        st = st
#        )

# To delete a phase shifter, we use circ.NxN_delete_ps(). For example:
# del_index = 0
# hi_or_lo = 'lo'
# ps_num, ps_all_ij, ps_num_hi, ps_num_lo, ps_all_io, ps_all_xs, ps_all_w1, ps_all_w2, ps_all_ang, ps_all_gap, ps_all_l, ps_all_os, all_layers_pad, all_xy_pad, ps_aoa_pi_ex, ps_aoa_po_ex = circ.NxN_delete_ps(del_index, hi_or_lo, ps_type, ps_num, ps_all_ij, ps_num_hi, ps_num_lo, ps_all_io, ps_all_xs, ps_all_w1, ps_all_w2, ps_all_ang, ps_all_gap, ps_all_l, ps_all_os, all_layers_pad, all_xy_pad, ps_aoa_pi_ex, ps_aoa_po_ex, st)
# Deletes the first phase shifter on the bottom
# If hi_or_lo = 'both', it deletes the phase shifters both above and below the circuit.


###############################################################################
###############################################################################
###############################################################################
##### aoa phase shifter variables, lists of lists for all ps ##################
#### Naming Scheme:
# 'all' means a list. 'oa' means 'one_all', 'aoa' means 'all_one_all'.
# A 'all_' list can be two things. It is either a list for one phase shifter,
# with many elements, or a list for all phase shifters, with one element per
# phase shifter.
# 'aoa_' lists are lists with one 'all_' list per phase shifter.
# For example, ps_oa_w is a list of all widths of one ps.
# ps_aoa_w is a list of all widths of all PSs. That is, a list of
# ps_oa_w lists, one 'ps_oa_w' for each PS.

def NxNm_get_ps_pts_data(N, ps_type, ps_w0, ps_w1, ps_w2, ps_r, ps_btype, ps_all_io,
        ps_aoa_pi_ex, ps_aoa_po_ex, all_pos, ps_all_ij, all_w, ps_all_w1,
        ps_all_w2, ps_all_gap, ps_all_l, ps_all_os, ps_y_hi, ps_y_lo, ps_x_hi,
        ps_x_lo, ps_dx_hi, ps_dx_lo, st=False):
    """
    Used in NxN unitary matrix circuit of Mach-Zehnder interferometer design.
    Gets ps_all_p0 and ps_aoa_pts1. Defines all origins in each ps, and defines
    the relative positions of each part of each ps. If eo, also gives
    ps_aoa_pts2.
    """
    ps_oa_w1 = [[],[]]
    ps_oa_w2 = [[],[]]     # Only for electro-optic
    ps_oa_r1 = [[],[]]
    ps_oa_r2 = [[],[]]     # Only for electro-optic
    ps_oa_btype1 = [[],[]]
    ps_oa_btype2 = [[],[]] # Only for electro-optic

    ps_all_p0 = [[],[]]    # Origins of all ps, one p0 for each ps
    ps_aoa_pts1 = [[],[]]  # Default points of one ps, where p0 is the origin
    ps_aoa_pts2 = [[],[]]  # Only for electro-optic ps (points of lower part)

    n1 = 1 + len(ps_aoa_pi_ex[0][0]) # Number of bent objects before active region
    n2 = 1 + len(ps_aoa_po_ex[0][0]) # Number of bent objects after active region

    if ps_type.lower() == 'to': # If thermo-optic
        ps_oa_w1 = n1*[ps_w0] + [ps_w1] + n2*[ps_w0] # Default widths of one ps
        ps_oa_r1 = (n1+n2+1)*[ps_r] # Default bending radii of one ps
        ps_oa_btype1 = (n1+n2+1)*[ps_btype] # Default bend types of one ps
        
        ps_all_p0, ps_aoa_pts1 = circ.NxN_ps_pts(N, ps_type, ps_all_io, ps_aoa_pi_ex,
            ps_aoa_po_ex, all_pos, ps_all_ij, all_w, ps_all_w1, ps_all_w2,
            ps_all_gap, ps_all_l, ps_all_os, ps_y_hi, ps_y_lo, ps_x_hi, ps_x_lo,
            ps_dx_hi, ps_dx_lo, st)
    elif ps_type.lower() == 'eo': # If electro-optic
        ps_oa_w1 = n1*[ps_w0] + [ps_w1] # Default widths of upper part of one ps
        ps_oa_w2 = [ps_w2] + n2*[ps_w0] # Default widths of lower part of one ps
        ps_oa_r1 = (n1+1)*[ps_r] # Default bending radii of upper part of one ps
        ps_oa_r2 = (n2+1)*[ps_r] # Default bending radii of lower part of one ps
        ps_oa_btype1 = (n1+1)*[ps_btype] # Default bend types of upper part of one ps
        ps_oa_btype2 = (n2+1)*[ps_btype] # Default bend types of lower part of one ps
        
        ps_all_p0, ps_aoa_pts1, ps_aoa_pts2 = circ.NxN_ps_pts(N, ps_type, ps_all_io,
            ps_aoa_pi_ex, ps_aoa_po_ex, all_pos, ps_all_ij, all_w, ps_all_w1,
            ps_all_w2, ps_all_gap, ps_all_l, ps_all_os, ps_y_hi, ps_y_lo, ps_x_hi,
            ps_x_lo, ps_dx_hi, ps_dx_lo, st)
    return ps_all_p0, ps_aoa_pts1, ps_aoa_pts2, ps_oa_w1, ps_oa_w2, ps_oa_r1, ps_oa_r2, ps_oa_btype1, ps_oa_btype2

# Calls NxNm_get_ps_pts_data()
ps_all_p0, ps_aoa_pts1, ps_aoa_pts2, ps_oa_w1, ps_oa_w2, ps_oa_r1, ps_oa_r2, ps_oa_btype1, ps_oa_btype2 = NxNm_get_ps_pts_data(N, ps_type, ps_w0, ps_w1, ps_w2, ps_r, ps_btype, ps_all_io, ps_aoa_pi_ex, ps_aoa_po_ex, all_pos, ps_all_ij, all_w, ps_all_w1, ps_all_w2, ps_all_gap, ps_all_l, ps_all_os, ps_y_hi, ps_y_lo, ps_x_hi, ps_x_lo, ps_dx_hi, ps_dx_lo, st)

def NxN_get_aoa_ps_data(ps_num_hi, ps_num_lo, ps_oa_w1, ps_oa_w2, ps_oa_r1,
                ps_oa_r2, ps_oa_btype1, ps_oa_btype2, st=False):
    """
    Used in NxN unitary matrix circuit of arbitrary design. Returns aoa lists
    with ps information. 'aoa' means 'all_one_all', meaning the list has
    the information for all ps in an entire circuit. An aoa lists is a list of
    two lists. List 1 is for ps that connect above the circuit, list 2 connects
    below the circuit.
    """
    ## Defines ps_aoa_w, the widths of all phase shifters
    ps_aoa_w1 = [ps_num_hi*[ps_oa_w1]] + [ps_num_lo*[ps_oa_w1]] # For thermo-optic: all parts, for electro-optic (eo): upper part
    ps_aoa_w2 = [ps_num_hi*[ps_oa_w2]] + [ps_num_lo*[ps_oa_w2]] # Lower part (for eo only)
    ## Defines ps_aoa_r, the bending radii of all phase shifters
    ps_aoa_r1 = [ps_num_hi*[ps_oa_r1]] + [ps_num_lo*[ps_oa_r1]]
    ps_aoa_r2 = [ps_num_hi*[ps_oa_r2]] + [ps_num_lo*[ps_oa_r2]] # For eo only
    ## Defines ps_aoa_btype, the bend type of all part of all phase shifters
    ps_aoa_btype1 = [ps_num_hi*[ps_oa_btype1]] + [ps_num_lo*[ps_oa_btype1]]
    ps_aoa_btype2 = [ps_num_hi*[ps_oa_btype2]] + [ps_num_lo*[ps_oa_btype2]] # For eo only
    ## Defines cross-sections and pad variables
    return ps_aoa_w1, ps_aoa_w2, ps_aoa_r1, ps_aoa_r2, ps_aoa_btype1, ps_aoa_btype2

# Calls NxN_get_aoa_ps_data()
ps_aoa_w1, ps_aoa_w2, ps_aoa_r1, ps_aoa_r2, ps_aoa_btype1, ps_aoa_btype2 = NxN_get_aoa_ps_data(ps_num_hi, ps_num_lo, ps_oa_w1, ps_oa_w2, ps_oa_r1, ps_oa_r2, ps_oa_btype1, ps_oa_btype2, st)

###############################################################################
###############################################################################
##################### NxN Circuit, Further Changes ############################
###############################################################################
###############################################################################
# Here you can manually change elements of any of the above lists
# The lists (that hold all of the information about the circuit) are:
# Photonic / wg and taper variables:
### all_pos (list of lists) - all wg positions (the wgs get connected automatically)
### cell_type (list of lists) - types of cells (wg or taper)
### all_w (list of lists) - all wg widths
### all_r (list of lists) - all wg bending radii
### all_btype (list of lists) - all wg bend types

# NOTE: All lists below are divided in two parts. They are lists of lists.
# For example, ps_all_xs has two elements, and both are lists. List 1 is for
# conductors that connect above the circuit, list 2 connects below the circuit.
# The conductors are ordered from left to right (in the physical circuit).
# Electronic / phase shifter variables:
### ps_all_xs (list) - all ps cross-sections
### ps_aoa_w1 (list of lists) - all ps widths (upper part)
### ps_aoa_w2 (list of lists) - all ps widths (lower part, only for eo ps)
### ps_aoa_r1 (list of lists) - all ps bending radii (upper part)
### ps_aoa_r2 (list of lists) - all ps bending radii (lower part, only for eo ps)
### ps_all_p0 (list) - all PS origin coordinates
### ps_aoa_pts1 (list of lists) - all the points of each ps, where ps_all_p0 is
# the origin (upper part) (the parts get connected automatically)
### ps_aoa_pts2 (list of lists) - all the points of each ps, where ps_all_p0 is
# the origin (lower part, only for eo ps)
### ps_aoa_btype1 (list of lists) - all ps bend types (upper part)
### ps_aoa_btype2 (list of lists) - all ps bend types (lower part, only for eo ps)
### all_layers_pad (list) - all ps pad layers
### all_xy_pad (list) - all ps pad dimensions

# You can experiment with them to see what they do.

# For example,
# If we have thermo-optic ps, ps_aoa_w1[0][2][4] is the width of the fifth bent
# object (5+1) of the third (2+1) hi ps (connecting above the circuit, 0).
# If we have electro-optic ps, ps_aoa_w2[1][0][0] is the width of the first bent
# object (0+1) of the first (0+1) lo ps (connecting below the circuit, 1).

###############################################################################
# The below function adds a point to phase shifter number 'ps_index'.
# The point is added at index 'index' in all_pts of the ps, at position 'point'.
# ori_type designates the origin of the coordinate system for position 'point'.
# w, r and btype describe the conductor right after the added point.

#ps_aoa_w1, ps_aoa_r1, ps_aoa_btype1, ps_aoa_pts1 = ncm.ps_add_point(
#    ps_aoa_w1, ps_aoa_r1, ps_aoa_btype1, ps_aoa_pts1, ps_all_p0,
#    index = 3,
#    w = ps_w0,
#    r = 10,
#    btype = 'bsb',
#    point = (-500,0,0),
#    ori_type = 'avg',
#    ori_ang = True,
#    ps_index = 0,
#    hi_or_lo = 'hi',
#    st = st
#    )

# If w, r, btype and point are NoneType, the added point merely splits the
# conductor in half, keeping all of its properties before and after the added point.
# Between 0 and 4 of these values can be NoneType.

#ps_aoa_w1, ps_aoa_r1, ps_aoa_btype1, ps_aoa_pts1 = ncm.ps_add_point(
#    ps_aoa_w1, ps_aoa_r1, ps_aoa_btype1, ps_aoa_pts1, ps_all_p0,
#    index = 3,
#    w = None,
#    r = None,
#    btype = None,
#    point = (1000,1000,0),
#    ori_type = 'avg',
#    ori_ang = True,
#    ps_index = 1,
#    hi_or_lo = 'lo',
#    st = st
#    )



###############################################################################
###############################################################################
###############################################################################
# Use the above lists in NxN_circuit(), creating the circuit
c = circ.NxN_circuit(all_pos, cell_type, all_w, all_r, all_btype, ps_type, ps_all_xs, ps_aoa_w1, ps_aoa_w2,
    ps_aoa_r1, ps_aoa_r2, ps_all_p0, ps_aoa_pts1, ps_aoa_pts2, ps_aoa_btype1, ps_aoa_btype2, all_layers_pad, all_xy_pad, pads, st)
c_ = c.put(0) # Put the circuit

nd.export_gds()

