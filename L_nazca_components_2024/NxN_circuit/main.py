# -*- coding: utf-8 -*-
"""
Created on Thu May  2 2024
Last edit  Thu Jun 27 2024

Allows for the creation of universal multiport interferometers on photonic
circuit platforms. Call NxN_circuit() to make a cell of the circuit. See
the template file 'NxN Unitary Circuit Template' for one proposed way to
use this module. Currently supports interferometer designs 'clem16' and
'bell21clem', from "Optimal design for universal multiport interferometers" by
Clements et al., contra "Further Compactifying Linear Optical Unitaries" by Bell
and Walmsley. The building block of the first is an MZI containing two 50:50
beam splitters, one phase shifter between them and one preceding both beam
splitters. The second improves said design such that phase shifters only appear
between the 50:50 beam splitters, making the circuit more compact. Feel free to
copy our code and implement more designs.

@author: LeoKrisK
"""
# All lengths and coordinates are in microns (μm).

import numpy as np
import nazca as nd
import nazca.demofab as demo
from nazca.interconnects import Interconnect
import random

import L_nazca_components_2024.main as ncm

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
###############################################################################
########### NxN Unitary Circuit, Abstract Waveguide Grid Functions ############
###############################################################################
###############################################################################

def NxNm_get_M(c_type, N, st=False):
    """
    Used in NxN unitary matrix circuit with Mach-Zehnder interferometer design
    for wgs. Gets the number of rows in all_pos, which is the number of wg
    'columns' in the circuit.
    """
    M = None
    M = 8*N + 5
    if N==2: # Special case, 2x2, single MZI
        M -= 8
    return M

def NxN_apy_make_column(N, y_init, dy, st=False):
    """
    Used in NxN unitary matrix circuit with arbitrary design for wgs. Used for
    getting all_pos. Creates a vector with N elements, representing the positions
    of one vertical layer of N modes in the circuit.
    """
    y_column = []

    dy_list = [dy, -dy] # Change in y from y_init
    for j in range(N): # All Mach-Zehnder interferometer y values
        # j is the row
        y_column.append( y_init[j] + dy_list[j%2] )
    return y_column

def NxNm_ap_y(c_type, N, w, gap_large, gap_small, st=False):
    """
    Used in NxN unitary matrix circuit with Mach-Zehnder interferometer design
    for wgs. Used for getting all_pos. Gets all y positions in the circuit.
    """
    dy = abs((gap_large - gap_small)/2)
    y_list0 = []
    for n in range(N):
        y_list0.append( -n*(gap_large+w ) )
    y_tot = [y_list0] # Start of matrix
    
    dy_bt0, dy_bt1 = None, None
    dy_bt0 = [-dy, 0, dy, 0, -dy, 0, dy, 0]
    dy_bt1 = [dy, 0, -dy, 0, dy, 0, -dy, 0]
        
    dy_bt01 = dy_bt0 + dy_bt1
    
    dy_all_list = [0] + dy_bt01 * int((N+1)/2)
    
    if N%2==1:
        # Erases last block of type 1, not used
        dy_all_list = dy_all_list[:-9]
    if N==2: # Special case, 2x2, single MZI
        dy_all_list = dy_all_list[:-9]
    
    for dy in dy_all_list:
        new_column = NxN_apy_make_column(N, y_tot[-1], dy, st)
        y_tot.append(new_column)
        
    y_tot = [y_tot[0]] + y_tot + [y_tot[-1]]
    if N==2: # Special case, 2x2, single MZI
        y_tot += [y_tot[-1]]*14
    return y_tot

def NxNm_ap_x(N, delta_x_list, st=False):
    """
    Used in NxN unitary matrix circuit with Mach-Zehnder interferometer design
    for wgs. Used for getting all_pos. Gets all x positions in the circuit.
    """
    M = len(delta_x_list)+1
    # This is 2 less then the total M, which includes the i/o taper (columns)
    
    # Gets the x values of one row
    x_one = [0]
    for i in range(M-1): # Columns of circuit
        x_one.append( x_one[i] + delta_x_list[i] )
    
    # Multiplies for all other rows
    x_tot = []
    for i in range(M): # Columns of circuit
        x_tot.append( N * [x_one[i]] ) # Multiplied by number of rows N
    return x_tot

def NxN_ap_tups(x_tot, y_tot, st=False):
    """
    Used in NxN unitary matrix circuit of arbitrary design for wgs. Used for
    getting all_pos. Compiles all position information into one matrix of tuples
    (x,y,a), dimensions MxN.
    """
    N = len(x_tot[0])
    M = len(x_tot)

    all_pos    = []
    pos_column = [] # Positions of one column
    for i in range(M):     # Columns of circuit
        for j in range(N): # Rows of circuit
            # Gets pos for one column
            pos_column.append( ( x_tot[i][j], y_tot[i][j], 0 ) )
        all_pos.append( pos_column ) # Adds column
        pos_column = []
    return all_pos

def NxNm_get_all_pos(c_type, N, w, gap_large, gap_small, gap_large_t1,
                     gap_small_t1, delta_x_list, st=False):
    """
    Used in NxN unitary matrix circuit with Mach-Zehnder interferometer design
    for wgs. Returns all_pos, a list of tuples for all of the waveguide i/o's
    in the circuit.
    """
    # 'ap' stands for 'all_pos'
    y_tot = NxNm_ap_y(c_type, N, w, gap_large, gap_small, st)
    x_tot = NxNm_ap_x(N, delta_x_list, st)
    
    all_pos = NxN_ap_tups(x_tot, y_tot, st)
    return all_pos
    
def NxN_ap_add_taper(all_pos, all_lt, st=False):
    """
    Used in NxN unitary matrix circuit of arbitrary design for wgs. Adds i/o
    tapers to all_pos. Used after calling NxNm_get_all_pos() (if using MZI
    design).
    """
    # 'ap' stands for 'all_pos'
    M = len(all_pos)
    N = len(all_pos[0])

    taper_pos1 = [[]]
    taper_pos2 = [[]]
    for j in range(N): # Rows of circuit
        taper_pos1[0] += [(   all_pos[0][j][0]-all_lt[0][j],   all_pos[0][j][1], 0 )]
        taper_pos2[0] += [( all_pos[M-1][j][0]+all_lt[M-1][j], all_pos[M-1][j][1], 0 )]
    all_pos = taper_pos1 + all_pos + taper_pos2
    return all_pos

###############################################################################
###############################################################################
############ NxN Unitary Circuit, Abstract Phase Shifter Functions ############
###############################################################################
###############################################################################

def NxNm_ps_num(c_type, N, init_layer=False, final_layer=False, st=False):
    """
    Used in NxN unitary matrix circuit with Mach-Zehnder interferometer design.
    Returns the number of phase shifters. The final layer adds N phase shifters.
    """
    ps_num = None
    if   c_type.lower() == 'clem16':
        ps_num = N*(N-1)
    elif c_type.lower() == 'bell21clem':
        ps_num = N**2

    if final_layer:
        ps_num += N
    if init_layer:
        ps_num += N
    
    return ps_num

def NxNm_all_ij_list(c_type, N, M, init_layer, final_layer, st=False):
    """
    Used in NxN unitary matrix circuit with Mach-Zehnder interferometer design
    for phase shifters. Returns all_ij_list, a list of (i,j) indices,
    representing where ps appear in the circuit. Called in NxNm_ps_all_ij().
    """
    ### Make i_list and j_list
    # i_list is a list of all columns of all_pos where ps occur
    # j_list is a list of 3 elements: the rows of ps in block_types 0, 1 and 2
    i_list = []
    j_list = [ [], [], [] ]
    if init_layer:
        i_list.append(1)
        j_list[2] = [i for i in range(N)] # j_list[block_type=2]
    
    if   c_type.lower() == 'clem16':
        for k in range(2*N):
            i_list.append(2+4*k)
    elif c_type.lower() == 'bell21clem':
        for k in range(N):
            i_list.append(6+8*k)
            
    # j_list[block_type=0] and j_list[block_type=1]
    if   c_type.lower() == 'clem16':
        for k in range(N-1):
            j_list[ k%2 ].append(k) 
    elif c_type.lower() == 'bell21clem':
        for k in range(N):
            j_list[0].append(k)
            j_list[1].append(k)
    
    if N==2: # Special case, 2x2, single MZI
        if   c_type.lower() == 'clem16':
            i_list = i_list[:-2]
        elif c_type.lower() == 'bell21clem':
            i_list = i_list[:-1]
    if final_layer:
        i_list.append(M-3)
        j_list[2] = [i for i in range(N)] # j_list[block_type=2]
        
    ### Make block_dict
    # Keys are i (in i_list) and values are block_type (0,1,2)
    block_dict = {}
    for i in i_list:
        if i != M-3 and i != 1:
            # Either type 0 (block_type=0) or 1 (block_type=1)
            block_type = int((i-2)/8) % 2
        else:
            # If initial/final layer of ps, block_type=2
            block_type = 2
        block_dict[i] = block_type
    
    ### Make all_ij_list
    all_ij_list = []
    for i in i_list:
        block_type = block_dict[i]
        for j in j_list[block_type]:
            all_ij_list.append((i,j)) # Nice list of all (i,j) pairs
    
    return all_ij_list

def NxN_get_ntype(N, j, j_edge, ignore=False, st=False):
    """
    Used in NxN unitary matrix circuit of arbitrary design for phase shifters.
    Returns the ntype (nearest (edge) type), with values None,0,1, representing
    undefined, closest to the top and closest to the bottom. If j <= j_edge then
    ntype is 0, if j >= N-2-j_edge, then ntype is 1. If ignore=True, ignores
    j_edge and instead gives ntype=0 if j < N/2, ntype=1 if j > N/2, ntype=None
    if j=N/2
    """
    ntype = None
    if not ignore:
        if j <= j_edge:
            ntype = 0
        elif j >= N-2 - j_edge:
            # -2 gives correct values if the ps is on the upper of the modes in all
            # beam splitters (the even modes).
            # -1 gives correct values when counting all modes, odd and even.
            # If the ps is on the lower (odd) modes, then the "j <= j_edge"
            # statement above should be changed, and the other criterion should
            # be "j >= N-j_edge".
            ntype = 1
    else:
        if j < (N-1)/2:
            ntype = 0
        elif j > (N-1)/2:
            ntype = 1
    return ntype

def NxN_all_ntype_list(c_type, N, M, all_ij_list, ps_mid_ntype, j_edge, st=False):
    """
    Used in NxN unitary matrix circuit of arbirary design for phase shifters.
    Gets the list all_ntype_list, used in NxNm_ps_all_ij() (in an MZI design).
    Element k is the ntype of element k in all_ij_list, coordinates (i,j). ntype
    is 0 (near top) or 1 (near bottom).
    """
    all_ntype_list = []
    k = 0
    final_layer, init_layer = False, False
    k_crit = None
    j_it_1, j_ot_1 = None, None
    for i, j in all_ij_list:
        ignore = False
        if c_type.lower() == 'bell21clem':
            ignore = True
        if i == M-3 or i == 1: # Initial and final layers
            ignore = True
            if i == M-3:
                if final_layer == False:
                    k_crit = k # Saves the k when final_layer starts
                ntype_temp1 = NxN_get_ntype(N, j, j_edge, ignore, st)
                if ntype_temp1==1 and j_ot_1 == None:
                    j_ot_1 = j # Saves first ntype=1 in the final layer
                final_layer = True
            if i == 1:
                ntype_temp2 = NxN_get_ntype(N, j, j_edge, ignore, st)
                if ntype_temp2==1 and j_it_1 == None:
                    j_it_1 = j # Saves first ntype=1 in the initial layer
                init_layer = True
        all_ntype_list.append( NxN_get_ntype(N, j, j_edge, ignore, st) )
        k += 1

    n_mid = 0
    for k in range(len(all_ntype_list)):
        # Defines ntype for the ones near the middle (in y)
        # They are undecided up until now (ntype=None). ps_mid_ntype decides
        # what the new ntype becomes
        if all_ntype_list[k] == None:
            all_ntype_list[k] = ps_mid_ntype[n_mid]
            n_mid += 1
    
    # Counts hi and lo lists
    ps_num_hi = 0
    ps_num_lo = 0
    for k in range(len(all_ntype_list)):
        if all_ntype_list[k] == 0:
            ps_num_hi += 1
        else:
            ps_num_lo += 1
    
    # If initial/final layers, balances hi and lo lists so they have about the
    # same amount of ps
    if abs(ps_num_hi - ps_num_lo) > 1 and (final_layer or init_layer):
        pass_final = False
        more_lo = False
        if ps_num_hi > ps_num_lo:
            more_lo = True
        for k in range(N):
            if more_lo != None:
                # Both initial and final layer at the same time
                if init_layer:
                    if j_it_1-1-k >= 0:
                        if     more_lo and all_ntype_list[j_it_1-1-k]==0:
                            all_ntype_list[j_it_1-1-k] = 1
                            ps_num_hi -= 1
                            ps_num_lo += 1
                    if not more_lo and all_ntype_list[k]==1:
                        all_ntype_list[k] = 0
                        ps_num_hi += 1
                        ps_num_lo -= 1
                if ps_num_hi == ps_num_lo:
                    pass_final = True # Skips final_layer
                
                if not pass_final and k_crit != None and final_layer:
                    if j_ot_1-1-k >= 0:
                        if     more_lo and all_ntype_list[k_crit+j_ot_1-1-k]==0:
                            all_ntype_list[k_crit+j_ot_1-1-k] = 1
                            ps_num_hi -= 1
                            ps_num_lo += 1
                    if not more_lo and all_ntype_list[k_crit+k]==1:
                        all_ntype_list[k_crit+k] = 0
                        ps_num_hi += 1
                        ps_num_lo -= 1
            pass_final = False
            if abs(ps_num_hi - ps_num_lo) < 2:
                more_lo = None
    
    return all_ntype_list

def NxN_format_aij_list(all_ij_list, all_ntype_list, flip=False, st=False):
    """
    Used in NxN unitary matrix circuit with arbitrary design for phase shifters.
    Formats all_ij_list so that phase shifters that connect to above the circuit
    (hi ps) and ps that connect to below the circuit (lo ps) are symmetric in
    left-right order.
    """
    all_ij_list_new = []
    ij_hi = []
    ij_lo = []
    k = 0
    i_old = None
    for i,j in all_ij_list:
        if i_old != i:
            if not flip:
                ij_tot = ij_hi + ij_lo[::-1]
            else:
                ij_tot = ij_hi[::-1] + ij_lo
            all_ij_list_new += ij_tot
            ij_hi = []
            ij_lo = []
        if all_ntype_list[k]==0:
                ij_hi.append((i,j))
        else:
                ij_lo.append((i,j))
        i_old = i
        k += 1
    # Final iteration outside of loop:
    if not flip:
        ij_tot = ij_hi + ij_lo[::-1]
    else:
        ij_tot = ij_hi[::-1] + ij_lo
    all_ij_list_new += ij_tot

    return all_ij_list_new

def NxNm_ps_all_ij(c_type, N, M, init_layer, final_layer, ps_mid_ntype, j_edge,
                   flip=False, st=False):
    """
    Used in NxN unitary matrix circuit with Mach-Zehnder interferometer design
    for phase shifters. Gets a list of two lists of (i,j) indices, list 1 and 2,
    representing the positions in all_pos where phase shifter active regions
    start. Each (i,j) element represents the 'start' of a phase shifter. List 1
    is for ps that connect to above the circuit, list 2 connects to below the
    circuit. Both lists go from left to right in the circuit. If flip=True, the
    order of j is changed such that the order of the ps connections within each
    all_pos column is flipped.
    """
    # Get all_ij_list
    all_ij_list = NxNm_all_ij_list(c_type, N, M, init_layer, final_layer, st)
    # Get all_ntype_list
    all_ntype_list = NxN_all_ntype_list(c_type, N, M, all_ij_list, ps_mid_ntype,
                                        j_edge, st)
    all_ij_list = NxN_format_aij_list(all_ij_list, all_ntype_list, flip, st)
    
    # Create final hi and lo lists
    all_ij_hi = []
    all_ij_lo = []
    n_hi = 0
    n_lo = 0
    k = 0
    for i, j in all_ij_list:
        if all_ntype_list[k] == 0: # If high
            all_ij_hi.append( (i,j) )
            n_hi += 1
        elif all_ntype_list[k] == 1: # If low
            all_ij_lo.append( (i,j) )
            n_lo += 1
        k += 1
    
    ps_all_ij = [all_ij_hi] + [all_ij_lo] # List of two lists
    return ps_all_ij

def NxN_ps_all_io(ps_all_ij, y_hi, y_lo, x_hi, x_lo, dx_hi, dx_lo,
                  ps_all_ang, st=False):
    """
    Used in NxN unitary matrix circuit of arbirary design for phase shifters.
    Gets ps_all_io, which is a list of two lists, with all i/o positions of
    the phase shifters. List 1 is for ps that connect to above the circuit,
    list 2 connects to below the circuit. Both lists go from left to right
    in the circuit.
    """
    all_io_hi = []
    all_io_lo = []
    k_hi_max = len(ps_all_ij[0])
    k_lo_max = len(ps_all_ij[1])
    for k in range(k_hi_max):  # Connections go up (hi ps)
        ang = ps_all_ang[0][k] # Angle at io
        ps_io = [ [dx_hi*k       + x_hi, y_hi, -ang],  # Input of ps
                  [dx_hi*(k+1/2) + x_hi, y_hi, +ang] ] # Output of ps
        all_io_hi.append( ps_io )

    for k in range(k_lo_max):  # Connections go down (lo ps)
        ang = ps_all_ang[1][k] # Angle at io
        ps_io = [ [dx_lo*k       + x_lo, y_lo, +ang],  # Input of ps
                  [dx_lo*(k+1/2) + x_lo, y_lo, -ang] ] # Output of ps
        all_io_lo.append( ps_io )
    
    ps_all_io = [all_io_hi] + [all_io_lo] # List of two lists
    return ps_all_io

def NxN_ps_ieo_get_ioco(is_eo, all_w, all_pos, a_io, a_gap, a_w1, a_w2, a_os, a_l,
                        i, j, st=False):
    """
    Used in NxN unitary matrix circuit of arbitrary design for phase shifters.
    For one ps, get the points at the input, output, start of active region
    and end of active region in the main loop of NxN_ps_ioe_pts(). For electro-
    optic ps, also has the start and end of the secondary active region. This
    refers to the conductor itself, which has two parts in an eo ps: one is
    above the wg and one is below. If the above active region is region 1,
    the region below is the secondary active region. Both are needed for the
    electric field.
    """
    #### When a variable ends with '_', the variable is for one ps, ###########
    #### as opposed to all ps.
    pos_i_ = a_io[0] # Absolute position of input
    pos_o_ = a_io[1] # Absolute position of output
    
    gap_, gap2_ = 0, 0 # 'Real' gap, the gap in the code
    if a_gap != None:
        if is_eo:
            # So the chosen gap is the distance between top and bottom active
            # regions.
            gap_  = a_gap/2 + a_w1/2
            gap2_ = a_gap/2 + a_w2/2
        else:
            # Gap from wg
            gap_ = a_gap + (all_w[i][j] + a_w1)/2

    # Starting point of the active region
    pos_ij_   = ncm.add_pos(all_pos[i][j],  (a_os,0,0))
    pos_co_i_ = ncm.add_pos(pos_ij_,   (0,gap_,0))
    pos_co_i_ = ncm.add_pos(pos_co_i_, pos_i_, minus=True)
    
    # Ending point of the active region
    pos_i2_j_ = ncm.add_pos(pos_ij_,   (a_l,0,0)) # a_l ahead of pos_ij_
    pos_co_o_ = ncm.add_pos(pos_i2_j_, (0,gap_,0))
    pos_co_o_ = ncm.add_pos(pos_co_o_, pos_i_, minus=True)
    
    pos_co2_i_, pos_co2_o_ = None, None
    if is_eo:
        # Starting point of the secondary active region
        pos_co2_i_ = ncm.add_pos(pos_ij_,    (0,-gap2_,0))
        pos_co2_i_ = ncm.add_pos(pos_co2_i_, pos_i_, minus=True)
        
        # Ending point of the secondary active region
        pos_co2_o_ = ncm.add_pos(pos_i2_j_,  (0,-gap2_,0))
        pos_co2_o_ = ncm.add_pos(pos_co2_o_, pos_i_, minus=True)
    
    # Final output point
    pos_o_ = ncm.add_pos(pos_o_, pos_i_, minus=True)
    return pos_i_, pos_o_, pos_co_i_, pos_co_o_, pos_co2_i_, pos_co2_o_

def NxN_ps_ioe_add_extra(is_hi, pos_i_ex_, pos_o_ex_inv_, pos_o_, st=False):
    """
    Used in NxN unitary matrix circuit of arbitrary design for phase shifters.
    Returns _pos_i_ex, _pos_o_ex in the main loop of NxN_ps_ioe_pts().
    These are the extra points at the i/o of one ps.
    """
    #### When an external variable ends with '_', the variable is for one ps, #
    #### as opposed to all ps.
    if pos_i_ex_ != None:
        if pos_i_ex_[0][1] > 0:
            for m in range(len(pos_i_ex_)):
                pi = pos_i_ex_[m] # Shorthand
                pos_i_ex_[m] = (pi[0], -pi[1], -pi[2])
        
        for m in range(len(pos_i_ex_)):
            pi = pos_i_ex_[m] # Shorthand
            if is_hi:
                pi = (pi[0],pi[1],pi[2])
            else:    # Inverts in y (lo)
                pi = (pi[0],-pi[1],-pi[2])
            pos_i_ex_[m] = pi
    else:
        pos_i_ex_ = [] # If no extra points
    
    if pos_o_ex_inv_ != None:
        pos_o_ex_inv_ = pos_o_ex_inv_[::-1] # Inverts it
        pos_o_ex_ = []
        for m in range(len(pos_o_ex_inv_)):
            po = pos_o_ex_inv_[m] # Shorthand
            po = ncm.add_pos(po, pos_o_) # Starts from the correct output
            pos_o_ex_.append( po )
    
        for m in range(len(pos_o_ex_)): # Same thing for o
            po = pos_o_ex_[m] # Shorthand
            if is_hi:
                po = (po[0],po[1],-po[2])
            else:    # Inverts in y (lo)
                po = (po[0],-po[1],po[2])
            pos_o_ex_[m] = po
    else:
        pos_o_ex_ = [] # If no extra points
    return pos_i_ex_, pos_o_ex_

def NxN_ps_ioe_pts(is_eo, is_hi, all_w, all_pos, all_ij, all_io, all_gap, all_w1,
                   all_w2, all_os, all_l, aoa_pi_ex, aoa_po_ex, st=False):
    """
    Used in NxN unitary matrix circuit of arbitrary design for phase shifters.
    Gets ps_all_p0, all phase shifter starting positions (one for each ps),
    and ps_aoa_pts, the relative positions of all objects within each phase
    shifter (one for each ps, where each element is a list of positions (x,y,a)).
    """
    ps_all_p0   = []
    ps_aoa_pts1 = [] # For thermo-optic (to) and electro-optic
    ps_aoa_pts2 = [] # Only for electro-optic (eo)
    for k in range(len(all_ij)):
        #### When a variable ends with '_', the variable is for one ps, #####
        #### as opposed to all ps.
        i, j = all_ij[k]
        
        # Gets all standard ps points
        pos_i_, pos_o_, pos_co_i_, pos_co_o_, pos_co2_i_, pos_co2_o_ = NxN_ps_ieo_get_ioco(
            is_eo, all_w, all_pos, all_io[k], all_gap[k], all_w1[k], all_w2[k],
            all_os[k], all_l[k], i, j, st)
        
        # Gets potential extra points
        if aoa_pi_ex != [None] or aoa_po_ex != [None]:
            pos_i_ex_, pos_o_ex_inv_ = None, None
            if aoa_pi_ex != [None]:
                pos_i_ex_ = aoa_pi_ex[k]
            if aoa_po_ex != [None]:
                pos_o_ex_inv_ = aoa_po_ex[k]
            
            pos_i_ex_, pos_o_ex_ = NxN_ps_ioe_add_extra(is_hi, pos_i_ex_,
                                        pos_o_ex_inv_, pos_o_, st)
        else:
            pos_i_ex_, pos_o_ex_ = [], []
        
        # Adds the origin p0 of the ps        
        ps_all_p0.append( pos_i_ )

        # Saves all relative points for one ps
        if is_eo:
            ps_oa_pts1 = [(0,0,0)] + pos_i_ex_ + [pos_co_i_, pos_co_o_]
            ps_aoa_pts1.append ( ps_oa_pts1 )

            ps_oa_pts2 = [pos_co2_i_, pos_co2_o_] + pos_o_ex_ + [pos_o_]
            ps_aoa_pts2.append ( ps_oa_pts2 )
        else:
            ps_oa_pts = [(0,0,0)] + pos_i_ex_ + [pos_co_i_, pos_co_o_] + pos_o_ex_ + [pos_o_]
            ps_aoa_pts1.append ( ps_oa_pts )
    return ps_all_p0, ps_aoa_pts1, ps_aoa_pts2

def NxN_ps_pts(N, ps_type, ps_all_io, ps_aoa_pi_ex, ps_aoa_po_ex, all_pos,
        ps_all_ij, all_w, ps_all_w1, ps_all_w2, ps_all_gap, ps_all_l,
        ps_all_os, y_hi, y_lo, x_hi, x_lo, dx_hi, dx_lo, st=False):
    """
    Used in NxN unitary matrix circuit of arbitrary design for phase shifters.
    Gets ps_all_p0 and ps_aoa_pts by calling NxN_ps_ioe_pts() for the upper
    ps's (hi) and lower ps's (lo).
    """
    
    is_eo = False
    if ps_type.lower() == 'eo': # Either eo or to
        is_eo = True
    
    ps_all_p0   = []
    ps_aoa_pts1 = [] # For electro-optic (eo)
    ps_aoa_pts2 = [] # For electro-optic
    
    for n in [0,1]: # Upper (0) and lower (1) phase shifters
        all_ij  = ps_all_ij[n] # Picks out the correct lists
        all_io  = ps_all_io[n]
        all_gap = ps_all_gap[n]
        all_w1  = ps_all_w1[n]
        all_w2  = ps_all_w2[n]
        all_os  = ps_all_os[n]
        all_l   = ps_all_l[n]
        aoa_pi_ex = ps_aoa_pi_ex[n]
        aoa_po_ex = ps_aoa_po_ex[n]
        
        if n==0:
            is_hi = True
        else:
            is_hi = False
        
        ap0, ap1, ap2 = NxN_ps_ioe_pts(is_eo, is_hi, all_w, all_pos, all_ij, all_io,
            all_gap, all_w1, all_w2, all_os, all_l, aoa_pi_ex, aoa_po_ex, st)

        ps_all_p0.append(ap0)
        ps_aoa_pts1.append(ap1)
        ps_aoa_pts2.append(ap2)
            
    if is_eo:
        return ps_all_p0, ps_aoa_pts1, ps_aoa_pts2
    else:
        return ps_all_p0, ps_aoa_pts1

def NxN_format_ps_aoa(ps_all_xs, ps_aoa_w1, ps_aoa_w2, ps_aoa_r1, ps_aoa_r2,
        ps_all_p0, ps_aoa_pts1, ps_aoa_pts2, ps_aoa_btype1, ps_aoa_btype2,
        all_layers_pad, all_xy_pad, st=False):
    """
    Used in NxN unitary matrix circuit of arbitrary design for phase shifters.
    Before this function, all of the ps lists are lists of two lists. List 1
    is for phase shifters that connect above the circuit, list 2 connects below.
    After this function, all of the lists are a concatenation of lists 1 and 2.
    """
    n_max = len(ps_all_xs)
    k_num_list = []
    for k in ps_all_xs:
        k_num = len(k)
        k_num_list.append(k_num)
    
    ps_all_xs_new = []
    ps_aoa_w1_new = []
    ps_aoa_w2_new = []
    ps_aoa_r1_new = []
    ps_aoa_r2_new = []
    ps_all_p0_new = []
    ps_aoa_pts1_new = []
    ps_aoa_pts2_new = []
    ps_aoa_btype1_new = []
    ps_aoa_btype2_new = []
    all_layers_pad_new = []
    all_xy_pad_new = []
    for n in range(n_max):
        ps_all_xs_new += ps_all_xs[n]
        ps_aoa_w1_new += ps_aoa_w1[n]
        ps_aoa_w2_new += ps_aoa_w2[n]
        ps_aoa_r1_new += ps_aoa_r1[n]
        ps_aoa_r2_new += ps_aoa_r2[n]
        ps_all_p0_new += ps_all_p0[n]
        ps_aoa_pts1_new += ps_aoa_pts1[n]
        ps_aoa_pts2_new += ps_aoa_pts2[n]
        ps_aoa_btype1_new += ps_aoa_btype1[n]
        ps_aoa_btype2_new += ps_aoa_btype2[n]
        all_layers_pad_new += all_layers_pad[n]
        all_xy_pad_new += all_xy_pad[n]
    return ps_all_xs_new, ps_aoa_w1_new, ps_aoa_w2_new, ps_aoa_r1_new, ps_aoa_r2_new, ps_all_p0_new, ps_aoa_pts1_new, ps_aoa_pts2_new, ps_aoa_btype1_new, ps_aoa_btype2_new, all_layers_pad_new, all_xy_pad_new

###############################################################################
###############################################################################
################### NxN Unitary Circuit, Nazca Functions ######################
###############################################################################
###############################################################################

def NxN_make_all_ps(ps_type, ps_all_xs, ps_aoa_w1, ps_aoa_w2, ps_aoa_r1, ps_aoa_r2,
        ps_all_p0, ps_aoa_pts1, ps_aoa_pts2, ps_aoa_btype1, ps_aoa_btype2,
        all_layers_pad=[], all_xy_pad=[], pads=False, st=False):
    """
    Makes all phase shifters in NxN matrix circuit of arbitrary design, given
    correct lists (aoa and others).
    """
    ps_num = len(ps_all_xs)
    elms = []

    for i in range(ps_num):
        # Gets the specific lists for one ps
        ps_xs = ps_all_xs[i]
        ps_oa_w1 = ps_aoa_w1[i]
        ps_oa_r1 = ps_aoa_r1[i]
        ps_oa_btype1 = ps_aoa_btype1[i]
        if ps_type.lower() == 'eo': 
            ps_oa_w2 = ps_aoa_w2[i]
            ps_oa_r2 = ps_aoa_r2[i]
            ps_oa_btype2 = ps_aoa_btype2[i]
        ps_p0 = ps_all_p0[i]
        ps_oa_pts1 = ps_aoa_pts1[i]
        ps_oa_pts2 = []
        if ps_type.lower() == 'eo':
            ps_oa_pts2 = ps_aoa_pts2[i]
        
        layers_pad = None
        xy_pad = None
        if pads:
            if all_layers_pad==[]:
                raise ValueError('Pad layers have not been defined.')
            if all_xy_pad==[]:
                raise ValueError('Pad sizes have not been defined.')
            layers_pad = all_layers_pad[i]
            xy_pad = all_xy_pad[i]
        
        ps_origin = ncm.add_pos(ps_p0, ps_oa_pts1[0])
        if ps_type.lower() == 'to': # Thermo-optic
            elms.append( ncm.make_to_ps_(ps_xs, ps_oa_w1, ps_oa_r1, ps_p0,
                ps_oa_btype1, ps_oa_pts1, layers_pad, xy_pad, pads, flip=False,
                st=st).put(ps_origin) )
        elif ps_type.lower() == 'eo': # Electro-optic
            elms.append( ncm.make_eo_ps_(ps_xs, ps_oa_w1, ps_oa_w2, ps_oa_r1,
                ps_oa_r2, ps_p0, ps_oa_btype1, ps_oa_btype2, ps_oa_pts1,
                ps_oa_pts2, layers_pad, xy_pad, pads, flip=False,
                st=st).put(ps_origin) )
    return elms

def NxN_circuit(all_pos, cell_type, all_w, all_r, all_btype, ps_type, ps_all_xs,
        ps_aoa_w1, ps_aoa_w2, ps_aoa_r1, ps_aoa_r2, ps_all_p0, ps_aoa_pts1,
        ps_aoa_pts2, ps_aoa_btype1, ps_aoa_btype2, all_layers_pad=[], all_xy_pad=[],
        pads=False, st=False):
    """
    Makes entire NxN matrix circuit of arbitrary design, given all_pos, wg data
    and ps data (aoa and other lists).
    """
    # Format the ps lists a little nicer
    ps_all_xs, ps_aoa_w1, ps_aoa_w2, ps_aoa_r1, ps_aoa_r2, ps_all_p0, ps_aoa_pts1, ps_aoa_pts2, ps_aoa_btype1, ps_aoa_btype2, all_layers_pad, all_xy_pad = NxN_format_ps_aoa(ps_all_xs, ps_aoa_w1, ps_aoa_w2, ps_aoa_r1, ps_aoa_r2, ps_all_p0, ps_aoa_pts1, ps_aoa_pts2, ps_aoa_btype1, ps_aoa_btype2, all_layers_pad, all_xy_pad, st)
        
    M = len(all_pos) - 1 # Note the minus 1
    N = len(all_pos[0])
    
    cname = str(N) + 'x'+str(N)+' Unitary Matrix Circuit'
    with nd.Cell(name=cname+ncm.get_n_str()) as C_cell:
        ## Creates all wg and taper cells
        wg_cells = []
        for i in range(M):     # Columns of circuit
            for j in range(N): # Rows of circuit
                if cell_type[i][j].lower() == '':
                    # Makes bent or strt wg
                    wg_cells.append( ncm.bent_wg(
                        w     = all_w[i][j],
                        x_tot = all_pos[i+1][j][0] - all_pos[i][j][0],
                        y_tot = all_pos[i+1][j][1] - all_pos[i][j][1],
                        r     = all_r[i][j],
                        btype = all_btype[i][j],
                        st    = st
                        ).put( all_pos[i][j] ) )
                    
                elif cell_type[i][j].lower() == 'taper':
                    # Makes taper
                    wg_cells.append( ncm.taper(
                        w_in  = all_w[i][j],
                        w_out = all_w[i+1][j],
                        l     = all_pos[i+1][j][0] - all_pos[i][j][0],
                        st    = st
                        ).put( all_pos[i][j] ) )

        ## Creates all ps cells
        PS_cells = NxN_make_all_ps(ps_type, ps_all_xs, ps_aoa_w1, ps_aoa_w2,
                ps_aoa_r1, ps_aoa_r2, ps_all_p0, ps_aoa_pts1, ps_aoa_pts2,
                ps_aoa_btype1, ps_aoa_btype2, all_layers_pad, all_xy_pad, pads, st)
        
    return C_cell

def NxN_swap_ps_io(ps_all_io, hi_or_lo, index1, index2, io1, io2, st=False):
    """
    Used in NxN unitary matrix circuit of arbitrary design for phase shifters.  
    Given ps_all_io, swaps two i/o positions of the conductors. Swaps
    [index1][io1] and [index2][io2], either of hi (ps_all_io[0]) or lo
    (ps_all_io[1]) phase shifters. hi_or_lo can be 'hi', 'lo' or 'both'. Both
    means that hi and lo both swap, keeping the circuit ps i/o connections
    symmetric. io1 an io2 have values 0 (input) or 1 (output) of the phase
    shifters at indeces index1 and index2.
    """
    hl_list = 0 # Hi/lo list
    if hi_or_lo.lower() == 'hi':
        hl_list = [0]
    elif hi_or_lo.lower() == 'lo':
        hl_list = [1]
    elif hi_or_lo.lower() == 'both':
        hl_list = [0,1]
    
    for hl in hl_list: # hl either 0 or 1
        pos_temp = ps_all_io[hl][index1][io1]
        if io1 != io2: # Flips by 180 deg
            pos_temp = (pos_temp[0], pos_temp[1], pos_temp[2]+180)
        pos_temp2 = ps_all_io[hl][index2][io2]
        if io1 != io2: # Flips by 180 deg
            pos_temp2 = (pos_temp2[0], pos_temp2[1], pos_temp2[2]+180)
        ps_all_io[hl][index1][io1] = pos_temp2
        ps_all_io[hl][index2][io2] = pos_temp
    return ps_all_io

def NxN_delete_ps_list(list0, index, hl, st=False):
    """
    Used in NxN unitary matrix circuit of arbitrary design for phase shifters.
    Deletes element 'index' of a generic phase shifter list.
    """
    list1 = list0[hl][:index]
    list2 = list0[hl][index+1:]
    list0[hl] = list1 + list2
    return list0

def NxN_delete_ps(index, hi_or_lo, ps_type, ps_num, ps_all_ij, ps_num_hi, ps_num_lo, ps_all_io, ps_all_xs,
        ps_all_w1, ps_all_w2, ps_all_ang, ps_all_gap, ps_all_l, ps_all_os,
        all_layers_pad, all_xy_pad, ps_aoa_pi_ex, ps_aoa_po_ex, st=False):
    """
    Used in NxN unitary matrix circuit of arbitrary design for phase shifters.
    Deletes all traces of one phase shifter, given index, hi_or_lo and the
    phase shifter lists and matrices.
    """
    hl_list = 0 # Hi/lo list
    if hi_or_lo.lower() == 'hi':
        hl_list = [0]
    elif hi_or_lo.lower() == 'lo':
        hl_list = [1]
    elif hi_or_lo.lower() == 'both':
        hl_list = [0,1]
    
    lists_of_lists = [ps_all_io, ps_all_xs, ps_all_w1, ps_all_ang, ps_all_gap, ps_all_l, ps_all_os, all_layers_pad, all_xy_pad, ps_aoa_pi_ex, ps_aoa_po_ex]
    if ps_type.lower() == 'eo':
        lists_of_lists = [ps_all_io, ps_all_xs, ps_all_w1, ps_all_w2, ps_all_ang, ps_all_gap, ps_all_l, ps_all_os, all_layers_pad, all_xy_pad, ps_aoa_pi_ex, ps_aoa_po_ex]
    
    for hl in hl_list: # hl either 0 or 1
        ps_num -= 1 # One less
        if   hl==0:
            ps_num_hi -= 1
        elif hl==1:
            ps_num_lo -= 1
        
        ps_all_ij = NxN_delete_ps_list(ps_all_ij, index, hl, st)
        ps_all_io = NxN_delete_ps_list(ps_all_io, index, hl, st)
        ps_all_xs = NxN_delete_ps_list(ps_all_xs, index, hl, st)
        ps_all_w1 = NxN_delete_ps_list(ps_all_w1, index, hl, st)
        if ps_type.lower() == 'eo':
            ps_all_w2 = NxN_delete_ps_list(ps_all_w2, index, hl, st)
        ps_all_ang = NxN_delete_ps_list(ps_all_ang, index, hl, st)
        ps_all_gap = NxN_delete_ps_list(ps_all_gap, index, hl, st)
        ps_all_l = NxN_delete_ps_list(ps_all_l, index, hl, st)
        ps_all_os = NxN_delete_ps_list(ps_all_os, index, hl, st)
        all_layers_pad = NxN_delete_ps_list(all_layers_pad, index, hl, st)
        all_xy_pad = NxN_delete_ps_list(all_xy_pad, index, hl, st)
        ps_aoa_pi_ex = NxN_delete_ps_list(ps_aoa_pi_ex, index, hl, st)
        ps_aoa_po_ex = NxN_delete_ps_list(ps_aoa_po_ex, index, hl, st)
    return ps_num, ps_all_ij, ps_num_hi, ps_num_lo, ps_all_io, ps_all_xs, ps_all_w1, ps_all_w2, ps_all_ang, ps_all_gap, ps_all_l, ps_all_os, all_layers_pad, all_xy_pad, ps_aoa_pi_ex, ps_aoa_po_ex