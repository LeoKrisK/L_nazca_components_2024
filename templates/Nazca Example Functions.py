# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 2024

Example functions from 'L_nazca_components_2024' module, templates. The components
need to be 'put', (.put(x,y,a)), and then the file needs to be generated, for
example using 'nd.export_gds()'.

@author: LeoKrisK
"""

import L_nazca_components_2024 as ncm

xs_ls = {}
xs1        = 'XS1'  # Name
xs_ls[xs1] = ['L1'] # List of the layers in XS1
xs2        = 'XS2'  # Name
xs_ls[xs2] = ['L2'] # List of the layers in XS2

############################# Default Values ##################################
st = False      # 'Show Text' - Print out errors, warnings etc.
xs = xs1        # Default cross-section, which is 'XS1'
w = 1.2         # Default width
btype = 'cobra' # Default bend type for bent objects (a cobra bend)


################################ Functions ####################################

# Plus marker
plus1 = ncm.plus(
    x_wide = 40,
    y_wide = 4,
    layers = xs_ls[xs],
    st = st
    )

# Taper
taper1 = ncm.taper(
    w_in  = 0,
    w_out = 5,
    l     = 100,
    layers = xs_ls[xs],
    st = st
    )

# Straight waveguide
strt_wg1 = ncm.strt_wg(
    w = w,
    l  = 100,
    xs = xs,
    st = st
    )

# Bent waveguide
bent_wg1 = ncm.bent_wg(
    w = w,
    x_tot = 150,
    y_tot = 120,
    r     = 70,
    btype = btype,
    xs = xs,
    st = st
    )

# Pair of general bent waveguides
gen_dbwg1 = ncm.generic_bent_wg_pair(
    w_hi = w,
    r_hi = 20,
    x_tot_b_hi = 25,
    y_tot_b_hi = 50,
    btype_hi = btype,

    w_lo = w,
    r_lo = 20,
    x_tot_b_lo = 25,
    y_tot_b_lo = -50,
    btype_lo = btype,
    
    gap_in  = 20,
    #gap_out = 20, # Optional
    xs = xs,
    st = st
    )
    
# Squircle / racetrack
sqc1 = ncm.squircle(
    w = w,
    x_os = 30,
    y_os = 30,
    r    = 20,
    xs = xs,
    st = st
    )
    
# Rectangle (approx.)
rect1 = ncm.squircle(
    w = w,
    x_os = 64,
    y_os = 64,
    r    = w,
    xs = xs,
    st = st
    )
    
# Circle
circ1 = ncm.circle(
    w = w,
    r  = 35,
    xs = xs,
    st = st
    )

# General beam splitter
gen_bs1 = ncm.generic_beamsplit(
    w_in = w,
    l_in = 50,
    
    w_hi = w,
    l_hi = 50/2,
    r_hi = 20/2,
    x_tot_b_hi = 25,
    y_tot_b_hi = 50/2,
    btype_hi = btype,
    
    w_lo = w,
    l_lo = 50,
    r_lo = 20,
    x_tot_b_lo = 100,
    y_tot_b_lo = -50,
    btype_lo = btype,
    gap  = -10,
    flip = False,
    xs = xs,
    st = st
    )
    
# General beam splitter 2
gen_bs2 = ncm.generic_beamsplit(
    w_in = 0,
    l_in = 0,
    
    w_hi = 10,
    l_hi = 50,
    r_hi = 20,
    x_tot_b_hi = 50,
    y_tot_b_hi = 50,
    btype_hi = 'cobra',
    
    w_lo = 5,
    l_lo = 50,
    r_lo = 20,
    x_tot_b_lo = 50,
    y_tot_b_lo = -50,
    btype_lo = 'bsb',
    gap  = None,
    flip = True,
    xs = xs,
    st = st
    )
    
# Beam splitter, bottom and top symmetric
sym_bs1 = ncm.sym_beamsplit(
    w = w,
    l_in  = 50,
    l_out = 50,
    x_tot_b = 0,
    y_tot_b = 50,
    r = 20,
    gap  = 5,
    flip = True,
    btype = btype,
    xs = xs,
    st = st
    )

# Waveguide with circular bump + potential straight part
bump_wg1 = ncm.bump_wg(
    w = w,
    l_bump = 20,
    l_tot  = 150,
    x_bump_middle = 75,
    bump_height   = 20,
    r_hi = 20,
    r_lo = 20,
    ang  = 60,
    btype_in  = 'bsb',
    btype_out = 'bsb',
    flip = False,
    xs = xs,
    st = st
    )

# Straight racetrack resonator
strt_rt1 = ncm.strt_racetrack_reson(
    gap  = 20,
    w_co = w,
    l_co = 100,
    
    w_tk = w,
    x_os_tk = 50,
    y_os_tk = 20,
    r = 50,
    flip = False,
    xs = xs,
    st = st
    )
    
# Straight ring resonator
strt_ri1 = ncm.strt_ring_reson(
    gap  = 20,
    w_co = w,
    l_co = 50,
    w_ring = w,
    r = 50,
    flip = False,
    xs = xs,
    st = st
    )
    
# Symmetric racetrack resonator
sym_rt1 = ncm.sym_racetrack_reson(
    gap  = 20,
    w_co = w,
    l_co = 50,
    r_co = 18,
        
    w_tk = w,
    x_os_tk = 50,
    y_os_tk = 20,
    r_tk = 50,

    l_tot = 200,
    bump_height = 20,
    ang = 45,
    btype_in  = 'bsb',
    btype_out = 'bsb',
    flip = False,
    xs = xs,
    st = st
    )

# Symmetric ring resonator
sym_ri1 = ncm.sym_ring_reson(
    gap  = 20,
    w_co = w,
    l_co = 0,
    r_co = 18,
        
    w_ring = w,
    r_ring = 50,

    l_tot = 200,
    bump_height = 20,
    ang = 45,
    btype_in  = 'bsb',
    btype_out = 'bsb',
    flip = False,
    xs = xs,
    st = st
    )
    
# Pulley racetrack resonator
pul_rt1 = ncm.pulley_racetrack_reson(
    gap  = 20,
    w_co = w,
    l_co_extra = 0,
    r_co = 18,
    
    w_tk = w,
    x_os_tk = 50,
    y_os_tk = 30,
    r_tk = 50,
    
    l_tot = 250,
    bump_height = -50,
    ang = 60,
    btype_in  = 'bsb',
    btype_out = 'bsb',
    flip = False,
    xs = xs,
    st = st
    )
    
# Pulley ring resonator
pul_ri1 = ncm.pulley_ring_reson(
    gap = 20,
    w_co = w,
    l_co_extra = 0,
    r_co = 18,
    
    w_ring = w,
    r_ring = 50,

    l_tot = 200,
    bump_height = -50,
    ang = 60,
    btype_in  = 'bsb',
    btype_out = 'bsb',
    flip = False,
    xs = xs,
    st = st
    )
    
# Rectangle spiral delay lines
rec_del1 = ncm.rect_spiral_delay(
    w = w,
    gap = 10,
    btype = 'bsb', 
    l_tot = 700,
    x_os = None,
    y_os = None,
    r_lo = None,
    r_b = 4,
    n_spirals = 4,
    
    dx_start = 10,
    dy_start = 10,
    drl_start = 10,
    max_steps = 100,
    l_error = 0.01,
    xy_same = False,
    
    flip = False,
    xs = xs,
    st = st
    )

rec_del2 = ncm.rect_spiral_delay(
    w = w,
    gap = 5,
    btype = 'cobra',
    l_tot = 1000,
    x_os = None,
    y_os = 0,
    r_lo = 12,
    r_b = None,
    n_spirals = 3,
    
    dx_start = 10,
    dy_start = 10,
    drl_start = 10,
    max_steps = 100,
    l_error = 0.01,
    xy_same = False,
    
    flip = True,
    xs = xs,
    st = st
    )

rec_del3 = ncm.rect_spiral_delay(
    w = w,
    gap = 5,
    btype = 'bsb',
    l_tot = None,
    x_os = 10,
    y_os = 10,
    r_lo = 12,
    r_b = 6,
    n_spirals = 6,
    # The optimization variables (see rec_del2) do not need to be defined
    flip = False,
    xs = xs,
    st = st
    )

# Circular spiral delay line
cir_del1 = ncm.circ_spiral_delay(
    w = w,
    gap = 5,
    btype = 'bsb',
    l_tot = 1200,
    r_lo = None,
    r_b = 20,
    n_spirals = 3,
    drl_start = 10,
    max_steps = 100,
    l_error = 0.01,
    
    flip = True,
    xs = xs,
    st = st
    )

# Thermo-optic phase shifter
to_ps1 = ncm.thermo_ps(
    all_w = [w,w,3*w,w,w],
    all_r = [10,10,10,10,10],
    p0 = (0,0,0),
    all_btype = ['bsb','bsb','bsb','bsb','bsb'],
    all_pts = [(0,0,0),(20,-100,-45),(40,-200,0),(100,-200,0),(120,-100,45),(140,0,0)],
    
    xy_pad = [10,10,10,10],
    pads = True,
    flip = False,
    xs     = xs2,
    layers_pad = xs_ls[xs2],
    st = st
    )

to_ps2 = ncm.thermo_ps(
    all_w = [w,2*w,3*w,4*w,5*w,6*w],
    all_r = [10,20,30,40,50,60],
    p0 = (0,0,0),
    all_btype = ['bsb','sbend','cobra','sin','bsb','cobra'],
    all_pts = [(0,0,0),(20,-70,-90),(40,-140,-90),(70,-240,0),(190,-200,0),(240,-20,60),(300,120,0)],
    
    xy_pad = [5,10,10,20],
    pads = True,
    flip = True,
    xs     = xs2,
    layers_pad = xs_ls[xs2],
    st = st
    )

# Electro-optic phase shifter
eo_ps1 = ncm.electro_ps(
    all_w1 = [w,w,3*w],
    all_w2 = [3*w,w,w],
    all_r1 = [10,10,10],
    all_r2 = [10,10,10],
    p0 = (0,0,0),
    all_btype1 = ['bsb','bsb','bsb'],
    all_btype2 = ['bsb','bsb','bsb'],
    all_pts1 = [(0,0,0),(20,-100,-45),(40,-190,0),(100,-190,0)],
    all_pts2 = [(40,-210,0),(100,-210,0),(120,-100,45),(140,0,0)],
    
    xy_pad = [10,10,10,10],
    pads = True,
    flip = False,
    xs     = xs2,
    layers_pad = xs_ls[xs2],
    st = st
    )

