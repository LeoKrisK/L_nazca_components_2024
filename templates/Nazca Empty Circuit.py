# -*- coding: utf-8 -*-
"""
"""

import numpy as np
import nazca as nd
import nazca.demofab as demo
from nazca.interconnects import Interconnect
import random

import L_nazca_components_2024 as ncm

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
###############################################################################
###############################################################################




# Here you can write your code




nd.export_gds()




