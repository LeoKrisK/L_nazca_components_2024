# Nazca Components and Photonic Universal Multiport Interferometer

This package contains (i) a library of useful photonic components and (ii) a straight-forward way to create a simple universal multiport interferometer on a photonic integrated circuit platform. The package relies on Nazca, which is free to download from https://nazca-design.org/.

## (i) Nazca Photonic Components
In (i), each component is a function that returns a Nazca cell (see https://nazca-design.org/manual/getting_started.html). There are two functions associated with each component: a user-friendly function ("A") and a "make" function ("B"). "B" might sometimes be less intuitive than "A", so we recommend using "A" in most cases. "A" formats the input variables and then calls "B". "B" (the "make" functions) all have "make" at the start of their function names. The library contains straight waveguides, four types of bent waveguides, beam splitters, many types of resonator cavities, spiral delay lines and thermo-optic and electro-optic phase shifters. More precisely, the functions in (i) are as follows.

| Component  | User-Friendly Function | "Make" Function |
| ------------- | ------------- | ------------- |
| Straight Taper  | taper()  | make_tp_() |
| Straight Waveguide  | strt_wg()  | make_strt_object_() |
| Bent Waveguide  | bent_wg()  | make_b_object_() |
| Beam Splitter | generic_beamsplit()  | make_bs_() |
| Symmetric Beam Splitter  | sym_beamsplit()  | make_sym_bs_() |
| Straight Racetrack Resonator | strt_racetrack_reson()  | make_rtr_() |
| Straight Ring Resonator | strt_ring_reson()  | make_rtr_() |
| Symmetric Racetrack Resonator | sym_racetrack_reson()  | make_rtr_() |
| Symmetric Ring Resonator | sym_ring_reson()  | make_rtr_() |
| Pulley Racetrack Resonator | pulley_racetrack_reson()  | make_rtr_() |
| Pulley Ring Resonator | pulley_ring_reson()  | make_rtr_() |
| Rectangle Spiral Delay Line  | rect_spiral_delay()  | make_delay_line_() |
| Circular Spiral Delay Line  | circ_spiral_delay()  | make_delay_line_() |
| Thermo-Optic Phase Shifter  | thermo_ps()  | make_to_ps_() |
| Electro-Optic Phase Shifter  | electro_ps()  | make_eo_ps_() |

There are also some miscellaneous functions, either used by the above functions or useful on their own. They are as follows.

| Component  | User-Friendly Function | "Make" Function |
| ------------- | ------------- | ------------- |
| Plus Marker  | plus()  | make_plus_() |
| Bent Waveguide Pair  | generic_bent_wg_pair()  | make_bo_pair_() |
| Ring Resonator Cavity  | circle()  | make_sqrc_() |
| Racetrack Resonator Cavity  | squircle()  | make_sqrc_() |
| Waveguide with Bump | bump_wg()  | make_bump_wg_() |
| Conductor Pad  | pad()  | make_pad_() |

In the "templates" folder, the file "Nazca Example Functions" shows one or more examples of each function. The functions in the code have comments and docstrings as well.

## (ii) Universal Multiport Interferometer
Given any arbitrary NxN unitary matrix U, (ii) creates a multiport interferometer performing the matrix operation. In (ii), the shape of the interferometer only depends on the dimension N, not the specific elements of U. *L_nazca_components_2024/NxN_circuit* contains all the universal multiport interferometer functions. The "templates" folder contains "NxN Unitary Circuit Template", a recommended template for using (ii), with comments. In the template, we use the functions NxNm_get_M, NxNm_get_all_pos, NxN_ap_add_taper, NxNm_ps_num, NxNm_ps_all_ij, NxN_ps_pts and NxN_circuit. The last function creates the circuit using Nazca, while the other functions create the lists / matrices of data neccessary in NxN_circuit, given the input variables in "NxN Unitary Circuit Template". There are also three extra functions that can change the layout and order of phase shifters and their conducting connections. These are NxN_swap_ps_io, which swaps inputs/outputs between two phase shifters, NxN_delete_ps, which deletes a phase shifter, and ps_add_point, which adds a point that a phase shifter must go through. These are used to avoid overlap between the conductors of the phase shifters.

There are currently two interferometer designs: "clem16" and "bell21clem", from "Optimal design for universal multiport interferometers" by Clements et al., contra "Further Compactifying Linear Optical Unitaries" by Bell and Walmsley. The building block of the first is an MZI containing two 50:50 beam splitters, one phase shifter between them and one preceding both beam splitters. The second improves said design such that phase shifters only appear between the 50:50 beam splitters, making the circuit more compact. The program assumes that the phase shifters have input and output conductors at regular intervals above and below the circuit.

## Getting Started
The "tutorial" folder contains some thorough tutorials and things to keep in mind. Files in the "templates" folder can be copied and used directly. Nazca is required for the package, so download it first at https://nazca-design.org/.

