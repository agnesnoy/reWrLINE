#! /usr/bin/env python3

"""
        __        __    _     ___ _   _ _____
 _ __ __\ \      / / __| |   |_ _| \ | | ____|
| '__/ _ \ \ /\ / / '__| |    | ||  \| |  _|
| | |  __/\ V  V /| |  | |___ | || |\  | |___
|_|  \___| \_/\_/ |_|  |_____|___|_| \_|_____|

reWrLINE: A reimplementation of WrLINE

(c) 2019 George D. Watson, University of York
https://georgewatson.me

Based on WrLINE
by Thana Sutthibutpong, Sarah Harris, and Agnes Noy.
Please cite
Sutthibutpong T, Harris S A and Noy A 2015 J. Chem. Theory Comput. 11 2768-75
https://doi.org/10.1021/acs.jctc.5b00035
"""

import sys
import os
import writhe
import caxislib

print(__doc__)
print("---\n")

name = sys.argv[1]
top = sys.argv[2]
traj = sys.argv[3]
num_bp = int(sys.argv[4])
num_steps = int(sys.argv[5])
try:
    linear = sys.argv[6] not in ('0', 'False')
except IndexError:
    linear = False

# Strip trajectory to get C1' coordinates
os.system(f'mkdir -p {name}')
os.system('\n'.join(['cpptraj <<EOF',
                     f'parm {top}',
                     f'trajin {traj}',
                     "strip !(@C1') outprefix C1",
                     f'trajout {name}/C.mdcrd',
                     'EOF']))

print(f"Processing {name}")
print(f"Treating system as {'linear' if linear else 'circular'}")

print("Reading files & initialising arrays")
strand_a, strand_b, midpoints = caxislib.read(name, num_bp, num_steps,
                                              linear=linear)

print("Calculating first-order helical axis")
helix_axis = caxislib.helix_axis(num_bp, num_steps, midpoints, strand_a,
                                 linear=linear)

print("Calculating twist")
twist = caxislib.full_twist(name, num_bp, num_steps, strand_a, strand_b,
                            helix_axis, linear=linear)

print("Calculating helical axis")
caxis = caxislib.caxis(name, num_bp, num_steps, midpoints, twist,
                       linear=linear)

print("Calculating register angles")
sinreg = caxislib.sinreg(name, num_bp, num_steps, midpoints, caxis)

print("Writing output .xyz and .3col files")
caxislib.make_files(name, num_bp, num_steps, midpoints, caxis)

print("Calculating writhe")
wr = writhe.main(name, num_bp, num_steps, linear)

print(f"Job {name} done!")
