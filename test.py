#! /usr/bin/env python3

# pylint: disable=anomalous-backslash-in-string

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

This is the test suite.
It is recommended that you run this before using this implementation.
Tests ensure consistency with WrLINE, not correctness.
"""  # noqa

import filecmp
import sys
import numpy as np
import caxislib
import writhe

try:
    from termcolor import colored
except ImportError:
    def colored(s, _):
        return s

print(__doc__)
print("---\n")

tolerance = 0.1

name = 'test'
num_bp = 336
num_steps = 8

a = np.array([[330.0, 330.0],
              [50.0,  40.0],
              [0.0,   10.0]])
b = np.array([[330.0, 335.0],
              [45.0,  40.0],
              [20.0,  15.0]])

a1 = a[:, 0]
a2 = a[:, 1]
b1 = b[:, 0]
b2 = b[:, 1]
z = np.array([5.0, 0.0, -2.0])


print(f"Processing {name}")
print("Reading files & initialising arrays")
strand_a, strand_b, midpoints = caxislib.read(name, num_bp, num_steps)

print("Calculating first-order helical axis")
helix_axis = caxislib.helix_axis(num_bp, num_steps, midpoints, strand_a)

print("Calculating twist")
twist = caxislib.full_twist(name, num_bp, num_steps, strand_a, strand_b,
                            helix_axis)

print("Calculating helical axis")
caxis = caxislib.caxis(name, num_bp, num_steps, midpoints, twist)

print("Calculating register angles")
sinreg = caxislib.sinreg(name, num_bp, num_steps, midpoints, caxis)

print("Writing output .xyz and .3col files")
caxislib.make_files(name, num_bp, num_steps, midpoints, caxis)

print("Reading coordinates")
read_coords = writhe.read_3col(name + '/C1.3col', num_bp, num_steps)

print("Calculating writhe")
wr = writhe.writhe(read_coords, 2, len(read_coords[0]))
full_writhe = writhe.main(name, num_bp, num_steps)

# Linear

print("Reading files & initialising arrays as if linear")
linear_strand_a, linear_strand_b, linear_midpoints = caxislib.read(name,
                                                                   num_bp,
                                                                   num_steps)

print("Calculating linear first-order helical axis")
linear_axis = caxislib.helix_axis(num_bp, num_steps, linear_midpoints,
                                  linear_strand_a, linear=True)

print("Calculating linear twist")
linear_twist = caxislib.full_twist(name, num_bp, num_steps, linear_strand_a,
                                   linear_strand_b, linear_axis, linear=True,
                                   write=False)

print("Calculating linear helical axis")
linear_caxis = caxislib.caxis(name, num_bp, num_steps, linear_midpoints,
                              linear_twist, linear=True)

tests = {
    "cross\t": [sum(sum(caxislib.cross(a, b))), -8850],
    "dot\t": [sum(caxislib.dot(a, b)), 223450],
    "norm 1\t": [sum(caxislib.norm(a)), 666.33216833189],
    "norm 2\t": [sum(caxislib.norm(b)), 671.36690822942],
    "rotate_to_x": [np.linalg.norm(caxislib.rotate_to_x(a1)), 1.7320508075689],
    "rotate_to_z": [np.linalg.norm(caxislib.rotate_to_z(a1)), 1.7320508075689],
    "twist 1\t": [caxislib.twist(a1, a2, b1, b2, z), 71.997556736987],
    "twist 2\t": [caxislib.twist(a1, b1, a2, b2, z), -15.06995574963],
    "twist 90deg": [caxislib.twist(np.array([0, 0, 0]),
                                   np.array([0, 1, 0]),
                                   np.array([0, 0, 0]),
                                   np.array([1, 0, 0]),
                                   np.array([0, 0, 1])), -90],
    "read strand_a": [sum(sum(sum(strand_a))), 1057248.34],
    "read strand_b": [sum(sum(sum(strand_b))), 1057277.65],
    "read midpoints": [sum(sum(sum(midpoints))), 1057262.995],
    "helix_axis": [sum(sum(sum(helix_axis))), 1057262.995],
    "helix_axis lin": [sum(sum(linear_axis[:, :, 150])),
                       sum(sum(helix_axis[:, :, 150]))],
    "full_twist": [sum(sum(twist)), 83325.202562320],
    "full_twist lin": [sum(linear_twist[:, 150]), sum(twist[:, 150])],
    "caxis\t": [sum(sum(sum(caxis))), 1057251.5419271],
    "caxis lin": [sum(sum(linear_caxis[:, :, 150])),
                  sum(sum(caxis[:, :, 150]))],
    "sinreg\t": [sum(sum(sinreg)), 46.523100971271],
    "read_3col": [sum(sum(sum(read_coords))), sum(sum(sum(caxis)))],
    "writhe.writhe": [wr, -1.013515045594],
    "writhe.main": [sum(sum(full_writhe)), 28.1093914578840],
    "C.3col\t": [filecmp.cmp(f'{name}/C.3col',
                             f'{name}/C.3col.original'), True],
    "C.xyz\t": [filecmp.cmp(f'{name}/C.xyz',
                            f'{name}/C.xyz.original'), True],
    "C1.3col\t": [filecmp.cmp(f'{name}/C1.3col',
                              f'{name}/C1.3col.original'), True],
    "C1.xyz\t": [filecmp.cmp(f'{name}/C1.xyz',
                             f'{name}/C1.xyz.original'), True],
    "tw.ser\t": [filecmp.cmp(f'{name}/tw.ser',
                             f'{name}/tw.ser.original'), True],
    "sinreg.ser": [filecmp.cmp(f'{name}/sinreg.ser',
                               f'{name}/sinreg.ser.original'), True],
    "writhe.ser": [filecmp.cmp(f'{name}/writhe.ser',
                               f'{name}/writhe.ser.original'), True],
}

pass_text = colored('[PASS]', 'green')
fail_text = colored('[FAIL]', 'red')

failures = []

print()
for test, assertion in tests.items():
    print(test, end="\t")
    try:
        if abs(assertion[1] - assertion[0]) < tolerance:
            print(pass_text)
        else:
            failures.append(test)
            print(fail_text)
            print(f"\tExpected\t{assertion[1]}")
            print(f"\tGot\t\t{assertion[0]}")
    except Exception as e:
        failures.append(test)
        print(fail_text)
        print("\tThe following exception was raised:")
        print(f"\t{e}")

print()
if failures:
    print(f"Total failures:\t{len(failures)}")
    for test in failures:
        print(f"\t{test}")
    sys.exit(len(failures))
else:
    print("All tests passed!")
