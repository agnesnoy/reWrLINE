# reWrLINE
**reWrLINE: A reimplementation of WrLINE**

A script to extract
the helical axis
and calculate writhe
from an AMBER trajectory
of DNA.

Based on
[WrLINE (agnesnoy/WrLINE)](https://github.com/agnesnoy/WrLINE)
by Thana Sutthibutpong, Sarah A Harris, and Agnes Noy (2015).

**Please cite**:
Sutthibutpong T, Harris S A, Noy A
2015
"Comparison of Molecular Contours for Measuring Writhe in Atomistic Supercoiled DNA"
*J. Chem. Theory. Comput.*
**11**
2768–2775
[https://doi.org/10.1021/acs.jctc.5b00035](https://doi.org/10.1021/acs.jctc.5b00035)

## Changes from WrLINE

* Fully reimplemented in Python 3
* Now supports linear DNA
  (although note that writhe is not a well defined quantity for open curves)
* A slightly friendlier interface

## Getting started

### Requirements

You will need the following things:

* [Python 3](https://www.python.org/)
* [NumPy](https://numpy.org/)

The following Python modules are necessary,
but are part of the standard library
so should be available unless your Python installation is broken:

  * `sys`
  * `os`

It is strongly recommended that you also install:

* CPPTRAJ
  (part of
  [AmberTools](https://ambermd.org/AmberTools.php))

CPPTRAJ is used
to extract the atoms of interest
and ensure that the trajectory is in the right format
for further analysis.
If reWrLINE cannot call CPPTRAJ,
you will need to do it yourself
by running
the following CPPTRAJ script:

```
  parm <topology_filename>
  trajin <trajectory_filename>
  strip !(@C1') outprefix C1
  trajout <path>/C.mdcrd
```

The test suite looks slightly nicer
if you also install
the
[`termcolor`](https://pypi.org/project/termcolor/)
package,
but this is optional
and is not used by the main software.

### Operating system

reWrLINE has only been tested on Linux systems.
It is expected that
it should behave as expected
on all systems
with a correctly configured
Python 3 environment,
but it may be necessary
on some systems
to perform the CPPTRAJ preprocessing
described above
manually.

If reWrLINE does not work correctly
on your system,
please
[open an issue](https://github.com/georgewatson/reWrLINE/issues)
on GitHub.

### Downloads

The easiest way to obtain
a copy of reWrLINE
is by cloning the
[`georgewatson/reWrLINE`](https://github.com/georgewatson/reWrLINE)
git repository
using one of the below methods:

```sh
git clone git@github.com:georgewatson/reWrLINE.git
```

```sh
git clone https://github.com/georgewatson/reWrLINE.git
```

By regularly checking the origin for new commits,
you can ensure that your copy is always up to date.

Alternatively,
stable releases are available for download
[via GitHub](https://github.com/georgewatson/reWrLINE/releases),
but these may lack the newest features
and must be upgraded manually.

### Tests

It is recommended that you run the test suite
after downloading reWrLINE.
You can do this
by entering the reWrLINE directory
and running the file
`test.py`
with no arguments.

The tests ensure that:

* The software is not corrupted and is able to run on your system
* The results of mathematical operations are correct
* The full output agrees with the original WrLINE for a circular system

Note that this last bullet point
is a test of
consistency with the original implementation,
not necessarily of correctness
(although
of course
there is no reason to suspect that the two are at odds).

If any tests fail,
ensure that:

* You have downloaded the newest version of reWrLINE
* Your system meets all the requirements listed above
* Your Python 3 environment is otherwise functioning normally

If the tests are still failing,
[open an issue](https://github.com/georgewatson/reWrLINE/issues)
on GitHub.

## Running reWrLINE

If you have used WrLINE before,
reWrLINE can be used as a drop-in replacement.
To analyse a linear system,
simply add a final argument with a truthy value
(such as `1`)
if your system is linear.

If you have not used WrLINE before,
please keep reading.

If your permissions are correctly configured,
you should be able to run:

```sh
/path/to/WrLINE.py <path> <topology_filename> <trajectory_filename> <num_bp> <num_steps> <is_linear>
```

If the above command does not work,
try:

```sh
python3 /path/to/WrLINE.py <path> <topology_filename> <trajectory_filename> <num_bp> <num_steps> <is_linear>
```

(If that doesn't work either,
something has broken.
Check the command line arguments
and your Python environment,
and
[open an issue](https://github.com/georgewatson/reWrLINE/issues)
on GitHub
if you get stuck.)

The command line arguments are:

1. `path`:
   The path to the directory in which to store the output
   (for example,
   `WrLINE`;
   to use the present working directory,
   type `.`).
   If this directory does not already exist,
   the software will attempt to create it.
2. `topology_filename`:
   The name of
   (and path to,
   if not in the present working directory)
   the AMBER topology
   (parm)
   file corresponding to the trajectory
3. `trajectory_filename`:
   The name of
   (and path to,
   if not in the present working directory)
   the MD trajectory
   to be analysed
4. `num_bp`:
   The number of DNA base pairs in the system
5. `num_steps`:
   The number of steps in the MD trajectory
6. `is_linear`:
   1 if the system is open,
   0 (or omitted entirely) if it is covalently closed into a loop

The topology and trajectory files
can be in any format supported by CPPTRAJ.
If you have already
run CPPTRAJ manually
as discussed above
and the file `<path>/C1.trj`
therefore already exists,
you can enter invalid filenames
in fields 2 and 3
to avoid running CPPTRAJ again.

Field 6 can be omitted entirely
if your system is circular,
allowing for backwards-compatibility with the original WrLINE.

### Example

Try it yourself by manually running the test analysis.

This is a
short (8 step) simulation of
a 336 bp piece of covalently closed circular DNA.
The topology file is `test.prmtop`
and the trajectory file is `test.mdcrd`.
We would like our output in the `test` subdirectory.

To do this,
navigate to your reWrLINE directory and run:

```sh
./WrLINE.py test test.prmtop test.mdcrd 336 8 0
```

(The final `0` is optional here.)

### Output

Output will be produced
in the directory
passed as the first argument to reWrLINE.

The most useful output files are:

* `C1.<topology_filename>`:
  An AMBER topology file
  containing only the C1 atoms,
  useful for visualising the output trajectories
* `C1.3col`:
  The trajectory of
  the helical axis
  in 3col format
* `C1.xyz`:
  The trajectory of
  the helical axis
  in xyz format
* `tw.ser`:
  2D array of
  twist
  at each base pair
  and time step
* `writhe.ser`:
  Time series
  of writhe
  at each time step

The following output will also be produced:

* `C.3col`:
  The trajectory of
  the midpoints of the DNA base pairs
  in 3col format
* `C.xyz`:
  The trajectory of
  the midpoints of the DNA base pairs
  in xyz format
* `sinreg.ser`:
  2D array of
  bending register angles
  at each base pair
  and time step

-- TRANSFORM FILES TO AMBER FORMAT --

To adapt WrLINE helical axis trajectory to AMBER format, we use vmd, although you can use other programs
>> vmd -e xtctomdcrd.tcl

Trick: leave 1st line of axis trajectory in AMBER format in blank for further analysis in ptraj
