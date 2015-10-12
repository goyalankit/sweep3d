The ASCI SWEEP3D README File
============================

------------------------------------------------------------------------

Table of contents
-----------------

-   [Code Description](#Code%20Description)
-   [Building Issues](#Building)
-   [Files in this distribution](#Files)
-   [Optimization Constraints](#Optimization)
-   [Execution Issues](#Execution)
-   [Timing Issues](#Timing)
-   [Storage Issues](#Storage)
-   [About the data](#About)
-   [Expected Results](#Expected)
-   [Modification Record](#Modification)
-   [Record of Formal Questions and Answers](#Record)

Code Description
----------------

SWEEP3D: 3D Discrete Ordinates Neutron Transport

A. General Description

The benchmark code SWEEP3D represents the heart of a real ASCI
application. It solves a 1-group time-independent discrete ordinates
(Sn) 3D cartesian (XYZ) geometry neutron transport problem. The XYZ
geometry is represented by an IJK logically rectangular grid of cells.
The angular dependence is handled by discrete angles with a spherical
harmonics treatment for the scattering source. The solution involves two
steps: the streaming operator is solved by sweeps for each angle and the
scattering operator is solved iteratively.

The two step solution coded in SWEEP3D is known as an inner iteration. A
realistic Sn code would solve a multi-group problem, which in simple
terms is nothing more than a group-ordered iterative solution on top of
what SWEEP3D does. Groups are solved sequentially since there is strong
coupling between groups due to downscattering. The multi-group
iterations are known as outer iterations. Finally, a realistic ASCI Sn
code would include time-dependence with 1000's of time steps on top of
what a multi-group code does. Work and memory increases substantially
for these real problems. The increased work is predominantly just more
inner iterations. Multi-group problems add another dimension to the flux
moments array while time-dependent problems are even worse since they
require saving the angle-dependent phi array for all grid points, all
angles, and all energy groups at both the old and new times.

An Sn sweep proceeds as follows. For a given angle, each grid cell has 4
equations with 7 unknowns (6 faces plus 1 central); boundary conditions
complete the system of equations. The solution is by a direct ordered
solve known as a sweep. Three known inflows allow the cell center and 3
outflows to be solved. Each cell's solution then provides inflows to 3
adjoining cells (1 each in the I, J, & K directions). This represents a
wavefront evaluation with a recursion dependence in all 3 grid
directions. For XYZ geometries, each octant of angles has a different
sweep direction through the mesh, but all angles in a given octant sweep
the same way.

This version of SWEEP3D is based on blocks of IJK & angles with a
stride-1 line-recursion in the I-direction as the innermost work unit.
The stride-1 I-line recursion is good for microprocessors and the
blocked domains provide parallelism opportunities.

B. Coding

This version of SWEEP3D is written entirely in Fortran77 except that it
requires automatic arrays and a C timer routine is used. All automatic
arrays are in inner\_auto(). The Sn solution is carried out by inner()
and sweep(): the source iterations are handled by inner() and the sweeps
are handled in sweep(). This solution process is all that is timed.

Explicit parallelism is by domain decomposition and message-passing.
This version of SWEEP3D supports both PVM and MPI message passing
libraries as well as a single processor version. All message library
dependencies are isolated to the msg\_stuff.cpp file which is
preprocessed to activate the appropriate code. Parallelism by
multitasking the loops in inner() and sweep() is also possible. Both
forms of parallelism can be exploited at the same time.

C. Parallelization

The only inherent parallelism is over the discrete angles: each sweep
for a given angle is independent. But possible reflective boundary
conditions restrict this form of parallelism to a single octant of
angles. Typical Sn orders of S4, S6, or S8 would provide only 3, 6, or
10-way parallelism so this doesn't go very far. Even for the case of all
vacuum BCs, there would only be 24, 48, or 80-way parallelism. This also
doesn't work that well across distributed memory systems such as MPPs or
clusters of machines as it requires multiple copies of the full spatial
grid.

SWEEP3D exploits parallelism via a wavefront process. First, a 2D
spatial domain decomposition onto a 2D array of processors in the I- and
J-directions is used. A single wavefront solve on these domains provides
limited parallelism. So to improve parallel efficiency, blocks of work
are pipelined through the domains. SWEEP3D is coded to pipeline blocks
of MK k-planes and MMI angles through this 2D processor array. Through a
specific ordering of octants, an upper and lower octant pair also gets
pipelined. And lastly, the sweeps of the next octant pair start before
the previous wavefront is completed; the octant order required for
reflective BCs limits this overlap to two octant pairs at a time. The
overall combination is sufficient to give good theoretical parallel
utilization. The resulting pipelined wavefronts are depicted below just
as the second wavefront gets started in the lower-right corner and heads
toward the upper-left corner (the first wavefront is still finishing up
its traversal from upper-right to lower-left).


                     Parallelism via Domain Decomposition

                               IJ (top) View



                               NPE_I domains

       -----------------------------------------------------------------

       |               |               |               |               |

       |               |               |               |               |

       |               |               |               |               |

       |               |               |               |               |

       |               |               |               |               |

       -----------------------------------------------------------------

    N  |               |               |               |               |

    P  |  octant 2     |               |               |               |

    E  |  m-block MMO  |               |               |               |

    |  |  k-block KB   |               |               |               |

    J  |               |               |               |               |

       -----------------------------------------------------------------

    d  |               |               |               |               |

    o  |  octant 2     |  octant 2     |               |               |

    m  |  m-block MMO  |  m-block MMO  |               |               |

    a  |  k-block KB-1 |  k-block KB   |               |               |

    i  |               |               |               |               |

    n  -----------------------------------------------------------------

    s  |               |               |               |               |

       |  octant 2     |  octant 2     |  octant 2     |  octant 3     |

       |  m-block MMO  |  m-block MMO  |  m-block MMO  |  m-block 1    |

       |  k-block KB-2 |  k-block KB-1 |  k-block KB   |  k-block 1    |

       |               |               |               |               |

       -----------------------------------------------------------------

The theoretical parallel processor utilization based on this 2D
decomposition with pipelining (assuming MK is an integral factor of KT
with KB=KT/MK and MMI is an integral factor of MM with MMO=MM/MMI) is
given by


                            8*MMO*KB * (NPE_I*NPE_J)

     -----------------------------------------------------------------------

      2*[2*MMO*KB+(NPE_J-1) + 2*MMO*KB+(NPE_I-1)+(NPE_J-1)] * (NPE_I*NPE_J)

Multitasking parallelism can also be exploited within each block of work
but requires a special ordering of that work. A JK-diagonal wavefront
yields independent I-line work units at each (J,K), thus enabling
multitasking work sharing. The amount of parallelism builds up and falls
off based on the length of the JK-diagonal, so parallel processor
utilization isn't perfect. Also, synchronization is required after each
diagonal. To help improve the situation, MMI angles are pipelined on
JK-diagonals to increase the amount of I-lines that can be solved in
parallel. This MMI-pipelined JK-diagonal wavefront is depicted below at
the 4th stage of a 3-deep wavefront that started in the upper right
corner.


              Parallelism via Diagonal Wavefront Multitasking

                             JK (side) View

     

                              JT grid points

         ---------------------------------------------------------

         |       |       |       |       |       |       |       |

         |       |       |       |  mi=1 |  mi=2 |  mi=3 |       |

         |       |       |       |       |       |       |       |

     M   ---------------------------------------------------------

     K   |       |       |       |       |       |       |       |

         |       |       |       |       |  mi=1 |  mi=2 |  mi=3 |

     g   |       |       |       |       |       |       |       |

     r   ---------------------------------------------------------

     i   |       |       |       |       |       |       |       |

     d   |       |       |       |       |       |  mi=1 |  mi=2 |

         |       |       |       |       |       |       |       |

     p   ---------------------------------------------------------

     o   |       |       |       |       |       |       |       |

     i   |       |       |       |       |       |       |  mi=1 |

     n   |       |       |       |       |       |       |       |

     t   ---------------------------------------------------------

     s   |       |       |       |       |       |       |       |

         |       |       |       |       |       |       |       |

         |       |       |       |       |       |       |       |

         ---------------------------------------------------------

The theoretical processor utilization for the MMI-pipelined diagonal
wavefront on NCPU processors (assuming MK is an integral factor of KT,
MMI is an integral factor of MM, and NPE\_J is an integral factor of
JT\_G with JT=JT\_G/NPE\_J) is given by


                                   MMI*JT*MK

    ----------------------------------------------------------------------

    SUM [

      (SUM [max(min(i-mi+1,JT,MK,JT+MK-i+mi-1),0);mi=1,mmi] + NCPU-1)/NCPU

        ;i=1,JT+MK-1+MMI-1 ] * NCPU

Multitasking the rest of sweep() is also possible. Multitasking on the
I, J, or K loops in loop nests is likely better than on the MI loop
because MMI is so small. Also, both face and leakage are reductions:
face on angles and leakage on angles and spatial boundary faces. The
MMI-pipelined JK-diagonals are such that these reductions are not a
problem inside of the main diagonal loop. But care must be taken when
multitasking other loops involving these reductions. Fusing loops where
possible might help with multitasking efficiency. Multitasking the
source() and flux\_err() routines should be easy; they can even be
rewritten in a single long 1D IJK-loop form. This form of multitasking
has been proven possible on a CRI YMP using their autotasking support.

As one might expect, there are several factors and tradeoffs in parallel
efficiencies. The most important is that the explicit domain
decomposition and multitasking efficiencies are inversely correlated in
MK and MMI. The best domain decomposition efficiency is when MK=1 and
MMI=1 and the best multitasking efficiency is when MK=KT and MMI=MM.
Another tradeoff is that adding another row or column of processors will
reduce efficiency in both cases. Non-integral blocking of k-planes can
lead to load imbalances. Message passing latencies are best hidden when
fewer larger messages are used, but this implies larger blocks via
larger MK & MMI, which in turn decreases the explicit parallel
efficiency. Multiprocessor synchronization overhead is best hidden with
larger blocks (ie. longer diagonals) as well. The size of a block can
also have cache effects. Selection of "best" values for MK and MMI are
likely system dependent, where the system in this case is the cluster of
N-way SMPs.

Asynchronous (or non-blocking) send/receives can also be used. A simple
scheme would simply overlap the pair of sends and the pair of receives.
One could also overlap communication for the next block with the work on
the current block. But this would require double buffering the Phiib,
Phiij, & Phikb arrays. It is not clear how much this will help since the
neighbor processors are likely to be being computing the needed messages
at the same time the current processor is busy anyway. Communication can
also be overlapped with the source, flux, and DSA current computations
within a block by pulling them out of the diagonal loop and computing
them either before the receives or after the sends. But this would
require promoting several arrays from I-line temporaries to full block
size and adding separate loop nests, and this is likely to slow things
down due to the increased memory traffic.

Building Issues
---------------

Type "make" at a shell prompt and follow the instructions. The vendor is
expected to provide the appropriate system calls in timers.f to return
local processor CPU time and global elapsed (wall clock) time. The
vendor is expected to modify the compiler related macros and PVM and MPI
macros in the makefile.

Three versions of the code can be built: single processor, PVM, & MPI.
This is all handled in msg\_stuff.cpp via the CPP C-preprocessor and
ifdef's. The executables are named sweep3d.single, sweep3d.pvm, &
sweep3d.mpi. The PVM version requires external setup of the PVM
configuration and uses a master/slave spawning model with the
appropriate number of slaves spawned with one call to pvmfspawn() using
the default round-robin scheme. The PVM version also has psend()/precv()
and reductions, but they have not been tested. The MPI version runs
under mpirun and uses the rank=0 member as the master.

Files in this distribution
--------------------------


          README.html             README in HTML format

          README                  README in ASCII text

          driver.f                main

          inner_auto.f            automatic arrays

          inner.f                 inner iteration loop

          sweep.f                 sweeper

          source.f                source moments

          flux_err.f              convergence test

          read_input.f            input

          decomp.f                domain decomposition

          initialize.f            problem setup

          octant.f                sweep directions

          timers.c                timers (Fortran callable)

          msg_stuff.cpp           message passing support library

          msg_stuff.h             message passing support

          makefile                make file

          input*                  input files

          output*                 sample output files

Optimization Constraints
------------------------

The entire code can be compiled with either F77 or F90 compilers, but
the selected compiler and its options must be the same for all source
files except inner\_auto.f as explained below.

This version of the SWEEP3D code requires automatic arrays which is a
typical extension of Fortran77 and a standard feature of Fortran90. All
automatics are isolated to the file inner\_auto.f. Vendors must retain
these automatic arrays, but may compile inner\_auto.f differently from
the rest of the code (eg. different compiler options or F90 for just
this file).

The cross sections, external source, and geometry arrays (Sigt, Sigs,
Srcx, hi, hj, hk, di, dj, & dk) are set to constants in the code to
provide a simple problem setup. In real practice, this won't be the
case, so optimizations that "demote" these arrays to constants or
scalars will not be allowed. Also, the code options ISCT, IBC, JBC, KBC,
& IDSA must remain viable; optimizing the code for the specific settings
will not be allowed.

The MPI message passing version must be used for cases involving domain
decomposition. The PVM version is simply an alternative for testing.

Vendors are permitted to make changes to the supplied Fortran codes to
optimize for data layout and alignment. Wholesale changes of the
parallel algorithms are also permitted as long as the full capabilities
of the code are maintained and a full description is provided.

The constraints for audit trails following code modification as noted in
the RFP apply.

Execution Issues
----------------

Benchmark code SWEEP3D must be run using 64-bit floating-point precision
for all non-integer data variables.

Virtually all of the time is spent in sweep() performing the sweeps for
all of the angles.

Timing Issues
-------------

All timings made by this code are via the timers() subroutine in
timers.f. The source should be modified by the vendor to provide both
CPU and wall (elapsed) times for the particular hardware with acceptable
accuracy and resolution.

Only the Sn solution process is timed via two timing calls in inner.f.
Only the master process does the timing, but all processes are barriered
before each timing call so as to measure the collective time of all
processors. The solution time is also expressed in terms of a "grind"
time, which is the effective time to "update" a spatial & angle
phase-space cell. The T3D PVM version can run the required 150-cubed
problem on 128 PEs with an elapsed grind time of 0.037 microseconds. The
best elapsed grind times we have ever seen are in the 0.008 microsecond
range for a 256-cubed version of the same problem on 512 PEs of a T3D on
a "cousin" version of SWEEP3D which uses SHMEM communications.

Storage Issues
--------------

Almost all of the storage requirements comes from the automatic arrays
in inner\_auto(). The printout gives an estimate of memory usage in
MB/domain based on these arrays and the decomposition and blocking.

The small 50-cubed problem only requires 16MB total. The 150-cubed
problem requires 434MB total. These are small because this is only a
1-group time-independent problem. A 20 group problem would be 2.5GB, and
a 20-group time-dependent problem would be 54GB and would run for a week
or so even at a 0.010 microsecond grind time.

The minumum memory configuration required to run the problems in each
configuration must be reported (OS + buffers + code + data + ...).

About the data
--------------

Two inputs are to be run and their outputs reported: input.50 and
input.150. These two input files correspond to 50-cubed and 150-cubed
spatial grid point problems. Each is expected to be run and the results
reported on

1.  a single processor
2.  on all processors of a single SMP
3.  on all processors in a cluster of SMPs comprising a 30 PEAK GFLOPS
    system

Case 2 is expected to use just domain decomposition & MPI within an SMP,
just multitasking, or both simultaneously. Case 3 is expected to use
just domain decomposition & MPI both within and across SMPs or domain
decomposition & MPI across SMPs and multitasking within SMPs. The
blocking factors MK and MMI can be chosen as desired. Alternative
parallel implementations are also acceptable as long as they rely on MPI
and/or multitasking.

The elapsed wall clock time and minimum memory configuration required to
run the 150-cubed problem will be scored.

The input file must be named "input" and be in the current working
directory. A problem description is printed which reflects the input
parameters. It also gives memory estimates, global messaging counts, and
theoretical domain decomposition and multitasking parallel efficiencies
based on the pipelined domain and diagonal JK-wavefronts described
above. The multitasking efficiency is based on NCPU read in from the
input file.

Only the first line of these files should be changed by the vendor. The
first line of input controls both the explicit domain decomposition and
multitasking parallelism. The blocking factors MK and MMI can be chosen
as desired. The NCPUS value is only used in calculating the estimated
theoretical multitasking efficiency. The rest of the input file controls
problem setup and code options and is not to be changed.

The format of the "input" file is as follows:


    Line #  Variables   Type        Explanation and Value(s)



      1       NPE_I     integer     # PEs in the I-direction

              NPE_J     integer     # PEs in the J-direction

              MK        integer     # K-planes for blocking

              MMI       integer     # angles for blocking (must factor MM)

              NCPU      integer     # CPUs per SMP to calc multitasking effic



      2       IT_G      integer     # grid points in the I-direction

              JT_G      integer     # grid points in the J-direction

              KT        integer     # grid points in the K-direction

              MM        integer     # angles per octant (3 or 6)

              ISCT      integer     Pn scattering order (0 or 1)



      3       DX        real        delta-x for I-direction

              DY        real        delta-y for J-direction

              DZ        real        delta-z for K-direction

              EPSI      real        convergence control

                                      < 0.0   =>  iteration count

                                      > 0.0   =>  epsilon criterion



      4       IBC       integer     BC flag for I-direction

              JBC       integer     BC flag for J-direction

              KBC       integer     BC flag for K-direction

                                       0 = vacuum

                                       1 = reflective



      5       IPRINT    integer     flux print flag

                                       0 = off

                                       1 = on

              IDSA      integer     DSA face-currents

                                       0 = off

                                       1 = on

              IFIXUPS   integer     flux fixup flag

                                       0 = off

                                       1 = on

                                     < 0 = on after -IFIXUPS iterations

Expected Results
----------------

Printed results for the two problem sizes are provided in the
output.50\* and output.150\* files. Vendor's results should match those
provided in these files to at least 10 significant digits for
floating-point quantities (iteration error and balance quantities).
Integer results should agree exactly. As a physical consistency check,
absorption plus leakages should balance the external source (accuracy
improves with the number of iterations).

Modification Record
-------------------

Version 2.2b (12/20/95)

Record of Formal Questions and Answers
--------------------------------------

For information about this page contact:\
 Harvey Wasserman -- <hjw@lanl.gov>

This is LANL code LA-CC-97-5.\

[![](../../bar-la3a.GIF)](http://ext.lanl.gov/)\
[LANL Disclaimers](http://www.lanl.gov/misc/disclaimer.html)\
 [Copyright © UC 1999](http://www.lanl.gov/Misc/copyright.html)
