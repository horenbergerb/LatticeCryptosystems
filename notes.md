# Thoughts

## Code

Currently working on my Z-Module class. My big to-dos include:

1) Make utility functions uniform in what arguments they take
2) Add Gramm matrix, LLL algorithm, BKZ-LLL, and a selection of good bases for Lattices/Z-Modules to Utilities.py
3) Implement more basic Z-Module tools (r_stars, covering radius, rank/span/basefield, isLattice)
4) Implement Lattice class which extends Z-Module
  - Voronoi cells
  - Area of fund. domain
  - Transformations of basis
  - SVP/CVP, apprSVP/apprCVP
  - Hermite's constant, Gaussian expected shortest length

Other goals:
* Create abstract crytographic system class
  - Create complexity analysis methods
  - Implement GGH using abstract class
  - Implement NTRU using abstract class
* Make 2DLattice class which extends Lattice and has visualizations
* Make 2DZModule class which extends ZModule and has visualizations
* Alternately make a visualization factory class or something?

## Theory

So what am I working towards here? I need to establish whether there are analogues of the CVP and SVP for Z-modules. I guess that's step one.
I'm worried about degenerate cases; can there be limits of a sequence of lattice pts which are not contained in the lattice?
Additionally, what happens with the matrix math for Z-modules? Seems like transformations are no longer square matrices, so we're gonna have some nastier linear algebra.

Gonna have to start crunching on that ASAP.

