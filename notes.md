# Thoughts

## Code

Currently working on my Z-Module class. My big to-dos include:

1) Migrate utility functions to their own file
2) Make utility functions which take instances of Z-Modules rather than raw vectors(such as hadamard, etc)
3) Generally think about when to make a function vs a method
4) Implement more basic Z-Module tools (r_stars, transformations, Voronoi)

I also want to make an implementation of NTRUE. Shouldn't be too hard; just matrix math.

## Theory

So what am I working towards here? I need to establish whether there are analogues of the CVP and SVP for Z-modules. I guess that's step one.
I'm worried about degenerate cases; can there be limits of a sequence of lattice pts which are not contained in the lattice?
Additionally, what happens with the matrix math for Z-modules? Seems like transformations are no longer square matrices, so we're gonna have some nastier linear algebra.

Gonna have to start crunching on that ASAP.

I also want to take some real notes on the quasicrystals book. I should do that soon.