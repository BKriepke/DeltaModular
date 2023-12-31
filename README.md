# Delta-modular

Sage code, C++-code and resulting data from implementing the algorithm 1 described in the paper [On the size of integer programs with bounded non-vanishing subdeterminants](https://arxiv.org/abs/2309.03772). This is joint work with Gohar Kyureghyan and Matthias Schymura.


## Content:

### Code:

* **inequiv_hnfs.sage**:            Computes one square matrix in Hermite Normal Form (HNF) of given size of every equivalence class under the equivalence relation described in the paper (see Definition 4.5).

* **singleExample.sage**:           Computes g(\Delta, 2) for all \Delta in a given range and gives a single matrix with that number of columns.

* **allExamples.sage**:             Finds at least one representative of every equivalence class of matrices with g(\Delta, 2) many columns. 

* **main**:                         C++ implementation of algorithm 1 of the paper. Computes g(\Delta, r) for given values of \Delta and r>=3 and gives one or more matrices with that number of columns. Compiled using g++ 7.5.0 on Ubuntu 18.04.


### Data:

This folder includes lists of generic/nongeneric \Delta-modular matrices with the maximal number of columns. It also includes lists of inequivalent matrices in HNF under Definition 4.5 in the paper.

* **hnfs/[\Delta]_[r].txt**:                                                A list of all inequivalent matrices in HNF with determinant \Delta and of size rxr.

* **[Generic/Nongeneric]/[r]/values.txt**:                                  A list of pairs \Delta, g/h(\Delta, r).

* **[Generic/Nongeneric]/[r]/[SingleExample]/Delta=[\Delta].txt**:          A single matrix which is generic/nongeneric and has g/h(\Delta, r) many columns.

* **Generic/[r]/AllExamples/Delta=[\Delta].txt**:                           A list containing at least one (in general more than one) representative of every equivalence class of generic \Delta-modular matrices with the maximal number of columns.
