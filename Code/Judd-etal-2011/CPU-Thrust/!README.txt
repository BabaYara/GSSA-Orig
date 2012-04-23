The folder 'Code/Judd-etal-2011/CPP' contains a C++ modification of
software for generalized stochastic-simulation algorithm (GSSA)
accompanying the article "Numerically Stable and Accurate Stochastic
Simulation Approaches for Solving Dynamic Economic Models" by Kenneth
L. Judd, Lilia Maliar and Serguei Maliar, published in Quantitative
Economics (2011), 2/2, 173-210. The original software is contained in
the directory 'Code/Judd-etal-2011/Original'. This software is
distributed under the original license provided by Lilia and Serguei
Maliar.
 
This version: 28 March 2012.

The following items are provided: 

I. LICENSE AGREEMENT.

II. Folder "GSSA_N_countries" contains C++ files that solve the N-country model

  1. "Main_GSSA_N.cpp"      computes GSSA solutions to the N-country model
  2. "Accuracy_Test_N.cpp"  computes approximation errors in the optimality 
                            conditions on a given set of points in the state 
                            space, for the N-country model 
  3. "Productivity.cpp"     generates random draws of the productivity shocks  
                            and simulates the corresponding series of the  
                            productivity levels
  4. "Num_Stab_Approx.cpp"  implements the numerically stable LS and LAD 
                            approximation methods
  5. "Ord_Polynomial_N.cpp" constructs the sets of basis functions for ordinary
                            polynomials of the degrees from one to five, for
                            the N-country model
  6. "Monomials_1.cpp"      constructs integration nodes and weights for an N-
                            dimensional monomial (non-product) integration rule 
                            with 2N nodes 
  7. "Monomials_2.cpp"      constructs integration nodes and weights for an N-
                            dimensional monomial (non-product) integration rule 
                            with 2N^2+1 nodes
  8. "GH_Quadrature.cpp"    constructs integration nodes and weights for the 
                            Gauss-Hermite rules with the number of nodes in 
                            each dimension ranging from one to ten                     
  9. "aT20200N10.dat"       contains the series of the productivity levels of 
                            length 20,200 for 10 countries that are used for 
                            computing solutions and for evaluating accuracy 
  10. "ran.h"               generates a normal random variable
  11. "cholesky.cpp"        computes a Cholesky decomposition of a matrix
  12. "qrdcmp.cpp"          computes the QR decomposition of square matrix and
                            provides a function to solve a linear system with
                            the resulting decomposition
  13. "simplx.cpp"          implements the simplex method to solve a linear
                            programming problem in standard form
  14. "svdcmp.cpp"          computes the singular value decomposition of a
                            rectangular matrix and provides a function to
                            solve a linear system with the resulting
                            decomposition 
  15. "auxfuncs.cpp"        various helper functions

  To solve the model, execute "Main_GSSA_N".

Author:  Eric M. Aldrich
Contact: ealdrich@gmail.com

-------------------------------------------------------------------------
Copyright © 2012 by Eric M. Aldrich. All rights reserved. The code may
be used, modified and redistributed under the terms provided in the file
"License_Agreement.txt".
-------------------------------------------------------------------------
