% Ord_Herm_Pol_1.m is a routine that constructs the basis functions of 
% complete ordinary and Hermite polynomials of the degrees from one to five
% for the two-dimensional (two-state variables) case as described in 
% "Numerically Stable and Accurate Stochastic Simulation Approaches for 
% Solving Dynamic Economic Models" by Kenneth L. Judd, Lilia Maliar and 
% Serguei Maliar, (2011), Quantitative Economics 2/2, 173–210 (henceforth, 
% JMM, 2011).
%
% This version: July 14, 2011. First version: August 27, 2009.
% -------------------------------------------------------------------------
% Inputs:  "z" is the data points on which the polynomial basis functions  
%          must be constructed; T-by-2; 
%          "D"  is the degree of the polynomial whose basis functions must 
%          be constructed; (can be 1,2,3,4 or 5)
%          "PF" is the polynomial family chosen; 0=Ordinary, 1=Hermite;
%          "zb" is Matrix of means and standard deviations of state
%          variables; it is used to normalize these variables in the 
%          Hermite polynomial; 2-by-2            

% Output:  "basis_fs" is the matrix of basis functions of a complete 
%          polynomial of the given degree 
% -------------------------------------------------------------------------
% Copyright © 2011 by Lilia Maliar and Serguei Maliar. All rights reserved. 
% The code may be used, modified and redistributed under the terms provided 
% in the file "License_Agreement.pdf".
% -------------------------------------------------------------------------

function basis_fs = Ord_Herm_Pol_1(z,D,PF,zb)

% A polynomial is given by the sum of polynomial basis functions, phi(i),
% multiplied by the coefficients; see condition (13) in JMM (2011). By 
% convention, the first basis function is one. 

n_rows = size(z,1);   % Infer the number of rows in the data z on which the 
                      % polynomial basis functions must be constructed

if PF == 1;      % If the polynomial family chosen is Hermite, ...
        zc1 = (z(:,1)-zb(1,1))/zb(2,1); 
                 % Normalize the variable z(:,1); n_rows-by-1
        zc2 = (z(:,2)-zb(1,2))/zb(2,2); 
                 % Normalize the variable z(:,2); n_rows-by-1
        p1 = zc1; p2 = zc1.^2-1; p3 = zc1.^3-3*zc1; p4 = zc1.^4-6*zc1.^2+3; p5 = zc1.^5-10*zc1.^3+15*zc1;
                 % p1,...,p5 are the vectors obtained by evaluating the Hermite 
                 % polynomial basis functions from the second to sixth 
                 % (excluding the first basis function equal to one) in all 
                 % values of zc1; each vector is n_rows-by-1
        q1 = zc2; q2 = zc2.^2-1; q3 = zc2.^3-3*zc2; q4 = zc2.^4-6*zc2.^2+3; q5 = zc2.^5-10*zc2.^3+15*zc2;
                 % q1,...,q5 are the vectors obtained by evaluating the Hermite 
                 % polynomial basis functions from the second to sixth 
                 % (excluding the first basis function equal to one) in all 
                 % values of zc2; each vector is n_rows-by-1
else             % If the polynomial family chosen is ordinary, ...
        zc1 = z(:,1); % No normalization for z(:,1); n_rows-by-1 
        zc2 = z(:,2); % No normalization for z(:,2); n_rows-by-1
        p1 = zc1; p2 = zc1.^2; p3 = zc1.^3; p4 = zc1.^4; p5 = zc1.^5;
                % p1,...,p5 are the vectors obtained by evaluating the ordinary 
                 % polynomial basis functions from the second to sixth 
                 % (excluding the first basis function equal to one) in all 
                 % values of zc1; each vector is n_rows-by-1

        q1 = zc2; q2 = zc2.^2; q3 = zc2.^3; q4 = zc2.^4; q5 = zc2.^5;
                 % q1,...,q5 are the vectors obtained by evaluating the ordinary 
                 % polynomial basis functions from the second to sixth 
                 % (excluding the first basis function equal to one) in all 
                 % values of zc2; each vector is n_rows-by-1

end

% Construct the matrix of the basis functions
%--------------------------------------------
if D == 1;        
    basis_fs = [ones(n_rows,1) p1 q1]; 
           % The matrix of basis functions of the first-degree polynomial
elseif D == 2;
    basis_fs = [ones(n_rows,1) p1 q1 p2 p1.*q1 q2];
           % The matrix of basis functions of the second-degree polynomial
elseif D == 3;
    basis_fs = [ones(n_rows,1) p1 q1 p2 p1.*q1 q2 p3 p2.*q1 p1.*q2 q3];
           % The matrix of basis functions of the third-degree polynomial               
elseif D == 4;
    basis_fs = [ones(n_rows,1) p1 q1 p2 p1.*q1 q2 p3 p2.*q1 p1.*q2 q3 p4 p3.*q1 p2.*q2 p1.*q3 q4];
           % The matrix of basis functions of the fourth-degree polynomial
elseif D == 5;
    basis_fs = [ones(n_rows,1) p1 q1 p2 p1.*q1 q2 p3 p2.*q1 p1.*q2 q3 p4 p3.*q1 p2.*q2 p1.*q3 q4 p5 p4.*q1 p3.*q2 p2.*q3 p1.*q4 q5];
           % The matrix of basis functions of the fifth-degree polynomial
end    
