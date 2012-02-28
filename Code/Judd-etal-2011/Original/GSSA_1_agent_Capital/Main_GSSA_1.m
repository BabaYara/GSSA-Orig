% MATLAB software that solves the standard one-agent growth model using the 
% generalized stochastic-simulation algorithm (GSSA), as described in the 
% article "Numerically Stable and Accurate Stochastic Simulation Approaches
% for Solving Dynamic Economic Models" by Kenneth L. Judd, Lilia Maliar and 
% Serguei Maliar, (2011), Quantitative Economics 2/2, 173–210 (henceforth, 
% JMM, 2011).
% 
% This version: July 14, 2011. First version: August 27, 2009.
% 
% ------------------------------------------------------------------------
% The software includes the following files: 
% ------------------------------------------------------------------------
% 1. "Main_GSSA_1.m"      computes GSSA solutions (in the form of capital 
%                         policy function) to the standard one-agent model
% 2. "Accuracy_Test_1.m"  computes Euler equation approximation errors on 
%                         a given set of points in the state space, for the
%                         standard one-agent model 
% 3. "Num_Stab_Approx.m"  implements the numerically stable LS and LAD 
%                         approximation methods 
% 4. "Ord_Herm_Pol_1.m"   constructs the sets of basis functions for ordinary
%                         and Hermite polynomials of the degrees from one 
%                         to five, for the one-agent model  
% 5. "GH_Quadrature.m"    constructs integration nodes and weights for 
%                         Gauss-Hermite quadrature integration 
% 6. "epsi10000.mat"      contains the series of the productivity shocks of 
%                         length 10,000 that are used for computing solutions
% 7. "epsi_test.mat"      contains the series of the productivity shocks of 
%                         length 10,200 that are used for evaluating
%                         accuracy 
% -------------------------------------------------------------------------
% Copyright © 2011 by Lilia Maliar and Serguei Maliar. All rights reserved. 
% The code may be used, modified and redistributed under the terms provided 
% in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

clc;
clear all;

% 1. Choose the simulation length 
% ------------------------------- 
T     = 10000;   % Choose the simulation length for the solution procedure,
                 % T<=10,000                    
                 
% To solve models with T>10,000, one needs to simulate new series of the 
% productivity levels by enabling the code in paragraph 6  

% 2. Model's parameters
% ---------------------
gam     = 1;         % Utility-function parameter
alpha   = 0.36;      % Capital share in output
beta    = 0.99;      % Discount factor
delta   = 0.02;      % Depreciation rate 
rho     = 0.95;      % Persistence of the log of the productivity level
sigma   = 0.01;      % Standard deviation of shocks to the log of the 
                     % productivity level

% 3. Steady state of capital
% --------------------------
ks = ( (1-beta+beta*delta) / (alpha*beta) )^(1/(alpha-1) );
                     % Use the Euler equation

% 4. Initial condition
% --------------------
k(1:T+1,1) = ks;  % Initial condition for capital (is equal to steady state)
a(1:T,1) = 1;     % Initial condition for the productivity level (is equal to 
                  % steady state)

% 5. Construct the productivity levels, a, for the solution procedure 
% -------------------------------------------------------------------
%epsi10000 = randn(10000,1);    % Generate a random draw of the productivity 
                                % shocks of length 10,000 
%save epsi10000 epsi10000;      % Save the series of the productivity shocks  
                                % into a file "epsi10000.mat" 
load epsi10000;                 % Load the previously saved series of the 
                                % productivity shocks of length 10,000 
epsi = epsi10000(1:T,1)*sigma;  % Compute the error terms in the process for 
                                % productivity level 
for t = 2:T; 
    a(t,1) = a(t-1,1)^rho*exp(epsi(t,1)); 
                            % Compute the next-period productivity levels
                            % using condition (4) in JMM (2011)
end;

% _________________________________________________________________________
%                               
% Compute a first-degree polynomial solution using the one-node Monte Carlo  
% integration method (this solution will be used as an initial guess for
% the other cases) 
% _________________________________________________________________________
%
tic;                  % Start counting time needed to compute the solution
                           
% 6. The GSSA parameters  
% ---------------------
kdamp     = 0.01;     % Damping parameter for (fixed-point) iteration on 
                      % the coefficients of the capital policy function
dif_GSSA_1d = 1e+10;  % Set the initial difference between the series from
                      % two iterations in the convergence criterion (condition
                      % (10) in JMM, 2011) to a very large number

% 7. Initialize the first-degree capital policy function 
%-------------------------------------------------------                         
bk_1d    	= [0; 0.95; ks*0.05];  	          
                      % Vector of polynomial coefficients; 3-by-1
                      
% 8. Initialize the capital series
% --------------------------------
k_old = ones(T+1,1)+ks;% Initialize the series of next-period capital; this
                       % series are used to check the convergence on the 
                       % subsequent iteration (initially, capital can take 
                       % any value); (T+1)-by-1

% 9. The main iterative cycle of GSSA
% -----------------------------------              
while dif_GSSA_1d > 1e-4*kdamp    % 1e-4*kdamp is a convergence parameter,
                                  % adjusted to the damping parameter; see 
                                  % JMM (2011) for a discussion
    
    % 9.1 Generate time series of capital
    % -----------------------------------
    for t = 1:T
      x(t,:) = [1 k(t,1) a(t,1)]; % The basis functions of the first-degree
                                  % polynomial at time t
      k(t+1,1) = x(t,:)*bk_1d;    % Compute next-period capital using bk_1d
    end;
   
    % Compute time series of consumption, c, and Monte Carlo realizations 
    % of the right side of the Euler equation, y, defined in condition (44) 
    % in JMM (2011)
    % -----------------------------------------------------------------------
    c    = k(1:T,1).^alpha.*a(1:T,1)  + (1-delta)*k(1:T,1)-k(2:T+1,1); 
                                 % T-by-1
    y    =  beta*c(2:T,1).^(-gam)./c(1:T-1,1).^(-gam).*(1-delta+alpha*k(2:T,1).^(alpha-1).*a(2:T,1)).*k(2:T,1);
                                 % (T-1)-by-1
                                 
    % 9.3 Evaluate the percentage (unit-free) difference between the series  
    % from the previous and current iterations
    % ---------------------------------------------------------------------
    dif_GSSA_1d = mean(abs(1-k./k_old))
                   % Compute a unit-free difference between the series from
                   % two iterations; see condition (10) in JMM (2011)
                                 
                                 
    % 9.4 Compute and update the coefficients of the capital policy function
    % ----------------------------------------------------------------------
    bk_hat_1d = inv(x(1:T-1,:)'*x(1:T-1,:))*x(1:T-1,:)'*y(1:T-1,:);
                                 % Compute new coefficients of the capital 
                                 % policy function using the OLS
    bk_1d = kdamp*bk_hat_1d + (1-kdamp)*bk_1d;  
                                 % Update the coefficients of the capital  
                                 % policy function using damping
                                     
    % 9.5 Store the capital series 
    %-----------------------------
    k_old = k;         % The stored capital series will be used for checking 
                       % the convergence on the subsequent iteration
    
end;

% 10. Time needed to compute the initial guess 
% --------------------------------------------
time_GSSA_1d     = toc; 
% _________________________________________________________________________                              
%
% Compute polynomial solutions of the degrees from one to D_max using either
% Monte Carlo or Gauss-Hermite quadrature integration methods
% _________________________________________________________________________

tic;                  % Start counting time needed to compute the solution

% 11. The GSSA parameters  
% -----------------------
kdamp = 0.1;          % Damping parameter for (fixed-point) iteration on 
                      % the coefficients of the capital policy function
dif_GSSA_D = 1e+10;   % Set the initial difference between the series from
                      % two iterations in the convergence criterion (condition
                      % (10) in JMM, 2011) to a very large number
                 
% 12. The matrices of the polynomial coefficients
% -----------------------------------------------                             
D_max  = 5;            % Maximum degree of a polynomial: the program computes
                       % polynomial solutions of the degrees from one to D_max;
                       % (D_max can be from 1 to 5) 
npol = [3 6 10 15 21]; % Number of coefficients in polynomials of the degrees                    
                       % from one to five
BK = zeros(npol(D_max),D_max); % Matrix of polynomial coefficients of the 
                               % capital policy function for the polynomial
                               % solutions of the degrees from one to D_max;
                               % npol(D_max)-by-D_max
                 
% 13. Choose an integration method for computing solutions  
% -------------------------------------------------------- 
IM  = 10;        % 0=a one-node Monte Carlo method (default);
                 % 1,2,..,10=Gauss-Hermite quadrature rules with 1,2,...,10 
                 % nodes, respectively
[n_nodes,epsi_nodes,weight_nodes] = GH_Quadrature(IM,1,sigma^2);
                 % n_nodes is the number of integration nodes, epsi_nodes  
                 % are integration nodes, and weight_nodes are integration 
                 % weights for Gauss-Hermite quadrature integration rule
                 % with IM nodes; for a one-dimensional integral, n_nodes  
                 % coincides with IM
a1 = a(1:T,:).^rho*exp(epsi_nodes');
                 % Compute the next-period productivity levels for each 
                 % integration node using condition (4) in JMM (2011); 
                 % T-by-n_nodes
                 
% 14. Choose a regression specification 
% ------------------------------------ 
RM    = 6;       % Choose a regression method: 
                 % 1=OLS,          2=LS-SVD,   3=LAD-PP,  4=LAD-DP, 
                 % 5=RLS-Tikhonov, 6=RLS-TSVD, 7=RLAD-PP, 8=RLAD-DP
normalize = 1;   % Option of normalizing the data; 0=unnormalized data; 
                 % 1=normalized data                    
penalty = 7;     % Degree of regularization for a regularization methods, 
                 % RM=5,6,7,8 (must be negative, e.g., -7 for RM=5,7,8 
                 % and must be positive, e.g., 7, for RM=6)
PF = 0;          % Choose a polynomial family; 0=Ordinary (default);  
                 % 1=Hermite
zb = [[mean(k(1:T,1)) mean(a(1:T,1))]; [std(k(1:T,1)) std(a(1:T,1))]];
                 % Matrix of means and standard deviations of state
                 % variables, k and a; it is used to normalize these 
                 % variables in Hermite polynomial; 2-by-2            

% 15. Initialize the capital series
% ---------------------------------
k_old = ones(T+1,1)+ks;% Initialize the series of next-period capital; this
                       % series are used to check the convergence on the 
                       % subsequent iteration (initially, this series can 
                       % take any value); (T+1)-by-1

% 16. Compute the polynomial solutions of the degrees from one to five
% --------------------------------------------------------------------
for D = 1:D_max;    
   
    % 16.1 Using the previously computed capital series, compute the initial 
    % guess on the coefficients under the  selected approximation method
    % ----------------------------------------------------------------------
    X = Ord_Herm_Pol_1([k(1:T,1) a(1:T,1)],D,PF,zb);
                 % Construct the matrix of explanatory variables X on the  
                 % series of state variables from the previously computed 
                 % time-series solution; columns of X are given by the 
                 % basis functions of the polynomial of degree D 
    bk_D = Num_Stab_Approx(X(1:T-1,:),y(1:T-1,:),RM,penalty,normalize);
                 % Compute the initial guess on the coefficients using the 
                 % chosen regression method
               
    % 16.2 The main iterative cycle of GSSA
    % ------------------------------------- 
    while dif_GSSA_D > 1e-4/10^D*kdamp;  
                 % Check convergence; the convergence parameter depends on 
                 % the polynomial degree D and the damping parameter kdamp;
                 % see the discussion in JMM (2011)
      
       % 16.2.1 Generate time series of capital
       %---------------------------------------
       for t = 1:T
           X(t,:) = Ord_Herm_Pol_1([k(t,1) a(t,1)],D,PF,zb);
                 % The basis functions of a polynomial of degree D at time t
           k(t+1) = X(t,:)*bk_D;
                 % Compute next-period capital using bk_D
       end;
   
        % 16.2.2 Compute time series of consumption, c
        %---------------------------------------------
        c    = k(1:T,1).^alpha.*a(1:T,1)  + (1-delta)*k(1:T,1)-k(2:T+1,1);
                 % T-by-1
                    
        % 16.2.3 Approximate conditional expectation, y, defined in (8) of
        % JMM (2011), for t=1,...T-1 using the integration method chosen
        %------------------------------------------------------------------
            % 16.2.3.1 The one-node Monte Carlo integration method approximates 
            % y with a realization of the integrand in the next period
            % ------------------------------------------------------------------
            if IM == 0;
                y    =  beta*c(2:T,1).^(-gam)./c(1:T-1,1).^(-gam).*(1-delta+alpha*k(2:T,1).^(alpha-1).*a(2:T,1)).*k(2:T,1);
                % Condition (44) in JMM (2011); (T-1)-by-1
            
            % 16.2.3.2 Deterministic integration methods approximate y as  
            % a weighted average of the values of the integrand in the given 
            % nodes with the given weights 
            % --------------------------------------------------------------
            else
                 k1  = k(2:T+1,:); % For each t, k1 is next period capital 
                                   % (chosen at t); T-by-1 
                 for i = 1:n_nodes
                     k2 = Ord_Herm_Pol_1([k1 a1(:,i)],D,PF,zb)*bk_D;
                     % For each t, k2 is capital of period t+2 (chosen at 
                     % t+1); T-by-1
                     c1(:,i) = k1.^alpha.*a1(1:T,i)  + (1-delta)*k1-k2;
                     % For each t, c1 is next-period consumption; T-by-1
                 end
                 k1_dupl = k1*ones(1,n_nodes);  
                     % Duplicate k1 n_nodes times to create a matrix with 
                     % n_nodes identical columns; T-by-n_nodes
                 c_dupl = c(1:T,:)*ones(1,n_nodes);
                     % Duplicate c n_nodes times to create a matrix with 
                     % n_nodes identical columns; T-by-n_nodes 
                 y  =  beta*c1.^(-gam)./c_dupl.^(-gam).*(1-delta+alpha*k1_dupl.^(alpha-1).*a1).*k1_dupl*weight_nodes;
                     % Condition (8) in JMM (2011); T-by-1
            end
     
         % 16.2.4 Evaluate the percentage (unit-free) difference between  
         % the capital series from the previous and current iterations
         % -------------------------------------------------------------
         D
         dif_GSSA_D = mean(abs(1-k./k_old))
                           % The convergence criterion is adjusted to the 
                           % damping parameter 
                           
         % 16.2.5 Compute and update the coefficients of the capital 
         % policy function 
         % ---------------------------------------------------------
         bk_hat_D = Num_Stab_Approx(X(1:T-1,:),y(1:T-1,:),RM,penalty,normalize);
                            % Compute new coefficients of the capital policy 
                            % function using the chosen approximation method
         bk_D = kdamp*bk_hat_D + (1-kdamp)*bk_D; 
                            % Update the coefficients of the capital policy  
                            % function using damping 
                                        
         % 16.2.6 Store the capital series 
         %--------------------------------
         k_old = k;         % The stored capital series will be used for checking 
                            % the convergence on the subsequent iteration                
             
   end;
   
   % 16.3 The GSSA output for the polynomial solution of degree D
   % ------------------------------------------------------------         
   BK(1:npol(D),D) = bk_D;  % Store the coefficients of the polynomial  
                            % of degree D that approximates capital policy 
                            % function 
   time_GSSA(D) = toc;      % Time needed to compute the polynomial solution 
                            % of degree D 
end

% 17. Accuracy test of the GSSA solutions: errors on a stochastic simulation 
% --------------------------------------------------------------------------

   % 17.1 Specify a set of points on which the accuracy is evaluated
   %----------------------------------------------------------------
   T_test = 10200;           % Choose the simulation length for the test 
                             % on a stochastic simulation, T_test<=10,200
                             
   %epsi_test = randn(10200,1); % Generate a random draw of the productivity 
                             % shocks of length 10,200 
   %save epsi_test epsi;     % Save the series of the productivity shocks  
                             % into a file "epsi_test.mat"                          
  
   load epsi_test;           % Load the previously saved series of the 
                             % productivity shocks of length 10,200 
                             
   epsi_test = epsi_test(1:T_test,1)*sigma;     
                             % Compute the error terms in the process for 
                             % productivity level 
                             
   a_test(1,1) = 1;          % Initial condition for the productivity level
  
   for t = 2:T_test; 
       a_test(t,1) = a_test(t-1,1)^rho*exp(epsi_test(t,1));
                             % Compute the next-period productivity levels
                             % using condition (4) in JMM (2011)
   end;
  
   k_test(1,1) = ks;          % Initial condition for capital (equal to steady state)
 
   % 17.2 Choose an integration method for evaluating accuracy of solutions
   %-----------------------------------------------------------------------
   IM_test = 10;                  % See paragraph 13 for the integration 
                                  % options

   % To implement the test on a stochastic simulation with T_test>10,200, one
   % needs to simulate new series of the productivity levels with larger T_test 
   % by enabling the code in paragraph 6.

   % 17.3 Compute errors on a stochastic simulation for the GSSA polynomial 
   % solution of the degrees from one to D_max
   % -----------------------------------------------------------------------
   for D = 1:D_max
      
      % 17.3.1 Simulate the time series solution under the given capital-
      % policy-function coefficients, BK(:,:,D) with D=1,...,D_max 
      %------------------------------------------------------------------
      for t = 1:T_test
        X_test = Ord_Herm_Pol_1([k_test(t,1) a_test(t,1)],D,PF,zb);
        % The basis functions of a polynomial of degree D at time t
        k_test(t+1,1) = X_test*BK(1:npol(D),D);
        % Compute next-period capital using BK(1:npol(D),D)
      end
      
      % 17.3.2 Errors across 10,200 points on a stochastic simulation
      %--------------------------------------------------------------
      discard = 200;  % Discard the first 200 observations to remove the 
                      % effect of the initial conditions 
      [Errors_mean(D) Errors_max(D) time_test(D)] = Accuracy_Test_1(sigma,rho,beta,gam,alpha,delta,k_test(1:T_test,1),a_test(1:T_test,1),BK(1:npol(D),D),D,IM_test,PF,zb,discard);
      % Errors_mean is the unit-free average absolute Euler equation error 
      % (in log10); 
      % Errors_max is the unit-free maximum absolute Euler equation error 
      % (in log10)
  
   end;

% 18. Display the results for the polynomial solutions of the degrees from 
% one to five  
% ------------------------------------------------------------------------
disp(' '); disp('           GSSA OUTPUT:'); disp(' '); 
disp('RUNNING TIME (in seconds):'); disp('');
disp('a) for computing the solution'); 
disp('      1        2         3         4        5');disp(time_GSSA);
disp('b) for implementing the accuracy test'); 
disp('      1        2         3         4        5'); disp(time_test);
disp('APPROXIMATION ERRORS (log10):'); disp(''); 
disp('a) mean error in the Euler equation'); 
disp('      1        2         3         4        5');disp(Errors_mean)
disp('b) max error in the Euler equation'); 
disp('      1        2         3         4        5');disp(Errors_max)

Table=[Errors_mean' Errors_max' time_GSSA']
