function [ff_IJ,A_MN,Phi_IN,Phi_JM] = FourierCosine_2D(x_I,y_J,f_IJ,N,M)
% FOURIERCOSINE_2D takes a 2D array F_IJ and computes the Fourier cosine
% coefficients. The domain of X_I and Y_J are [0,1]
% INPUTS:
% X_I: Column vector of "x" values with "I" elements, [0,1]
% Y_J: Column vector of "y" values with "J" elements, [0,1]
% F_IJ: 2D matrix of function to Fouerier decompose
% N: Number of terms in the p(x) eigenfunction sum
% M: Number of terms in the p(y) eigenfunction sum
% OUTPUTS:
% A_MN: 2D array of Fourier Cosine Coefficients, MxN elements
% Phi_IN: 2D array of "cos(npix)" eigenfunctions for all "n" 
% Phi_JM: 2D array of "cos(mpiy)" eigenfunctions for all "m" 
% ff_IJ : Approximated version of f_IJ

% J.F. Caneses Marin
% Created 2019_03_10


% Check size of independent variables
% -------------------------------------------------------------------------
I = length(x_I);
J = length(y_J);

dx = mean(diff(x_I));
dy = mean(diff(y_J));

% Check dimensions of inputs:
% -------------------------------------------------------------------------
if ~iscolumn(x_I)
    x_I = x_I';
end
if ~iscolumn(y_J)
    y_J = y_J';
end

if size(f_IJ) ~= [I,J]
   f_IJ = f_IJ';
end

% Define "m" and "n" numbers for the infinite sum
% -------------------------------------------------------------------------
n_N = [0:(N-1)]';
m_M = [0:(M-1)]';


% Define Kronecker delta functions
% -------------------------------------------------------------------------
d_mm = 0.5*ones(size(m_M))'; d_mm(1) = 1;
d_nn = 0.5*ones(size(n_N))'; d_nn(1) = 1;
D_MM = diag(1./d_mm);
D_NN = diag(1./d_nn);

% Define eigenfunction matrices
% -------------------------------------------------------------------------
Phi_IN = cos(pi*x_I*n_N');
Psi_IN = Phi_IN*D_NN;
Phi_JM = cos(pi*y_J*m_M');
Psi_JM = Phi_JM*D_MM;

% Compute the Fourier cosine components A_MN
% -------------------------------------------------------------------------
A_MN = Psi_JM'*f_IJ'*Psi_IN*dx*dy;

% Compute the approximation of f_IJ
% -------------------------------------------------------------------------
ff_IJ = Phi_IN*A_MN'*Phi_JM';
end

