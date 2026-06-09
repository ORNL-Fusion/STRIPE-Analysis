function [fs_ij,phi_in,F_nj] = FourierCosineSeries1D_v2(x_i,f_ij,N)
% This function computes the 1D fourier cosine series of a data set
% V2 is more systematic with matrix multiplication and explicit
% on the dimensions of the elements with better naming convention
% INPUT: f_ij means f(x_i,t_j)
% x_i : Column vector of independent variable of the series, interval [0,1]
% f_ij: Matrix of data to approximate with a cosine series
% N: number of terms in the cosine series
% OUTPUT: F_nj means Fn(t_j)
% fs_ij: column vector Fourier cosine series approximation of fdata
% phi_in: Matrix containing the Eigenfunctions
% F_nj: column vector of fourier coefficients

% Define "i" as index for vector "x" and "f"

% Define the number of terms to include in the series:
n = [0:N]'; % Column vector
dx = mean(diff(x_i)); 

% Assemble eigenfunctions
% Apply outer product to produce matrix r
r_in = pi*x_i*n';
phi_in = cos(r_in);
phi_ni = phi_in';

% Define Kronecker delta function
d_nn = 0.5*ones(size(n)); 
d_nn(1) = 1;
M_nn = diag(1./d_nn);

% Compute series coefficients:
% phi*fdata leads to a column vector
F_nj = M_nn*phi_ni*f_ij*dx;

% Assemble Fourier cosine series:
fs_ij = phi_in*F_nj;
end

