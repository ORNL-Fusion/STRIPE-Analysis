function [fs,phi,a] = FourierCosineSeries1D(x,fdata,N)
% This function computes the 1D fourier cosine series of a data set
% INPUT:
% x: Column vector of Independent variable of the series, interval [0,1]
% fdata: Column vector of data to approximate with a cosine series
% N: number of terms in the cosine series
% OUTPUT:
% fs: column vector Fourier cosine series approximation of fdata
% phi: Matrix containing the Eigenfunctions
% a: column vector of fourier coefficients

% Define the number of terms to include in the series:
n = [0:N]'; % Column vector
dx = mean(diff(x)); 

% Assemble eigenfunctions
% Apply outer product to produce matrix r
r = pi*n*x';
phi = cos(r);

% Define Kronecker delta function
d = 0.5*ones(size(n)); 
d(1) = 1;

% Compute series coefficients:
% phi*fdata leads to a column vector
a = phi*fdata*dx./d;

% Assemble Fourier cosine series:
fs = phi'*a;
end

