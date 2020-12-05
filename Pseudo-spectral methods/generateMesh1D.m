function [x, k] = generateMesh1D(a, b, N)
%GENERATEMESH1D  One-dimensional mesh and vector of wave numbers.
%   [x, k] = generateMesh1D(a, b, N) creates a one-dimensional mesh of
%   N evenly distributed grid points within the inverval [a, b] and the
%   corresponding vector of wavenumbers in the Fourier space. It returns
%   the coordinates of nodes in x and the wavenumbers in k as column
%   vectors.
% 
%   This helper function is used in several scripts such as kdv.m, bo.m and
%   so on.
%
%   See also kdv

%       Author: Nguyen Lu Trong Khiem
%       Email:  herokhiem@yahoo.com

L = b - a;
h = L/N;
scale = L/(2*pi);
if mod(N, 2) == 0
    x = (a : h : b-h)';
    k = [0:N/2-1, 0, -N/2+1 : -1]'/scale;
else
    x = (a+h/2 : h : b-h/2)';
    k = [0:(N-1)/2, -(N-1)/2: -1]'/scale;
end