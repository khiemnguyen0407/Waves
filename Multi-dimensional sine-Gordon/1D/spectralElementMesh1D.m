function Mesh = spectralElementMesh1D(a, b, nElements, nElemNodes, varargin)
%SPECTRALELEMENTMESH1D  Mesh for one interval [a, b] for spectral elements.
%
% mesh = spectralElementMesh1D(a, b, nElements, order, varargin) generates a mesh
% consisting of nElements spectral elements, each of which contains "nElemNodes" nodes.
% The mesh is created for the interval [a,b].
%
% This function allows only one type of element in the mesh. That is, every
% element has the same number nodes and they are positioned using the same
% distribution rule xi = -cos(linspace(0, 1, N)*pi) within one reference
% interval [-1, 1]. The coordinates of nodes in the reference interval will
% be mapped to the physical coordinates by a linear mapping and
% constitutive a grid of nodes.
%
%   Author: Khiem Nguyen  -- Email: herokhiem@yahoo.com


x = linspace(a, b, nElements+1);
xi = -cos(linspace(0, 1, nElemNodes)*pi);
eta = xi(1:end-1);
nNodes = (nElemNodes - 1)*nElements + 1;     % number of nodes
Mesh.Nodes = zeros(1, nNodes);
left  = 0.5*(1-eta).*x(1:end-1)';
right = 0.5*(1+eta).*x(2:end)';
Mesh.Nodes(1:end-1) = reshape((left + right)', 1, []);
% Last node is assigned seperately as it doesn't follow the above rule.
Mesh.Nodes(end) = b;

Mesh.Elements = zeros(nElemNodes, nElements);     % connectivity of nodes within finite elements
Mesh.Elements(1,:) = 1 : (nElemNodes-1) : (nElemNodes-1)*(nElements - 1) + 1;
for i = 2 : nElemNodes
    Mesh.Elements(i,:) = i : (nElemNodes-1) : (nElemNodes-1)*nElements + 1;
end
