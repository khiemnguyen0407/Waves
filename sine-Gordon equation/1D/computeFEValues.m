function [B, WxJ, varargout] = computeFEValues(Mesh, nGaussPoints)

% Each shape function is characterized by the polynomial coefficients.
elements = Mesh.Elements;      % extract element info. from mesh
nodes    = Mesh.Nodes;
elemDoFs = size(elements, 1);       % number of DoFs per element.
nElements = size(elements, 2);      % number of elements.


[gPoint, gWeight] = legpts(nGaussPoints, 'fastsmall'); % Legendre-Gauss quadrature points.
gPoint = gPoint(:).';
gWeight = gWeight(:).';

N = size(elements, 1);      % number of nodes per element
shape = zeros(N, N);        % shape functions in reference domain [-1,1].
dshape = zeros(N, N-1);     % derivative of shape functions
natuaralDerivatives = zeros(N, nGaussPoints);  % derivatives at quadrature points
xi = -cos(linspace(0, 1, N)*pi);    % nodes for each element
for i = 1:N
    pp = poly( xi( (1:N) ~= i ) );
    shape(i,:) = pp./polyval(pp, xi(i));    % derive shape functions
    dshape(i,:) = polyder(shape(i,:));      % derivative of shape functions
    natuaralDerivatives(i,:) = polyval(dshape(i,:), gPoint);
end

% Although B-matrix at each quadrature point can be stored as row vector,
% this function creates a list of matrices of the size (1 x n). Each matrix
% is can be accessed by the third dimension indexing.
B = zeros(1, elemDoFs, nElements * nGaussPoints);
WxJ = zeros(1, nElements * nGaussPoints);
for e = 1:nElements
    for g = 1:nGaussPoints
        Jac = nodes(elements(:,e)) * natuaralDerivatives(:,g);  % Jacobian matrix
        physicalDerivatives = natuaralDerivatives(:,g)/Jac;
        
        gIndex = nGaussPoints*(e-1)+g;      % indexing into gauss points
        WxJ(gIndex) = gWeight(g)*det(Jac);  % weight x det(J)
        B(:, :, gIndex) = physicalDerivatives';   % assemble components of B-matrix
    end
end

nargout_default = 2;
nargout_option = nargout - nargout_default;
if nargout_option > 0
    A = zeros(1, N, size(gPoint, 2));
    for i = 1:N
%         for g = 1:nGaussPoints
%             A(:,i,g) = polyval(shape(i,:), gPoint(g));
%         end
        func = @(x) polyval(shape(i,:), x);
        A(:,i,:) = arrayfun(func, gPoint);
    end
    varargout{1} = A;
end
