%SINEGORDON1D  Numerical scheme for one-dimensional sine-Gordon equation
%
% This script implements a numerical method for solving the
% multi-dimensional sine-Gordon equation. We use the spectral element
% method for spatial discretization and the trapezoidal rule for temporal
% discretization.
%
% The solver uses the Newton-Raphson method to deal with nonlinearity.
% After linearization, we obtain the so-called effective stiffness matrix
% and effective residual force vector. They are combinations of consistent
% tangent stiffness corresponding to the solution $U^{(i+1),t+\Delta t}$
% and the standard mass matrix arising from the second-order time
% derivative term.
%
% Sine-Gordon equation accepts exact solutions such as kink-type waves,
% breather solutions. These solutions will be used as sanity check of the
% proposed numerical scheme.

%% Geometry of the problem and time steps
a = -20;    b = 30;         % problem domain is [a, b]
nElements = 100;             % number of elements
nElemNodes = 4;             % number of nodes per element
Mesh = spectralElementMesh1D(a, b, nElements, nElemNodes); % mesh generation
gDoF = length(Mesh.Nodes);  % total number of global DoFs.

tmax = 20;                  % last time of simulation
dt = 1e-2;                  % time step
nTimeSteps = fix(tmax/dt);  % number of time steps

% Boundary conditions: The Neumann boundary conditions on both the left-
% and right-most nodes of problem domain are applied. There
% is no Dirichlet (essential) boundary condition and the solutions at all
% nodes are considered as free DoFs.

%% Computation of iteration-independent variables.
% Number of Gauss points used for integrating the consistent tangent
% stiffness matrix and mass matrix must be equal to number of nodes per
% elements.

nGaussPoints = nElemNodes;
[gPoint, gWeight] = legpts(nGaussPoints);   % Legendre-Gauss quadrature points on interval [-1,1]

% B is the standard B-matrix in FEM. WxJ is the Jacobian determinant
% multiplied by the corresponding weight at all quadrature points. A is the
% values of shape functions at the quadrature points. The ndarray A is used
% to compute the mass matrix.
[B, WxJ, A] = computeFEValues(Mesh, nGaussPoints);
AA = zeros(nElemNodes, nElemNodes, nGaussPoints);
for g = 1:nGaussPoints, AA(:,:,g) = A(:,:,g)' * A(:,:,g); end

%% Initial condition.
initialCondition = 'kink-antikink';
[initialSolution, initialVelocity, initialAcceleration] = getInitialCondition(initialCondition);

u = zeros(gDoF, nTimeSteps+1);      % solution vector
% Interpolate the initial conditions at the nodes.
u(:,1) = initialSolution(Mesh.Nodes)';
v_old = initialVelocity(Mesh.Nodes)';
a_old = initialAcceleration(Mesh.Nodes)';

% Project the initial condition solution onto the finite element space.
% to-do

%% Compute constant matrices
TOL = 1e-8;             % Tolerance for stopping the Newton-Raphson.
R = zeros(gDoF, 1);     % Effective load vector.
kt = zeros(gDoF, gDoF); % Consistent ttangent operator.
T = zeros(gDoF, gDoF);  % Effective tangent stiffness matrix (taking into account mass matrix)
F = zeros(gDoF, 1);     % Internal force vector.

a0 = 4/(dt*dt);         % Algorithmic constants.
a1 = 4/dt;

% Compute the constant MASS matrix and constant STIFFNESS matrix.
M = zeros(gDoF, gDoF);
K = zeros(gDoF, gDoF);
for e = 1:nElements
    eDoF = Mesh.Elements(:,e);     % DoFs associated with the element
    for g = 1:nGaussPoints
        gIndex = nGaussPoints*(e-1) + g;
        M(eDoF, eDoF) = M(eDoF, eDoF) + AA(:,:,g)' * WxJ(gIndex);
        K(eDoF, eDoF) = K(eDoF, eDoF) + B(:,:,gIndex)' * B(:,:,gIndex) * WxJ(gIndex);
    end
end
M = sparse(M); K = sparse(K);

%% Time-stepping integration.
t = zeros(1, nTimeSteps + 1);
maxIter = 100;
% Waitbar for showing the process of solving the equation.
processing_bar = waitbar(0, 'Please wait');
activeDoF = 1:gDoF;
for it = 1:nTimeSteps
    t(it+1) = it * dt;
    u(:,it+1) = u(:,it);        % initialize solution for the current time step
    kt(:) = 0;                  % reset tangent stiffness
    % Assemble global tangent stiffness and internal force.
    for e = 1:nElements
        eDoF = Mesh.Elements(:, e);  % DoFs associated with the element
        for g = 1:nGaussPoints
            gIndex = nGaussPoints*(e - 1) + g;
            u_current = A(:,:,g)*u(eDoF, it+1); % solution at quadrature points.
            % Compute consistent tangent operator.
            kt(eDoF, eDoF) = kt(eDoF, eDoF) + cos(u_current)*AA(:,:,g)*WxJ(gIndex);
        end
    end
    kt = sparse(kt);
    kt = kt + K;            % complete the computation of "kt".
    T = kt + a0 * M;        % Compute effective tangent stiffness matrix.
    for iter = 1:maxIter        
        F(:) = 0;       % reset internal force vector
        % Assemble global tangent stiffness and internal force.
        for e = 1:nElements
            eDoF = Mesh.Elements(:, e);    
            for g = 1:nGaussPoints
                gIndex = nGaussPoints*(e - 1) + g;
                u_current = A(:, :, g) * u(eDoF, it+1); % solution at quadrature points.
                F(eDoF) = F(eDoF) + A(:, :, g)'*sin(u_current)*WxJ(gIndex);
            end
        end
        F = F + K * u(:,it+1);  % Complete the computation of F.
        % Compute effective residual force vector.
        R = -F - M * (a0*(u(:,it+1) - u(:,it)) - a1*v_old - a_old); 
        
        % If Dirichlet boundary conditions are applied, further processing
        % of the effeective matrix and residual force vector is required
        % before computing the incremetal solution. Since only the Neumann
        % boundary conditions are used in this script, we don't do such
        % processing and thus activeDoF consists of the degrees of
        % freedom associated with the solution evaluated at all nodes.
        du = zeros(gDoF, 1);    % Incremental solution.
        du(activeDoF) = T(activeDoF, activeDoF) \ R(activeDoF);
        u(activeDoF, it+1) = u(activeDoF, it+1) + du;
        
        % Convergence test.
        residual = norm(R(activeDoF), 2) / norm(u(:,it+1), 2);
        if residual < TOL
            waitbar(it/nTimeSteps, processing_bar, ...
                ['Process: ' num2str(floor(it/nTimeSteps*100)), '%']);
            % Update solution.
            a_new = a0 * ( u(:,it+1) - u(:,it) ) - a1*v_old - a_old;
            v_old = v_old + 0.5*dt*(a_old + a_new);
            a_old = a_new;
            break
        end
    end
end
delete(processing_bar);
