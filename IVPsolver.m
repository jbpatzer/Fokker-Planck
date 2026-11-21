function IVPsolver(v,epsilon,Nx,Ny,t_end,dt,dtout,zmax,moviefile,movietime)

% Domain and resolution
% v, epsilon     - parameters
% Nx = 40;            % Number of Chebyshev points in x
% Ny = 30;            % Number of Fourier cosine modes in y
% t_end = 1.0;        % Final time
% dt = timestep
% dtout = output time step

figure(1),clf
subplot(1,2,2)
findevs(epsilon,v,Nx,Ny,1,false,true,false);
drawnow

if ~isempty(moviefile)
    set(gca,'nextplot','replacechildren');
    mv = VideoWriter(moviefile,'MPEG-4');
    Nout = ceil(t_end/dtout);
    mv.FrameRate = min(60,ceil(Nout/movietime));
    open(mv);
else
    mv = [];
end

% Initial conditions
IC_func = @(x, y) 0.1*exp(-100*((x - 0.5).^2+(y/2/pi-0.5).^2));

% =========================================================================
% 1. SPATIAL DISCRETIZATION (x-direction: Chebyshev)
% =========================================================================
[D_cheb, x_cheb] = cheb(Nx - 1);

% Map from standard domain [-1, 1] to [0, 1]
x = (x_cheb + 1) / 2;
Dx  = 2 * D_cheb;          % 1st derivative in x
Dxx = Dx^2;                % 2nd derivative in x

Ix = speye(Nx);

% =========================================================================
% 2. SPATIAL DISCRETIZATION (y-direction: Fourier Cosine)
% =========================================================================
% Modes: k = 0, 1, ..., Ny-1
k = (0:Ny-1)';
Dyy = spdiags(-k.^2, 0, Ny, Ny);

% Matrix C represents multiplication by cos(y) in the cosine basis.
% cos(y)*cos(ky) = 0.5 * (cos((k+1)y) + cos(|k-1|y))
% Pre-build diagonals vectorized
diags = 0.5 * ones(Ny, 2);
C = spdiags(diags, [-1, 1], Ny, Ny);

% Correct the k=0 boundary case (O(1) assignment is negligible)
C(2, 1) = 1.0;

Iy = speye(Ny);

% =========================================================================
% 3. TENSOR PRODUCT OPERATOR (L)
% =========================================================================
% The PDE is u_t = -v*cos(y)*u_x + u_yy + epsilon*u_xx
% We use kron(A, B) where A acts on y-modes and B acts on x-points.
L = -v * kron(C, Dx) + kron(Dyy, Ix) + epsilon * kron(Iy, Dxx);

% =========================================================================
% 4. BOUNDARY CONDITIONS (x = 0 and x = 1)
% =========================================================================
% Indices for boundaries in the flattened (Nx*Ny, 1) state vector:
% x(1)  = 1  -> rows 1, Nx+1, 2Nx+1, ...
% x(Nx) = 0  -> rows Nx, 2Nx, 3Nx, ...
idx_x1 = 1 : Nx : (Nx*Ny);
idx_x0 = Nx : Nx : (Nx*Ny);

% BC operator: epsilon u_x - v cos(y) u = 0
B_x1 = epsilon * kron(Iy, Dx(1, :))  - v * kron(C, Ix(1, :));
B_x0 = epsilon * kron(Iy, Dx(Nx, :)) - v * kron(C, Ix(Nx, :));

% =========================================================================
% 5. CRANK-NICOLSON MATRICES
% =========================================================================
% (I - dt/2 L) U^{n+1} = (I + dt/2 L) U^n
M_impl = speye(Nx*Ny) - (dt/2) * L;
M_expl = speye(Nx*Ny) + (dt/2) * L;

% Enforce BCs perfectly implicitly on the LHS matrix
% 1. Create a diagonal matrix mask to quickly zero out the boundary rows
mask = ones(Nx*Ny, 1);
mask([idx_x1, idx_x0]) = 0;
P = spdiags(mask, 0, Nx*Ny, Nx*Ny);

% Zero out the boundary rows via fast sparse matrix multiplication
M_impl = P * M_impl;

% 2. Extract (row, col, value) triplets from the BC matrices
[i1, j1, v1] = find(B_x1);
[i0, j0, v0] = find(B_x0);

% 3. Map the local row indices (1 to Ny) to their global indices in the flattened system
I_bc = [idx_x1(i1)'; idx_x0(i0)'];
J_bc = [j1; j0];
V_bc = [v1; v0];

% 4. Build the full BC matrix from triplets and add it into the system
B_full = sparse(I_bc, J_bc, V_bc, Nx*Ny, Nx*Ny);
M_impl = M_impl + B_full;

% =========================================================================
% 6. INITIAL CONDITIONS & TRANSFORMATION
% =========================================================================
% Uniform y-grid for evaluation (represents domain [0, pi])
y = pi * (0:Ny-1)' / (Ny-1);

[X_grid, Y_grid] = ndgrid(x, y);
Xplot = [X_grid X_grid(:,Ny-1:-1:1)];
Yplot = [Y_grid 2*pi-Y_grid(:,Ny-1:-1:1)];
U_phys = IC_func(X_grid, Y_grid); % physical values [Nx by Ny]

% Transform physical values to Fourier Cosine coefficients
% T relates coefficients to physical space: U_phys = U_mat * T'
T = cos(y * k');
U_mat = U_phys / T.';

% Vectorize for time-stepping
U = U_mat(:);

figure(1)
subplot(1,2,1)
U0 = reshape(U, Nx, Ny);
U0 = U0 * T.';
output(Xplot,Yplot,U0,Ny,zmax,0)

if ~isempty(moviefile)
    writeVideo(mv,getframe(gcf));
end

% =========================================================================
% 7. TIME STEPPING (Crank-Nicolson)
% =========================================================================
t = 0; tout = 0;
while t<=t_end
    % Compute explicit right-hand side
    rhs = M_expl * U;

    % Overwrite boundary equations with the homogeneous BC value (0)
    rhs(idx_x1) = 0;
    rhs(idx_x0) = 0;

    % Solve linear system for next time step
    U = M_impl \ rhs;

    t = t+dt; tout = tout+dt;

    if tout>=dtout-1e-12

        % =========================================================================
        % 8. RECONSTRUCT AND VISUALIZE
        % =========================================================================
        U_mat = reshape(U, Nx, Ny);
        U_final_phys = U_mat * T.';
        figure(1)
        subplot(1,2,1)
        output(Xplot,Yplot,U_final_phys,Ny,zmax,t)
        tout = tout-dtout;

        if ~isempty(moviefile)
            writeVideo(mv,getframe(gcf));
        end
    end
end

if ~isempty(moviefile)
    close(mv);
end

%%
function output(X,Y,U,Ny,zmax,t)
surf(X,Y/2/pi,[U U(:,Ny-1:-1:1)], 'EdgeColor', 'interp');
xlim([0 1])
ylim([0 1])
zlim([0 zmax])
xlabel('x'); ylabel('\phi/2\pi'); zlabel('p(x,\phi,t)');
title(sprintf('t = %.2f', t));
set(gca, 'FontSize', 18);
view([-45, 30]);
drawnow
