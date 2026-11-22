function [lambda, x_grid, y_grid, Vout] = findevs(epsilon, v, Nx, Ny, num_modes, adjointflag, plotflag, fullflag)
% SOLVE_SPECTRAL_EIGENVALUE_PROBLEM Solves the 2D PDE eigenvalue problem.
%
% Inputs:
%   epsilon      - Diffusion parameter
%   v            - Advection parameter
%   Nx           - Number of points in x-direction
%   Ny           - Number of Fourier modes in y-direction (must be even)
%   num_modes    - Number of eigenvalues to plot/extract
%   adjointflag  - solve adjoint problem
%   plotflag     - Plot num_modes eigensolutions
%   fullflag     - Solve the eigenvalue problem using eig (eigs if false)

% ==========================================
% 1. 1D Operators
% ==========================================

% --- X-direction (Chebyshev physical space) ---
[Dx_cheb, X] = cheb(Nx - 1);
% Map from Chebyshev domain X in [1, -1] to x in [0, 1]
% x = (1 - X)/2 => dx = -0.5 dX => d/dx = -2 d/dX
Dx = -2 * Dx_cheb;
Dxx = Dx^2;
x_grid = (1 - X) / 2; % Physical x grid for plotting

% Cast to sparse to ensure kron products remain memory efficient
sDx = sparse(Dx); sDxx = sparse(Dxx); Ix = speye(Nx);

% --- Y-direction (Cosine mode space for Even Solutions) ---
% Mode indices for cosine basis: cos(0*y), cos(1*y), ..., cos((Ny-1)*y)
k = (0 : Ny - 1)';
Dyy = spdiags(-k.^2, 0, Ny, Ny);

% cos(y) * cos(ky) = 0.5 * cos((k-1)y) + 0.5 * cos((k+1)y)
Cy = spdiags(0.5 * ones(Ny, 2), [-1, 1], Ny, Ny);
% Fix the k=0 boundary condition: cos(y) * 1 = 1 * cos(y)
Cy(2, 1) = 1.0;

Iy = speye(Ny);

% ==========================================
% 2. 2D Tensor Product Operators
% ==========================================
% Vector u is stacked column by column (x varies fastest)

if adjointflag
    L = epsilon * kron(Iy, sDxx) + kron(Dyy, Ix) + v * kron(Cy, sDx);
else
    L = epsilon * kron(Iy, sDxx) + kron(Dyy, Ix) - v * kron(Cy, sDx);
end

% The field equation is L*u = -lambda*u => -L*u = lambda*u. Let A = -L.
A = -L;

% Boundary operator: B_eq = epsilon * u_x - v * cos(y) * u = 0 or u_x = 0.
if adjointflag
    B_op = kron(Iy, sDx);
else
    B_op = epsilon * kron(Iy, sDx) - v * kron(Cy, Ix);
end

% ==========================================
% 3. Extract Submatrices for Schur Complement
% ==========================================
% Indices for boundaries x=0 (index 1) and x=1 (index Nx)
idx_B0 = 1 : Nx : Nx*Ny;       % x = 0
idx_B1 = Nx : Nx : Nx*Ny;      % x = 1
idx_B  = [idx_B0, idx_B1];     % All boundary indices
idx_I  = setdiff(1:Nx*Ny, idx_B); % Interior indices

% Partition the system matrices
A_II = A(idx_I, idx_I);
A_IB = A(idx_I, idx_B);

% Extract boundary condition rows
B_matrix = B_op(idx_B, :);
B_BI = B_matrix(:, idx_I);
B_BB = B_matrix(:, idx_B);

% ================================
% 4. Solve the Eigenvalue Problem 
% ================================
A_Schur = A_II - A_IB * (B_BB \ B_BI);

if fullflag
    % Convert to a full matrix for eig()
    [V_all, E_all] = eig(full(A_Schur));
else
    [V_all, E_all] = eigs(A_Schur,num_modes,0.01);
end

% Extract eigenvalues to a vector
lambda_all = diag(E_all);

% Sort eigenvalues by real part
[~, sort_idx] = sort(real(lambda_all));

% Keep only the requested number of eigenvalues and eigenvectors
lambda = lambda_all(sort_idx(1:num_modes));

if plotflag
    Vout = cell(1,num_modes);
    V_I = V_all(:, sort_idx(1:num_modes));

    % Reconstruct full eigenvectors (including boundaries)
    V_full = zeros(Nx*Ny, num_modes);
    V_full(idx_I, :) = V_I;
    V_full(idx_B, :) = - (B_BB \ (B_BI * V_I));

    % ==========================================
    % 5. Transform and Plot the Eigensolutions
    % ==========================================

    % Determine optimal grid for subplots (closest to a square)
    num_cols = ceil(sqrt(num_modes));
    num_rows = ceil(num_modes / num_cols);

    % Setup physical y-grid matching the spectral resolution (Ny elements)
    y_grid = linspace(0, 2*pi, Ny)';
    C_trans = cos(y_grid * k'); % Size: Ny x Ny
    figure(1),clf

    for i = 1:num_modes

        % Reshape to Nx x Ny (rows = x-direction, cols = y-direction modes) and normalise
        V_hat = reshape(V_full(:,i), Nx, Ny);

        % Transform Fourier modes back to physical space in y
        V_phys = real(V_hat * C_trans.');
        % % Normalize so the maximum amplitude is real and positive
        V_phys = V_phys/l2_norm_eigenfunction(V_hat, Nx, Ny, x_grid, idx_I, idx_B, B_BB, B_BI);
        [~, max_idx] = max(abs(V_phys(:)));
        phase_angle = angle(V_phys(max_idx));
        V_phys = V_phys * exp(-1i * phase_angle);

        Vout{i} = V_phys;

        % Plotting
        subplot(num_rows, num_cols, i);
        % subplot(1,2,2)
        surf(x_grid, y_grid/pi, real(V_phys)','EdgeColor','none','FaceColor','interp');
        % contourf(x_grid, y_grid/2/pi, V_phys', 20, 'LineColor', 'none');
        % colorbar
        view(-45,30)

        % Formatting
        ylabel('\phi/\pi');
        xlabel('x');
        zlabel('p(x,\phi)')
        % axis image
        set(gca, 'FontSize', 18);

        % Clean formatting for eigenvalues (handling purely real vs complex)
        lam = lambda(i);
        if abs(imag(lam)) < 1e-10
            title(sprintf('\\lambda = %.4f', real(lam)));
        else
            title(sprintf('\\lambda = %.4f %+.4fi', real(lam), imag(lam)));
        end
    end
    sgtitle(sprintf('\\epsilon = %.1e, v = %.1g',epsilon,v))
    set(gca, 'FontSize', 18);
end
end