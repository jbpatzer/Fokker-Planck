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
        ylim([0 1])

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

%%
function l2_val = l2_norm_eigenfunction(v_vec, Nx, Ny, x_grid, idx_I, idx_B, B_BB, B_BI)
% L2_NORM_EIGENFUNCTION Computes the L2 norm of an eigensolution over
% [0, 1] x [0, pi) using Parseval's identity for cosine modes in y and
% Clenshaw-Curtis quadrature weights in x.
%
% Inputs:
%   v_vec      - Eigenvector (length Nx*Ny for full vector, or length(idx_I) for interior)
%   Nx, Ny     - Number of x points and cosine modes in y
%   x_grid     - Physical x grid on [0, 1] (corresponding to x = (1 - X)/2)
%   idx_I      - Interior indices (1..Nx*Ny)
%   idx_B      - Boundary indices (1..Nx*Ny)
%   B_BB, B_BI - Boundary submatrices from the Schur complement setup
%
% Output:
%   l2_val     - Scalar L2 norm of the eigensolution over [0, 1] x [0, 2*pi)

% Reconstruct full vector if only interior values were passed
if numel(v_vec) == numel(idx_I)
    v_full = zeros(Nx * Ny, 1);
    v_full(idx_I) = v_vec(:);
    v_full(idx_B) = -(B_BB \ (B_BI * v_vec(:)));
elseif numel(v_vec) == Nx * Ny
    v_full = v_vec(:);
else
    error('v_vec must have length Nx*Ny or length(idx_I).');
end

% Reshape into Nx x Ny modal coefficients (columns correspond to k = 0, 1, ..., Ny-1)
V_hat = reshape(v_full, Nx, Ny);

% 1. Exact L2 inner product over y in [0, 2*pi) for each x:
% \int_0^{2*pi} cos(m*y)*cos(n*y) dy = 2*pi for m=n=0, pi for m=n>=1, and 0 otherwise
wy_modal = ones(Ny, 1) * pi;
wy_modal(1) = 2 * pi;
iy_sq = (abs(V_hat).^2) * wy_modal;

% 2. Clenshaw-Curtis quadrature weights on [0, 1]
wx = clenshaw_curtis_weights_01(Nx, x_grid);

% 3. Combine and take square root for the L2 norm
l2_val = sqrt(max(0, real(wx.' * iy_sq))/2);
end

%%
function w01 = clenshaw_curtis_weights_01(N, x_grid)
% Computes Clenshaw-Curtis weights on [0, 1] for N Chebyshev-Gauss-Lobatto nodes.
if N == 1
    w01 = 1.0;
    return;
end
n = N - 1;
theta = (0:n)' * (pi / n);
w_cheb = zeros(n + 1, 1);

for j = 0:n
    s = 0;
    for m = 1:floor(n / 2)
        if 2 * m == n
            b_m = 1;
        else
            b_m = 2;
        end
        s = s + (b_m / (4 * m^2 - 1)) * cos(2 * m * theta(j + 1));
    end
    if j == 0 || j == n
        c_j = 0.5;
    else
        c_j = 1.0;
    end
    w_cheb(j + 1) = (2 / n) * c_j * (1 - s);
end

x_cheb_01 = (1 - cos(theta)) / 2;
w_01_cheb_order = 0.5 * w_cheb;

w01 = zeros(N, 1);
for idx = 1:N
    [~, loc] = min(abs(x_cheb_01 - x_grid(idx)));
    w01(idx) = w_01_cheb_order(loc);
end
end