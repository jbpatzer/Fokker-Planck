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