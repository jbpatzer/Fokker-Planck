v = 1; eps = 0.01;
Nphi = 10:5:50; Nx = 5:5:75;
N = length(Nphi); M = length(Nx);
error = zeros(N,M);
format compact
for i = 1:N
    for j = 1:M
        l0 = findevs(eps, v, Nx(j), Nphi(i) , 12, false, true, false);
        error(i,j) = abs(l0(8)-4);
        disp([Nphi(i) Nx(j) log10(error(i,j))])
    end
end
figure(1),clf
imagesc(Nx,Nphi,log10(error)')
set(gca,'YDir','normal')
grid on
ylabel('$N_{\phi}$','interpreter','latex')
xlabel('$N_x$','interpreter','latex')
c = colorbar;
c.Label.String = 'log_{10}(error)';
set(gca, 'FontSize', 18);