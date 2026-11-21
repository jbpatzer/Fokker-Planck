l = zeros(1,41); Ny0 = 3;
format compact
for Ny = Ny0:50
    l0 = findevs(0.001, 0.1, 2*Ny, Ny , 9, false, true, false);
    l(Ny) = l0(9);
    disp([Ny l0(9)-1])
end
figure(1),clf
semilogy(Ny0:50,abs(l(Ny0:50)-1),'ok')
grid on
xlabel('$N_{\phi}$','interpreter','latex')
ylabel('$|\lambda_8-1|$','interpreter','latex')
set(gca, 'FontSize', 18);