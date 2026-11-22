l = zeros(1,100); Nphi0 = 3;
format compact
Nphi = Nphi0:50; N = length(Nphi);
for i = 1:N
    l0 = findevs(0.001, 0.1, 2*Nphi(i), Nphi(i) , 11, false, true, false);
    l(i) = l0(11);
    disp([Nphi(i) l0(11)-1])
end
figure(1),clf
semilogy(Nphi,abs(l(1:N)-1),'ok')
grid on
xlabel('$N_{\phi}$','interpreter','latex')
ylabel('$|\lambda_{10}-1|$','interpreter','latex')
set(gca, 'FontSize', 18);