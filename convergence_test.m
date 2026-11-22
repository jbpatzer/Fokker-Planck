l = zeros(1,100); Nx0 = 5;
format compact
Nx = Nx0:5:100; N = length(Nx);
for i = 1:N
    l0 = findevs(0.01, 1, Nx(i), 25 , 12, false, true, false);
    l(i) = l0(8);
    disp([Nx(i) l0(8)-4])
end
figure(1),clf
semilogy(Nx,abs(l(1:N)-4),'ok')
grid on
xlabel('$N_x$','interpreter','latex')
ylabel('$|\lambda_8-1|$','interpreter','latex')
set(gca, 'FontSize', 18);