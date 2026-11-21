epsilon = 1;
v = [0.1 1 10]; 
figure(2),clf
for k = 1:3
    [lambda, x_grid, y_grid, u] = findevs(epsilon, v(k), 250, 50, 2, false, true, false);
    figure(2)
    subplot(2,3,k)
    surf(x_grid, y_grid/pi, u{1}','EdgeColor','none','FaceColor','interp');
    view(-45,30)
    ylim([0 1])
    if v(k)==1
        title('$v = 1$','interpreter','latex')
    elseif v(k)==10
        title('$v = 10$','interpreter','latex')
    else
        title(sprintf('$v = 10^{%i}$',log10(v(k))),'interpreter','latex')
    end
    set(gca,'YDir','normal')
    ylabel('$\phi/\pi$','interpreter','latex');
    xlabel('$x$','interpreter','latex');
    zlabel('$\Psi_0$','interpreter','latex')
    set(gca, 'FontSize', 18);
    ylim([0 1])

    subplot(2,3,k+3)
    surf(x_grid, y_grid/pi, u{2}','EdgeColor','none','FaceColor','interp');
    view(-45,30)
    ylim([0 1])
    set(gca,'YDir','normal')
    ylabel('$\phi/\pi$','interpreter','latex');
    xlabel('$x$','interpreter','latex');
    zlabel('$\Psi_1$','interpreter','latex')
    set(gca, 'FontSize', 18);
    drawnow
end