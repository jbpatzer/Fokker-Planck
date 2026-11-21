format compact
format short e
tic
N = 250; M = 250;
epsilon = logspace(-3,0,N);
v = linspace(0,0.5,M);
imax = zeros(N,M);
Ny = 50; Nx = 100;
Nev = 50;
parfor j = 1:M
    imax1 = zeros(N,1);
    for k = 1:N
        l = findevs(epsilon(k),v(j),Nx,Ny,Nev,true,false,false);   %#ok<PFBNS>
        imax1(k) = sum(abs(imag(l)));
    end
    imax(:,j) = imax1;
    disp([j M])
end


figure(1),clf
imagesc(v,log10(epsilon),log10(1+imax));
set(gca,'YDir','normal')
colormap default
ylabel('$\log_{10} \epsilon$','interpreter','latex')
xlabel('$v$','interpreter','latex')
set(gca,'FontSize',18)
c = colorbar;
c.Label.String = '$\log_{10}(1+S_I(v,\epsilon))$';
c.Label.Interpreter = 'latex';
set(gca,'FontSize',18)

save imag_data2 epsilon v imax Nev v1 eps1
toc