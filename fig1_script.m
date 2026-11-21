format compact
format short e
tic
N = 100; M = 100;
epsilon = logspace(-3,-0.8,N);
v = linspace(0,1,M);
l1 = ones(N,M);
parfor j = 1:M
    k = 1; l2 = ones(N,1);
    while ((k==1)||(abs(l2(k-1)-1)>1e-5))&&(k<=N)
        l = findevs(epsilon(k),v(j),100,50,2,true,false,false);     %#ok<PFBNS>
        l2(k) = l(2);
        k = k+1;
    end
    l1(:,j) = l2;
    disp([j k M N])
end

eps1 = logspace(-1.01,-3,250);
v1 = findl1(eps1,[]);

figure(1),clf
imagesc(v,log10(epsilon),l1);
set(gca,'YDir','normal')
colormap default
ylabel('$\log_{10} \epsilon$','interpreter','latex')
xlabel('$v$','interpreter','latex')
set(gca,'FontSize',18)
hold on
[C,h] = contour(v,log10(epsilon),l1,linspace(0,0.9,10),'k',"ShowText",true,"LabelFormat","%0.1f",'LabelColor','k','Linewidth',1.4);
clabel(C,h,"FontSize",16)
plot(v1,log10(eps1),'k--','Linewidth',1.4)
set(gca,'FontSize',18)

save lambda1_data epsilon v l1 eps1 v1
toc