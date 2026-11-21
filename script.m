format compact
format short e
s = 0.01;
N = length(s);
l1 = zeros(N,101);
for k = 1:N
    th = linspace(0,pi/2,101);
    for j = 1:101
        epsilon = s(k)*cos(th(j));
        l = findevs(epsilon,s(k)*sin(th(j)),min(250,max(25,ceil(2.5/sqrt(epsilon)))),50,2,true,false,false);
        l1(k,j) = l(2);
        disp([k j epsilon l(1) l(2)])
    end
end
semilogy(th,l1,'linewidth',1.4)