function v = findl1(epsilon,logv0)
N = length(epsilon);
v = zeros(size(epsilon)); dlogv = 0.25;
if isempty(logv0), logv0 = [-10 0]; end
options = optimset('TolX',1e-4,'display','off');
for k = 1:N
    fval = 1;
    while fval>1e-8
        [logv, fval] = fminbnd(@(logv) f(logv,epsilon(k)),logv0(1),logv0(2),options);
        v(k) = exp(logv);
        if fval>1e-8
            disp('fail')
        else
            disp([epsilon(k) v(k) fval])
        end
        logv0 = [logv0(1)+dlogv logv0(2)+dlogv];
    end
    logv0 = [logv logv+dlogv];
end

%%
function f = f(logv,epsilon)
v = exp(logv);
l = findevs(epsilon,v,min(250,max(25,ceil(7.5/sqrt(epsilon)))),32,3,true,false,false);
f = real(l(3)-l(2))^2;