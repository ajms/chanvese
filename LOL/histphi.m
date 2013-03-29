function [ histp, histm ] = histphi(usmooth,phi,sigma,beta)

[M, N] = size(usmooth);
umin = floor(min(usmooth(:)));
umax = ceil(max(usmooth(:)));
uinc = (umax-umin)/500;
mphi = mean(usmooth(abs(phi)<=0.5));
stepsm = mphi:uinc:umax;
stepsp = umin:uinc:mphi-uinc;
histm = zeros(1,size(stepsm,2));
histp= zeros(1,size(stepsp,2));

for j=1:size(stepsm,2)
    histm(1,j) = sum(reshape(PW(stepsm(1,j),beta,usmooth),M*N,1));
end

for j=1:size(stepsp,2)
    histp(1,j) = sum(reshape(PW(stepsp(1,j),beta,usmooth),M*N,1));
end
histsum = sum([histp,histm]);
histm = histm/histsum;
histp = histp/histsum;
hold on;
plot(stepsm,histm,'r');
plot(stepsp,histp,'b');
hold off;
end