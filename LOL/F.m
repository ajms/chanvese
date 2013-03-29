function [ Fu ] = F( u,phi,nu,mu,sigma,beta,lambda1,lambda2 )
%F Summary of this function goes here
%   Detailed explanation goes here

histo = zeros(101,101);
for i = min(min(usmooth)):(max(max(usmooth))-min(min(usmooth)))/100: ...
      max(max(usmooth))
    for j = j = min(min(phi)):(max(max(phi))-min(min(phi)))/100: ...
            max(max(phi))
        histo(i,j) = histiphi(i(k),-20,usmooth,phi,beta)*(meanm+i(k));
Fu = lambda1*sum(sum(abs(histp-mean(histp))) + lambda2*sum(abs(histm-mean(histm)));
hist3([histp,histm],[hist2p,hist2m]);
end

