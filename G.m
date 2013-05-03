function [ absgrad ] = G( phi, h )
    a = (phi(2:end-1,2:end-1)-phi(1:end-2,2:end-1))/h;
    b = (phi(3:end,2:end-1)-phi(2:end-1,2:end-1))/h;
    c = (phi(2:end-1,2:end-1)-phi(2:end-1,1:end-2))/h;
    d = (phi(2:end-1,3:end)-phi(2:end-1,2:end-1))/h;
    absgrad = zeros(size(phi));
    absgrad(2:end-1,2:end-1) = (phi(2:end-1,2:end-1)>0).*(sqrt(max(((a>0).*a).^2,((b<0).*b).^2) + max(((c>0).* c).^2,((d<0).*d).^2))-1)+ (phi(2:end-1,2:end-1)<0).*(sqrt(max(((a<0).*a).^2,((b>0).*b).^2) + max(((c<0).*c).^2,((d>0).*d).^2))-1);
end