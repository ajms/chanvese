function [ absgrad ] = G( phi_n, dtau, e )

a = phi_n(2:end-1,2:end-1) - phi_n(1:end-2,2:end-1)
b = phi_n(3:end,2:end-1) - phi_n(2:end-1,2:end-1)
c = phi_n(2:end-1,2:end-1) - phi_n(2:end-1,1:end-2)
d = phi_n(2:end-1,3:end) - phi_n(2:end-1,2:end-1)
absgrad = sqrt(max(ab

end