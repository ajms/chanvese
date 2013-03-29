function [ sign ] = S( phi, e )
    sign = phi./sqrt(phi.^2+e)   
end