function [ sign ] = S( phi, epsilon )
    sign = phi./sqrt(phi.^2*epsilon);
end