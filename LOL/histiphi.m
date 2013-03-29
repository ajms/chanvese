function [ histo ] = histiphi(i,j,usmooth,phi,beta)

histo = sum(sum(PW(i,0,usmooth).*PW(j,beta,phi)));

end