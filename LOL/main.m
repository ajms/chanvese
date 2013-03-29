clear all;
close all;
%{
[FileName, PathName, FilterIndex] = uigetfile('*', 'Select image');
filename = strcat(PathName, FileName);
[pathstr, name, ext] = fileparts(filename);

if strcmpi(ext, '.mat')
    S = load(filename);
    u0 = S.I/255;
else
    RGB = imread(filename);
    u0 = double(rgb2gray(RGB))/255;
end
%}
RGB = imread('test_images/test2.jpg');
u0 = double(rgb2gray(RGB))/255;
[M, N] = size(u0);

%figure('Position', [100 150 300 500])

% parameters
circle = 1;
h = 1.0;
dt = 0.1;
di = 0.5;
lambda1 = 1;
lambda2 = 1;
mu = 0;%0.01*255^2;
nu = 200;
sigma = 0.5;
beta = 0.01;
doreinit = 0;

usmooth = real(ifftn(scalen(fftn(u0),[2,2],[0,0])));
umin = floor(min(usmooth(:)));
umax = ceil(max(usmooth(:)));
uinc = (umax-umin)/100;
ulev = umin:uinc:umax;
umean = mean(reshape(usmooth,M*N,1));
isocurves = contourc(usmooth,ulev);
%plot(ulev,histc(reshape(usmooth,1,M*N),ulev));
%meanpw = PW(umean,0.005,usmooth);
%level0 = sum(reshape(meanpw,M*N,1));

phi = -ones(M,N);
if circle == 1
    [X Y] = meshgrid(1:M);
    phip = (X-10).^2 + (Y-10).^2; 
    phi(phip <= 10^2) = 1;
elseif circle == 0
    phip = zeros(M,N);
    for i=umin:0.005:umean
        phip = phip + PW(i,0.005,usmooth);
    end
    phi(phip >= 0.5) = 1;
end    
phi = init(phi);

meanp = mean(mean(usmooth(phi>=0)));
meanm = mean(mean(usmooth(phi<0)));
histo = zeros(101,101);
i = min(min(usmooth)):(max(max(usmooth))-min(min(usmooth)))/100: ...
      max(max(usmooth));
j = min(min(phi)):(max(max(phi))-min(min(phi)))/100: ...
            max(max(phi));
a = cputime();
for k=1:101
    for l=1:101
        histo(k,l) = histiphi(i(k),j(l),usmooth,phi,beta);
    end
end
cputime()-a

surf(j,i,histo);




%surf(1:M, 1:N, phi, 'EdgeColor', 'none');
%contour(phi,[0 0])
%histphi(usmooth,phi,sigma,beta);



% Size of bins 
% for-loop through the positive / negative levelset
% do i need a signed distance map? Yes, I suppose so.