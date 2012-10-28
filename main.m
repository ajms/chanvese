%clear all;
close all;

%[FileName, PathName, FilterIndex] = uigetfile('*', 'Select image');
%filename = strcat(PathName, FileName);
%[pathstr, name, ext] = fileparts(filename);

%if strcmpi(ext, '.mat')
%    S = load(filename);
%    u0 = S.I;
%else
%    RGB = imread(filename);
%    u0 = double(rgb2gray(RGB));
%end

RGB = imread('/home/albert/Documents/Uni/ProjNat/LevelSet/test_images/test2.jpg');
u0 = double(rgb2gray(RGB));

%u0 = 255*ones(100, 100);
%u0(30:50, 40:50) = 0;

figure('Position', [100 150 300 500])

h = 1.0;
dt = 0.05;
lambda1 = 1;
lambda2 = 15;
mu = 0.5*255^2;
nu = 0;
doreinit = 0;
[M, N] = size(u0);

phi = -ones(M, N);

% Circle
[X Y] = meshgrid(1:M);
Z = (X-50).^2 + (Y-50).^2; 
phi(Z <= 10^2) = 1;

% Rectangle
%phi(45:55, 45:55) = 1;

% Initialize to a signed distance function
phi = reinit(phi);
tempdphi = 5000;
for i=1:200
    clf;
    subplot(2,1,1);
    imagesc(u0);
    hold on;
    [C H] = contour(phi, [0 0], 'r');
    set(H, 'LineWidth', 3);
    axis tight;
    colormap gray;
    xlabel('x');
    ylabel('y');
    hold off;
    
    subplot(2,1,2);
    surf(1:M, 1:N, phi, 'EdgeColor', 'none');
    axis tight;
    xlabel('x');
    ylabel('y');
    zlabel('\phi^n(x, y)');
    
    pause(0.1);
    
    phi_n = levelset(phi, u0, lambda1, lambda2, mu, nu, dt, h, doreinit);
    
    % Stopping condition
    dphi = norm(phi_n - phi, Inf);
    if dphi < h || abs(dphi-tempdphi) < 0.1
        fprintf('Stopping @ iteration %d\n', i);
        break;
    end
    
    phi = phi_n;
    tempdphi = dphi;
end

set(gcf, 'PaperPositionMode', 'auto');
print('-f1', '-djpeg', strcat(name, '-contour.jpg'));
