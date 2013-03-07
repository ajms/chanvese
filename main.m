clear all;
close all;

[FileName, PathName, FilterIndex] = uigetfile('*', 'Select image');
filename = strcat(PathName, FileName);
[pathstr, name, ext] = fileparts(filename);

if strcmpi(ext, '.mat')
    S = load(filename);
    u0 = S.I;
else
    RGB = imread(filename);
    u0 = double(rgb2gray(RGB));
end

%RGB = imread('/home/albert/Dropbox/Uni/ProjNat/LevelSet/test_images/test2.jpg');
%u0 = double(rgb2gray(RGB));

figure('Position', [100 150 300 500])

% parameters
h = 1.0;
dt = 0.1;
lambda1 = 1;
lambda2 = 1;
mu = 1;%0.05*255^2;
nu = 0;
doreinit = 0;
[M, N] = size(u0);

phi = -ones(M, N);

% Initialise phi_0 as a circle
[X Y] = meshgrid(1:M);
Z = (X-50).^2 + (Y-50).^2; 
phi(Z <= 10^2) = 1;

% Initialise phi_0 as a rectangle
%phi(5:15, 5:15) = 1;

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
    fprintf('Iteration: %d\n',i);
    phi_n = chlevelset(phi, u0, lambda1, lambda2, mu, nu, dt, h);
    
    % Stopping condition
    dphi = norm(phi_n - phi, Inf);
    if dphi < 10*h || abs(dphi-tempdphi) < 0.1
        fprintf('Stopping @ iteration %d\n', i);
        break;
    end
    phi = phi_n;
    tempdphi = dphi;
end

set(gcf, 'PaperPositionMode', 'auto');
print('-f1', '-djpeg', strcat(name, '-contour.jpg'));
