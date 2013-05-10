%clear all;
close all;
%{
[FileName, PathName, FilterIndex] = uigetfile('*', 'Select image');
filename = strcat(PathName, FileName);
[pathstr, name, ext] = fileparts(filename);
if strcmpi(ext, '.mat')
    S = load(filename);
    I = S.I;
else
    RGB = imread(filename);
    I = double(rgb2gray(RGB));
end

I = real(ifftn(scalen(fftn(I),[10,10],[0,0])));
%}
figure('Position', [100 150 300 500])

% parameters
h = 1.0;
dt = 0.1;
lambda1 = 1;
lambda2 = 1;
mu = 0.4*255^2;
nu = 0;
[M, N] = size(I);

phi = -ones(M+2, N+2);

% Initialise phi_0 as a circle
[X Y] = meshgrid(1:M+2);
Z = (X-floor(M/2)).^2 + (Y-floor(N/2)).^2; 
phi(Z <= 20^2) = 1;

% Initialize to a signed distance function
phi = init(phi);
tempdphi = 5000;
tempphi = zeros(M+2,N+2);
dphi = zeros(100,1);
ddphi = zeros(100,1);
for i=1:100
    fprintf('Iteration: %d\n',i);
    phi_n = chlevelset(phi, I, lambda1, lambda2, mu, nu, dt, h);
    
    subplot(2,1,1);
    imagesc(I);
    hold on;
    [C H] = contour(phi_n(2:end-1,2:end-1), [0 0], 'r');
    set(H, 'LineWidth', 3);
    axis tight;
    xlabel('x');
    ylabel('y');
    hold off;
    
    subplot(2,1,2);
    surf(phi_n(2:end-1,2:end-1), 'EdgeColor', 'none');
    axis tight;
    xlabel('x');
    ylabel('y');
    zlabel('\phi^n(x, y)');
    pause(0.1);
    
    % Stopping condition, print timesteps. Procent of max value
    dphi(i+1) = norm(phi_n(2:end-1,2:end-1) - phi(2:end-1,2:end-1));
    ddphi(i+1) = abs(dphi(i+1)-dphi(i));
    fprintf('Difference in phi: %f, %f\n',dphi(i+1),ddphi(i+1));
    if dphi(i+1) < 0.1*h ||  ddphi(i+1) < 0.01*h
        fprintf('Stopping @ iteration %d\n', i-1);
        break;
    end
    phi = phi_n;
end

%set(gcf, 'PaperPositionMode', 'auto');
%print('-f1', '-djpeg', strcat(name, '-contour.jpg'));
