clear all;
close all;



%{
filename = '/home/albert/Dropbox/Uni/Bachelorprojekt/LOL/medical_images/images1.mat';

S = load(filename);
I = S.I;
%}
filename = '/home/albert/Dropbox/Uni/Bachelorprojekt/LOL/test_images/test2.jpg';
RGB = imread(filename);
I = double(rgb2gray(RGB));

figure('Position', [100 150 300 500])

% parameters
h = 1.0;
dt = 0.5;
lambda1 = 1;
lambda2 = 1;
mu = 0.5*255^2;
nu = 0;%200;
doreinit = 0;
[M, N] = size(I);

phi = -ones(M+2, N+2);

% Initialise phi_0 as a circle
[X Y] = meshgrid(1:M+2);
Z = (X-floor(M/2)).^2 + (Y-floor(N/2)).^2; 
phi(Z <= 40^2) = 1;

% Initialise phi_0 as a rectangle
%phi(5:15, 5:15) = 1;

% Initialize to a signed distance function
phi = init(phi);
tempdphi = 5000;
tempphi = zeros(M+2,N+2);
for i=1:600
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
    surf(1:M, 1:N, phi_n(2:end-1,2:end-1), 'EdgeColor', 'none');
    axis tight;
    xlabel('x');
    ylabel('y');
    zlabel('\phi^n(x, y)');
    pause(0.1);
    
    % Stopping condition, print timesteps. Procent of max value
    dphi = norm(phi_n - tempphi);
    ddphi = abs(dphi-tempdphi);
    fprintf('Difference in phi: %f, %f\n',dphi,ddphi);
    if dphi < 0.1 ||  ddphi< 0.1*dt
        fprintf('Stopping @ iteration %d\n', i);
        break;
    end
    tempphi = phi_n;
    if i>=1
        phi = init(phi_n);
        re = phi_n;
    else
        phi = phi_n;
    end
    tempdphi = dphi;
end

%set(gcf, 'PaperPositionMode', 'auto');
%print('-f1', '-djpeg', strcat(name, '-contour.jpg'));
