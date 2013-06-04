close all;

% Import image
I = rgb2gray(imread('/path/to/image.jpg'));

% Parameters
no = '0'; % number of image for figures
maxit = 60; % maximal number of iterations
h = 1; % space step size
dt = 0.5; % time step size
lambda1 = 1; % parameter weigting the inner segment
lambda2 = 1; % parameter weighting the outer segment
mu = 0.1*255^2; % parameter weighting the length regularization
nu = 0; % parameter weighting the area regularization
k = 100; % stopping criterion for difference in F
l = 0.5*10^5; % stopping criterion for difference in phi
m = 30; % size of circle for initial phi
[M, N] = size(I);

% Initialize phi as a circle
phi = -ones(M+2, N+2);
[X Y] = meshgrid(1:M+2);
Z = (X-50).^2 + (Y-50).^2; 
phi(Z <= m^2) = 1;
phi = init(phi);

% Initial plot
initial = figure;
hold on;
imagesc(I/255);
colorbar();
colormap('gray');
contour(phi(2:end-1,2:end-1), [0 0], 'Color', [1 0 0],'LineWidth',3);
axis tight;
print(initial,'-dpsc',strcat('I',no,'initcv.eps'));

dphi = zeros(maxit,1);
F = zeros(maxit,1);

% Iteration over time steps
tic;
for i=1:maxit
    fprintf('Iteration: %d\n',i);
    [phi_n, F(i+1)] = chlevelset(phi, I, lambda1, lambda2, mu, nu, dt, h);
    dphi(i) = norm(phi_n(2:end-1,2:end-1) - phi(2:end-1,2:end-1));
    dF = abs(F(i+1)-F(i));
    fprintf('Difference in F, difference in phi: %f, %f\n',dF,dphi(i));
    if dphi(i) < k && dF < l
        fprintf('Stopping @ iteration %d\n', i);
        break;
    end
    phi = phi_n;
end
toc;

% Plot of final segmentation
seg = figure;
hold on;
imagesc(I/255);
colorbar();
colormap('gray');
contour(phi(2:end-1,2:end-1), [0 0], 'Color', [1 0 0],'LineWidth',3);
axis tight;
hold off;
print(seg,'-dpsc',strcat('I',no,'segcv.eps'));

% Convergence plot over F
Fplot = figure;
plot(1:sum(F>0),F(F>0));
set(gca,'xtick',0:1:maxit);
title('Convergence plot','interpreter','latex','FontSize',15);
xlabel('Iteration no.','interpreter','latex','FontSize',15);
ylabel('Value of $F$','interpreter','latex','FontSize',15);
print(Fplot,'-dpsc',strcat('I',no,'concv.eps'));

% Convergence plot over phi
phiplot = figure;
plot(1:sum(dphi>0),dphi(dphi>0));
set(gca,'xtick',0:1:maxit);
title('Convergence plot','interpreter','latex','FontSize',15);
xlabel('Iteration no.','interpreter','latex','FontSize',15);
ylabel('Difference in $\varphi$','interpreter','latex','FontSize',15);
print(phiplot,'-dpsc',strcat('I',no,'con2cv.eps'));

% For syntetic images: difference to reference image
inseg = (phi(2:end-1,2:end-1)>=0);
dev = sum(abs(inseg(:)-Iref(:)))/(M*N);
fprintf('Deviation in pct. %f\n',dev);