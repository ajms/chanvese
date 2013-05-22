close all;

% parameters
no = '8';
maxit = 60;
h = 1;
dt = 0.5;
lambda1 = 1;
lambda2 = 1;
mu = 0.1*255^2;
nu = 0;
k = 100;
l = 0.5*10^5;
m = 30;
[M, N] = size(I);

phi = -ones(M+2, N+2);

% Initialise phi_0 as a circle
[X Y] = meshgrid(1:M+2);
Z = (X-50).^2 + (Y-0).^2; 
phi(Z <= m^2) = 1;


% Initialize to a signed distance function
phi = init(phi);

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

tic;
for i=1:maxit
    fprintf('Iteration: %d\n',i);
    [phi_n, F(i+1)] = chlevelset(phi, I, lambda1, lambda2, mu, nu, dt, h);

    % Stopping condition, print timesteps. Procent of max value
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


seg = figure;
hold on;
imagesc(I/255);
colorbar();
colormap('gray');
contour(phi(2:end-1,2:end-1), [0 0], 'Color', [1 0 0],'LineWidth',3);
axis tight;
hold off;
print(seg,'-dpsc',strcat('I',no,'segcv.eps'));

Fplot = figure;
plot(1:sum(F>0),F(F>0));
set(gca,'xtick',0:1:maxit);
title('Convergence plot','interpreter','latex','FontSize',15);
xlabel('Iteration no.','interpreter','latex','FontSize',15);
ylabel('Value of $F$','interpreter','latex','FontSize',15);
print(Fplot,'-dpsc',strcat('I',no,'concv.eps'));

phiplot = figure;
plot(1:sum(dphi>0),dphi(dphi>0));
set(gca,'xtick',0:1:maxit);
title('Convergence plot','interpreter','latex','FontSize',15);
xlabel('Iteration no.','interpreter','latex','FontSize',15);
ylabel('Difference in $\varphi$','interpreter','latex','FontSize',15);
print(phiplot,'-dpsc',strcat('I',no,'con2cv.eps'));

inseg = (phi(2:end-1,2:end-1)>=0);
dev = sum(abs(inseg(:)-Iref(:)))/(M*N);
fprintf('Deviation in pct. %f\n',dev);