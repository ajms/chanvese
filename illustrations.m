close all;
clear all;

S = load('testimages.mat');
I = S.I1;

x = -5:0.1:5;
Hside = figure;
plot(x,hside(x,0.1),x,hside(x,0.5),x,hside(x,1));
title('Regularized Heavyside function $H_\epsilon$','interpreter', 'latex','FontSize', 15);
xlabel('$x$', 'interpreter', 'latex', 'FontSize', 15);
ylabel('$y$', 'interpreter', 'latex','FontSize', 15);
Hleg = legend('$H_{0.1}$', '$H_{0.5}$', '$H_1$');
set(Hleg,'Interpreter', 'latex', 'location', 'SouthEast', 'FontSize', 15);
print(Hside,'-dpsc','Hside.eps');
uiwait;

Dirac = figure;
plot(x,drac(x,0.1),x,drac(x,0.5),x,drac(x,1));
title('Regularized Dirac measure $\delta_\epsilon$','interpreter', 'latex','FontSize', 15);
xlabel('$x$', 'interpreter', 'latex','FontSize', 15);
ylabel('$y$', 'interpreter', 'latex','FontSize', 15);
dleg = legend('$\delta_{0.1}$', '$\delta_{0.5}$', '$\delta_1$');
set(dleg,'Interpreter', 'latex', 'location', 'SouthEast', 'FontSize', 15);
print(Dirac,'-dpsc','Dirac.eps');

uiwait;

RandIm = figure;
imagesc(I);
title('Image with Uniformly distributed noise','interpreter', 'latex','FontSize', 15);
colorbar();
colormap('gray');
print(RandIm,'-deps','RandIm.eps');

uiwait;
i = min(I(:)):(max(I(:))-min(I(:)))/99:max(I(:));
histd = histc(I(:),i);
HistRandIm = figure;
hold on;
stairs(i,histd/sum(histd));

hists = zeros(length(i),1);
for k=1:length(i)
    hists(k) = sum(sum(PW(i(k),0.02,I)));
end

plot(i,hists/sum(hists),'Color',[1 0 0]);
title('Histogram of random image','interpreter', 'latex','FontSize', 15);
xlabel('Image intensity', 'interpreter', 'latex','FontSize', 15);
ylabel('Frequency', 'interpreter', 'latex','FontSize', 15);
hleg = legend('Normal histogram, 100 bins', 'Smooth histogram, $\beta = 0.01$');
set(hleg,'Interpreter', 'latex', 'location', 'NorthEast', 'FontSize', 15);
hold off;
print(HistRandIm,'-dpsc','HistRandIm.eps');

uiwait;

Ismooth = real(ifftn(scalen(fftn(I),[1 1],[0,0])));
Ismooth = Ismooth;

SmoothRandIm = figure;
imagesc(Ismooth);
title('Smoothed image, $\sigma=4$','interpreter', 'latex','FontSize', 15);
colorbar();
colormap('gray');
print(SmoothRandIm,'-dpsc','SmoothRandIm.eps');

uiwait;

i = min(Ismooth(:)):(max(Ismooth(:))-min(Ismooth(:)))/99:max(Ismooth(:));
histd = histc(Ismooth(:),i);

HistSmoothRandIm = figure;
hold on;
stairs(i,histd/sum(histd));

hists = zeros(length(i),1);
for k=1:length(i)
    hists(k) = sum(sum(PW(i(k),0.001,Ismooth)));
end

plot(i,hists/sum(hists),'Color',[1 0 0]);
title('Histogram of smoothed image','interpreter','latex','FontSize',15);
xlabel('Image intensity', 'interpreter', 'latex','FontSize', 15);
ylabel('Frequency', 'interpreter', 'latex','FontSize', 15);
hleg = legend('Normal histogram, 100 bins', 'Smooth histogram, $\beta = 0.001$');
set(hleg,'Interpreter', 'latex', 'location', 'NorthEast', 'FontSize', 15);
hold off;
print(HistSmoothRandIm,'-dpsc','HistSmoothRandIm.eps');

uiwait;

isophote1 = PW(0.7,0.1,Ismooth);
isophote2 = PW(0.5,0.1,Ismooth);

SoftIso1 = figure;
imagesc(isophote1);
title('Soft isophote, $\beta = 0.1$, $i=0.8$','interpreter','latex','FontSize',15);
colorbar();
colormap('gray');
print(SoftIso1,'-dpsc','SoftIso1.eps');

uiwait;

SoftIso2 = figure;
imagesc(isophote2);
title('Soft isophote, $\beta = 0.1$, $i=0.5$','interpreter','latex','FontSize',15);
colorbar();
colormap('gray');
print(SoftIso2,'-dpsc','SoftIso2.eps');