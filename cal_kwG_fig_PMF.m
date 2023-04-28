clearvars;
close all;
clc;

i = 1;       % salt concentration order
cc = [0.05 0.15 0.3 0.5 1];   % salt concentration
B = [];
C = [];
G = [];
X = [];
Y = [];
num1 = 11;
num2 = 11;
wrange = [28.76 33.08; 29.39 33.48; 29.76 33.57; 30.04 33.71; 30.31 34.02]; % twist range for different salt concentrations
Grange = [0.38 1.35; 0.31 1.18; 0.29 1.08; 0.25 1.02; 0.21 0.95];    % groove range for different salt concentrations
w1 = wrange(i,1);
w2 = wrange(i,2);
G1 = Grange(i,1);
G2 = Grange(i,2);
edgesw = linspace(w1,w2,num1);
edgesG = linspace(G1,G2,num2);

A1   = importdata(['C:\Users\Tian\Desktop\Flexibility\RNA_PO\Na\inter_bp\' int2str(i) '.txt']);  % inter-bp parameters path
A2 = importdata(['C:\Users\Tian\Desktop\Flexibility\RNA_PO\Na\groove\' int2str(i) '.txt']);   % groove parameters path
A = [A1(:,8) A2(:,1)];    % [twist, major groove width]
N = histcounts2(A(10000:60000,1),A(10000:60000,2),edgesw,edgesG,'Normalization','probability');
    for j = 1:num1-1
        for k = 1:num2-1
            if N(j,k) == 0
                B(j,k) = -log(0.0000001);
            else
                B(j,k) = -log(N(j,k));
            end
            if B(j,k) >6
                continue
            end
            X = [X; (edgesw(j)+edgesw(j+1))/2];
            Y = [Y; (edgesG(k)+edgesG(k+1))/2];
            C = [C; -log(N(j,k))];
        end
    end
    
B = B - min(min(B));           % for 2D fitting
[xang,ygroove] = find(B==min(min(B)));
mintwist = (edgesw(xang)+edgesw(xang+1))/2;
mingroove = (edgesG(ygroove)+edgesG(ygroove+1))/2;
X0 = X(:)-mintwist;
Y0 = Y(:)-mingroove;
Z0 = C - min(min(C));          % to obtain 2D fitting deviation

dw1 = (w2-w1)/(num1-1);
dG1 = (G2-G1)/(num2-1);
xx0 = linspace(G2-dG1/2,G1+dG1/2,num2-1);
yy0 = linspace(w1+dw1/2,w2-dw1/2,num1-1);
B1 = fliplr(B);
B2 = B1';

ix   = 0;
for x = xx0
    ix = ix + 1;    
    iy = 0; 
    for y = yy0  
        iy  = iy + 1;
        xd  = x - mingroove;
        yd  = y - mintwist;            
        if( B2(ix,iy)<6 ) 
           G = [G; xd,yd,B2(ix,iy)];
        end
    end
end
kws  = 1:0.1:4;
kGs  = 30:1:80;
kwGs = 5:0.1:10;
Gmin = 10000;
NG = length(G);
for kw = kws
    for kG = kGs
        for kwG = kwGs
            Gtheory        = 0.5*kG*G(:,1).^2 + 0.5*kw*G(:,2).^2 + kwG*G(:,1).*G(:,2);
            Devsum         = sum(abs(Gtheory - G(:,3)))/NG;
            if(Devsum < Gmin)
                kwmin  = kw;
                kGmin  = kG;
                kwGmin = kwG;
                Gmin   = Devsum;
            end
        end
    end
end
Gmin            % 2d fitting minimum deviation
kwmin          % kw
kGmin          % kG
kwGmin        % kwG

xx1 = linspace(w1+dw1/2,w2-dw1/2,num1-1);
yy1 = linspace(G1+dG1/2,G2-dG1/2,num2-1);
[X1,Y1] = meshgrid(xx1,yy1);
num3 = 201;
num4 = 201;
dw2 = (w2-w1)/(num3-1);
dG2 = (G2-G1)/(num4-1);
xx2 = linspace(w1+dw1/2,w2-dw1/2,num3-1);
yy2 = linspace(G1+dG1/2,G2-dG1/2,num4-1);
[X2,Y2] = meshgrid(xx2,yy2);
B3 = interp2(X1,Y1,B,X2,Y2);
x0 = [w1+dw2/2 w2-dw2/2];
y0 = [G1+dG2/2 G2-dG2/2];
clims0 = [0 6];

FZ = 18;
LZ = 4;
SZ = 16;

figure(1)
imagesc(x0,y0,B3',clims0)
set(gca,'YDir','normal')
colormap(jet)
c = colorbar;
set(get(c,'label'),'string','Potential of Mean Force (k_BT)','fontsize',FZ,'fontname','Arial')
xlabel('Twist, \omega (degree)','fontsize',FZ,'fontname','Arial')
ylabel('Major groove width, {\itG} (nm)','fontsize',FZ,'fontname','Arial')
ylim([G1 G2])
xlim([w1 w2])
set(gca,'fontsize',FZ,'linewidth',LZ,'ticklength',[0.01 0.01],'xtick',[w1 w1+2*dw1 w1+4*dw1  w1+6*dw1 w1+8*dw1 w2],'ytick',[G1 G1+2*dG1 G1+4*dG1  G1+6*dG1 G1+8*dG1 G2],'plotboxaspectratio',[4 2.6 1],'fontname','Arial','Position',[0.1 0.2 0.8 0.6])
set(gcf,'unit','normalized','position',[0.2 0.1 0.5 0.7])
saveas(1,'fig_PMF_vs_major_groove_width_twist','jpg')