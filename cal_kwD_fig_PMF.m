clearvars;
close all;
clc;

i = 1;       % salt concentration order
cc = [0.05 0.15 0.3 0.5 1];   % salt concentration
B = [];
C = [];
D = [];
X = [];
Y = [];
num1 = 11;
num2 = 11;
wrange = [28.76 33.08; 29.39 33.48; 29.76 33.57; 30.04 33.71; 30.31 34.02]; % twist range for different salt concentrations
Drange = [1.948 2.027; 1.946 2.026; 1.945 2.025; 1.944 2.024; 1.941 2.021];   % groove range for different salt concentrations
w1 = wrange(i,1);
w2 = wrange(i,2);
D1 = Drange(i,1);
D2 = Drange(i,2);
edgesw = linspace(w1,w2,num1);
edgesD = linspace(D1,D2,num2);

A1   = importdata(['C:\Users\Tian\Desktop\Flexibility\DNA_PO\Na\inter_bp\' int2str(i) '.txt']);  % inter-bp parameters path
A2 = importdata(['C:\Users\Tian\Desktop\Flexibility\DNA_PO\Na\diameter\' int2str(i) '.txt']);   % diameter path
A = [A1(:,8) A2(:,1)];    % [twist, diameter]
N = histcounts2(A(10000:60000,1),A(10000:60000,2),edgesw,edgesD,'Normalization','probability');
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
            Y = [Y; (edgesD(k)+edgesD(k+1))/2];
            C = [C; -log(N(j,k))];
        end
    end
    
B = B - min(min(B));           % for 2D fitting
[xang,ygroove] = find(B==min(min(B)));
mintwist = (edgesw(xang)+edgesw(xang+1))/2;
mingroove = (edgesD(ygroove)+edgesD(ygroove+1))/2;
X0 = X(:)-mintwist;
Y0 = Y(:)-mingroove;
Z0 = C - min(min(C));          % to obtain 2D fitting deviation

dw1 = (w2-w1)/(num1-1);
dD1 = (D2-D1)/(num2-1);
xx0 = linspace(D2-dD1/2,D1+dD1/2,num2-1);
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
           D = [D; xd,yd,B2(ix,iy)];
        end
    end
end
kws  = 1:0.1:4;
kDs  = 30:1:80;
kwDs = 5:0.1:10;
Dmin = 10000;
ND = length(D);
for kw = kws
    for kD = kDs
        for kwD = kwDs
            Dtheory        = 0.5*kD*D(:,1).^2 + 0.5*kw*D(:,2).^2 + kwD*D(:,1).*D(:,2);
            Devsum         = sum(abs(Dtheory - D(:,3)))/ND;
            if(Devsum < Dmin)
                kwmin  = kw;
                kDmin  = kD;
                kwDmin = kwD;
                Dmin   = Devsum;
            end
        end
    end
end
Dmin            % 2d fitting minimum deviation
kwmin          % kw
kDmin          % kD
kwDmin        % kwD

xx1 = linspace(w1+dw1/2,w2-dw1/2,num1-1);
yy1 = linspace(D1+dD1/2,D2-dD1/2,num2-1);
[X1,Y1] = meshgrid(xx1,yy1);
num3 = 201;
num4 = 201;
dw2 = (w2-w1)/(num3-1);
dD2 = (D2-D1)/(num4-1);
xx2 = linspace(w1+dw1/2,w2-dw1/2,num3-1);
yy2 = linspace(D1+dD1/2,D2-dD1/2,num4-1);
[X2,Y2] = meshgrid(xx2,yy2);
B3 = interp2(X1,Y1,B,X2,Y2);
x0 = [w1+dw2/2 w2-dw2/2];
y0 = [D1+dD2/2 D2-dD2/2];
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
ylabel('Diameter, {\itD} (nm)','fontsize',FZ,'fontname','Arial')
ylim([D1 D2])
xlim([w1 w2])
set(gca,'fontsize',FZ,'linewidth',LZ,'ticklength',[0.01 0.01],'xtick',[w1 w1+2*dw1 w1+4*dw1  w1+6*dw1 w1+8*dw1 w2],'ytick',[D1 D1+2*dD1 D1+4*dD1  D1+6*dD1 D1+8*dD1 D2],'plotboxaspectratio',[4 2.6 1],'fontname','Arial','Position',[0.1 0.2 0.8 0.6])
set(gcf,'unit','normalized','position',[0.2 0.1 0.5 0.7])
saveas(1,'fig_PMF_vs_diameter_twist','jpg')