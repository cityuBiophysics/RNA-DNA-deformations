clearvars;
close all;
clc;

E = cell(5,1);

for i = 1:5
A = importdata(['C:\Users\Tian\Desktop\Flexibility\RNA_PO\Rb\inter_bp\' int2str(i) '.txt']);   % inter-bp parameters
B = importdata(['C:\Users\Tian\Desktop\Flexibility\RNA_PO\Rb\groove\' int2str(i) '.txt']);     % groove parameters
C = importdata(['C:\Users\Tian\Desktop\Flexibility\RNA_PO\Rb\axis\' int2str(i) '.txt']);          % axis parameters
D = importdata(['C:\Users\Tian\Desktop\Flexibility\RNA_PO\Rb\diameter\' int2str(i) '.txt']);  % diameter
G = [A, B, C, D];
G = G(10001:60000,:);
Z = zscore(G);
rho1 = corr(G);
rho2 = corr(Z);
E{i,1} = rho1;
[r,lags] = xcorr(G(:,8),G(:,9));
stem(lags,r)
end





