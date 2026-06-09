% Read data test:
clear all
% close all
clc

% Select case:
% =========================================================================
caseData = 8;

% Read data:
% =========================================================================
d{1} = readcell(['T2D_BackSurface_case_',num2str(caseData),'.txt']);
d{2} = readcell(['q2D_FrontSurface_case_',num2str(caseData),'.txt']);

% Spatial coordinates:
% =========================================================================
xx = str2num(cell2mat(d{1}(2)));
zz = str2num(cell2mat(d{1}(3)));
NS = 150;
tt = linspace(0,1.5,NS);

% Temperature surface:
% =========================================================================
kk = 5;
for jj = 1:NS
    for ii = 1:numel(xx)
        T(ii,:,jj) = str2num(cell2mat(d{1}(kk)))';
        q(ii,:,jj) = str2num(cell2mat(d{2}(kk)))';
        kk = kk + 1;
    end
    kk = kk + 1;
end

% Plot data:
% =========================================================================
if 0
    figure('color','w'); 
    for jj = 1:1:150
        mesh(T(:,:,jj));
        caxis([30,1000])
        zlim([30,1000])
        view([30,10])
        drawnow
    end
end

if 0
    figure('color','w'); 
    for jj = 1:1:150
        mesh(q(:,:,jj));
        caxis([0,5]*1E6)
        zlim([0,5]*1E6)
        view([30,10])
        drawnow
    end
end

figure('color','w'); 
plot(tt,permute(T(50,50,:),[3,1,2]))
grid on
ylim([300,1000])

figure('color','w'); 
plot(tt,permute(q(50,50,:)*1E-6,[3,1,2]))
grid on
ylabel('[MWm$^{-2}$]','interpreter','latex')
xlabel('[s]','interpreter','latex')