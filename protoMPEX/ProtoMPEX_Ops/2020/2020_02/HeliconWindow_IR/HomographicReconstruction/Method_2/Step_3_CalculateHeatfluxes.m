% Objective:
% Load the Composite data and calculate the time dependent heat flux on
% each pixel

clc
close all
clear all

saveData = 1;
saveFig  = 1;

for kk = [2:8]
    kk
clearvars -except kk saveData saveFig
% Load composite data:
% =========================================================================
fileName = ['CompositeData_ShotSeries_',num2str(kk),'.mat'];
load(fileName);    

% Material properties:
% =========================================================================
rho = 3300; 
kt   = 180 ;
cp  = 740 ;
a = kt/(rho*cp); % thermal diffusivity

% Geometry
% =========================================================================
Lz = 6.3500/1000; % ALN wall thickness in ProtoMPEX

% Normalize data:
% =========================================================================
% Order of magnitude heat flux:
q0 = 1e6; 
% Characteristic time scale:
t_star = Lz*Lz/a;
% Characteristic temperature:
T_star = q0*Lz/kt;
% Initial temperature distribution:
T0 = f.dT(:,:,1);
% Define size of data:
[I,J,R] = size(f.dT);
% Normalized data:
for r = 1:R
    u(:,:,r) = (f.dT(:,:,r) - T0)/T_star;
end
t_u = f.t_dT/t_star;
dt_u = t_u(2)-t_u(1);
    
% Permute data:
% =========================================================================
z    = permute(u,[3,1,2]);
t_z  = t_u;
dt_z = dt_u;
    
% Compute Toeplitz matrix:
% =========================================================================
xx0 = 0.0005;
% To capture the details of the front surface we need to use 1e4
% partial sum terms and use xx0 = 0.01 to 0.05
Ns = 1e4;
for r = 1:R
   K0(r) = G_Impulse_1D(xx0,t_z(r),Ns);
end
P = toeplitz(K0, zeros( size(K0) ) );

%  Conjugate gradient method:
% =========================================================================
% Input data:
% P: convolution operator
% q{1}: initial guess at heat flux, column vector
% z: experimental data, column vector
% Ni: number of iteration
% Output data:
% u: calculated temperature, column vector
% J: residual, column vector
% q: minimized heat flux, structure, column vector
disp('Applying inverse method...')
initType = 3;
tic
for ii = 1:I
    for jj = 1:J
%         tic
        zdata = z(:,ii,jj);
        % Initial guess:
        switch initType
            case 1
                qq{1} = zeros(size(t_z))';
                Ni = 350;
            case 2
                qq{1} = diff(zdata)/dt_z;
                qq{1} = sgolay_t([qq{1};qq{1}(end)],3,5);
                Ni = 150;
            case 3
                qq{1} = diff(zdata)/dt_z;
                qq{1} = [qq{1};qq{1}(end)];
                Ni = 200;
        end
        % Inverse solution:
        [a1,a2,a3] = IHCP_ConjugateGradient(P,qq,zdata,dt_z,Ni,3,[]);
        % Normalized heat flux:
        q(:,ii,jj) = a1{end};
        % Reconstructed normalized temperature:
        v(:,ii,jj) = a2{end};
        % Residual:
        res(:,ii,jj) = a3;
%         toc
    end
end
toc
disp('Computation completed!') 
figure; plot(z(:,24:29,66));  hold on; plot(v(:,24:29,66),'r')
hold on; plot(q(:,24:29,66),'r')

%%
figure('color','w')
bb = 0.7;
La = 30/100;
qq = permute(q,[2,3,1]);
mean_RF = round(mean(f.rfPwr));
for fr = 5:1:40
        heatflux = (qq(:,:,fr) + qq(:,:,fr+1) + qq(:,:,fr+2))/3;
        surf(f.phi_2D,(f.s_2D/bb)*1e2,heatflux*1e3,'LineStyle','none')
        set(gca,'XTick',[0:45:360],'XDir','reverse')
        set(gca,'YTick',[0:5:(La/bb)*1e2])
        xlabel('Angle [deg]','Interpreter','Latex','FontSize',14)
        ylabel('z [cm]','Interpreter','Latex','FontSize',14)
        view([0,90])
        axis tight
        colormap('bone')
        colormap('hot')
        xlim([0,360])
        caxis([0,800])
        colorbar
        title(['Heat flux [kW/m^2] , ',cell2mat(f.limitMode),' at ',num2str(mean_RF),' kW , frame: ',num2str(fr)])

        % Draw antenna:
        lw = 13;
        % Transverse straps:
        line([360,000],[26,26],[500,500],'color','k','LineWidth',lw)
        line([360,000],[04,04],[500,500],'color','k','LineWidth',lw)
        % Bottom helical strap:
        line([225,135],[04,26],[500,500],'color','k','LineWidth',lw)
        % HV top side helical strap:
        line([045,000],[04,15],[500,500],'color','k','LineWidth',lw)
        % GND top side helical strap:
        line([360,315],[15,26],[500,500],'color','k','LineWidth',lw)
        drawnow
        pause(0.1)
end

% Organizing output data into a structure:
% =========================================================================
% Heat flux data:
w.qnorm   = qq;
w.q0      = q0;
w.t_qnorm = t_z;
w.t_star  = t_star; 
w.phi_2D  = f.phi_2D;
w.z_2D    = f.s_2D/bb;
w.bb      = bb;

% Metadata: 
w.limitMode    = f.limitMode;
w.shots        = f.shots;
w.rfPwr        = f.rfPwr;
w.thermalParam = f.thermalParam;
w.comment{1}   = 'Heat flux in [W] is obtained by qnorm*q0';
w.comment{2}   = 'Time in [s] is obtained by t_qnorm*t_star, where t_star is the characteristic time constant';

% Window properties in SI units:
w.windowData.rho  = rho; 
w.windowData.kt   = kt;
w.windowData.cp   = cp;
w.windowData.a    = a ; % thermal diffusivity
w.windowThickness = Lz; % Thickness in [m]

% Save data in .mat files:
% =========================================================================
if saveData
    t1 = tic;
    disp('Saving data ...')
    fileName = ['HeatfluxData_ShotSeries_',num2str(kk),'.mat'];
    save(fileName,'w')
    t1 = toc(t1);
    disp(['Data Saved!! took ',num2str(t1),' seconds'])
    beep
end

% Save figures:
% =========================================================================
if saveFig
    figureName = ['HeatFlux_',cell2mat(f.limitMode),'_RF_',num2str(mean_RF),'kW'];
    saveas(gcf,figureName,'tiffn')
end

end
disp('End of script!!')