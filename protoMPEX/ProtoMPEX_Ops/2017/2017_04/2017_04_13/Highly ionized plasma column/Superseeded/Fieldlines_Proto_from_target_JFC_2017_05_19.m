close all
clear all
verbose = 1;
% #########################################################################
Config = '4 kA';
switch Config
    case '6 kA'
        helicon_current = 260;
        current_A = 5900;
        current_B = 5900;
        current_C = 0;
        config = 'flat';
        ConfigName = 'TR1 0A, TR2 260 A, PS1/2 5900 A, Downstream fueling';
    case '4 kA'
        helicon_current = 160;
        current_A = 4000;
        current_B = 4000;
        current_C = 600;
        config = 'flat';
        ConfigName = 'TR1 600A, TR2 160 A, PS1/2 4000 A, Upstream fueling';
end
skimmer = 1;
shot = 0;
target_position = 2; %=1 puts it at 7.5, =2 puts it at 11.5
sleeve = 1;    
shot = 7333;    % MS edit - move up the entry for the shot number to be found or conditions simulated 
%% Build coils
[coil,current] = build_Proto_coils(helicon_current,current_A,current_B,config,verbose,current_C);
geo = get_Proto_geometry(0,0,skimmer,target_position,sleeve);
% plotit,newfig,add_skimmer,target_position,add_sleeve
%%  Calculate field lines
close all
plasma_radius_cm = 1.54;
% plasma_radius_cm = 1.52; % MAB conditions

bfield.coil = coil;
bfield.current = current;
bfield.type = 'just_coils';
bfield.vessel_clip_r   = geo.vessel_clip_r;
bfield.vessel_clip_z   = geo.vessel_clip_z;
bfield.vessel_clip_phi = 0;
bfield.stop_at_vessel = 1;
bfield.nsym = 1;

num_lines = 8; 
rr = linspace(1e-3,plasma_radius_cm/100,num_lines);  % Initial R positions %%MS - adjusted so not multiplying plasma_radius by 2
zz = geo.target.z*ones(size(rr));                        % Initial Z positions
L = 4;                                                 % Length to follow field lines; %% MS- adjusted L to 3 not 2.5 - map back further
dl = -0.01;                                              % Step size (m)
nsteps = abs(L/dl);
phistart = zeros(size(rr));
tic;
f2 = follow_fieldlines_rzphi_dz(bfield,rr,zz(1),phistart,dl,nsteps);
for i =1:length(zz)
    f.r = f2.r(:,i); f.z = f2.z(:,i); f.phi = f2.phi(:,i);
    fsave{i} = f;
end
toc
% #########################################################################
% Preview the field lines
if 0
    figure; hold on
    for i = 1:length(zz)
        fprintf('Line %d of %d\n',i,num_lines)
        f = follow_fieldlines_rzphi_dz(bfield,rr(i),zz(i),phistart(i),dl,nsteps);
        plot(f.z,f.r,'b','linewidth',2)
        fsave{i} = f;
    end
end
% #########################################################################
% Plot Field lines

% ########################################################################
% FIGURE 1B
% ########################################################################

figure; hold on
hf = gcf;
% Define the window white in color and with no menubar. the normalized
% units are relative to the computer screen
set(hf,'Menubar','figure','color','w','Units','normalized')
% Set the figure to cover the Right half of the screen
set(hf,'Position',[0.3 0.1 0.5 0.6]); % [Left Bottom Width height]hold on; box on;
for i =1:length(zz)
    h(i) = plot(fsave{i}.z,fsave{i}.r,'b','linewidth',0.5);
end
h(end).Color = 'r'; 
h(end).LineWidth = 1;
tt = xlabel('z [m]','fontsize',10); tt.HorizontalAlignment = 'right';
ylabel('r [m]','fontsize',10)
% title(['Shot ',num2str(shot)])
geo = get_Proto_geometry(1,0,skimmer,target_position,sleeve);
% To get the mean location of the coil use coil_zmean(i) = mean(zcoil(i,:))
for i = 1:size(geo.cmax,1)
    zmean_coil(i) =  mean([geo.cmin(i),geo.cmax(i)]);
end
for i = 1:(size(geo.cmax,1)-1)
    Spool(i) = 0.5*(zmean_coil(i+1)+zmean_coil(i));
end
set(gca,'fontsize',10)

% set(gcf,'Position',[0.3 0.3 0.6 0.6],'Units','normalized')
scale = 2;
switch scale
    case 1
axis([0,5,0,0.15]); box on
set(gca,'Position',[0.1 0.3 0.8 0.4],'Units','normalized','FontName','Times')
case 2
axis([0,5,0,0.15]); box on
set(gca,'Position',[0.1 0.3 0.8 0.3],'Units','normalized','FontName','Times')
end

% title(sprintf('Shot %d, I_H^*=%3.0f A',shot,helicon_current*(3300/current_A)))
hold on
for i = [1,4,6,10]
    hsp(i) = plot(Spool(i),0,'ksq');
    hsp(i).MarkerFaceColor = 'k';
%     line([Spool(i),Spool(i)],[-0.1,0])
end
set(gca,'XTick',[0 0.5 1 1.5 2.5 3 3.5 4 4.5 5])
% set(gca,'XTick',[0 0.5 1 2 3 3.5 4 4.5 5])

RootAddress = cd;
print(figure(1),[RootAddress,'\Figure1B'],'-dpng')

% #########################################################################
% Plot mod B
% JF Caneses edit 2017_05_18
zzB = linspace(0,5,600);
for s = 1:length(zzB);[Bx(s),By(s),Bz(s),Btot(s)]=bfield_bs_jdl(0,0,zzB(s),bfield.coil,bfield.current);end

% ########################################################################
% FIGURE 1A
% ########################################################################

figure; hold on
% hf = gcf;
set(gcf,'Menubar','figure','color','w','Units','normalized')

switch scale
    case 1
axis([0,5,0,1]); box on
set(gca,'Position',[0.1 0.3 0.8 0.4],'Units','normalized','XTickLabel',[],'YTick',[0:0.2:1],'FontName','Times')
% set(gcf,'Position',[0.3 0.1 0.5 0.6]); % [Left Bottom Width height]hold on; box on;
case 2
axis([0,5,0,1]); box on
set(gca,'Position',[0.1 0.3 0.8 0.2],'Units','normalized','XTickLabel',[],'YTick',[0:0.2:1],'FontName','Times')
end

plot(zzB,Btot,'k')
grid on
RootAddress = cd;
print(figure(2),[RootAddress,'\Figure1A'],'-dpng')

%%
% ########################################################################
% FIGURE 1
% ########################################################################
% full image:
close all
figure
fmag = gcf;
set(fmag,'Menubar','figure','color','w','Units','normalized')

% Part (b): Geometry and field lines
geo = get_Proto_geometry(1,0,skimmer,target_position,sleeve);
axmag(2) = gca;
set(axmag(2),'Units','Normalized','box','on','NextPlot','add')
% Draw the field lines
for i =1:length(zz)
    h(i) = plot(fsave{i}.z,fsave{i}.r,'b','linewidth',0.5,'Parent',axmag(2));
end
h(end).Color = 'r'; 
h(end).LineWidth = 1;
ylabel('$r$ $[m]$','fontsize',12,'interpreter','Latex')

% To get the mean location of the coil use coil_zmean(i) = mean(zcoil(i,:))
for i = 1:size(geo.cmax,1)
    zmean_coil(i) =  mean([geo.cmin(i),geo.cmax(i)]);
end
for i = 1:(size(geo.cmax,1)-1)
    Spool(i) = 0.5*(zmean_coil(i+1)+zmean_coil(i));
end
% Draw the location of the probes
for i = [1,4,6,10]
    hsp(i) = plot(Spool(i),0,'ksq','Parent',axmag(2));
    hsp(i).MarkerFaceColor = 'k';
%     line([Spool(i),Spool(i)],[-0.1,0])
end
line([1.596,1.896],0.06287*[1,1],'color','k','LineWidth',3,'Parent',axmag(2))
line([0.5, 0.5],[0,0.15],'color','k','LineWidth',3,'Parent',axmag(2))
line([4.33, 4.33],[0,0.04],'color','k','LineWidth',3.5,'Parent',axmag(2))
set(axmag(2),'XTick',[0 0.5 1 1.5 2.5 3 3.5 4 4.5 5])


% Part (a): On-axis magnetic field
axmag(1) = axes;
set(axmag(1),'Units','Normalized','box','on','NextPlot','add')
plot(zzB,Btot,'k','Parent',axmag(1)); 
set(axmag(1),'XTicklabel',[],'YTick',[0:0.2:1]); grid(axmag(1),'off') 
ylabel('$ B_0$ $[mT]$','fontsize',12,'interpreter','Latex')

axis(axmag(2),[0,5,0,0.15]); box on
set(axmag(1),'Position',[0.14 0.6 0.8 0.2],'Units','normalized','FontName','Times','fontsize',10)
set(axmag(2),'Position',[0.14 0.25 0.8 0.27],'Units','normalized','FontName','Times','fontsize',10)
text(axmag(1),-0.6,1,'a)','FontSize',14,'FontName','Times')
text(axmag(2),-0.6,0.15,'b)','FontSize',14,'FontName','Times')

