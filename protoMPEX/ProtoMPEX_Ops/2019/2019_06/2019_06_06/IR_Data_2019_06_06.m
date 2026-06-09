close all 
clear all
 
ViewType = 3;
switch ViewType
    case 1
    % Middle-High view, no HV observed
    shot(1) = 26569; % ?
    case 2
    % higher view, HV side hot spot can be seen
    shot(1) = 26570; % 300 + 2 notches
    case 3
    % Middle-center view, no HV view
    shot(1) = 26572; % 300 + 1 notches
    case 4
    % Lower view
    shot(1) = 26574; % 320 + 0 notches
end
        
a{1} = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2019_06_06';

mov = 1;

for s = 1:length(shot)
PATHNAME = [a{s},'\'];
try
FILENAME = ['Shot ',num2str(shot(s)),'.seq'];
catch
end
% FILENAME = ['Shot-0',num2str(shot(s)),'.seq'];    
videoFileName=[PATHNAME FILENAME];


% Load the Atlats SDK
atPath = getenv('FLIR_Atlas_MATLAB');
atImage = strcat(atPath,'Flir.Atlas.Image.dll');
asmInfo = NET.addAssembly(atImage);
%open the IR-file'
file = Flir.Atlas.Image.ThermalImageFile(videoFileName);
seq = file.ThermalSequencePlayer();
%Get the pixels
img = seq.ThermalImage.ImageProcessing.GetPixelsArray;
im = double(img);
    
Data{s}(:,:,1) = im;
fr = 1;
if(seq.Count > 1)
    while(seq.Next())
        img = seq.ThermalImage.ImageProcessing.GetPixelsArray;
        im = double(img);
%         Data(:,:,fr) = im;
        Data{s}(:,:,fr) = im(end:-1:1,:);         
        fr = fr + 1;

    end
end
end

% Gather FP data (window temperature)
b = 0;

% To convert voltage to temperature, multiply by 20
Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];

% Gather data
DA{1} = [RootAddress,'FLUOROPT_1']; 
DA{2} = [RootAddress,'FLUOROPT_2']; 
DA{3} = [RootAddress,'FLUOROPT_3']; 
DA{4} = [RootAddress,'FLUOROPT_4']; 

[f1,t_f1]   = my_mdsvalue_v2(shot(1),DA(1));
[f2,t_f2]   = my_mdsvalue_v2(shot(1),DA(2));
[f3,t_f3]   = my_mdsvalue_v2(shot(1),DA(3));
[f4,t_f4]   = my_mdsvalue_v2(shot(1),DA(4));

%%
for s = 1:size(Data{1},3)
   SS(s) = sum(sum(Data{1}(:,:,s))); 
end

figure; 
plot(SS,'k.')

[~,n1] = max(diff(SS));
[~,n2] = min(diff(SS));

rng = [n1-10:n2+40];
D0 = Data{1}(:,:,n1-10);
% D0 = Data{1}(:,:,b+20);

figure; 
hold on
plot(SS,'k.')
plot(rng,SS(rng),'g')

figure; 
hold on
plot(diff(SS),'k.')
plot(rng(1:end-1),diff(SS(rng)),'g')

%% Image prior to RF
figure
 surf(D0,'LineStyle','none')
 view([0,90])
%  caxis([-200,1800])
 colormap('hot')
  title(['shot: ',num2str(shot(1))])
  set(gca,'PlotBoxAspectRatio',[1 0.5 1])
  
C = {'k','bl','r','g'};

figure; hold on
plot(t_f1{1}(1:end-1),(f1{1}-b*min(f1{1}))*20,C{1},'LineWidth',2)
plot(t_f2{1}(1:end-1),(f2{1}-b*min(f2{1}))*20,C{2},'LineWidth',2) % Ground Side
plot(t_f3{1}(1:end-1),(f3{1}-b*min(f3{1}))*20,C{3},'LineWidth',2) % High Voltage side
plot(t_f4{1}(1:end-1),(f4{1}-b*min(f4{1}))*20,C{4},'LineWidth',2)

ylim([0,100])
ylabel('${\Delta}T$ $[C]$','Interpreter','latex','FontSize',14)
xlabel('$t$ $[sec]$','Interpreter','latex','FontSize',14)
box on
set(gcf,'color','w')
grid on


%% Window heating with RF
figure
for s = rng
     surf(Data{1}(:,:,s)-1*D0,'LineStyle','none')
      view([0,90])
      caxis([0,1200])
      colormap('hot')
      title(['shot: ',num2str(shot(1)),', frame: ',num2str(s)])
      set(gca,'PlotBoxAspectRatio',[1 0.5 1])
%      drawnow
     pause(0.1)
end

%%
close all
figure
for s = rng(1:end-10)
     surf(Data{1}(:,:,s+1)-Data{1}(:,:,s),'LineStyle','none')
      view([0,90])
      caxis([0,100])
      colormap('hot')
      title(['shot: ',num2str(shot(1)),', frame: ',num2str(s)])
      set(gca,'PlotBoxAspectRatio',[1 0.5 1])
%      drawnow
     pause(0.01)
end
