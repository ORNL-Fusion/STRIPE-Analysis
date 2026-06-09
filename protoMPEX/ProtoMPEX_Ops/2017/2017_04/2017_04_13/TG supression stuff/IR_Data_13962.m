% IR data analysis April 13th 2017, shot 13962
close all 
clear all

shot = 13962;
mov = 1;
UploadType = 2;

switch UploadType
    case 1
        a = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2017_04_13';
    case 2
        a = cd;
end

PATHNAME = [a,'\'];
FILENAME = ['shot ',num2str(shot),'.seq'];
videoFileName=[PATHNAME FILENAME];

% Load the Atlats SDK
atPath = getenv('FLIR_Atlas_MATLAB');
atImage = strcat(atPath,'Flir.Atlas.Image.dll');
asmInfo = NET.addAssembly(atImage);
%open the IR-file
file = Flir.Atlas.Image.ThermalImageFile(videoFileName);
seq = file.ThermalSequencePlayer();
%Get the pixels
img = seq.ThermalImage.ImageProcessing.GetPixelsArray;
im = double(img);
    
Data(:,:,1) = im;
fr = 1;
if(seq.Count > 1)
    while(seq.Next())
        img = seq.ThermalImage.ImageProcessing.GetPixelsArray;
        im = double(img);
%         Data(:,:,fr) = im;
        Data(:,:,fr) = im(end:-1:1,:);         
        fr = fr + 1;

    end
end

%% Plot raw data:
if mov
figure(1)
for s = 1:60
   imshow(Data(:,:,s),[]);
   drawnow;
   title(num2str(s))
   pause(0.1)
end
end

%% Determine the first RF frame
close all
for s = 1:59
    dIR{s+1}(:,:) = Data(:,:,s+1)-Data(:,:,s);
    SdIR(s+1) = sum(sum(abs(dIR{s+1}(:,:).^2),1));
end

figure; hold on
m = diff(SdIR - SdIR(1));
nm = find(m>max(m)/30);
plot(m)
plot(nm(1),m(nm(1)),'ro')
ylim([0,max(m)])

figure;
subplot(1,2,1);
contourf(dIR{nm(1)+0},20,'lineStyle','none');
caxis([0,200])
axis('square')
title(['frame: ',num2str(nm(1))])

subplot(1,2,2);
contourf(dIR{nm(1)+1},20,'lineStyle','none');
caxis([0,200])
axis('square')
title(['frame: ',num2str(nm(1)+1)])

% determine the location of profile
figure; hold on
Sx = abs(sum(dIR{nm(1)+1},1));
Sx_s = sgolay_t(Sx,3,21);
plot(Sx,'k')
plot(Sx_s,'k')
[nx,x] = peakseek(Sx_s,0.5*max(Sx));
plot(nx,x,'ko')

Sy = abs(sum(dIR{nm(1)+1},2));
Sy_s = sgolay_t(Sy,3,21);
plot(Sy,'r')
plot(Sy_s,'r')
[ny,y] = peakseek(Sy_s,30,0.4*max(Sy));
plot(ny,y,'ro')

ny_mean = mean(ny);
dny     = diff(ny);
ny1 = ny_mean-dny;
ny2 = ny_mean+dny;

nx_mean = nx;
dnx     = dny;
nx1 = nx_mean-dnx;
nx2 = nx_mean+dnx;

figure; 
contourf(dIR{nm(1)+1},20,'lineStyle','none');
ylim([ny1,ny2])
xlim([nx1,nx2])
caxis([0,200])
axis('square')

%%
close all
% figure
% for s = nm(1):60
%    RadialHeatProfile(:,s) = Data(:,359,s)-Data(:,359,20); 
%    RadialHeatProfile2D(:,:,s) = Data(ny1:ny2,nx1:nx2,s)-Data(ny1:ny2,nx1:nx2,nm(1)); 
%    plot(RadialHeatProfile(:,s));
%    ylim([0,8000])
%    drawnow;
%    title(num2str(s))
%    pause(0.15)
% end

figure;
for s = nm(1)-4:38
   surf(dIR{s}(ny1:ny2,nx1:nx2),'lineStyle','none');
   view([60,60])
   set(gca,'DataAspectRatio',[1 1 5])
   zlim([0,1000])
%    xlim([0,500])
%    ylim([0,500])
   drawnow;
   title(num2str(s))
   pause(0.05)
end
%%
close all
figure;
subplot(1,2,1);
dIR_t_4_2s = dIR{nm(1)+2}(ny1:ny2,nx1:nx2);
contourf(dIR_t_4_2s,20,'lineStyle','none');
caxis([0,1000])
axis('square')
% title(['frame: ',num2str(nm(1)+2)])
text(20,190,'t = 4.2 s','Color','w','FontName','times','FontSize',12)
colormap('hot')
ylim([30,200])
xlim([10,190])
set(gca,'YTick',[],'XTick',[])
ylabel('y','FontName','times','FontSize',12)
xlabel('x','FontName','times','FontSize',12)

subplot(1,2,2);
dIR_t_4_43s = dIR{38}(ny1:ny2,nx1:nx2);
contourf(dIR_t_4_43s,20,'lineStyle','none');
caxis([0,1000])
axis('square')
% title(['frame: ',num2str(38)])
text(20,190,'t = 4.43 s','Color','w','FontName','times','FontSize',12)
colormap('hot')
ylim([30,200])
xlim([10,190])
set(gca,'YTick',[],'XTick',[])
ylabel('y','FontName','times','FontSize',12)
xlabel('x','FontName','times','FontSize',12)

set(gcf,'color','w')
box on
set(gca,'FontName','times','FontSize',11)

figure;
dt = (seq.Duration.TotalSeconds)/double(seq.Count);
for s = nm(1):42
    RadialHeatFlux(:,s) = dIR{s}(ny1:ny2,nx_mean);
end
t_rhf = [ [1:length(RadialHeatFlux(1,:))].*dt - 0.47  + 4.15 ]';
r_rhf = [ 1:length(RadialHeatFlux(:,1)) ]';
contourf(t_rhf,r_rhf,RadialHeatFlux,20,'lineStyle','none');
line([4.2,4.2],[0,400],'color','w','linewidth',2,'linestyle','--')
line([4.434,4.434],[0,400],'color','w','linewidth',2,'linestyle','--')
colormap('hot')
caxis([0,1000])
xlim([4.12,4.50])
ylim([25,200])
box on
xlabel('t [s]','Interpreter','latex','FontSize',12)
ylabel('y [a.u]')
set(gca,'YTick',[])
set(gca,'FontName','times','FontSize',11,'PlotBoxAspectRatio',[2 1 1])
set(gcf,'color','w')

if 0
save('IR_data_13962_2017_04_13','dIR_t_4_2s','dIR_t_4_43s','t_rhf','r_rhf','RadialHeatFlux')
end


%% 2019_01_10
% here we are using the semi-infinite slab solution with surface heat flux
% to figure out what the heat flux is

T0 = Data(:,:,24);
for s = 24:size(Data,3)
     dT(:,:,s-23) = (Data(:,:,s) - T0)*1e-3;
     dT2(:,:,s-23) = dT(:,:,s-23).^2;
end

dt = 1/60;

for s = 1:size(dT2,3)-1
    ddT2(:,:,s) = (dT2(:,:,s+1) - dT2(:,:,s))*dt;
    SqddT2(:,:,s) = sqrt(ddT2(:,:,s));
end

 figure; 
%surf(ddT2(:,:,3),'LineStyle','none')
surf(real(SqddT2(:,:,24)),'LineStyle','none')