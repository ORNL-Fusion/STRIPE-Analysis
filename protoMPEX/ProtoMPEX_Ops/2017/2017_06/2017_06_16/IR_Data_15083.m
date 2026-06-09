close all 
clear all

shot = 15083;
mov = 1;

UploadType = 2;
switch UploadType
    case 1
        a = '\\mpexserver\ProtoMPEX_Data\IR_Camera\2017_06_16';
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
%%
if mov
figure(1)
for s = 1:40
   imshow(Data(:,:,s),[]);
   drawnow;
   title(num2str(s))
   pause(0.5)
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
nm = find(m>max(m)/180);
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

Rloc = 352;
figure; hold on
for s = 2:60
   S_RadialHeatProfile(:,s) = Data(:,Rloc,s)-Data(:,Rloc,12); 
   RadialHeatProfile2D(:,:,s) = Data(:,:,s)-Data(:,:,s-1); 
   plot(S_RadialHeatProfile(:,s))
   ylim([0,1e4])
   drawnow;
   pause(0.2)
end

figure
surf(S_RadialHeatProfile,'LineStyle','none')

figure
for s = 1:30 % use 13 to 22
   surf(RadialHeatProfile2D(200:300,310:410,s),'lineStyle','none');
    zlim([0,9000])
    xlim([0,100])
    ylim([0,100])
    drawnow;
   title(num2str(s))
   pause(0.2)
end

figure;
subplot(1,2,1);
IR_15083_t_4_43 = RadialHeatProfile2D(153:381,252:480,21);
contourf(IR_15083_t_4_43,20,'LineStyle','none');

caxis([0,3000])
axis('square')
% title(['frame: ',num2str(nm(1)+2)])
text(20+10,190-10,'t = 4.43 s','Color','w','FontName','times','FontSize',12)
colormap('hot')
ylim([30,200]-10)
xlim([10,190]+10)
set(gca,'YTick',[],'XTick',[])
ylabel('y','FontName','times','FontSize',12)
xlabel('x','FontName','times','FontSize',12)
set(gcf,'Color','w')
box on

if 0
   save('IR_hollow_15083_2017_06_16','IR_15083_t_4_43')
end

%======================
return
%=====================

%%
close all
figure; 
surf(S_RadialHeatProfile(:,1:60),'lineStyle','none')

% Data(:,359,25)-Data(:,359,24)
% Data(:,:,25)-Data(:,:,15)
figure; hold on
for s = 1:45
    S_RPD(:,s) = sgolay_t(S_RadialHeatProfile(:,s),3,7);
    RPD2D(:,:,s) = sgolay_t(RadialHeatProfile2D(:,:,s),3,27);
    plot(S_RPD(:,s))
end
% xlim([150,350])

%%
figure; hold on
dt = 19.8543/1000;
t = [1:size(S_RPD,2)-1]*dt -0.45 + 4.15 ;
plot(t,diff(S_RPD(255,:)),'k.-')
plot(t,diff(S_RPD(228,:)),'r.-')

for s = 1:44
    ff(:,s) = S_RPD(:,s+1)-S_RPD(:,s);
    ff2D(:,:,s) = RPD2D(:,:,s+1)-RPD2D(:,:,s);
end

figure;
R = 1:length(ff(:,1));
%contourf(t(20:40),R,ff(:,20:40),20,'lineStyle','none')
%caxis([0,800])
% surf(t(20:40),R,ff(:,20:40),'lineStyle','none')
surf(ff(:,:),'lineStyle','none')
colormap('hot')
%%
figure
for s = 1:24;
   surf(ff2D(:,:,s),'lineStyle','none');
    view([-75,50])
%     view([0,90])
    zlim([0,1000])
    caxis([0,1000])
    colormap('hot')
    ylim([100,400])
    xlim([200,500])
    axis('square')
   drawnow;
   title(num2str(s))
   pause(0.1)
end

%%
% close all

f0 = 20;
f1 = 40;
% 26 begins, 28 looks good, 30 helicon has developed
% 36 is fully developed

% for 300 ms use frame 42 for fully developed heat,

rngy = [280:480] - 30;
rngx = [100:310]  + 60;
IR_f10 = Data(rngx,rngy,f1)-Data(rngx,rngy,f0);
IR_f10n = IR_f10/max(max(IR_f10));

a = 30e4;

% Delta emission
% absolute
figure;
contourf(IR_f10,80,'LineStyle','none')
caxis_max = 300;
caxis([0,caxis_max]*a*1e-4)
axis(gca,'equal')
colormap('jet')
th = title([FILENAME,', caxis: ',num2str(caxis_max),', f = ',num2str(f1)])
set(th,'FontSize',14)

% normalized
% figure;
% contourf(IR_f10n,80,'LineStyle','none')
% zlim([0,1])
% caxis([0,0.9])
% axis(gca,'equal')
% colormap('jet')
% title(FILENAME)


