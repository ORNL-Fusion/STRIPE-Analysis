% Created 2017_08_14, JF CANESES
close all
clear all

% 16 seconds to load movies
% > 60 seconds to load and analyse data: 520 MB
% 20 seconds to load data: movies and analyzed data
%        6 seconds to load two variables

Retrieve = 0;

if Retrieve
shotlist = 13900 + [62,66];
% #########################################################################
% Loading the visible camera video
% s = 1; 13962: ECH OFF
% s = 2; 13966: ECH ON
% #########################################################################
%  Address = 'C:\Users\nfc\Documents\Proto-MPEX Data Analysis\2017\2017_04\2017_04_13\Highly ionized plasma column\';
Address = cd;
for s = 1:length(shotlist)
    FileName{s} = ['\slomo_c_',num2str(shotlist(s))];
    f{s} = VideoReader([Address,FileName{s},'.mov']);
    vidHeight = f{s}.Height;
    vidWidth  = f{s}.Width;
    S{s} = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap',[]);
    n = 1;
    while hasFrame(f{s})
        S{s}(n).cdata = readFrame(f{s});
        n = n+1;
    end
end

%##########################################################################
% Analyse Video: OnAxisEmission and RadialEmission
% #########################################################################
VidFrameOffset = 163 - 4165;
t = [[1:1:500]-VidFrameOffset]; % 1 ms per frame, since frame rate is 500 fr/s and 500 ms movie
for k = 1:length(shotlist)
for s = 1:1:500
       R{k}(:,:,s) = S{k}(s).cdata(:,:,1)/3;
       G{k}(:,:,s) = S{k}(s).cdata(:,:,2)/3;
       B{k}(:,:,s) = S{k}(s).cdata(:,:,3)/3;
       fr{k}(:,:,s) = (R{k}(:,:,s)+ G{k}(:,:,s) + B{k}(:,:,s));
   IntensityIntegrated{k}(s) = sum(sum(fr{k}(:,:,s),1),2);
   x1 = 440; x2 = 460;
   y1 = 5  ; y2 = 688-y1;
   rp{k}(:,:,s) = fr{k}(y1:y2,x1:x2,s); 
   RadialEmission{k}(:,s) = mean(rp{k}(:,:,s),2);
   OnAxisEmission{k}(:,s) = mean(mean(rp{k}(280:320,:,s),2),1);
end
end
    save('Figure6_data')
else
    load('Figure6_data','S','RadialEmission','t')
    x1 = 440; x2 = 460;
   y1 = 5  ; y2 = 688-y1;
end

%% Figure 6
close all
t1 = 4.177;
t2 = 4.43;
rng = find(t*1e-3 >= t1 & t*1e-3<=t2);

figure
fvid = gcf;
set(fvid,'Menubar','figure','color','w','Units','normalized')
for nax = 1:4
    axvid(nax) = axes;
    set(axvid(nax),'Units','Normalized','box','on')%,'NextPlot','add')
end
% a), b), c) , d) positions
x_abcd = -0.15;
y_abcd = 0.96;
abcd.Color = 'k';
abcd.FontSize = 15;
FigTitle.Size = 12;
FigTitle.x0 = 0;
FigTitle.y0 = 1.08

s = 1; % Shot
nax = 1; % axis number
Intensity_HeliconOnly_start = S{s}(rng(1)).cdata; %rgb2gray
imh(nax) = image(axvid(nax),Intensity_HeliconOnly_start);   %axvid(1) = gca;
text(axvid(nax),x_abcd,y_abcd,'a)','Units','normalized','color',abcd.Color,'FontSize',abcd.FontSize,'FontName','Times')
text(axvid(nax),FigTitle.x0,FigTitle.y0,'Start helicon pulse','Units','normalized','color','k','FontSize',FigTitle.Size,'FontName','Times')
text(axvid(nax),20,630,'t = 4.17 s','Color','w','FontName','times','FontSize',11,'FontWeight','bold')
line(axvid(nax),[x1,x1],[y1,y2],'Color','w','LineWidth',1,'LineStyle',':')
xlabel(axvid(nax),'z','interpreter','Latex','FontSize',13)
ylabel(axvid(nax),'y','interpreter','Latex','FontSize',13)
s = 1;
nax = 2;
Intensity_HeliconOnly_end = S{s}(rng(end)).cdata;
imh(nax) = image(axvid(nax),Intensity_HeliconOnly_end); hold on %axvid(2) = gca;
text(axvid(nax),x_abcd,y_abcd,'b)','Units','normalized','color',abcd.Color,'FontSize',abcd.FontSize,'FontName','Times')
text(axvid(nax),FigTitle.x0,FigTitle.y0,'End helicon pulse','Units','normalized','color','k','FontSize',FigTitle.Size,'FontName','Times')
text(axvid(nax),20,630,'t = 4.43 s','Color','w','FontName','times','FontSize',11,'FontWeight','bold')
line(axvid(nax),[x1,x1],[y1,y2],'Color','w','LineWidth',1,'LineStyle',':')
% xlabel(axvid(nax),'z','interpreter','Latex','FontSize',13)
ylabel(axvid(nax),'y','interpreter','Latex','FontSize',13)
s = 2;
nax = 3;
Intensity_HeliconOnly_end = S{s}(rng(end)).cdata;
imh(nax) = image(axvid(nax),Intensity_HeliconOnly_end); hold on %axvid(2) = gca;
text(axvid(nax),x_abcd,y_abcd,'c)','Units','normalized','color',abcd.Color,'FontSize',abcd.FontSize,'FontName','Times')
text(axvid(nax),FigTitle.x0,FigTitle.y0,'End helicon pulse + 28 GHz','Units','normalized','color','k','FontSize',FigTitle.Size,'FontName','Times')
text(axvid(nax),20,630,'t = 4.43 s','Color','w','FontName','times','FontSize',11,'FontWeight','bold')
line(axvid(nax),[x1,x1],[y1,y2],'Color','w','LineWidth',1,'LineStyle',':')
xlabel(axvid(nax),'z','interpreter','Latex','FontSize',13)
ylabel(axvid(nax),'y','interpreter','Latex','FontSize',13)

s = 1;
nax = 4;
dn = 10;
RE1 = mean(RadialEmission{s}(:,rng(1:dn)),2);
plot(RE1(end:-1:1),1:length(RE1),'Parent',axvid(nax),'color','k','LineWidth',2);
text(axvid(nax),0.645,0.5,'a)','Units','normalized','color','k','FontSize',12,'FontName','Times')

n2 = length(rng);n1 = n2-dn;
RE2 = mean(RadialEmission{s}(:,rng(n1:n2)),2);
plot(RE2(end:-1:1),1:length(RE2),'Parent',axvid(nax),'color','k','LineWidth',1,'LineStyle','-');
text(axvid(nax),0.38,0.5,'b)','Units','normalized','color','k','FontSize',10,'FontName','Times')

RE2_ECH = mean(RadialEmission{s+1}(:,rng(n1:n2)),2);
text(axvid(nax),0.1,0.5,'c)','Units','normalized','color','k','FontSize',12,'FontName','Times')
plot(RE2_ECH(end:-1:1),1:length(RE2_ECH),'Parent',axvid(nax),'color','k','LineWidth',2,'LineStyle',':');
ylim([-10,688-10])
ylabel(axvid(nax),'y','interpreter','Latex','FontSize',13)
xlabel(axvid(nax),'Line-integrated Intensity','FontWeight','bold','interpreter','Latex','FontSize',12)
text(axvid(nax),x_abcd,y_abcd,'d)','Units','normalized','color',abcd.Color,'FontSize',abcd.FontSize,'FontName','Times')

dx1 = 0.15;dx2 = 0.000015;
dy1 = 0.10;dy2 = 0.09;
nCol = 2;
nRow = 2;
w = (1- (2*dx1) - dx2)/nCol;
h = (1- (2*dy1) - dy2)/nRow;


axvid(4).Position = [dx1       dy1       w h];
axvid(3).Position = [dx1+w+dx2 dy1       w h];
axvid(2).Position = [dx1+w+dx2 dy1+h+dy2 w h];
axvid(1).Position = [dx1       dy1+h+dy2 w h];

axis(axvid,'square')
set(axvid,'XTick',[],'YTick',[])
