%##### Load image #####
%close all 
clear all
tic
[FILENAME, PATHNAME, FILTERINDEX] = uigetfile('*.jpg;*.seq;*.ats', 'Choose IR file (jpg) or radiometric sequence (seq)');
videoFileName=[PATHNAME FILENAME];

shotnumber = FILENAME(1:end-4);
shot=double(string(FILENAME(6:end-4)));
%if sandia image: shot = double(string(FILENAME(7:end-4);

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
Counts=seq.ThermalImage.Count;
frames=double(Counts);
%e1_list=zero(frames,1);
if(seq.Count > 1)
    while(seq.Next())
        img = seq.ThermalImage.ImageProcessing.GetPixelsArray;
        im = double(img);
        Data(:,:,fr) = im;
        %e1_list(fr)=ellipseAvg;
        fr = fr + 1;
    end
end

delta_t= input('what is the frame rate of the shot [in seconds]  ');%pulse length (sec)
% target_position = input('how far behind is the plate from 11.5 (in cm)? '); %moving target plate
% if target_position == 6
%     px_per_cm = 22.55;
% end
% if target_position == 4
%     px_per_cm = 23.38;
% end
% if target_position == 3
%     px_per_cm = 23.48;
% end
% if target_position == 2
%     px_per_cm = 23.54;
% end
% if target_position == 1.75
%     px_per_cm = 23.58;
% end
% if target_position == 1.5
%     px_per_cm = 23.58;
% end

px_per_cm = input('what is the px_per_cm conversion ?  ');

powersource = input('What power sources were applied? helicon only=0; ech+hel=1; ich+hel=2; all=3   ');
fullpulse = input('Did additional power source(s) have full pulse? no=0; yes=1; ');
%frame_analysis = input('Frame subtraction from first frame analysis [0], frame subtraction from previous frame analysis [1], or both [2]? ');
%First - find the approximate pixel location of the hottest spot
%guess max frame - just want to be able to narrow in on the right max area
%to find the correct max frame and max delta T
f0 = 1;
firstplasma=0;
f1 = 55;

%zoom in and shorten the frames - all interesting stuff happens within the
%first 40 @ 50 Hz, 80 @ 100 Hz
zoomframes=65;

Delta_Data_zoom1 = Data(:,:,f1)- Data(:,:,f0);
figure;
image(Delta_Data_zoom1,'CDataMapping','scaled');
colormap jet
c=colorbar;
grid on

title('Click center of plasma profile')
[x_c,y_c] = ginput(1);
%ellipse [Position] = [Xmin Ymin Width Height] = [columnMin rowMin columns
%rows]
TC_diameter = 18;
TC_radius = TC_diameter/2;
center_plate=[x_c,y_c];

%zoom_rect = imrect(gca,[x_c-100 y_c 200 200]);
zoom_x=double(int16(x_c));
zoom_y = double(int16(y_c));

%try to find the edge peak
average_data=zeros(zoomframes,1);
for i=1:zoomframes
    average_data(i)= mean2(Data(zoom_y-80:zoom_y+80,zoom_x-80:zoom_x+80,i));
    i=i+1;
end

average_data_delta=zeros(zoomframes,1);
for i=2:zoomframes
    average_data_delta(i)= average_data(i)-average_data(i-1);
    i=i+1;
end

average_data_delta_list=zeros(zoomframes,1);
j=1;
for i=1:zoomframes-5
   if average_data_delta(i)>3
       average_data_delta_list(j) = i;
       i=i+1;
       j=j+1;
   end
end

[max5, plasma5] = max(average_data_delta);

firstplasma = average_data_delta_list(1);
centerplasma=firstplasma+6;

Delta_Data_zoom = Data(zoom_y-80:zoom_y+80,zoom_x-80:zoom_x+80,firstplasma+3) - Data(zoom_y-80:zoom_y+80,zoom_x-80:zoom_x+80,f0);
figure;
image(Delta_Data_zoom,'CDataMapping','scaled');
TC_circle = viscircles([x_c y_c],TC_radius);
centerline_horiz = imline(gca, [x_c-5 y_c; x_c+5 y_c]);
centerline_vert = imline(gca, [x_c y_c+5;x_c y_c-5]);
colormap jet
c=colorbar;
grid on

title('Click on approximate edge hot spot position on SiC plate')
[column,row] = ginput(1);
drawnow;

Delta_Data_zoom2 = Data(zoom_y-80:zoom_y+80,zoom_x-80:zoom_x+80,firstplasma+10) - Data(zoom_y-80:zoom_y+80,zoom_x-80:zoom_x+80,f0);
figure;
image(Delta_Data_zoom2,'CDataMapping','scaled');
TC_circle = viscircles([x_c y_c],TC_radius);
centerline_horiz = imline(gca, [x_c-5 y_c; x_c+5 y_c]);
centerline_vert = imline(gca, [x_c y_c+5;x_c y_c-5]);
colormap jet
c=colorbar;
grid on

title('Click on approximate center hot spot position on SiC plate')
[center1,center2] = ginput(1);
drawnow;

title('Click approximate diameter width')
[x1,y1] = ginput(1);
[x2,y2] = ginput(1);
diameter_raw = sqrt((x2-x1)^2+(y2-y1)^2);
drawnow;

column = double(int16(column));
row = double(int16(row));
diameter = double(int16(diameter_raw));
center1 = double(int16(center1));
center2 = double(int16(center2));

%zoom in and shorten the frames - all interesting stuff happens within the
%first 60 frames
%zoomframes=40;
templist=zeros(zoomframes,1);
templist2=zeros(zoomframes,1);
TG_templist=zeros(zoomframes,1);
hel_templist=zeros(zoomframes,1);
Data_zoom = Data(zoom_y-80:zoom_y+80,zoom_x-80:zoom_x+80,1:zoomframes);
temp=zeros(zoomframes,1);
tempA=zeros(zoomframes,1);
tempB=zeros(zoomframes,1);
tempC=zeros(zoomframes,1);
tempD=zeros(zoomframes,1);
tempE=zeros(zoomframes,1);
tempF=zeros(zoomframes,1);
tempG=zeros(zoomframes,1);
tempH=zeros(zoomframes,1);
temp2=zeros(zoomframes,1);
Temperature = Data_zoom;
DeltaTMatrix=Data_zoom;
TG_templist_step=zeros(zoomframes,1);
hel_templist_step=zeros(zoomframes,1);
DeltaTMatrix_step = Data_zoom(:,:,1:zoomframes);
DeltaData = Data_zoom;
DeltaTemperature = Temperature;

% for i = 1:zoomframes
%     DeltaData(:,:,i)=Data_zoom(:,:,i)-Data_zoom(:,:,1);
%     i=i+1;
% end
% 
% for i = 1:zoomframes
%     DeltaTemperature(:,:,i) = arrayfun(@(DeltaData) seq.ThermalImage.GetValueFromEmissivity(0.69,DeltaData),DeltaData(:,:,i));
%     i=i+1;
% end

for i=1:zoomframes
    %grabs temperature and delta T from a zoomed in location on the
        %image - takes too long if keep full image - need to keep the delta
        %in the y-range 190, and the x-range 160 pixels. - Data zoom in @
        %center of the image
     
     Temperature(:,:,i)=arrayfun(@(Data_zoom) seq.ThermalImage.GetValueFromEmissivity(0.7,Data_zoom),Data_zoom(:,:,i));
     temp(i)=Temperature(row,column,i);
     tempA(i)=Temperature(row,column-2,i);
     tempB(i)=Temperature(row,column+2,i);
     tempC(i)=Temperature(row-2,column,i);
     tempD(i)=Temperature(row+2,column,i);
     
     temp2(i)=Temperature(center2,center1,i);
     tempE(i)=Temperature(center2,center1-2,i);
     tempF(i)=Temperature(center2,center1+2,i);
     tempG(i)=Temperature(center2+2,center1,i);
     tempH(i)=Temperature(center2-2,center1,i);
     
     temparray=[temp(i);tempA(i);tempB(i);tempC(i);tempD(i)];
     avgtemp=mean(temparray);
     templist(i)=avgtemp;
     TG_templist(i)=templist(i)-templist(1);
     
     temparray2=[temp2(i);tempE(i);tempF(i);tempG(i);tempH(i)];
     avgtemp2=mean(temparray2);
     templist2(i)=avgtemp2;
     hel_templist(i)=templist2(i)-templist2(1);
        %delta T general is by subtracting first frame from remaining frames. This
        %can be modified if needed for other image subtraction
     DeltaTMatrix(:,:,i)=Temperature(:,:,i)-Temperature(:,:,1);
     
     i=i+1;
end

%add analysis for the adjacent frame subtraction. 
TG_templist_step(1)=0;
hel_templist_step(1)=0;
DeltaTMatrix_step(:,:,1)=0;
for i=2:zoomframes
     TG_templist_step(i)=templist(i)-templist(i-1);
     hel_templist_step(i)=templist2(i)-templist2(i-1);
     DeltaTMatrix_step(:,:,i) = Temperature(:,:,i)-Temperature(:,:,i-1);
     
     i=i+1;
end

%%
%close all

%% Find the hottest frame ; or can also default for frame 36-37 
MaxDelta=DeltaTMatrix(:,:,1);
MaxDelta_step=DeltaTMatrix_step(:,:,1);

%shifted center because of zooming
x_c2 = 80;
y_c2 = 80;

%centerline_horiz coordinates: [x_c2-5 y_c2-87; x_c2+5 y_c2-87]);
%centerline_vert coordinates: [x_c2 y_c2-82;x_c2 y_c2-92]);

%create list of averages per frame:
meanlist=zeros(zoomframes,1);
meanlist_step=zeros(zoomframes,1);

for i=1:zoomframes
    meanlist(i)= mean2(DeltaTMatrix(y_c2-58:y_c2+68,x_c2-63:x_c2+63,i));
    meanlist_step(i) = mean2(DeltaTMatrix_step(y_c2-58:y_c2+68,x_c2-63:x_c2+63,i));
end

%x-range (num columns - second), y-range (num rows -first)
%first attempt at narrowing into max frame: meanlist frame average
% [max_delta_step,framenumber_step] = max(meanlist_step);
% framenumber1=framenumber_step+1; %max frame generally comes one frame after the max delta_step frame
% framenumber=framenumber1;
% %check that it's close enough to the max.- should be overall max - avg max 
% for i=framenumber1:zoomframes-1
%     if meanlist_step(i)>0.1
%         framenumber=framenumber+1;
%     end
%     i=i+1;
% end

[max_delta_avg,framenumber] = max(meanlist);

%find max TG frame
[TG_max,TG_frame]=max(TG_templist);
% [TG_max1,TG_frame1]=max(TG_templist_step);
% TG_frame_check = TG_frame1+1; % should be within one of each other.
% if abs(TG_frame_check-TG_frame)>3
%     if powersource==0 || fullpulse==0
%         TG_frame_check=TG_frame;
%     else
%     disp,'gap between max frames is large: manual analysis suggested' 
%     end
% end

%find helicon max = should match avg max. 
[hel_max,hel_frame]=max(hel_templist);
[hel_max1,hel_frame1]=max(hel_templist_step);
hel_frame_check = hel_frame1+1; % should be within one of each other.
% if hel_frame_check-hel_frame>2
%     if powersource~=0 %addition of the ECH causes a dip in the center heating as it ramps up
%         hel_frame_check=TG_frame;
%         if abs(hel_frame_check-hel_frame)<=2
%             if hel_templist_step(hel_frame)<0 %if the deriv of hel_templist is less than zero, then the heat is decreasing
%                 hel_frame=hel_frame-1;
%             end
%         else
%             disp,'gap between max frames is large: manual analysis suggested'
%         end
%     else
%         disp,'gap between max frames is large: manual analysis suggested'
%     end
% end

hel_max=hel_templist(hel_frame); %update max as necessary

%plot both possibilities for troubleshooting if necessary + general
%analysis 
figure;
plot(1:zoomframes, TG_templist_step,'-b.',1:zoomframes, TG_templist,'-b.',1:zoomframes,hel_templist_step,'-k.', 1:zoomframes,hel_templist,'-k.',1:zoomframes,meanlist,'-m.',1:zoomframes,meanlist_step,'-m.');
grid on
legend('TG Frame from Frame', 'TG Frame from Frame 1','Helicon Frame from Frame','Helicon Frame from Frame 1','Average Frame from Frame1','Average Frame from Frame'); 
legend('Location','Northwest');

finalframenumber=0;

if framenumber-hel_frame>=3
    if abs(meanlist(hel_frame)-meanlist(framenumber))<0.2
        finalframenumber=hel_frame;
    else
    %'error: gap between two possible max delta T frames is large - manual analysis suggested'
    finalframenumber=framenumber;
    end
end

if hel_frame-framenumber>=3
    if abs(meanlist(hel_frame)-meanlist(framenumber))<0.2
        finalframenumber=hel_frame;
    else
    %'error: gap between two possible max delta T frames is large - manual analysis suggested'
    finalframenumber=framenumber;
    end
end

if abs(hel_frame-framenumber)==2
    finalframenumber=(hel_frame+framenumber)/2;
end

if abs(hel_frame-framenumber)<2
    %if abs(hel_templist(hel_frame)-hel_templist(framenumber))>abs(meanlist(hel_frame)-meanlist(framenumber))
    if meanlist_step(framenumber)<hel_templist_step(hel_frame)
        finalframenumber=hel_frame;
    else finalframenumber=framenumber;
    end
end

if TG_templist(TG_frame)<TG_templist(finalframenumber)
    TG_frame=finalframenumber;
end
    
MaxDelta=DeltaTMatrix(:,:,finalframenumber);


deltamax = double(max(max(int16(round(MaxDelta,2,'significant')))));
deltamin = double(min(min(int8(round(MaxDelta,1,'significant')))));
scalestep = deltamax/10;

%% Heat Flux Measurement
%Calculates the heat flux of a image 

% Cp=500; %J/kg*C
% density=8030; %kg/m3
% thickness= .06*(2.54/100); %inches to m

delta_tframe = delta_t;%gap between frames


%Plot MaxDelta
figure;
figDeltaT=imagesc(MaxDelta, 'CDataMapping','scaled');
caxis([0 deltamax])
colormap jet
c=colorbar;
ylabel(c, 'Delta T [C]')
ax.FontSize = 13;
title([shotnumber,'; Delta Temperature vs. Camera Pixels'],'FontSize',13);
xlabel('Pixels','FontSize',13);
ylabel('Pixels','FontSize',13);

TC_circle = viscircles([x_c2 y_c2],TC_radius);
centerline_horiz = imline(gca, [x_c2-5 y_c2; x_c2+5 y_c2]);
centerline_vert = imline(gca, [x_c2 y_c2-5;x_c2 y_c2+5]);

%find average ellipse delta temp
%diameter = input('what is the diameter of the plasma? [in pixels]  ');
large_diameter = round(diameter*1.07,3,'sig');
r1=diameter/2;
r2=large_diameter/2;
%ellipse [Position] = [Xmin Ymin Width Height] = [columnMin rowMin columns
%rows] 
%Approximate 'E1' from ResearchIR
ellipse=imellipse(gca,[80-r2 80-r1 diameter large_diameter]);
mask=createMask(ellipse,figDeltaT);
avg_all=mean2(mask.*MaxDelta);
avg_mask=mean2(mask);
avg_delta=avg_all/avg_mask; %this is hand-wavy? 

%Approximate 'E2' from ResearchIR
ellipse2=imellipse(gca,[80-63 80-63 63*2 63*2]);
mask2=createMask(ellipse2,figDeltaT);
avg_all2=mean2(mask.*MaxDelta);
avg_mask2=mean2(mask2);
avg_delta2=avg_all2/avg_mask2; %this is hand-wavy?

%Draw hot spot for deltaT
ellipse3=imellipse(gca,[column-2 row-2 4 4]);

%Draw center spot ellipse deltaT
ellipse4 = imellipse(gca,[center1-2 center2-2 4 4]);
rect_mean = imrect(gca,[x_c2-63 y_c2-63 63*2 63*2]);


%Start with ECH/ICH intra-shot analysis
%find shot where non-preionization plasma first appears on the target plate 
plasmaframe_list=zeros(finalframenumber,1);
j=1;
for i=2:finalframenumber-2 %looks for first helicon/ECH heating
    if meanlist_step(i)>0.1
        plasmaframe_list(j)=i;
        j=j+1;
        i=i+1;
    end
end
%plasma frame where plasma first appears
plasmaframe = plasmaframe_list(1);

%in case range error - plasma too cool 
if plasmaframe==0
    plasmaframe=firstplasma;
end

% 
%Plot Frame with first helicon plasma
figure;
imagesc(DeltaTMatrix(:,:,plasmaframe), 'CDataMapping','scaled')
caxis([0 5])
colormap jet
c=colorbar;
ylabel(c, 'Delta T [C]')
ax.FontSize = 13;
title('First Helicon Plasma Frame','FontSize',13);
xlabel('Pixels','FontSize',13);
ylabel('Pixels','FontSize',13);


%ECH and ICH analysis
%Analyze the Temperature matrix to do the image subtraction to determine
%ECH/ICH - use input power source
%to skip this analysis - end the power source as helicon only [0]
%helicon_start = input('when did the helicon turn on after t=+4 sec [sec]?  '); %t+4 sec pulse start for helicon. - almost always constant around 4.157

%This portion of the code will be modified once the matlab code pulls from
%the MPEX tree/jscope
if powersource == 0
    helicon_start = heliconsourcetimes(shot);
end
if powersource == 1
    [helicon_start, ech_start_time,ech_end_time] = powersourcetimes(shot);
%     helicon_start = input('when did the helicon turn on after t=+4 sec [sec]?  '); %t+4 sec pulse start for helicon. - almost always constant around 4.157
%     ech_start_time = input('when did the ech turn on after t=+4 sec [sec]?  ');
%     ech_end_time = input('when did the ech turn off after t=+4 sec [sec]?  ');
    
    ech_start_frame = int16(((ech_start_time-helicon_start)/0.02))+plasmaframe;
    ech_end_frame = int16(((ech_end_time-helicon_start)/0.02))+plasmaframe;
    
    %delta of the ECH pulse
    figure;
    Delta_ECH = DeltaTMatrix(:,:,ech_end_frame) - DeltaTMatrix(:,:,ech_start_frame);
    delta_ech = double(max(max(int16(round(Delta_ECH,2,'significant')))));
    figDelta_ECH = imagesc(Delta_ECH, 'CDataMapping','scaled');
    caxis([0 delta_ech])
    colormap jet
    c=colorbar;
    ylabel(c, 'Delta T [C]')
    ax.FontSize = 13;
    title([shotnumber,'; Delta T from ECH pulse'],'FontSize',13);
    xlabel('Pixels','FontSize',13);
    ylabel('Pixels','FontSize',13);
        
    %ellipse [Position] = [Xmin Ymin Width Height] = [columnMin rowMin columns
    %rows] 
    %Approximate 'E1' from ResearchIR
    ellipseECH=imellipse(gca,[98-r2 120-r1 diameter large_diameter]);
    maskECH=createMask(ellipseECH,figDelta_ECH);
    avg_all_ECH=mean2(maskECH.*Delta_ECH);
    avg_mask_ECH=mean2(maskECH);
    avg_delta_ECH=avg_all_ECH/avg_mask_ECH;

    %Approximate 'E2' from ResearchIR
    ellipseECH2=imellipse(gca,[100-68 120-(68) 68*2 68*2]);
    maskECH2=createMask(ellipseECH2,figDelta_ECH);
    avg_all_ECH2=mean2(maskECH2.*Delta_ECH);
    avg_mask_ECH2=mean2(maskECH2);
    avg_delta_ECH2=avg_all_ECH2/avg_mask_ECH2; 
    
    mask_ellipse3=createMask(ellipse3,figDelta_ECH);
    avg_ellipse3=mean2(mask_ellipse3.*Delta_ECH);
    avg_mask_ellipse3=mean2(mask_ellipse3);
    avg_delta_ellipse3=avg_ellipse3/avg_mask_ellipse3;
    
    mask_ellipse4=createMask(ellipse4,figDelta_ECH);
    avg_ellipse4=mean2(mask_ellipse4.*Delta_ECH);
    avg_mask_ellipse4=mean2(mask_ellipse4);
    avg_delta_ellipse4=avg_ellipse4/avg_mask_ellipse4;
    
    %draw figures
    %Draw hot spot for deltaT
    ellipse3=imellipse(gca,[column-2 row-2 4 4]);
    %Draw center spot ellipse deltaT
    ellipse4 = imellipse(gca,[center1-2 center2-2 4 4]);
    %Draw center grids
    TC_circle = viscircles([x_c2 y_c2],TC_radius);
    centerline_horiz = imline(gca, [x_c2-5 y_c2; x_c2+5 y_c2]);
    centerline_vert = imline(gca, [x_c2 y_c2+5;x_c2 y_c2-5]);
end



%PLOTS!
%Delta T plot
figure;
plot(1:zoomframes,TG_templist_step,'-k.',1:zoomframes,hel_templist_step,'-b.',1:zoomframes,meanlist_step,'-m.')
ax.FontSize = 13;
title([shotnumber,'; Delta T between Frames'],'FontSize',13);
xlabel('Frames','FontSize',13);
ylabel('Delta T (deg C)','FontSize',13);
legend('TG','Helicon','Average');  
legend('Location','Northwest')

% %heat density plots
% TG_powerdensity_step = ((Cp*TG_templist_step*thickness*density)/delta_tframe)/1E6; %MW/m2
% hel_powerdensity_step = ((Cp*hel_templist_step*thickness*density)/delta_tframe)/1E6; %MW/m2
% avg_powerdensity_step = ((Cp*meanlist_step*thickness*density)/delta_tframe)/1E6; %MW/m2

%areas
px_per_m=px_per_cm*100;
hel_radius = 2/px_per_m;%pixels per m.
area_hel = pi*(hel_radius^2); %same as area_TG;
TG_radius = 2/px_per_m;
area_TG = pi*(TG_radius^2);
avg_side = (68*2)/px_per_m;
avg_area = avg_side^2;

% %power plots
% TG_power_step = TG_powerdensity_step*area_TG*1000;  %kW
% hel_power_step = hel_powerdensity_step*area_hel*1000; %kW
% avg_power_step = avg_powerdensity_step*avg_area*1000; %kW
% 
% %find max instantaneous heat flux
% [hel_max_step, hel_frame1] = max(hel_powerdensity_step);
% [TG_max_step, TG_frame1] = max(TG_powerdensity_step);
% [avg_max_step, avg_frame] = max(avg_powerdensity_step);
% 
% max_all = max([hel_max_step,hel_frame1;TG_max_step,TG_frame1;avg_max_step, avg_frame]);
% 
% DeltaMax_step_peak = DeltaTMatrix_step(:,:,max_all(2));
% PowerDensity_max = (Cp*DeltaMax_step_peak*thickness*density)/delta_tframe/1E6;


%find time step
frame2time = zeros(zoomframes,1);
    for i = 1:plasmaframe
        frame2time(i) = helicon_start-(plasmaframe-i)*delta_t;
        i=i+1;
    end
    for i = plasmaframe:zoomframes
        frame2time(i) = helicon_start+(i-plasmaframe)*delta_t;
        i=i+1;
    end

%outputs of valueable numbers
'Frame where helicon plasma first appears: ', disp(plasmaframe)
'Frame where overall deltaT is greatest', disp(finalframenumber)
'Average edge hot spot deltaT', disp(TG_templist(TG_frame))
'Average center hot spot deltaT', disp(hel_templist(hel_frame))
'Average deltaT (E1)', disp(avg_delta)
'Expanded average deltaT (E2)', disp(avg_delta2)

if powersource ==1
    'ECH edge hot spot deltaT', disp(avg_delta_ellipse3)
    'ECH center hot spot deltaT', disp(avg_delta_ellipse4)
    'ECH Average deltaT (E1)', disp(avg_delta_ECH)
    'ECH Expanded average deltaT (E2)', disp(avg_delta_ECH2)
    'ECH start and stop frames', disp(double(ech_start_frame)), disp(double(ech_end_frame))  
end
if powersource ==2
    'ICH edge hot spot deltaT', disp(avg_delta_ellipse3)
    'ICH center hot spot deltaT', disp(avg_delta_ellipse4)
    'ICH Average deltaT (E1)', disp(avg_delta_ICH)
    'ICH Expanded average deltaT (E2)', disp(avg_delta_ICH2)
    'ICH start and stop frames', disp(double(ich_start_frame)), disp(double(ich_end_frame))  
end
if powersource ==3
    'ECH+ICH edge hot spot deltaT', disp(avg_delta_ellipse3)
    'ECH+ICH center hot spot deltaT', disp(avg_delta_ellipse4)
    'ECH+ICH Average deltaT (E1)', disp(avg_delta_ECH_ICH)
    'ECH+ICH Expanded average deltaT (E2)', disp(avg_delta_ECH_ICH2)
    'ECH+ICH start and stop frames', disp(double(power_start_frame)), disp(double(power_end_frame))  
end

%save raw data of hottest frame as .mat file so can use it for the field line
%mapping
filename=strcat(FILENAME(6:end-4),'.mat');
% rawdata=Data(:,:,finalframenumber);
% rawdata_matrix=cell2mat(num2cell(rawdata));
%get the full image of the hottest frame in terms of delta T.
% emax=cell2mat(num2cell(EMax));
Frame=cell2mat(num2cell(MaxDelta));
save(filename,'DeltaTemperature','frame2time','Frame','delta_t','shotnumber','zoomframes','TG_frame','hel_frame','finalframenumber','TC_radius','px_per_cm','helicon_start','center_plate','TG_templist','hel_templist','TG_templist_step','hel_templist_step','DeltaTMatrix_step','meanlist_step','meanlist','area_TG','avg_area','x_c2','y_c2','diameter', 'plasmaframe','Temperature');

timeall=(toc)/60  