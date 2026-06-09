% Programming the Wavetek SG
clear all
close all

Ymax = 127;
Ymin = 0  ;

Xmin = 0  ;
Xmax = 255;

Ical = 382.6;

CurrentProgram = 5;

switch CurrentProgram
    case 1

    case 2
        I = [0,280 ,280 ,340  ,340  ,0];
        t = [0,0.5 ,4.2 ,4.3  ,5.5  ,6];
    case 3
        I = [0,240 ,240 ,400  ,400  ,0];
        t = [0,0.5 ,4.2 ,4.4  ,5.5  ,6];
    case 4
        I = [0,240 ,240 ,450  ,450  ,0];
        t = [0,0.5 ,4.2 ,4.4  ,5.5  ,6];
    case 5
        I = [0,240 ,240 ,370  ,370  ,0];
        t = [0,0.5 ,4.2 ,4.4  ,5.5  ,6];
end

% Convert the desired current to voltge request level:
V = I/Ical;

% Convert request voltage and time to bits:
Y = V*(Ymax/max(V));
X = t*(Xmax/max(t));

% Plot data
figure
subplot(2,2,1)
plot(t,I,'ko-')
box on
title('Current vs time')
grid on

subplot(2,2,2)
plot(t,V,'ko-')
box on
title('Request voltage vs time')
grid on

subplot(2,2,3)
plot(X,Y,'ko-')
title('7 bit amplitude vs address')
xlim([Xmin,Xmax])
grid on

set(gcf,'color','w','Position',[600  60  500  550])

