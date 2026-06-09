% definition of the RF power trace during 2017_12_03 experiments:

clear all 
close all

tn = [0,  3,102,104,153,155,204,206,254,255];
An = [0,127,127,102,102,76 ,76 ,51 , 51,  0];

tn = [0,  3,102,104,183,185,204,254,255];
An = [0,127,127,102,102,76 ,76 ,76 ,  0];

figure; 
h(1) = line(tn,An);
h(1).Color = 'red';
h(1).LineWidth = 2;
xlim([0,255])
set(gcf,'color','w')
box on

if 0
    figure; 
    h(1) = line(tn*500/255,An*100/127)
    xlim([0,500])
end