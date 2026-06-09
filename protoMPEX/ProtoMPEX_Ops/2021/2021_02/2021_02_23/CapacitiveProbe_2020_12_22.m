% Reading the capacitive probe data:
clear all
close all

% shotlist = 30000 + [224 ,225,227,229,230,231,232,233,234,235,238,236,239 ,240 ,241,244,245,246,247 ,248 ,249 ,250 ];
% rp       =         [10.0,5.0,4.5,4.0,3.5,3.0,2.5,2.0,1.5,1.0,0.5,0.0,3.75,3.25,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0];

shotlist = 31000 + [392,393];
rp       =         [3.5,3.5];

shotlist = 31000 + [437,439];
rp       =         [3.5,3.5];

dataAddress{1} = ['\MPEX::TOP.MACHOPS1:BDOT_10']; % 
dataAddress{2} = ['\MPEX::TOP.MACHOPS1:RF_FWD_PWR']; % 

for ii = 1:numel(shotlist)
    [f,t_f] = my_mdsvalue_v2(shotlist(ii),dataAddress(1));
    v{ii} = f{1};
    vrms{ii} = 40.95.*(-v{ii}).^(0.503);
    t_v{ii} = t_f{1};
    
    [f,t_f] = my_mdsvalue_v2(shotlist(ii),dataAddress(2));
    rf{ii} = f{1};
    t_rf{ii} = t_f{1};
end

figure('color','w')
hold on
for ii = 1:numel(shotlist)
   y = sgolay_t(vrms{ii},3,201);
   plot(t_v{ii}(1:end-1),y);
   hf(ii) = plot(t_v{ii}(1:end-1),vrms{ii});
   legendText{ii} = num2str(shotlist((ii)));
end
hrf = plot(t_rf{1}(1:end-1),40*rf{1},'k');
legend(hf,legendText)
xlim([4,5])
% ylim([-1,2])
box on
ylabel('RF potential [Vrms]')
xlabel('time [s]')

figure('color','w')
hold on
for ii = 1:numel(shotlist)
   y = sgolay_t(vrms{ii},3,201);
   hf(ii) = plot3(t_v{ii}(1:end-1),rp(ii)*ones(size(y)),y);
   legendText{ii} = num2str(shotlist((ii)));
end
hrf = plot3(t_rf{1}(1:end-1),rp(1)*ones(size(rf{1})),40*rf{1},'k');
legend(hf,legendText)
xlim([4,5])
zlim([-10,60])
box on
zlabel('RF potential [Vrms]')
xlabel('time [s]')
view([30,30])