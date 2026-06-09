% Get IR data
clear all
close all

DownloadData = 1;
SaveData = 0;
DataName = 'GetRawData_IRheatFlux';

switch DownloadData
    case 1
 
    IR_ShotList = 24000 + [526,528,529,530,532,533];
    % Data is not available so will need to analyze myself with analytical
    % code
    PS1      =         [2.0,3.0,4.0,5.0,6.0,7.0];
    PS2      =         [2.0,2.0,2.0,2.0,2.0,2.0];
    
    mdsconnect('mpexserver')
for s = 1:length(IR_ShotList)
   [heat_flux{s},heat_peak{s},heat_peak_x{s}...
       ,heat_peak_y{s},heat_power{s}...
       ,heat_radius{s},T{s},X{s},Y{s}...
       ,comments{s}] = get_IR(IR_ShotList(s));
   
   RawData_IR.heat_flux{s} = heat_flux{s};
   RawData_IR.heat_peak{s} = heat_peak{s};
   RawData_IR.heat_peak_x{s} = heat_peak_x{s};
   RawData_IR.heat_peak_y{s} = heat_peak_y{s};
   RawData_IR.heat_power{s}  = heat_power{s};
   RawData_IR.T{s}  = T{s};
   RawData_IR.X{s}  = X{s};
   RawData_IR.Y{s}  = Y{s};
end
    
RawData_IR

if SaveData == 1
save(DataName,'RawData_IR')
end

case 0

end
% #########################################################################
% #########################################################################
% End of Gather Raw Data
% #########################################################################
% #########################################################################