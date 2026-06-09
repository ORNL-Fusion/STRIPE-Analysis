% This script searches through the folders contained in "ProtoMPEX Data
% analysis" and produces an excell spreadsheet with dates and topics for
% all the experiments performed

close all
clear all

a = dir(cd);
x = 0;
p = 0;
for sx = 3:length(a);
    if strcmp(a(sx).name(1:3),'201') | strcmp(a(sx).name(1:3),'202')
        x = x + 1;
        Year{x} = a(sx).name;
        b = dir([cd,'\',Year{x}]);
        y = 0;
        for sy = 3:length(b);
           if strcmp(b(sy).name(1:3),'201') | strcmp(b(sy).name(1:3),'202')
              y = y + 1;
              Month{x}{y} = b(sy).name;
                c = dir([cd,'\',Year{x},'\',Month{x}{y}]);
                z = 0;
                for sz = 3:length(c);
                    if strcmp(c(sz).name(1:3),'201') | strcmp(c(sz).name(1:3),'202')
                       if  length(c(sz).name)>10
                           continue
                       end
                       p = p + 1;
                       z = z + 1;
                       Day{x}{y}{z} = c(sz).name;
                       Address{x}{y}{z} = [cd,'\',Year{x},'\',Month{x}{y},'\',Day{x}{y}{z}];
                       Date{p} = Day{x}{y}{z};
                       try
                           Topic{p} = cell2mat(importdata([Address{x}{y}{z},'\Topic_',Day{x}{y}{z},'.txt'],'',3));
                       catch
                           Topic{p} = [];
                       end
                       if isstruct(Topic{p}) | isempty(Topic{p})
                           p
                       end
                    end
                end
           end
        end   
     end
end

% Create table
xlswrite('ProtoMPEX_XP_list',[Date',Topic'])


