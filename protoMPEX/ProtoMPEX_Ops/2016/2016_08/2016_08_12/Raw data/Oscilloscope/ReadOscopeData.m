% This script reads out the Oscope data taken during the Bdot probe XP on
% August 12th 2016
% Writte by Juan F Caneses

clear all
close all

N = [0:5,10,11,21:55];
f = @(x,t) x(1).*sin(x(2).*(t-x(3)) ) + x(4);

for s = 1:length(N)
    if length(num2str(N(s))) == 1
        FileName{s} = ['tek00','0',num2str(N(s)),'ALL.csv'];
    else
        FileName{s} = ['tek00',num2str(N(s)),'ALL.csv'];
    end
    D{s} = importdata(FileName{s},',',21);
    % time signal
    if s == 38
        t{s} = D{s}.data(:,1)'; % time clipped
    else           
        t_c{s} = D{s}.data(:,1); % time clipped
        locs{s} = peakseek(diff(t_c{s}),9);
        locs0{s} = locs{s} - locs{s}(1) + 1;
        t_sparse{s} = t_c{s}(locs0{s});
        x0 = [locs0{s}(1),locs0{s}(end)];
        v0 = [t_sparse{s}(1),t_sparse{s}(end)];
        t{s} = interp1(x0,v0,1:length(t_c{s}),'linear','extrap');
    end
    
    fwd{s} = D{s}.data(:,2);
    bwd{s} = D{s}.data(:,3);
    Ch3{s} = D{s}.data(:,4); 
    Ch4{s} = D{s}.data(:,5);
    
    % fitting RF sines to data
    x0(1) = max(fwd{s});
    x0(2) = 2*pi*13.56e6;
    x0(3) = 0;
    x0(4) = 0*x0(1);
    
    res = @(x)  f(x,t{s}') - fwd{s};
    [x1{s}(:,s),ssq,cnt] = LMFnlsq(res,x0,'MaxIter',10);
end
%%
close all

Scenario = 2;
switch Scenario
    case 1
        L = [40,41]; % both at 90 degrees
    case 2
        L = [41,42]; % 90 and 270 degrees
    case 3
        L = [42,43]; % both at 270 degrees
end

% shot = [53];
% if length(shot) == 1
%     L = find(N == shot)*[1,1];
% elseif length (shot) == 2
%     L = [find(N == shot(1)),find(N == shot(2))];
% end

C = {'k','r','k:','r:'};
tStartPlot = 0.1342050;
tEndPlot   = 0.1342072;

%figure; 
 %s = L(2); plot(t{s},fwd{s}); hold on; plot(t{s},f(x1{s},t{s}),'r')

figure; hold on
s = L(1); plot(t{s},fwd{s}); hold on; plot(t{s},Ch4{s},'r');
s = L(2); plot(t{s},fwd{s},'bl:','LineWidth',1); hold on; plot(t{s},Ch4{s},'r:','LineWidth',1);
%xlim([tStartPlot,tEndPlot])
ylim([-1,1])

figure; hold on
s = L(1); plot(t{s},fwd{s}); plot(t{s},Ch3{s},'r');
s = L(2); plot(t{s},fwd{s},'bl:','LineWidth',1); plot(t{s},Ch3{s},'r:','LineWidth',1);
%xlim([tStartPlot,tEndPlot])
ylim([-1,1])
