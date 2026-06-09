% This script reads the data acquired with the RF phase detector installed
% in ProtoMPEX on Sept 16th 2016.
% Looking at the 180 degree phase shift of the data:

close all
clear all

shot = 1e4 + 300 + [70,72,73,78];
% 70: +90 deg
% 72: +90 deg
% 73: -90 deg
% 78: -90 deg

% Address of data:
Data{1} = ['\MPEX::TOP.MACHOPS1:INT_4MM_1'];
Data{2} = ['\MPEX::TOP.MACHOPS1:INT_4MM_2'];
Data{3} = ['\MPEX::TOP.MACHOPS1:INT_2MM_1'];
Data{4} = ['\MPEX::TOP.MACHOPS1:INT_2MM_2'];

% name of plots:
Title{1} = 'Vp0';
Title{2} = 'Vp90';
Title{3} = 'Vmag0';
Title{4} = 'None';
m_v = 1.8/60; % V/dB slope of the amplitude ratio signal
tStart = 4.23;
tEnd   = 4.33;
ClipInTime = 0;

%% Import and Plot data:
% in this section we import the data from the Server and then plot it to
% inspect it before any analysis can begin
% Two types of data come out of this: (1) Magnitude ratio and (2) Phase
% information
% The only analysis performed is to convert the dB ratio into linear ratio

figure;
for ch = 1:4 % For all channels on digitizer, 4 in this case
    [f{ch},t{ch}] = my_mdsvalue_v2(shot,Data(ch));
    % Remember the following: f{channel}{shot}
    
    subplot(2,2,ch);hold on
    
    for s = 1:length(shot)
        if ClipInTime == 1;
            rng = find(t{ch}{s}>= tStart & t{ch}{s}<= tEnd);
            t{ch}{s} = t{ch}{s}(rng);
            f{ch}{s} = f{ch}{s}(rng);
        end
            if ch == 3; % if ch == Vmag0
                % Convert Vmag0 to linear scale
                VmagRatio{ch}{s} = 10.^( ((f{ch}{s}-0.9)./m_v )/20);
            end
            if ch~=4 
                h(s) = plot(t{ch}{s},f{ch}{s});
                ylim([0,2]);xlim([4.25,4.33])
                title(Title{ch})
                ylabel('[V]')
                xlabel('t [s]')
            else
                plot(t{3}{s},VmagRatio{3}{s});
                ylim([0,10]); xlim([4.25,4.33])
                title('Amplitude Ratio (Linear)')
                ylabel('[Input/reference]')
                xlabel('t [s]')                        
            end                
            box on
    end
end
legend(h,num2str(shot'),'location','SouthWest')
set(gcf,'color','w')
 
%% Convert phase signals into a single angle

% Phase angle ANALYSIS starts in this section

% We assume that we are given the signals Vp0 and Vp90 as f{ch}{shot}
% ch = 1 is Vp0
% ch = 2 is Vp90
% ch = 3 is Vmag0

Vp0  = f{1}(:);  % Vp0{shots}, in this case we have 4 shots
Vp90 = f{2}(:);  % Vp90{shots}
N = length(Vp0); % N = number of shots 
M = length(Vp0{1}); % M = number of points in trace of single shot
m = 1.8/180; % [V/deg] slope of phase signal

% For all shots:
for s = 1:N % For all shots, 4 in this case
    for p = 1:M % For all data points in the trace of sth shot
        f0  = Vp0{s}(p);  % In phase trace
        f90 = Vp90{s}(p); % Quadrature phase
        %==================================================================
        % Determine which half of the complex plane the angle is: Upper or
        % Lower half
        if f90 > 0.9
            % QI or QII
            if f0 > 0.9;Q = 1;
            else Q = 2;
            end
        elseif f90 <= 0.9
            % QIII or QIV
            if f0 > 0.9;Q = 4;
            else Q = 3;
            end
        end
        %==================================================================
        
        %==================================================================
        % Choose the part of the trace less prone to error:
        if     f0 > 1.72; K = 2;
        elseif f0 <= 0.1; K = 2;
        else              K = 1;
        end
        %==================================================================
        % Now we have obtained "Q" and "K".
        % Q: defines the quadrant the angle is on
        % K: defines if we use Vph0 or Vph90, K:1 is Vph0 
        
        %==================================================================
        % Convert voltages to angles using the appropriate expression
        if K == 1; A = (f0 - 1.8)/m;
            if Q == 1 || Q == 2;     theta = -A;
            else                     theta = +A;
            end
        elseif K == 2;                A = (f90 - 0.9)/m;
            if Q == 2;                theta = +180 - A;
            elseif Q == 3;            theta = -180 - A;                
            elseif Q == 1 || Q == 4;  theta = A;  
            end            
        end
        P{s}(p) = theta;
        %==================================================================
    end
    % Once we have analized all the signals, we unwrap the angle
    P_u{s}= unwrap(P{s}*pi/180)*180/pi;
    PhaseEnd = mean(P_u{s}((M-100):M));
    P_u{s} = P_u{s} - PhaseEnd;
end
%%
figure; 
subplot(2,1,1); hold on
Offset = [360,0,0,360] + 360 + 80;
for s = 1:N
h(s) = plot(t{1}{1},P_u{s}-Offset(s));
end
legend(h,num2str(shot'),'location','West')
set(findobj('Type','line'),'LineWidth',1)
ylim([0,400]); xlim([tStart,tEnd])
set(gca,'YTick',[0:60:360]); grid on
ylabel('Angle [deg]')
title('Unwrapped phase'); box on

subplot(2,1,2); hold on
for s = 1:N
h(s) = plot(t{1}{1},P{s});
end
xlim([tStart,tEnd])
set(findobj('Type','line'),'LineWidth',1)
title('Phase before unwrapping'); box on
ylabel('Angle [deg]')
set(gca,'YTick',[-180:60:180]); grid on
set(gcf,'color','w')

% Results:
% the data clearly indicates a 180 degree phase shift due to the rotation
% of the probe.
% Clearly the unwrapped data is the best method to use