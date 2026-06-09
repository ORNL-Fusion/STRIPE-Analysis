% This script reads out the Oscope data taken during the Bdot probe XP on
% August 12th 2016

% Writte by Juan F Caneses on Octber 14th 2016

clear all
close all

% N = [0:5,10,11,21:55]; 
N = [0:5,10,11,21:55]; 


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
    
    % Calculating the hilbert transform of the data:
    yp = @(y) 1i*paddedhilbert(y) + y;
    H_fwd{s} = yp(fwd{s});
    H_bwd{s} = yp(bwd{s});
    H_Ch3{s} = yp(Ch3{s});
    H_Ch4{s} = yp(Ch4{s});
    
    % Phase between signal and reference:
    p3{s} = mean(unwrap( angle(H_Ch3{s}) - angle(H_fwd{s}) ));
    p4{s} = mean(unwrap( angle(H_Ch4{s}) - angle(H_fwd{s}) ));
    
    % Amplitude:
    A3{s} = mean(abs(H_Ch3{s}));
    A4{s} = mean(abs(H_Ch4{s}));
    
end


%% Plot all data N = 1 to N = 

close all

% Look at data from x' = 0 to 4.5
% rngName = [0,1,5,11,12,15:19];
% R = 0:0.5:4.5;
rng = [21 ,  30,  31,  34,  35,  40,  42   44,  45,  46];
R =   [4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0];

figure; hold on;
for s = 1:length(rng)
    k = find(N == rng(s));
    %plot(p4{k})
    pR(s) = p4{k};
    aR(s) = A4{k};
    plot(s,pR(s),'ro')
end
xlim([0,15])
    
figure; plot(-R,unwrap(pR))
figure; plot(-R,aR)

% Looking at 180 degree phase shift
% files 53 (theta) and 54(theta + 180 deg)

figure; hold on
s = find(N == 35);
hf(1) = plot(t{s},fwd{s},'k:');
plot(t{s},Ch4{s},'k')

H_fwd{s} = 1i*paddedhilbert(fwd{s}) + fwd{s};
fwd_angle{s} = angle(H_fwd{s});
H_Ch4{s} = 1i*paddedhilbert(Ch4{s}) + Ch4{s};
Ch4_angle{s} = angle(H_Ch4{s});
p{s} = unwrap(fwd_angle{s} - Ch4_angle{s});

s = find(N == 40);
hf(2) = plot(t{s},fwd{s},'r:');
plot(t{s},Ch4{s},'r')
legend(hf,'fwd','fwd')

H_fwd{s} = 1i*paddedhilbert(fwd{s}) + fwd{s};
fwd_angle{s} = angle(H_fwd{s});
H_Ch4{s} = 1i*paddedhilbert(Ch4{s}) + Ch4{s};
Ch4_angle{s} = angle(H_Ch4{s});
p{s} = unwrap(fwd_angle{s} - Ch4_angle{s});


figure; hold on
s = 23;
plot(t{s},p{s}/pi,'k')

s = 28;
plot(t{s},p{s}/pi,'r')
ylim([-2*pi,2*pi])

plot(t{s},(p{28}-p{23})/pi,'m')


% NOTE: the 180 degree phase shift is evident when one looks at the
% location of the reference first. notice that the fwd signal between shots
% has a phase shift, what matter is the phase difference between fwd and
% signal


