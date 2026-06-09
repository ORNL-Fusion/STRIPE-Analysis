% Script to convert Vp0 and Vp90 to angles:

% We assume that we are given the signals Vp0 and Vp90 as f{ch}{shot}
% ch = 1 is Vp0
% ch = 2 is Vp90
% ch = 3 is Vmag0

Vp0  = f{1}(:); % Vp0{shots}
Vp90 = f{2}(:); % Vp90{shots}
N = length(Vp0);
M = length(Vp0{1});
m = 1.8/180;

% For all shots:
for s = 1:N
    % For all data points in trace
    for p = 1:M
        % scan over all points
        f0  = Vp0{s}(p);
        f90 = Vp90{s}(p);
        
        if f90 > 0.9
            % QI or QII
            if f0 > 0.9
                Q = 1;
            else
                Q = 2;
            end
        elseif f90 <= 0.9
            % QIII or QIV
            if f0 > 0.9
                Q = 4;
            else
                Q = 3;
            end
        end
        
        if f0 > 1.72;
            K = 2;
        elseif f0 <= 0.1
            K = 2;
        else 
            K = 1;
        end
        % Now we have obtained "Q" and "K".
        % convert voltages to angles using the appropriate expression
        
        if K == 1
            A = (f0 - 1.8)/m;
            if Q == 1 || Q == 2;
                theta = -A;
            else 
                theta = +A;
            end
        elseif K == 2
            A = (f90 - 0.9)/m;
            if Q == 2
                theta = +180 - A;
            elseif Q == 3
                theta = -180 - A;                
            elseif Q == 1 || Q == 4
                theta = A;  
            end            
        end
        P{s}(p) = theta;
    end
    P_u{s}= unwrap(P{s}*pi/180)*180/pi;
    PhaseEnd = mean(P_u{s}((M-100):M));
    P_u{s} = P_u{s} - PhaseEnd;
end
%%
figure; 
hold on
for s = 1:N
h(s) = plot(P_u{s});
end
legend(h,num2str(shot'))
set(findobj('Type','line'),'LineWidth',1)

figure; 
hold on
for s = 1:N
h(s) = plot(P{s});
end
legend(h,num2str(shot'))
set(findobj('Type','line'),'LineWidth',1)