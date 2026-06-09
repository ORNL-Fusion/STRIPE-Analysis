close all
clear all
verbose = 0;
% #########################################################################
helicon_current = [160,180,200,200,220,220,240,240,240,260,260,260,280,280,300,300,330,330,160,160,160,140,140,120,120,100, 80, 80, 40, 40, 20, 20,  0,  0] ; %
%  Total time is 17 sec

current_A = 4000;
current_B = 4000;
current_C = 600;
config = 'flat';
skimmer = 1;
target_position = 2; %=1 puts it at 7.5, =2 puts it at 11.5
sleeve = 1; 

% Source location:
zz = linspace(0,5,600);
rng = find(zz>=1.72 & zz<=1.76);

for n = 1:length(helicon_current)
% Build coils
    [coil,current] = build_Proto_coils(helicon_current(n),current_A,current_B,config,verbose,current_C);
%     geo = get_Proto_geometry(0,0,skimmer,target_position,sleeve);
    
    bfield.coil = coil;
    bfield.current = current;

%% Mod B
    for s = 1:length(zz);
        [Bx{n}(s),By{n}(s),Bz{n}(s),Btot{n}(s)]=bfield_bs_jdl(0,0,zz(s),bfield.coil,bfield.current);
    end
    % Extract magnetic field at the source:
    B0Data.B0(n) = mean(Bz{n}(rng)); 
    B0Data.TR2(n) = helicon_current(n); 
end
save('SourceMagneticField_2017_12_01','B0Data')
%%
close all
figure; plot(B0Data.TR2,B0Data.B0,'ko')
