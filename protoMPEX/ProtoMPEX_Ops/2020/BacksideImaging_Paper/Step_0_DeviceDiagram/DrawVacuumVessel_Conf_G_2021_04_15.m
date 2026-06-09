% Define main dimensions of vacuum vessel and components:
% =========================================================================
% Inch to meter conversion factor:
in2m = 0.0254;

% Dump plate geometry:
rDump = 15.75*in2m/2;
zDump = 0.5;

% Central chamber:
zCC1 = coil{6}.z + 0.5*coil{6}.dz;
zCC2 = coil{7}.z - 0.5*coil{7}.dz;
rCC  = 24*in2m/2;

% Vaccum vessel:
rVac1 = 5.834*in2m/2; % Vacuum vessel rad that is on either side of ALN window
rVac2 = 19.25*in2m/2; % Wide sections dowmstream of central chamber to end
rVac3 = 4.272*2*in2m/2;  % Narrow sections are just inside of coil
rVac4 = 10*in2m/2; % 13.25 inch 5 way cross ID
rVac5 = 14.5*in2m/2; % 16.5 inch cross ID that connects to pump 
zVac1 = zCC1;
zVac2 = zCC2;
zVac3 = coil{end}.z + 0.5*coil{end}.dz + 3*in2m;
zVac4 = zVac3 + 20*in2m; % Length of 13.25 inch 5-way cross
zVac5 = zVac4 + 24.5*in2m; % Length of 16.5 inch target pump cross

% Pump duct Central chamber:
ID = 10*in2m;
L  = 10*in2m;
z1 = mean([coil{6}.z,coil{7}.z]) - ID/2;
z2 = z1 + ID;
r1 = rCC;
r2 = rCC + ID;
pumpDuct1.r = [r1,r2,r2,r1];
pumpDuct1.z = [z1,z1,z2,z2];

% Pump duct Target:
ID = 10*in2m;
L  = 10*in2m;
z1 = mean([zVac4,zVac5]) - ID/2;
z2 = z1 + ID;
r1 = rVac5;
r2 = rVac5 + ID;
pumpDuct2.r = [r1,r2,r2,r1];
pumpDuct2.z = [z1,z1,z2,z2];

% Helicon window:
L =  11.8*in2m; % L is 11.8" from Meitner
r1 = 4.95*in2m/2; % diameter is 5.1 inches --> update 4/21/16.  Nominal is ~4.95;
r1 = (5.426 - (2*0.25) )*in2m/2; % According to CAD dwg given to me by R. Goulding, OD is 5.426 inches and wall thickness is 0.25 inches
% r1 = 4.735*in2m/2; % New SiN window
r2 = rVac1;
z1 = 0.5*(coil{3}.z + coil{4}.z) - L/2;
z2 = 0.5*(coil{3}.z + coil{4}.z) + L/2;
heliconWindow.r = [r2,r1,r1,r2];
heliconWindow.z = [z1,z1,z2,z2];

% MPEX-like limiter:
limiterLength = 30e-2;
limiterWidth  = 3e-3;
limiterWidth  = 2.55e-3; % Email from John 2019-11-07: ID 12.09 cm and OD 12.6 cm (4.963 in)
z1 = heliconWindow.z(end) + 1e-2;
z2 = z1 + limiterLength;
r1 = 2.5*in2m - limiterWidth;
% r1 = 11.6/200; % New limiter that goes with SiN window
r2 = rVac1;
limiter.r = [r2,r1,r1,r2];
limiter.z = [z1,z1,z2,z2];

% Skimmer 1:
r1 = (7/100)/2;
r2 = rVac1;
z1 = coil{5}.z + 0.0702 - 0.5e-2;
z2 = coil{5}.z + 0.0702 + 0.5e-2;
skimmer1.r = [r2,r1,r1,r2];
skimmer1.z = [z1,z1,z2,z2];

% Skimmer 2:
r1 = (8.6/100)/2;
r2 = rVac3;
z1 = coil{7}.z  - 0.5e-2 - 4.25e-2;
z2 = coil{7}.z  + 0.5e-2 - 4.25e-2;
skimmer2.r = [r2,r1,r1,r2];
skimmer2.z = [z1,z1,z2,z2];

% ECH heating region:
r1 = rVac3;
r2 = rVac2;
z1 = coil{8}.z + 0.5*coil{8}.dz;
z2 = coil{9}.z - 0.5*coil{9}.dz;
echSection.r = [r1,r2,r2,r1];
echSection.z = [z1,z1,z2,z2];

% ICH Sleeve:
r1 = rVac3;
r2 = 0.08/2; % ID = 80 mm
z1 = coil{09}.z - 0.5*coil{09}.dz;
z2 = coil{10}.z + 0.5*coil{10}.dz;
ichSleeve.r = [r1,r2,r2,r1];
ichSleeve.z = [z1,z1,z2,z2];

% Space between coil 10-11
r1 = rVac3;
r2 = rVac2;
z1 = coil{10}.z + 0.5*coil{10}.dz;
z2 = coil{11}.z - 0.5*coil{11}.dz;
space1.r = [r1,r2,r2,r1];
space1.z = [z1,z1,z2,z2];

% Space between coil 11-12
r1 = rVac3;
r2 = rVac2;
z1 = coil{11}.z + 0.5*coil{11}.dz;
z2 = coil{12}.z - 0.5*coil{12}.dz;
space2.r = [r1,r2,r2,r1];
space2.z = [z1,z1,z2,z2];

% Space between coil 12-13
r1 = rVac3;
r2 = rVac2;
z1 = coil{end-1}.z + 0.5*coil{end-1}.dz;
z2 = coil{end}.z - 0.5*coil{end}.dz;
space3.r = [r1,r2,r2,r1];
space3.z = [z1,z1,z2,z2];

% Define coordinates of "base" Proto-MPEX vacuum vessel:
% =========================================================================
% 1- Define main features of vesse:
vessel_0.r = [0    ,rDump,rDump    ,rVac1     ,rVac1,rCC  ,rCC  ,rVac3,rVac3,rVac4,rVac4,rVac5,rVac5,0];
vessel_0.z = [zDump,zDump,zDump+0.2,zDump+0.2 ,zVac1,zVac1,zVac2,zVac2,zVac3,zVac3,zVac4,zVac4,zVac5,zVac5];

% 2- Add finer details to vaccum vessel:
vessel_0 = AddComponent(vessel_0,echSection);
vessel_0 = AddComponent(vessel_0,ichSleeve);
vessel_0 = AddComponent(vessel_0,space1);
vessel_0 = AddComponent(vessel_0,space2);
vessel_0 = AddComponent(vessel_0,space3);
vessel_0 = AddComponent(vessel_0,heliconWindow);
vessel_0 = AddComponent(vessel_0,skimmer1);
vessel_0 = AddComponent(vessel_0,skimmer2);

% Create upper and lower boundaries:
vessel_0_U = vessel_0;
vessel_0_L = vessel_0;

% Add pump duct to lower boundary:
vessel_0_L = AddComponent(vessel_0,pumpDuct1);
vessel_0_L = AddComponent(vessel_0_L,pumpDuct2);

% Proto-MPEX vacuum vessel with MPEX-like limiter:
vessel_1_U = AddComponent(vessel_0_U,limiter);
vessel_1_L = AddComponent(vessel_0_L,limiter);

% Segment the vacuum boundary:
vessel_0_U = SegmentBoundary(vessel_0_U,0.02);
vessel_0_L = SegmentBoundary(vessel_0_L,0.02);

% Segment Vessel with limiter:
vessel_1_U = SegmentBoundary(vessel_1_U,0.02);
vessel_1_L = SegmentBoundary(vessel_1_L,0.02);

% Plot "base" vessel:
% =========================================================================
testVesselBoundary = 1;
if testVesselBoundary
    figure;
    hold on

    plot(vessel_1_U.z,+vessel_1_U.r,'r.-') 
    plot(vessel_1_L.z,-vessel_1_L.r,'r.-') 
    plot(vessel_0_U.z,+vessel_0_U.r,'k.-')
    plot(vessel_0_L.z,-vessel_0_L.r,'k.-') 
    
    xlim([0,6])
    ylim([-1,1])
end