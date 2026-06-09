clc; clear all; close all;

% ======================= FILE DEFINITIONS ==========================
% Include ALL Neon charge states (Ne1+ … Ne10+) and D+
yields_D    = readmatrix("../ieads_D+/yields_D+.csv");
yields_Ne1  = readmatrix("../ieads_Ne1+/yields_Ne1+.csv");
yields_Ne2  = readmatrix("../ieads_Ne2+/yields_Ne2+.csv");
yields_Ne3  = readmatrix("../ieads_Ne3+/yields_Ne3+.csv");
yields_Ne4  = readmatrix("../ieads_Ne4+/yields_Ne4+.csv");
yields_Ne5  = readmatrix("../ieads_Ne5+/yields_Ne5+.csv");
yields_Ne6  = readmatrix("../ieads_Ne6+/yields_Ne6+.csv");
yields_Ne7  = readmatrix("../ieads_Ne7+/yields_Ne7+.csv");
yields_Ne8  = readmatrix("../ieads_Ne8+/yields_Ne8+.csv");
yields_Ne9  = readmatrix("../ieads_Ne9+/yields_Ne9+.csv");
yields_Ne10 = readmatrix("../ieads_Ne10+/yields_Ne10+.csv");

data_D    = readmatrix("../ieads_D+/Targets_D+.txt");
data_Ne1  = readmatrix("../ieads_Ne1+/Targets_Ne1+.txt");
data_Ne2  = readmatrix("../ieads_Ne2+/Targets_Ne2+.txt");
data_Ne3  = readmatrix("../ieads_Ne3+/Targets_Ne3+.txt");
data_Ne4  = readmatrix("../ieads_Ne4+/Targets_Ne4+.txt");
data_Ne5  = readmatrix("../ieads_Ne5+/Targets_Ne5+.txt");
data_Ne6  = readmatrix("../ieads_Ne6+/Targets_Ne6+.txt");
data_Ne7  = readmatrix("../ieads_Ne7+/Targets_Ne7+.txt");
data_Ne8  = readmatrix("../ieads_Ne8+/Targets_Ne8+.txt");
data_Ne9  = readmatrix("../ieads_Ne9+/Targets_Ne9+.txt");
data_Ne10 = readmatrix("../ieads_Ne10+/Targets_Ne10+.txt");

density_D = data_D(:,2);
flow_D = data_D(:,5);
ero_data_D = [0;yields_D(:,1)].*density_D.*flow_D;

density_Ne1 = data_Ne1(:,2);
flow_Ne1 = data_Ne1(:,5);
ero_data_D = [0;yields_Ne1(:,1)].*density_Ne1.*flow_Ne1;

density_Ne2 = data_Ne2(:,2);
flow_Ne2 = data_Ne2(:,5);
ero_data_D = [0;yields_Ne2(:,2)].*density_Ne2.*flow_Ne2;


% % ======================= TOTALS ==========================
% % Sum all Neon (1..10) + D+ (use col 1 for yields, col 4 for density, col 6 for flow)
% 
% yields_total_Ne = yields_Ne1(:,1)  + yields_Ne2(:,2)  + yields_Ne3(:,3)  + ...
%                yields_Ne4(:,4)+ yields_Ne5(:,5)  + yields_Ne6(:,6)  + yields_Ne7(:,7)  + ...
%                yields_Ne8(:,8)+ yields_Ne9(:,9)  + yields_Ne10(:,10);
% 
% yields_total_D = yields_D;
% 
% density_total_Ne =  data_Ne1(:,11)    + data_Ne2(:,11)    + data_Ne3(:,11)    + ...
%                 data_Ne4(:,11) + data_Ne5(:,11)    + data_Ne6(:,11)    + data_Ne7(:,11)    + ...
%                 data_Ne8(:,11) + data_Ne9(:,11)    + data_Ne10(:,11);
% 
% density_total_D = data_D(:,2);
% 
% % Keep same format (simple average of species flows). Now 11 species total.
% flow_total_Ne = ( data_Ne1(:,5)    + data_Ne2(:,5)    + data_Ne3(:,5)    + ...
%               data_Ne4(:,5)   + data_Ne5(:,5)    + data_Ne6(:,5)    + data_Ne7(:,5)    + ...
%               data_Ne8(:,5)   + data_Ne9(:,5)    + data_Ne10(:,5)) ./ 10;
% 
% flow_total_D = data_D(:,5);
% 
% % Shot number tag
% shot_tag = "IMAS123361";
% 
% writematrix(yields_total_Ne,  "yields_total_Ne_" + shot_tag + ".txt");
% writematrix(density_total_Ne, "density_total_Ne_" + shot_tag + ".txt");
% writematrix(flow_total_Ne,    "flow_total_Ne" + shot_tag + ".txt");
% 
% % % ======================= LOAD BASE SCALARS ==========================
% % data = readmatrix('Targets_Ne1.txt');
% % potential_data = data(:,1);
% % ne_data        = data(:,2);
% % te_data        = data(:,3);
% % ni_data        = data(:,4);
% % v_data         = data(:,6);
% 
% % ======================= EROSION FLUX ==========================
% ero_data = [0;yields_total_Ne] .* density_total_Ne .* flow_total_Ne + ...
%     [0;yields_total_D] .* density_total_D .* flow_total_D;
% 
% ero_data_Ne = [0;yields_total_Ne] .* density_total_Ne .* flow_total_Ne;
% ero_data_D =  [0;yields_total_D] .* density_total_D .* flow_total_D;



% ======================= GEOMETRY ==========================
if (exist('x1') == 0)
    fid = fopen(strcat('gitrGeometryPointPlane3d.cfg'));
    tline = fgetl(fid);
    tline = fgetl(fid);
    for i = 1:18
        tline = fgetl(fid);
        evalc(tline);
    end
    Zsurface = Z;
end

subset = 1:length(x1);

% ======================= FIGURE 1 ==========================

potential_data=data_Ne10(:,1);
figure(1)
X = [transpose(x1(subset)), transpose(x2(subset)), transpose(x3(subset))];
Y = [transpose(y1(subset)), transpose(y2(subset)), transpose(y3(subset))];
Z = [transpose(z1(subset)), transpose(z2(subset)), transpose(z3(subset))];
patch(transpose(X), transpose(Y), transpose(Z), potential_data, 'FaceAlpha', 1, 'EdgeAlpha', 0.3);
title('Yields')
colorbar('eastoutside')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')

% ======================= FIGURE 2 ==========================
figure(2)
X = [transpose(x1(subset)), transpose(x2(subset)), transpose(x3(subset))];
Y = [transpose(y1(subset)), transpose(y2(subset)), transpose(y3(subset))];
Z = [transpose(z1(subset)), transpose(z2(subset)), transpose(z3(subset))];
patch(transpose(X), transpose(Y), transpose(Z), density_total_D, 'FaceAlpha', 1, 'EdgeAlpha', 0.3);
title('Density')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
colorbar('eastoutside')

% ======================= FIGURE 3 ==========================
te_data=data_Ne10(:,3);
figure(3)
X = [transpose(x1(subset)), transpose(x2(subset)), transpose(x3(subset))];
Y = [transpose(y1(subset)), transpose(y2(subset)), transpose(y3(subset))];
Z = [transpose(z1(subset)), transpose(z2(subset)), transpose(z3(subset))];
patch(transpose(X), transpose(Y), transpose(Z), te_data, 'FaceAlpha', 1, 'EdgeAlpha', 0.3);
title('Temperature')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
colorbar('eastoutside')

% ======================= FIGURE 4 ==========================
figure(4)
X = [transpose(x1(subset)), transpose(x2(subset)), transpose(x3(subset))];
Y = [transpose(y1(subset)), transpose(y2(subset)), transpose(y3(subset))];
Z = [transpose(z1(subset)), transpose(z2(subset)), transpose(z3(subset))];
patch(transpose(X), transpose(Y), transpose(Z), ero_data_Ne, 'FaceAlpha', 1, 'EdgeAlpha', 0.3);
title('Flux')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
colorbar('eastoutside')

% ======================= FIGURE 5 ==========================
figure(5)
X = [transpose(x1(subset)), transpose(x2(subset)), transpose(x3(subset))];
Y = [transpose(y1(subset)), transpose(y2(subset)), transpose(y3(subset))];
Z = [transpose(z1(subset)), transpose(z2(subset)), transpose(z3(subset))];
patch(transpose(X), transpose(Y), transpose(Z), flow_total_D, 'FaceAlpha', 1, 'EdgeAlpha', 0.3);
title('Flow')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
colorbar('eastoutside')

% ======================= OPTIONAL EROSION RATE ==========================
% figure(6)
% erosion_rate = flux_data .* area;
% writematrix(erosion_rate,'erosion_rate.csv')
% imagesc(erosion_rate)
% title('Flux')
% xlabel('X [m]')
% ylabel('Y [m]')
% zlabel('Z [m]')
% colorbar('eastoutside')