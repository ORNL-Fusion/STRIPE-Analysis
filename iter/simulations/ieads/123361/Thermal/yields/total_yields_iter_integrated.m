clc; clear all; close all;

% ======================= FILE DEFINITIONS ==========================
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

% ======================= GEOMETRY ==========================
if (exist('x1','var') == 0)
    fid = fopen('gitrGeometryPointPlane3d.cfg');
    if fid == -1
        error('Failed to open gitrGeometryPointPlane3d.cfg');
    end
    fgetl(fid); fgetl(fid); % two header lines
    for i = 1:18
        tline = fgetl(fid);
        evalc(tline);
    end
    fclose(fid);
end

subset = 1:length(x1);

% Build face vertices arrays (Nfaces x 3)
Xtri = [transpose(x1(subset)), transpose(x2(subset)), transpose(x3(subset))];
Ytri = [transpose(y1(subset)), transpose(y2(subset)), transpose(y3(subset))];
Ztri = [transpose(z1(subset)), transpose(z2(subset)), transpose(z3(subset))];
Nfaces = size(Xtri,1);

% Triangle areas (m^2)
p1 = [Xtri(:,1), Ytri(:,1), Ztri(:,1)];
p2 = [Xtri(:,2), Ytri(:,2), Ztri(:,2)];
p3 = [Xtri(:,3), Ytri(:,3), Ztri(:,3)];
v1 = p2 - p1;
v2 = p3 - p1;
cp = cross(v1, v2, 2);
faceArea = 0.5 * sqrt(sum(cp.^2,2));
faceArea(faceArea==0) = eps;

Atot = sum(faceArea,'omitnan');
if ~isfinite(Atot) || Atot <= 0
    error('Total area non-positive. Check geometry.');
end

% ======================= PER-CHARGE DATA ==========================
labels = {'D+','Ne1+','Ne2+','Ne3+','Ne4+','Ne5+','Ne6+','Ne7+','Ne8+','Ne9+','Ne10+'};

% Yields per face (pad with leading 0 like you do)
Y_D    = [0; yields_D(:,1)];
Y_Ne1  = [0; yields_Ne1(:,1)];
Y_Ne2  = [0; yields_Ne2(:,2)];
Y_Ne3  = [0; yields_Ne3(:,3)];
Y_Ne4  = [0; yields_Ne4(:,4)];
Y_Ne5  = [0; yields_Ne5(:,5)];
Y_Ne6  = [0; yields_Ne6(:,6)];
Y_Ne7  = [0; yields_Ne7(:,7)];
Y_Ne8  = [0; yields_Ne8(:,8)];
Y_Ne9  = [0; yields_Ne9(:,9)];
Y_Ne10 = [0; yields_Ne10(:,10)];

Ylist = {Y_D,Y_Ne1,Y_Ne2,Y_Ne3,Y_Ne4,Y_Ne5,Y_Ne6,Y_Ne7,Y_Ne8,Y_Ne9,Y_Ne10};
for k = 1:numel(Ylist)
    yk = Ylist{k};
    if numel(yk) < Nfaces
        error('Yield vector %s shorter than Nfaces (len=%d, Nfaces=%d).', labels{k}, numel(yk), Nfaces);
    end
    Ylist{k} = yk(1:Nfaces);
end

% density and flow
n_D    = data_D(1:Nfaces,2);     v_D    = data_D(1:Nfaces,5);
n_Ne1  = data_Ne1(1:Nfaces,11);  v_Ne1  = data_Ne1(1:Nfaces,5);
n_Ne2  = data_Ne2(1:Nfaces,11);  v_Ne2  = data_Ne2(1:Nfaces,5);
n_Ne3  = data_Ne3(1:Nfaces,11);  v_Ne3  = data_Ne3(1:Nfaces,5);
n_Ne4  = data_Ne4(1:Nfaces,11);  v_Ne4  = data_Ne4(1:Nfaces,5);
n_Ne5  = data_Ne5(1:Nfaces,11);  v_Ne5  = data_Ne5(1:Nfaces,5);
n_Ne6  = data_Ne6(1:Nfaces,11);  v_Ne6  = data_Ne6(1:Nfaces,5);
n_Ne7  = data_Ne7(1:Nfaces,11);  v_Ne7  = data_Ne7(1:Nfaces,5);
n_Ne8  = data_Ne8(1:Nfaces,11);  v_Ne8  = data_Ne8(1:Nfaces,5);
n_Ne9  = data_Ne9(1:Nfaces,11);  v_Ne9  = data_Ne9(1:Nfaces,5);
n_Ne10 = data_Ne10(1:Nfaces,11); v_Ne10 = data_Ne10(1:Nfaces,5);

nList = {n_D,n_Ne1,n_Ne2,n_Ne3,n_Ne4,n_Ne5,n_Ne6,n_Ne7,n_Ne8,n_Ne9,n_Ne10};
vList = {v_D,v_Ne1,v_Ne2,v_Ne3,v_Ne4,v_Ne5,v_Ne6,v_Ne7,v_Ne8,v_Ne9,v_Ne10};

% ======================= TOTAL-AREA AVERAGES (NO WETTED MASK) ==========================
Yield_areaAvg_totalA   = nan(numel(labels),1);   % <Y>_A
EroFlux_areaAvg_totalA = nan(numel(labels),1);   % <Gamma>_A

for q = 1:numel(labels)
    y = Ylist{q}(:);
    n = nList{q}(:);
    v = vList{q}(:);

    ion_flux = n .* abs(v);       % #/m^2/s
    ero_flux = y .* ion_flux;     % #/m^2/s

    valid = isfinite(y) & isfinite(ion_flux) & isfinite(faceArea);

    % <Y>_A = sum(Y*A)/Atot
    Yield_areaAvg_totalA(q) = sum( y(valid).*faceArea(valid), 'omitnan') / Atot;

    % <Gamma>_A = sum(Gamma*A)/Atot, Gamma=Y*n|v|
    EroFlux_areaAvg_totalA(q) = sum( ero_flux(valid).*faceArea(valid), 'omitnan') / Atot;
end

% ======================= PRINT RESULTS ==========================
fprintf('\n=== TOTAL-AREA averages (NO wetted mask) ===\n');
fprintf('Atot = %.6e m^2\n\n', Atot);

for q = 1:numel(labels)
    fprintf('%-5s  <Y>_A=%.6e   <Gamma>_A=%.6e [#/m^2/s]\n', ...
        labels{q}, Yield_areaAvg_totalA(q), EroFlux_areaAvg_totalA(q));
end

% save to CSV
T = table(labels(:), Yield_areaAvg_totalA, EroFlux_areaAvg_totalA, ...
    'VariableNames', {'Species','Yield_areaAvg_totalA','EroFlux_areaAvg_totalA_per_m2_s'});
writetable(T, 'totalArea_averages_only.csv');

% ======================= TWO SUMMARY BAR PLOTS ==========================
figure('Name','Total-area averaged yield and erosion flux','Color','w','Position',[100 100 1100 420]);

subplot(1,2,1);
bar(Yield_areaAvg_totalA);
set(gca,'XTick',1:numel(labels),'XTickLabel',labels,'XTickLabelRotation',45);
ylabel('<Y>_A'); title('Area-averaged yield (total area)');
grid on;

subplot(1,2,2);
bar(EroFlux_areaAvg_totalA);
set(gca,'XTick',1:numel(labels),'XTickLabel',labels,'XTickLabelRotation',45);
ylabel('<\Gamma>_A  [#/m^2/s]'); title('Area-averaged erosion flux (total area)');
grid on;