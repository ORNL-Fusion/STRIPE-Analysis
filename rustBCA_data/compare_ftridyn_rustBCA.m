% compare_ftridyn_rustBCA_noInterp.m
% Compare variables from ftridynSelf.nc and rustbca_as_ftridynSelf_WonW_EXACT.nc
% WITHOUT interpolation. Each dataset plotted on its own native (A,E) grid.
%
% Canonical internal order:
%   2D -> (nE,nA)
%   3D -> (nE,nA,nBins)
%
% Then for 3D we also plot:
%   sum over bins map (E,A)
%   weighted-mean bin-coordinate map (E,A)
%   example PDFs at nearest (E,A) in each dataset

clear; close all; clc;

%% -------- inputs --------
f1 = "ftridynSelf.nc";
f2 = "ftridynSelf_new_higResol.nc";

outDir = "compare_plots_noInterp";
if ~exist(outDir,"dir"), mkdir(outDir); end
doSavePNGs = true;
dpi = 200;

% choose some reference energies/angles (in physical units)
% We will pick nearest indices in each dataset separately.
targetE_vals = [];     % e.g. [50 200 1000];
targetA_vals = [];     % e.g. [5 30 60];

%% -------- read grids --------
E1 = ncread(f1,"E");  A1 = ncread(f1,"A");
E2 = ncread(f2,"E");  A2 = ncread(f2,"A");
E1 = E1(:); A1 = A1(:);
E2 = E2(:); A2 = A2(:);

% bin axes
theta1 = readVarIf(f1,"thetaGrid");
theta2 = readVarIf(f2,"thetaGrid");
phi1   = readVarIf(f1,"phiGrid");
phi2   = readVarIf(f2,"phiGrid");
egrid1    = readVarIf(f1,"eDistEgrid");
egrid2    = readVarIf(f2,"eDistEgrid");
egridRef1 = readVarIf(f1,"eDistEgridRef");
egridRef2 = readVarIf(f2,"eDistEgridRef");

%% -------- pick overlay points if empty --------
if isempty(targetE_vals), targetE_vals = pickVals(E1,3); end
if isempty(targetA_vals), targetA_vals = pickVals(A1,3); end

idxE1 = arrayfun(@(x) nearestIndex(E1,x), targetE_vals);
idxA1 = arrayfun(@(x) nearestIndex(A1,x), targetA_vals);
idxE2 = arrayfun(@(x) nearestIndex(E2,x), targetE_vals);
idxA2 = arrayfun(@(x) nearestIndex(A2,x), targetA_vals);

%% ===================== 0) grid sanity =====================
fig = figure('Name','Grids','Color','w'); tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile; plot(E1,'-'); grid on; title('FTRI E'); xlabel('index'); ylabel('E');
nexttile; plot(E2,'-'); grid on; title('RUST E'); xlabel('index'); ylabel('E');
nexttile; plot(A1,'-'); grid on; title('FTRI A'); xlabel('index'); ylabel('A');
nexttile; plot(A2,'-'); grid on; title('RUST A'); xlabel('index'); ylabel('A');

saveFig(fig,outDir,"00_grids",doSavePNGs,dpi);

%% ===================== 1) yields =====================
vars2D = ["spyld","rfyld"];

for v = vars2D
    if ~(hasVar(f1,char(v)) && hasVar(f2,char(v))), continue; end

    V1 = readEA_canonical(f1,char(v));   % (nE1,nA1)
    V2 = readEA_canonical(f2,char(v));   % (nE2,nA2)

    fig = figure('Name',char(v),'Color','w'); tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

    % --- maps on native grids ---
    nexttile;
    imagesc(A1, E1, V1); axis xy; colorbar; grid on;
    xlabel('A'); ylabel('E'); title(v + " (FTRI native)");

    nexttile;
    imagesc(A2, E2, V2); axis xy; colorbar; grid on;
    xlabel('A'); ylabel('E'); title(v + " (RUST native)");

    % --- cuts vs E at selected A (nearest in each dataset) ---
    nexttile; hold on; grid on;
    for k=1:numel(idxA1)
        plot(E1, V1(:,idxA1(k)),'-','DisplayName',sprintf('FTRI A=%.3g',A1(idxA1(k))));
        plot(E2, V2(:,idxA2(k)),'--','DisplayName',sprintf('RUST A=%.3g',A2(idxA2(k))));
    end
    set(gca,'XScale','log');  % energies often log-spaced
    xlabel('E'); ylabel(char(v)); title('Cuts vs E (nearest A in each)'); legend('Location','best');

    % --- cuts vs A at selected E (nearest in each dataset) ---
    nexttile; hold on; grid on;
    for k=1:numel(idxE1)
        plot(A1, V1(idxE1(k),:),'-','DisplayName',sprintf('FTRI E=%.3g',E1(idxE1(k))));
        plot(A2, V2(idxE2(k),:),'--','DisplayName',sprintf('RUST E=%.3g',E2(idxE2(k))));
    end
    xlabel('A'); ylabel(char(v)); title('Cuts vs A (nearest E in each)'); legend('Location','best');

    saveFig(fig,outDir,"01_yield_"+v,doSavePNGs,dpi);
end

%% ===================== 2) angular distributions =====================
cosVars = ["cosXDist","cosYDist","cosZDist","cosXDistRef","cosYDistRef","cosZDistRef"];

for v = cosVars
    if ~(hasVar(f1,char(v)) && hasVar(f2,char(v))), continue; end

    C1 = readEAB_canonical(f1,char(v));  % (nE1,nA1,nB1)
    C2 = readEAB_canonical(f2,char(v));  % (nE2,nA2,nB2)

    sum1 = squeeze(sum(C1,3,'omitnan'));
    sum2 = squeeze(sum(C2,3,'omitnan'));

    x1 = pickBinAxis(theta1, phi1, size(C1,3));
    x2 = pickBinAxis(theta2, phi2, size(C2,3));

    mean1 = weightedMean3D(C1,x1);
    mean2 = weightedMean3D(C2,x2);

    fig = figure('Name',char(v),'Color','w'); tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

    nexttile; imagesc(A1,E1,sum1); axis xy; colorbar; grid on;
    xlabel('A'); ylabel('E'); title(v + " sumBins (FTRI native)");

    nexttile; imagesc(A2,E2,sum2); axis xy; colorbar; grid on;
    xlabel('A'); ylabel('E'); title(v + " sumBins (RUST native)");

    nexttile; imagesc(A1,E1,mean1); axis xy; colorbar; grid on;
    xlabel('A'); ylabel('E'); title(v + " meanBinCoord (FTRI native)");

    nexttile; imagesc(A2,E2,mean2); axis xy; colorbar; grid on;
    xlabel('A'); ylabel('E'); title(v + " meanBinCoord (RUST native)");

    saveFig(fig,outDir,"02_cos_"+v,doSavePNGs,dpi);

    % PDF overlay at reference point (nearest in each)
    e1 = idxE1(min(2,end)); a1 = idxA1(min(2,end));
    e2 = idxE2(min(2,end)); a2 = idxA2(min(2,end));

    pdf1 = squeeze(C1(e1,a1,:));
    pdf2 = squeeze(C2(e2,a2,:));

    fig2 = figure('Name',char(v)+"_PDF",'Color','w'); hold on; grid on;
    plot(x1, pdf1./max(pdf1+eps),'-','DisplayName',sprintf('FTRI @E=%.3g A=%.3g',E1(e1),A1(a1)));
    plot(x2, pdf2./max(pdf2+eps),'--','DisplayName',sprintf('RUST @E=%.3g A=%.3g',E2(e2),A2(a2)));
    xlabel('bin coordinate'); ylabel('PDF (norm by max)');
    title(v + " PDF overlay (nearest point in each)");
    legend('Location','best');
    saveFig(fig2,outDir,"02_cosPDF_"+v,doSavePNGs,dpi);
end

%% ===================== 3) energy distributions =====================
eVars = ["energyDist","energyDistRef"];

for v = eVars
    if ~(hasVar(f1,char(v)) && hasVar(f2,char(v))), continue; end

    D1 = readEAB_canonical(f1,char(v));  % (nE1,nA1,nB1)
    D2 = readEAB_canonical(f2,char(v));  % (nE2,nA2,nB2)

    if v=="energyDistRef"
        x1 = chooseOrIndex(egridRef1, size(D1,3));
        x2 = chooseOrIndex(egridRef2, size(D2,3));
    else
        x1 = chooseOrIndex(egrid1, size(D1,3));
        x2 = chooseOrIndex(egrid2, size(D2,3));
    end

    sum1 = squeeze(sum(D1,3,'omitnan'));
    sum2 = squeeze(sum(D2,3,'omitnan'));
    mean1 = weightedMean3D(D1,x1);
    mean2 = weightedMean3D(D2,x2);

    fig = figure('Name',char(v),'Color','w'); tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

    nexttile; imagesc(A1,E1,sum1); axis xy; colorbar; grid on;
    xlabel('A'); ylabel('E'); title(v + " sumBins (FTRI native)");

    nexttile; imagesc(A2,E2,sum2); axis xy; colorbar; grid on;
    xlabel('A'); ylabel('E'); title(v + " sumBins (RUST native)");

    nexttile; imagesc(A1,E1,mean1); axis xy; colorbar; grid on;
    xlabel('A'); ylabel('E'); title(v + " meanEnergy (FTRI native)");

    nexttile; imagesc(A2,E2,mean2); axis xy; colorbar; grid on;
    xlabel('A'); ylabel('E'); title(v + " meanEnergy (RUST native)");

    saveFig(fig,outDir,"03_energy_"+v,doSavePNGs,dpi);

    % PDF overlay at reference point (nearest in each)
    e1 = idxE1(min(2,end)); a1 = idxA1(min(2,end));
    e2 = idxE2(min(2,end)); a2 = idxA2(min(2,end));

    pdf1 = squeeze(D1(e1,a1,:));
    pdf2 = squeeze(D2(e2,a2,:));

    fig2 = figure('Name',char(v)+"_PDF",'Color','w'); hold on; grid on;
    plot(x1, pdf1./max(pdf1+eps),'-','DisplayName',sprintf('FTRI @E=%.3g A=%.3g',E1(e1),A1(a1)));
    plot(x2, pdf2./max(pdf2+eps),'--','DisplayName',sprintf('RUST @E=%.3g A=%.3g',E2(e2),A2(a2)));
    xlabel('sputtered energy'); ylabel('PDF (norm by max)');
    title(v + " PDF overlay (nearest point in each)");
    legend('Location','best');
    saveFig(fig2,outDir,"03_energyPDF_"+v,doSavePNGs,dpi);
end

disp("Done. Plots saved in: " + outDir);

%% ===================== Local helper functions =====================

function tf = hasVar(ncfile, varname)
info = ncinfo(ncfile);
names = string({info.Variables.Name});
tf = any(names == string(varname));
end

function v = readVarIf(ncfile, varname)
if hasVar(ncfile,varname)
    v = ncread(ncfile,varname);
    v = v(:);
else
    v = [];
end
end

function idx = nearestIndex(vec, val)
[~,idx] = min(abs(vec - val));
end

function vals = pickVals(vec, n)
vec = vec(:);
if numel(vec) < n, vals = vec'; return; end
q = linspace(0.2,0.8,n);
idx = max(1, min(numel(vec), round(q*(numel(vec)-1)+1)));
vals = vec(idx)';
end

function x = chooseOrIndex(xin, nbins)
if isempty(xin) || numel(xin) ~= nbins, x = (1:nbins)'; else x = xin(:); end
end

function x = pickBinAxis(theta, phi, nbins)
if ~isempty(theta) && numel(theta) == nbins
    x = theta(:);
elseif ~isempty(phi) && numel(phi) == nbins
    x = phi(:);
else
    x = (1:nbins)';
end
end

function m = weightedMean3D(arr, x)
x = x(:);
den = squeeze(sum(arr,3,'omitnan'));
num = squeeze(sum(arr .* reshape(x,1,1,[]), 3, 'omitnan'));
m = num ./ (den + eps);
end

function saveFig(fig, outDir, baseName, doSave, dpi)
if ~doSave, return; end
fn = fullfile(outDir, baseName + ".png");
set(fig,'InvertHardcopy','off');
exportgraphics(fig, fn, 'Resolution', dpi);
end

%% -------- canonical readers for your actual netCDF dim orders --------

function V = readEA_canonical(ncfile, varname)
% Yield vars in your files are stored as (nA,nE) (per ncdisp).
info = ncinfo(ncfile, varname);
dn = string({info.Dimensions.Name});
raw = ncread(ncfile, varname);

if isequal(dn, ["nA","nE"])
    V = raw.';   % -> (nE,nA)
elseif isequal(dn, ["nE","nA"])
    V = raw;     % already (nE,nA)
else
    error("%s %s: unsupported 2D dims %s", ncfile, varname, strjoin(dn,","));
end
end

function C = readEAB_canonical(ncfile, varname)
% 3D vars in your files are stored as (bin,nA,nE) (per ncdisp).
info = ncinfo(ncfile, varname);
dn = string({info.Dimensions.Name});
raw = ncread(ncfile, varname);

if numel(dn) ~= 3
    error("%s %s: expected 3D variable, got %dD", ncfile, varname, numel(dn));
end

% Expect (..., nA, nE) with bin first
if dn(2)=="nA" && dn(3)=="nE"
    C = permute(raw, [3 2 1]);   % (bin,A,E) -> (E,A,bin)
elseif dn(1)=="nE" && dn(2)=="nA"
    C = raw;                     % (E,A,bin) already
else
    error("%s %s: unsupported 3D dims %s", ncfile, varname, strjoin(dn,","));
end
end