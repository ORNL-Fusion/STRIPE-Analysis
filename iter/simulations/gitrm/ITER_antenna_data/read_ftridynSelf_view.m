function read_ftridynSelf_view(file)
% read_ftridynSelf_view
%
% Reader + quick plots for ftridynSelf.nc (per your ncdisp)
%
% Variables/dims:
%   E [nE=50], A [nA=40]
%   spyld, rfyld           [nA x nE]
%   energyDist             [nEdistBins=100 x nA x nE]
%   energyDistRef          [nEdistBinsRef=500 x nA x nE]
%   eDistEgrid, eDistEgridRef
%   cosX/Y/ZDist, cosX/Y/ZDistRef [nAdistBins=50 x nA x nE]
%   phiGrid, thetaGrid     [nAdistBins=50]
%
% Also computes:
%   RE_energy(angle,energy) = rfyld * <E_ref>/E_inc
% from energyDistRef + eDistEgridRef.
%
% Usage:
%   read_ftridynSelf_view('ftridynSelf.nc');
%   % or simply:
%   read_ftridynSelf_view;

    if nargin < 1 || isempty(file)
        file = 'ftridynSelf.nc';
    end

    %% ---- Read core grids/yields ----
    E = ncread(file,'E'); E = E(:);        % [50]
    A = ncread(file,'A'); A = A(:);        % [40]
    spyld = ncread(file,'spyld');          % [nA x nE]
    rfyld = ncread(file,'rfyld');          % [nA x nE]

    % sanity
    if ~isequal(size(spyld), [numel(A) numel(E)])
        error('spyld size %s does not match [numel(A) numel(E)] = [%d %d].', ...
            mat2str(size(spyld)), numel(A), numel(E));
    end

    %% ---- Read distributions + grids (optional but available) ----
    eDistEgrid    = ncread(file,'eDistEgrid');    eDistEgrid = eDistEgrid(:);       % [100]
    eDistEgridRef = ncread(file,'eDistEgridRef'); eDistEgridRef = eDistEgridRef(:); % [500]

    energyDist    = ncread(file,'energyDist');      % [100 x nA x nE]
    energyDistRef = ncread(file,'energyDistRef');   % [500 x nA x nE]

    phiGrid   = ncread(file,'phiGrid');   phiGrid = phiGrid(:);     % [50]
    thetaGrid = ncread(file,'thetaGrid'); thetaGrid = thetaGrid(:); % [50]

    cosXDist    = ncread(file,'cosXDist');    % [50 x nA x nE]
    cosYDist    = ncread(file,'cosYDist');
    cosZDist    = ncread(file,'cosZDist');
    cosXDistRef = ncread(file,'cosXDistRef');
    cosYDistRef = ncread(file,'cosYDistRef');
    cosZDistRef = ncread(file,'cosZDistRef');

    %% ---- Compute energy reflection yield RE(angle,energy) ----
    % RE = rfyld * (<E_ref>/E_inc)
    RE = zeros(numel(A), numel(E));
    for ia = 1:numel(A)
        for ie = 1:numel(E)
            p = energyDistRef(:, ia, ie);
            s = sum(p);
            if s > 0
                p = p./s;
                meanEref = sum(p .* eDistEgridRef);
            else
                meanEref = 0;
            end
            if E(ie) > 0
                RE(ia, ie) = rfyld(ia, ie) * (meanEref / E(ie));
            else
                RE(ia, ie) = 0;
            end
        end
    end
    RE = max(min(RE,1),0);

    %% ===================== Plots =====================

    % 1) Sputtering yield map
    figure('Name','spyld (map)','Color','w');
    h = pcolor(E, A, spyld); h.EdgeColor = 'none';
    set(gca,'XScale','log','ColorScale','log','YDir','normal');
    xlabel('Energy [eV]'); ylabel('Angle [deg]');
    title('Sputtering yield: spyld(A,E)');
    colorbar; grid on; box on;

    % 2) Reflection yield map
    figure('Name','rfyld (map)','Color','w');
    h = pcolor(E, A, rfyld); h.EdgeColor = 'none';
    set(gca,'XScale','log','YDir','normal');
    xlabel('Energy [eV]'); ylabel('Angle [deg]');
    title('Reflection yield: rfyld(A,E)');
    colorbar; grid on; box on;

    % 3) Energy reflection yield map (computed)
    figure('Name','RE_energy (computed map)','Color','w');
    h = pcolor(E, A, RE); h.EdgeColor = 'none';
    set(gca,'XScale','log','YDir','normal');
    xlabel('Energy [eV]'); ylabel('Angle [deg]');
    title('Energy reflection yield: RE = rfyld * <E_{ref}>/E_{inc}');
    colorbar; grid on; box on;

    % 4) Linecuts vs Energy at a few angles
    angles_to_plot = [0 15 30 45 60 75];
    aIdx = nearest_idx(A, angles_to_plot);

    figure('Name','spyld vs Energy (selected angles)','Color','w'); hold on;
    for k = 1:numel(aIdx)
        plot(E, spyld(aIdx(k),:), 'LineWidth', 1.6);
    end
    set(gca,'XScale','log','YScale','log');
    xlabel('Energy [eV]'); ylabel('spyld');
    title('spyld vs Energy at selected angles');
    legend(compose('A=%.1f°', A(aIdx)), 'Location','best');
    grid on; box on;

    % 5) Linecuts vs Angle at a few energies
    energies_to_plot = [20 50 100 500 1000 5000];
    eIdx = nearest_idx(E, energies_to_plot);

    figure('Name','spyld vs Angle (selected energies)','Color','w'); hold on;
    for k = 1:numel(eIdx)
        plot(A, spyld(:,eIdx(k)), 'LineWidth', 1.6);
    end
    set(gca,'YScale','log');
    xlabel('Angle [deg]'); ylabel('spyld');
    title('spyld vs Angle at selected energies');
    legend(compose('E=%.0f eV', E(eIdx)), 'Location','best');
    grid on; box on;

    %% ---- Optional: visualize one reflected-energy distribution ----
    % pick representative indices
    [~, ie1] = min(abs(E - 1000));
    [~, ia1] = min(abs(A - 45));

    pRef = energyDistRef(:, ia1, ie1);
    if sum(pRef) > 0, pRef = pRef./sum(pRef); end

    figure('Name','energyDistRef example','Color','w');
    plot(eDistEgridRef, pRef, 'LineWidth', 1.8);
    set(gca,'XScale','log');
    xlabel('Reflected energy bin [eV]');
    ylabel('PDF (normalized)');
    title(sprintf('energyDistRef (A=%.1f°, E=%.0f eV)', A(ia1), E(ie1)));
    grid on; box on;

    %% ---- Dump a struct to base workspace (optional) ----
    S = struct();
    S.file = file;
    S.E = E; S.A = A;
    S.spyld = spyld; S.rfyld = rfyld; S.RE = RE;
    S.eDistEgrid = eDistEgrid; S.eDistEgridRef = eDistEgridRef;
    S.energyDist = energyDist; S.energyDistRef = energyDistRef;
    S.phiGrid = phiGrid; S.thetaGrid = thetaGrid;
    S.cosXDist = cosXDist; S.cosYDist = cosYDist; S.cosZDist = cosZDist;
    S.cosXDistRef = cosXDistRef; S.cosYDistRef = cosYDistRef; S.cosZDistRef = cosZDistRef;

    assignin('base','ftridynSelf',S);
    fprintf('Loaded file and assigned struct "ftridynSelf" in base workspace.\n');

end

function idx = nearest_idx(grid, values)
    idx = zeros(size(values));
    for k=1:numel(values)
        [~, idx(k)] = min(abs(grid - values(k)));
    end
    idx = unique(idx);
end