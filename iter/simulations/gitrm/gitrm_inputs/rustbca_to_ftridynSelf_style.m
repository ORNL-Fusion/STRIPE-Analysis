function write_ftridynSelf_header_exact(infile, outfile)
% write_ftridynSelf_header_exact
% Create a new NetCDF whose *ncdump -h* header matches ftridynSelf:
%   spyld(nE,nA), rfyld(nE,nA)
%   cos*Dist(nE,nA,nAdistBins)
%   energyDist(nE,nA,nEdistBins)
%   energyDistRef(nE,nA,nEdistBinsRef)
% etc.
%
% IMPORTANT: Because MATLAB reverses dimension order relative to ncdump,
% we intentionally define variables in MATLAB as:
%   spyld(nA,nE)  -> ncdump prints spyld(nE,nA)
%   cosXDist(nAdistBins,nA,nE) -> ncdump prints cosXDist(nE,nA,nAdistBins)
%   energyDist(nEdistBins,nA,nE) -> ncdump prints energyDist(nE,nA,nEdistBins)
%
% Usage:
%   write_ftridynSelf_header_exact("rustbca_as_ftridynSelf_WonW_FIXED.nc", ...
%                                 "rustbca_as_ftridynSelf_WonW_EXACT.nc")

if nargin<1 || isempty(infile), infile="rustbca_as_ftridynSelf_WonW_FIXED.nc"; end
if nargin<2 || isempty(outfile), outfile="rustbca_as_ftridynSelf_WonW_EXACT.nc"; end
infile = char(infile); outfile = char(outfile);

if exist(infile,'file')~=2, error("Input not found: %s", infile); end
if exist(outfile,'file')==2, delete(outfile); end

% --- Read coordinate variables (trust these) ---
E = ncread(infile,"E"); E = E(:); nE = numel(E);
A = ncread(infile,"A"); A = A(:); nA = numel(A);

phiGrid = ncread(infile,"phiGrid");   phiGrid = phiGrid(:);   nAdistBins = numel(phiGrid);
thetaGrid = ncread(infile,"thetaGrid"); thetaGrid = thetaGrid(:);

eDistEgrid = ncread(infile,"eDistEgrid"); eDistEgrid = eDistEgrid(:); nEdistBins = numel(eDistEgrid);
eDistEgridRef = ncread(infile,"eDistEgridRef"); eDistEgridRef = eDistEgridRef(:); nEdistBinsRef = numel(eDistEgridRef);

% --- Read vars and FORCE into canonical MATLAB shapes first:
% canonical we want in MATLAB workspace:
%   spyld_c : (nE,nA)
%   cosX_c  : (nE,nA,nAdistBins)
%   eDist_c : (nE,nA,nEdistBins)
spyld_c = forceEA_2D(ncread(infile,"spyld"), nE, nA, "spyld");
rfyld_c = forceEA_2D(ncread(infile,"rfyld"), nE, nA, "rfyld");

cosX_c = forceEAB_to_shape(ncread(infile,"cosXDist"),    nE, nA, nAdistBins, "cosXDist");
cosY_c = forceEAB_to_shape(ncread(infile,"cosYDist"),    nE, nA, nAdistBins, "cosYDist");
cosZ_c = forceEAB_to_shape(ncread(infile,"cosZDist"),    nE, nA, nAdistBins, "cosZDist");
cosXR_c= forceEAB_to_shape(ncread(infile,"cosXDistRef"), nE, nA, nAdistBins, "cosXDistRef");
cosYR_c= forceEAB_to_shape(ncread(infile,"cosYDistRef"), nE, nA, nAdistBins, "cosYDistRef");
cosZR_c= forceEAB_to_shape(ncread(infile,"cosZDistRef"), nE, nA, nAdistBins, "cosZDistRef");

eDist_c   = forceEAB_to_shape(ncread(infile,"energyDist"),    nE, nA, nEdistBins,    "energyDist");
eDistR_c  = forceEAB_to_shape(ncread(infile,"energyDistRef"), nE, nA, nEdistBinsRef, "energyDistRef");

% --- Now convert to MATLAB-definition order that yields correct ncdump header:
% Define variables in MATLAB as:
%   spyld: (nA,nE)  so ncdump prints (nE,nA)
%   cosX:  (nAdistBins,nA,nE) so ncdump prints (nE,nA,nAdistBins)
%   eDist: (nEdistBins,nA,nE) so ncdump prints (nE,nA,nEdistBins)
spyld_w = spyld_c.';   % (nA,nE)
rfyld_w = rfyld_c.';   % (nA,nE)

cosX_w  = permute(cosX_c, [3 2 1]);   % (nAdistBins,nA,nE)
cosY_w  = permute(cosY_c, [3 2 1]);
cosZ_w  = permute(cosZ_c, [3 2 1]);
cosXR_w = permute(cosXR_c,[3 2 1]);
cosYR_w = permute(cosYR_c,[3 2 1]);
cosZR_w = permute(cosZR_c,[3 2 1]);

eDist_w  = permute(eDist_c, [3 2 1]);   % (nEdistBins,nA,nE)
eDistR_w = permute(eDistR_c,[3 2 1]);   % (nEdistBinsRef,nA,nE)

% --- Create file using HIGH-LEVEL nccreate with reversed dim order ---
% (This is the trick that makes ncdump show ftridynSelf ordering.)
fprintf("Writing: %s\n", outfile);

% coords (1D) are fine
nccreate(outfile,"E","Dimensions",{"nE",nE},"Datatype","double");
nccreate(outfile,"A","Dimensions",{"nA",nA},"Datatype","double");
nccreate(outfile,"eDistEgrid","Dimensions",{"nEdistBins",nEdistBins},"Datatype","double");
nccreate(outfile,"eDistEgridRef","Dimensions",{"nEdistBinsRef",nEdistBinsRef},"Datatype","double");
nccreate(outfile,"phiGrid","Dimensions",{"nAdistBins",nAdistBins},"Datatype","double");
nccreate(outfile,"thetaGrid","Dimensions",{"nAdistBins",nAdistBins},"Datatype","double");

% yields: define as (nA,nE) so ncdump prints (nE,nA)
nccreate(outfile,"spyld","Dimensions",{"nA",nA,"nE",nE},"Datatype","double");
nccreate(outfile,"rfyld","Dimensions",{"nA",nA,"nE",nE},"Datatype","double");

% cos*: define as (nAdistBins,nA,nE) so ncdump prints (nE,nA,nAdistBins)
nccreate(outfile,"cosXDist","Dimensions",{"nAdistBins",nAdistBins,"nA",nA,"nE",nE},"Datatype","double");
nccreate(outfile,"cosYDist","Dimensions",{"nAdistBins",nAdistBins,"nA",nA,"nE",nE},"Datatype","double");
nccreate(outfile,"cosZDist","Dimensions",{"nAdistBins",nAdistBins,"nA",nA,"nE",nE},"Datatype","double");
nccreate(outfile,"cosXDistRef","Dimensions",{"nAdistBins",nAdistBins,"nA",nA,"nE",nE},"Datatype","double");
nccreate(outfile,"cosYDistRef","Dimensions",{"nAdistBins",nAdistBins,"nA",nA,"nE",nE},"Datatype","double");
nccreate(outfile,"cosZDistRef","Dimensions",{"nAdistBins",nAdistBins,"nA",nA,"nE",nE},"Datatype","double");

% energyDist: define as (nEdistBins,nA,nE) so ncdump prints (nE,nA,nEdistBins)
nccreate(outfile,"energyDist","Dimensions",{"nEdistBins",nEdistBins,"nA",nA,"nE",nE},"Datatype","double");
nccreate(outfile,"energyDistRef","Dimensions",{"nEdistBinsRef",nEdistBinsRef,"nA",nA,"nE",nE},"Datatype","double");

% --- Write data ---
ncwrite(outfile,"E",E);
ncwrite(outfile,"A",A);
ncwrite(outfile,"eDistEgrid",eDistEgrid);
ncwrite(outfile,"eDistEgridRef",eDistEgridRef);
ncwrite(outfile,"phiGrid",phiGrid);
ncwrite(outfile,"thetaGrid",thetaGrid);

ncwrite(outfile,"spyld",spyld_w);
ncwrite(outfile,"rfyld",rfyld_w);

ncwrite(outfile,"cosXDist",cosX_w);
ncwrite(outfile,"cosYDist",cosY_w);
ncwrite(outfile,"cosZDist",cosZ_w);
ncwrite(outfile,"cosXDistRef",cosXR_w);
ncwrite(outfile,"cosYDistRef",cosYR_w);
ncwrite(outfile,"cosZDistRef",cosZR_w);

ncwrite(outfile,"energyDist",eDist_w);
ncwrite(outfile,"energyDistRef",eDistR_w);

fprintf("Done.\nNow run: system('ncdump -h %s')\n", outfile);

end

%% ---------- helpers ----------
function V = forceEA_2D(Vin, nE, nA, name)
s = size(Vin);
if isequal(s,[nE nA])
    V = Vin;
elseif isequal(s,[nA nE])
    V = Vin.';    % transpose
else
    error("%s: cannot map size %s to (nE,nA)=(%d,%d)", name, mat2str(s), nE, nA);
end
end

function C = forceEAB_to_shape(Cin, nE, nA, nB, name)
% Force 3D Cin into (nE,nA,nB) via permutation
if isempty(Cin)
    C = zeros(nE,nA,nB);
    return;
end
s = size(Cin);
if numel(s)~=3
    error("%s: expected 3D, got %dD", name, numel(s));
end
perms = {[1 2 3],[1 3 2],[2 1 3],[2 3 1],[3 1 2],[3 2 1]};
for k=1:numel(perms)
    cand = permute(Cin, perms{k});
    if size(cand,1)==nE && size(cand,2)==nA && size(cand,3)==nB
        C = cand;
        return;
    end
end
error("%s: cannot permute %s to (%d,%d,%d)", name, mat2str(s), nE, nA, nB);
end