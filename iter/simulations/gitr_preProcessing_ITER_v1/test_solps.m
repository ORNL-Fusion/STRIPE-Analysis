%% SOLPS Plasma Profile Exporter (D + Neon) - Robust Version
clear; clc;

% === Input files ===
fstate = 'b2fstate.mat';
fstati = 'b2fstati.mat';
fgmtry = 'b2fgmtry.mat';
outfile = 'solps_profiles.nc';

load(fstate);
load(fstati);
load(fgmtry);

%% --- Geometry extraction ---
if exist('b2fgmtry','var')
    geom = b2fgmtry;
elseif exist('gmtry','var')
    geom = gmtry;
else
    geom = struct();
end

R = []; Z = [];
if isfield(geom,'rm'), R = geom.rm; end
if isfield(geom,'zm'), Z = geom.zm; end
if isempty(R) && isfield(geom,'crx'), R = mean(geom.crx,3); end
if isempty(Z) && isfield(geom,'cry'), Z = mean(geom.cry,3); end

nx = size(R,1); ny = size(R,2);

%% --- Load plasma state ---
if exist('b2fstate','var')
    st = b2fstate;
else
    st = struct();
end

Te = []; Ti = []; ne = []; ni = []; upar = [];
if isfield(st,'te'), Te = st.te; end
if isfield(st,'ti'), Ti = st.ti; end
if isfield(st,'ne'), ne = st.ne; end
if isfield(st,'ni'), ni = st.ni; end
if isfield(st,'upar'), upar = st.upar; end

%% --- Identify species ---
names = {};
if isfield(st,'iznam')
    raw = st.iznam;
    if iscell(raw), names = raw;
    elseif isstring(raw) || ischar(raw), names = cellstr(raw); end
elseif exist('b2fstati','var') && isfield(b2fstati,'iznam')
    raw = b2fstati.iznam;
    if iscell(raw), names = raw;
    elseif isstring(raw) || ischar(raw), names = cellstr(raw); end
end

if isempty(names) && ndims(ni)==3 && size(ni,1)==2
    names = {'D','Ne'};
end

d_idx  = find(contains(upper(string(names)),"D"));
ne_idx = find(contains(upper(string(names)),"NE"));

%% --- Extract per-species data ---
if ndims(ni)==3
    if size(ni,1) < size(ni,3)
        nD  = squeeze(sum(ni(d_idx,:,:),1));
        nNe = squeeze(sum(ni(ne_idx,:,:),1));
        if ~isempty(upar)
            uD  = squeeze(mean(upar(d_idx,:,:),1));
            uNe = squeeze(mean(upar(ne_idx,:,:),1));
        else
            uD = []; uNe = [];
        end
    else
        nD  = squeeze(sum(ni(:,:,d_idx),3));
        nNe = squeeze(sum(ni(:,:,ne_idx),3));
        if ~isempty(upar)
            uD  = squeeze(mean(upar(:,:,d_idx),3));
            uNe = squeeze(mean(upar(:,:,ne_idx),3));
        else
            uD = []; uNe = [];
        end
    end
else
    nD = []; nNe = []; uD = []; uNe = [];
end

%% --- Create NetCDF using low-level interface ---
if exist(outfile,'file'), delete(outfile); end

ncid = netcdf.create(outfile,'NETCDF4');

% Define dimensions once
dimid_x = netcdf.defDim(ncid,'x',nx);
dimid_y = netcdf.defDim(ncid,'y',ny);

% Define variables
var_R      = netcdf.defVar(ncid,'R','double',[dimid_x dimid_y]);
var_Z      = netcdf.defVar(ncid,'Z','double',[dimid_x dimid_y]);
var_ne     = netcdf.defVar(ncid,'ne','double',[dimid_x dimid_y]);
var_Te     = netcdf.defVar(ncid,'Te','double',[dimid_x dimid_y]);
var_Ti     = netcdf.defVar(ncid,'Ti','double',[dimid_x dimid_y]);
var_nD     = netcdf.defVar(ncid,'n_D','double',[dimid_x dimid_y]);
var_nNe    = netcdf.defVar(ncid,'n_Ne','double',[dimid_x dimid_y]);
var_uD     = netcdf.defVar(ncid,'upar_D','double',[dimid_x dimid_y]);
var_uNe    = netcdf.defVar(ncid,'upar_Ne','double',[dimid_x dimid_y]);

% Global attributes
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'title','SOLPS Plasma Profiles (D + Neon)');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'history',['Created ' datestr(now) ' by MATLAB']);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'institution','ORNL');
if ~isempty(names)
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'species_names',strjoin(names,', '));
end

% End define mode
netcdf.endDef(ncid);

% Write data
netcdf.putVar(ncid,var_R,R);
netcdf.putVar(ncid,var_Z,Z);
if ~isempty(ne),  netcdf.putVar(ncid,var_ne,ne); end
if ~isempty(Te),  netcdf.putVar(ncid,var_Te,Te); end
if ~isempty(Ti),  netcdf.putVar(ncid,var_Ti,Ti); end
if ~isempty(nD),  netcdf.putVar(ncid,var_nD,nD); end
if ~isempty(nNe), netcdf.putVar(ncid,var_nNe,nNe); end
if ~isempty(uD),  netcdf.putVar(ncid,var_uD,uD); end
if ~isempty(uNe), netcdf.putVar(ncid,var_uNe,uNe); end

% Close file
netcdf.close(ncid);

disp('✅ NetCDF successfully created: solps_profiles.nc');