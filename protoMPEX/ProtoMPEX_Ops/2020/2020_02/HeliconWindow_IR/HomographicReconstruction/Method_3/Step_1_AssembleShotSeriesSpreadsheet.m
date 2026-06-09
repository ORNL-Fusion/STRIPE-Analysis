% Objective:
% =========================================================================
% Assemble the shotSeries spreadsheet to be used for extracting the seq
% files and ultimately producing the dT .mat files

% Context:
% =========================================================================
% The entire data set for the helicon window IR inferred surface heat
% fluxes are to be found in two experimental campaigns:
% - 2020_02_11: Pit side view
% - 2020_02_14: Non-pit side view

% Each of the folders for these experiments have a dedicated spreadsheet
% that provides all the metadata associated with each shot:
% Spreadsheet: HeliconWindowIR_2020_02_11 ~ 44 shots
% Spreadsheet: HeliconWindowIR_2020_02_14 ~ 38 shots

% These two spreadsheets have been put together into a single large
% spreadsheet called the PARENT spreadsheet:
% HeliconWindowIR_2020_02_XPs, which has 68 shots in total

% Not all these shots are required for the heat flux reconstruction; hence,
% we need to assemble a spreadsheet with just the right information we
% need. Hence, we need to create a CHILD spreadsheet. This is precisely
% what this script will do.

clear all
close all

t1 = tic;
% Read child spreadsheet:
% =========================================================================
spreadsheetNameChild = 'Step_1_ShotSeries_HeliconWindowIR_2020_02.xlsx';
Tchild = readtable(spreadsheetNameChild,'Sheet',1);

% Read parent spreadsheet:
% =========================================================================
homeAddress = cd;
cd ..
cd ..
spreadsheetNameParent = 'HeliconWindowIR_2020_02_XPs.xlsx';
Tparent = readtable(spreadsheetNameParent,'Sheet',1);
cd(homeAddress)
t1 = toc(t1);
disp(['Time to read spreadsheets: ',num2str(t1),' [s]'])

% Extract information from Tparent and assemble Tchild:
% =========================================================================
% Find the indices in Tparent of shots requested by Tchild:
 [shotList,nC,nP] = intersect(Tchild.shot,Tparent.shot,'rows','stable');

% Get all the fields Tparent:
fieldNames = fieldnames(Tparent);
 
% Assemble Tchild table:
for ii = 1:numel(fieldNames)
    % All the fields that we do not want:
    if strcmpi('group',fieldNames(ii)) | strcmpi('shot',fieldNames(ii)) ...
            | strcmpi('dataset',fieldNames(ii)) ...
            | strcmpi('properties',fieldNames(ii)) ...
            | strcmpi('viewAngle',fieldNames(ii))...
            | strcmpi('row',fieldNames(ii))...
            | strcmpi('Variables',fieldNames(ii));
        continue
    end
    % Add the fields from parent to child that we want:
    Tchild.(fieldNames{ii}) = Tparent.(ii)(nP);
end

writetable(Tchild,spreadsheetNameChild,'FileType','spreadsheet')