% Template

% =========================================================================
% Entry commands:
clc
close all
disp('Start of task X ##################################################');

%% SECTION 1
% =========================================================================
% Compute prompt:
InputStructure.prompt = {['Put question here? Yes [1], No [0]']};
cmpt = GetUserInput(InputStructure);

% =========================================================================
% Clear memory
if cmpt
    clearvars('-except','cmpt');
end

% #########################################################################
% SECTION 1:
sectionName = 'Load .mat file';
disp(['Section: ',sectionName])
% #########################################################################

% =========================================================================
% Compute
if cmpt
    % Add process here
end

% =========================================================================
% Define variables specific to this section:
if cmpt
    vars.section1 = '';
    vars.section1 = who;
    vars.section1 = setdiff(vars.section1,[]);
end

%% SECTION 2
% =========================================================================
% Clear memory and figures:
try
    clearvars('-except',vars.section1{:});
    close(findobj('Tag',sectionName))
end

% #########################################################################
% SECTION 2:
sectionName = 'Name';
disp(['Section: ',sectionName])
% #########################################################################

% =========================================================================
% Compute prompt:
InputStructure.prompt = {['Question? Yes [1], No [0]']};
cmpt = GetUserInput(InputStructure);

% =========================================================================
% Compute
if cmpt
    % ---------------------------------------------------------------------
    disp('Computing...') 
    tic
    % Add compute process here
    toc
else
    % ---------------------------------------------------------------------
    disp('Loading precalculated data...')
    fileName = ['Put file name here'];
    load(fileName)
    disp('Data loaded succesfully')
end

% =========================================================================
% Some calculation
% Add calculation here

% =========================================================================
% Define variables specific to this section:
if cmpt
    vars.section2 = '';
    vars.section2 = who;
    vars.section2 = setdiff(vars.section2,vars.section1);
end

% =========================================================================
% Save prompt:
if cmpt
    InputStructure.prompt = {['Would you like to save? Yes [1], No [0]']};
    InputStructure.option.WindowStyle = 'normal';
    svdt = GetUserInput(InputStructure);
    if svdt        
        variableNames = vars.section2;
        fileName = ['Put file name here'];
        SaveData
    end
end

%% SECTION 3
% =========================================================================
% Clear memory and figures:
try
    clearvars('-except',vars.section1{:},vars.section2{:});
    close(findobj('Tag',sectionName))
end

% #########################################################################
% SECTION 3:
sectionName = 'Put section name here';
disp(['Section: ',sectionName])
% #########################################################################

% =========================================================================
% Compute prompt:
InputStructure.prompt = {['Put Question here?, Yes [1], No [0]']};
cmpt = GetUserInput(InputStructure);

% =========================================================================
% Compute:
if cmpt
    % ---------------------------------------------------------------------
    disp('Computing...')
    tic
    % Add process here
    toc
end
 
% =========================================================================
% Compute:
if cmpt
    % ---------------------------------------------------------------------
    disp('Computing...')
    tic
    % Add compute process
    toc
    
    % ---------------------------------------------------------------------
    % Define variables specific to this section:
    vars.section3 = '';
    vars.section3 = who;
    vars.section3 = setdiff(vars.section3,[vars.section2;vars.section1]);
    
    % ---------------------------------------------------------------------
    % Save prompt:
    disp('Save prompt...')
    InputStructure.prompt = {['Would you like to save? Yes [1], No [0]']};
    InputStructure.option.WindowStyle = 'normal';
    svdt = GetUserInput(InputStructure);
    if svdt        
        variableNames = vars.section3;
        fileName = ['Put file name here'];
        SaveData
    else
        disp('User choose not to save data')
    end
else
    % ---------------------------------------------------------------------
    % Load existing data
    disp('Loading existing data...')
    fileName = ['Put file name here'];
    load(fileName)
    disp('Data loaded succesfully')
end

% =========================================================================
% Create figure:
figure('Tag',sectionName)
figure('Tag','section2','color','w')
ylabel('$\Delta T$ [C]','Interpreter','latex','FontSize',12)
set(gca,'FontName','times','FontSize',11)

% =========================================================================
% Save figure
saveFig = 1;

if saveFig
    disp('Saving figures...')
    hdum = findobj('Tag','GDT B profile');
    saveas(hdum,hdum.Tag,'tiffn')
end


% =========================================================================
% End of script
disp('Task X complete!!')