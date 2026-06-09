% Step 1, Read fast visible camera videos

% Strategy:
% Every shot consumes a large amount of memory and thus we need to perform
% the analysis shot per shot and clearing memory between successive shots.
% All we really need is statistics on the number of speckles

% use shot numners ot dfine the time interval to pliot and suit the figure

clear all
close all

% Define shot:
% =========================================================================
shot = 30000 + [290];
% shot = 13840;
% shot = 13809;

% For shot 13840, we need to crop the area of integration since the
% increased intensity caused by the probe messes up the calculation

CMPT = 1;

t0 = tic;
if CMPT
    % Import data:
    % =========================================================================
    disp('Importing data...')
    tic
    % Takes about 9 sec per shot
    % Each video is 1.3 GB
    % Deleting 2 of the 3 "color" lead to 0.45 GB per video
    d1 = importdata(['fvc3_',num2str(shot),'.mov']);
    f  = permute(d1,[1 2 4 3]);
    clear d1
    f(:,:,:,2:3) = [];
    toc
    disp('Importing data completed!')

    % Convert to double:
    % =========================================================================
    disp('convert to intergers')
    tic
    % Takes 1.5 sec per video
    % File size increases to 1.78 GB each 
    v = cast(f,'single');
    toc

    disp('Compute intensity')
    tic
    % Compute intensity:
    % =========================================================================
    % Takes 2.2 sec per video
    for ff = 1:size(f,3)
        inten(ff) = sum(sum(v(:,:,ff)));
    end

    toc

    % Compute gradient of intensity:
    % =========================================================================
    disp('Compute gradient')
    tic
    % Takes 33 sec
    % 1.2 GB
    for ff = 1:size(f,3)
        [fx(:,:,ff),~] = gradient(v(:,:,ff));
    end
    toc
    
    % Importing data from server:
    % =========================================================================
    disp('Get data from server')
    tic
    % Takes 7 sec
    Stem = '\MPEX::TOP.';
    Branch = 'MACHOPS1:';
    RootAddress = [Stem,Branch];
    DA{1} = [RootAddress,'PWR_28GHZ']; % Isx
    [ECH,t_ECH]   = my_mdsvalue_v2(shot,DA(1));
    DA{1} = [RootAddress,'RF_FWD_PWR'];
    [RF,t_RF]   = my_mdsvalue_v2(shot,DA(1));
    DA{1} = [RootAddress,'PS1_I'];
    [PS1,t_PS1]   = my_mdsvalue_v2(shot,DA(1));
    address{1} = '\MPEX::TOP.MACHOPS1:PG3';
    [PG6,t_PG6]   = my_mdsvalue_v2(shot,address(1));
    address{1} = '\MPEX::TOP.MACHOPS1:SXR_01';
    [SXR1,t_SXR1]   = my_mdsvalue_v2(shot,address(1));
    address{1} = '\MPEX::TOP.MACHOPS1:SXR_02';
    [SXR2,t_SXR2]   = my_mdsvalue_v2(shot,address(1));    
    toc

    % Count speckles:
    % =========================================================================
    % On the fvC data, RF starts on frame:
    fr_RF = find(inten > 1e6,1);

    % RF start time [s]:
    y = abs(RF{1})./max(abs(RF{1}));
    n_rfStart = find(y > 0.05,1);

    t_rfStart = t_RF{1}(n_rfStart);

    % Frame time resolution:
    dt = 1e-3;
    t_nSpeck = ([1:1:size(v,3)]')*dt;

    disp('Count speckles')
    tic
    % Calculate speckle count
    % Takes 4 sec
    try
        clear aa
    end

    nSpeck = zeros(size(v,3),1);
    for ff = 1:size(v,3)
        [xx{ff},~,aa{ff}] = find(fx(:,:,ff) > 35);

    %     if ff > 2
    %         [~,nn{ff}] = intersect(xx{ff-2},xx{ff});      
    %         if ~isempty(nn{ff})
    %             aa{ff} = [];
    %         end
    %     end

        nSpeck(ff) = numel(aa{ff});
    end
    toc
    
    varList2Keep = {'shot','inten','nSpeck','t_nSpeck',...
        'ECH','t_ECH','RF','t_RF','PS1','t_PS1','PG6','t_PG6','t_rfStart'};

    disp('Saving data ...')
    tic
    fileName = ['fvc3_Data_shot_',num2str(shot),'.mat'];
    save(fileName,varList2Keep{:});
    toc
    disp('Data saved!')
else
    disp('Loading data ...')
    tic
    fileName = ['fvc3_Data_shot_',num2str(shot),'.mat'];
    load(fileName)
    toc
    disp('Data loaded!')
end
disp(['Total time taken to CMPT/LOAD data:'])
t0 = toc(t0)

%%

figure; 
plot(t_nSpeck,nSpeck)
% Plot data:
% =========================================================================
figure
hold on
plot(t_RF{1}(1:end-1),abs(RF{1}))
plot(t_ECH{1}(1:end-1),sgolay_t(ECH{1},3,21))
plot(t_nSpeck + t_rfStart - t_nSpeck(fr_RF) + 0.0135,inten/max(inten))
plot(t_nSpeck + t_rfStart - t_nSpeck(fr_RF) + 0.0135,2*nSpeck/max(nSpeck))
plot(t_PG6{1}(1:end-1),PG6{1})
plot(t_PS1{1}(1:end-1),PS1{1})
xlim([4.11,4.8])
title(['shot: ',num2str(shot)])