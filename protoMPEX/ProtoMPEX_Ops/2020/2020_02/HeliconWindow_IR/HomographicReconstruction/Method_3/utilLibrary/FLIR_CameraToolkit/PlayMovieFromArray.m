function [flag] = PlayMovieFromArray(A,options)
% This function is used to provide a quick preview of the IR data during 
% experiments or data analysis.

% A: array with 3 dimensions
% options:
% .frames = n1:n2 or "all"
% .magnitudePlotMode = 1 or 2, 1: delta, 2: absolute
% .colorbar = 0 or 1
% .removeAxesTicks = 1 or 0
% .mirrorImage == 1 or 0
% .shot
% .zlim = [0,12];

% Example for options:
% options.frames = 1:3:size(A,3);
% options.colorbar = 1;
% options.magnitudePlotMode = 1;
% options.removeAxesTicks = 0;
% options.shot = 0000;
% options.zlim = [0,12];
% options.mirrorImage = 1;

% Check data size:
% -------------------------------------------------------------------------
if numel(size(A))<3
    flag = 0;
    disp('Cannot create movie, incorrect size of input Array')
    return
end

% Check mirroring option:
% -------------------------------------------------------------------------
if isfield(options,'mirrorImage')
    if options.mirrorImage == 1
        for ii = 1:size(A,3)
            A(:,:,ii) = A(:,end:-1:1,ii);
        end
    end
end

% Define size of array to plot:
% -------------------------------------------------------------------------
ny = size(A,1);
nx = size(A,2);
if isfield(options,'zlim')
    maxA = options.zlim(2);
    minA = options.zlim(1);
else
    maxA = max(max(max(A)));
    minA = min(min(min(A(:,:,1:end-1))));
end

% Check frame option:
% -------------------------------------------------------------------------
if isfield(options,'frames')
    if strcmpi(options.frames,'all')
        frameRng = 1:size(A,3);
    else
        frameRng = options.frames;
    end
else
    frameRng = 1:size(A,3);
end

% Plot data:
% -------------------------------------------------------------------------
figure('color','w')
for ii = frameRng
    switch options.magnitudePlotMode
        case 1
        surf(A(:,:,ii)-A(:,:,1),'LineStyle','none')
        caxis([0,maxA-minA])
        case 2
        surf(A(:,:,ii),'LineStyle','none')
        caxis([minA,maxA])
        case 3
        surf(A(:,:,ii+1)-A(:,:,ii),'LineStyle','none')
        caxis([minA,maxA])
    end
    
    % Set default plot options:
    axis('equal')
    view([0,90])
    xlim([0,nx])
    ylim([0,ny])
    colormap(flipud(hot))
        
    % Modify output based on user-defined contents of "options":   
    if isfield(options,'removeAxesTicks')
        if options.removeAxesTicks == 1
            set(gca,'XTick',[],'YTick',[]) 
        end
    end
    
    % Colobar option:
    if isfield(options,'colorbar')
       if options.colorbar == 1
            colorbar
       end
    end
             
    % title option:
    if isfield(options,'shot')
        title(['Frame: ',num2str(ii),' , shot: ',num2str(options.shot)])
    else
        title(['Frame: ',num2str(ii)])
    end
        
    % Plot frame:
    drawnow
end


end

