function [intensity,t_intensity,seq] = ExtractDataFromSeqFile(fileName,pathName,thermalParam,extractOptions)
% ExtractDataFromSeqFile:
% =========================================================================
% The objective of this function is to extract the data from a FLIR
% radiometric file for the purpose of numerical analysis.
% Given a fileName and a pathName, this function will extract 3D array of
% double precision floating numbers from the associated .seq file. In
% addition, the radiometric data will be
% extracted based on the contents of given the structure thermalParam.
% The effects of mirrors and associated mirror images are not taken into
% acccount in this function. The effect of mirrors needs to be incorporated
% outside of this function.
% INPUTS:
% =========================================================================
%  - thermaParam:
% thermalParam.ExternalOpticsTransmission = 0.7;
% thermalParam.AtmosphericTemperature = 24;
% thermalParam.Distance = 1;
% thermalParam.ExternalOpticsTemperature = 24;
% thermalParam.ReferenceTemperature = 24;
% thermalParam.Transmission = 1;
% thermalParam.RelativeHumidity = 0;
% thermalParam.ReflectedTemperature = 24;
%  - fileName: name of the .seq file to use such as "shot_23465.seq"
%  - pathName: location of file, such as C:\Users\nfc\Documents\ProtoMPEX_Ops\2020\2020_02\2020_02_05\IR_Data
% OUTPUTS:
%  - intensityRaw: 3D array of double precision numbers representing the
%  raw output of the .seq file. The output has a 16 bit resolution so each
%  pixel can take values from 0 to 65535
%  - seq: FLIR object or class that contains all the information on the
%  radiometric file. This file is required to computing the surface
%  temperature given an emissivity.

% Load the Atlats SDK:
% =========================================================================
atPath = getenv('FLIR_Atlas_MATLAB');
atImage = strcat(atPath,'Flir.Atlas.Image.dll');
asmInfo = NET.addAssembly(atImage);

% Get file:
% =========================================================================
videoFileName=[pathName fileName];
file = Flir.Atlas.Image.ThermalImageFile(videoFileName);

% Define the seq file and associated thermal parameters:
% =========================================================================
try
    seq = file.ThermalSequencePlayer();
    seq.ThermalImage.ThermalParameters.ExternalOpticsTransmission = thermalParam.ExternalOpticsTransmission;
    seq.ThermalImage.ThermalParameters.AtmosphericTemperature = thermalParam.AtmosphericTemperature;
    seq.ThermalImage.ThermalParameters.Distance = thermalParam.Distance;
    seq.ThermalImage.ThermalParameters.ExternalOpticsTemperature = thermalParam.ExternalOpticsTemperature;
    seq.ThermalImage.ThermalParameters.ReferenceTemperature = thermalParam.ReferenceTemperature;
    seq.ThermalImage.ThermalParameters.Transmission = thermalParam.Transmission;
    seq.ThermalImage.ThermalParameters.RelativeHumidity = thermalParam.RelativeHumidity;
    seq.ThermalImage.ThermalParameters.ReflectedTemperature = thermalParam.ReflectedTemperature;
    disp('Using user defined Thermal Parameters...')
catch
    disp('Warning: Default Thermal Parameters selected')
end

% Get the pixels and assign to 3D array of double precision numbers:
% =========================================================================
if isfield(extractOptions,'frames')
    frames = extractOptions.frames;
else
    frames = 1:1:500;
end

frameRate = extractOptions.frameRate;
t_intensity = (frames-1)/frameRate;

% #######
% Change made on 2020_12_06 while analysing data from 2017
% intensity = zeros(240,640,numel(frames));
intensity = zeros(seq.ThermalImage.Size.Height,seq.ThermalImage.Size.Width,numel(frames));
% #######

fr = 1;
ii = 1;
if(seq.Count > 1)
    while(seq.Next())
        img = seq.ThermalImage.ImageProcessing.GetPixelsArray;
        im = img.double;
        
        if fr>=frames(1) && fr<=frames(end)
            % Image needs to be flipped row-wise to allow correct plotting
            % of data when using surf. This process has nothing to do with
            % the use of mirrors or their effects on images:
            intensity(:,:,ii) = im(end:-1:1,:);
            ii = ii + 1;
            if fr == frames(end)
                break
            end
        end
        fr = fr + 1;
    end
end

% Note:
% The processes that takes the longest in this script are the following: 
% 1- intensityRaw(:,:,fr) = im(end:-1:1,:);  (2.75 sec)
% 2- img = seq.ThermalImage.ImageProcessing.GetPixelsArray; (1.70 sec)
% 3- im = img.double; (0.91 sec)

end

