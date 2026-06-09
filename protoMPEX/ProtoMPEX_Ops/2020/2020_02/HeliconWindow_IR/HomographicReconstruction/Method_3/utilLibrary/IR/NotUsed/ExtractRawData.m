function [RawData,frames,seq, frame_rate] = ExtractRawData(PATHNAME,FILENAME,shot)
%%%This function reads in the IR camera data and returns RawData which is a
%%%3D array containing a 2D image from the IR camera over frames amount of
%%%pictures

%%%%%%%%%%%%%%%%%%%READ RAW INTENSITY DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Read in video File
videoFileName=[PATHNAME FILENAME];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Get path info for FLIR
asmInfo = NET.addAssembly([getenv('FLIR_Atlas_MATLAB') 'Flir.Atlas.Image.dll']);
%%%Get the sequence file
seq = Flir.Atlas.Image.ThermalImageFile(videoFileName).ThermalSequencePlayer();
%%%Get Pixel Size
Pixels = size(double(seq.ThermalImage.ImageProcessing.GetPixelsArray));
%%%Preallocate size of the data for a little faster processing 
RawData = zeros(Pixels(1),Pixels(2),seq.ThermalImage.Count);
frame_rate = round( (numel(RawData(1,1,:)) - 1) / seq.Duration.TotalSeconds);
if frame_rate > 55
    frame_rate = 100;
end
RawData(:,:,1) = double(seq.ThermalImage.ImageProcessing.GetPixelsArray);
%%Loop through frame counting
frames = 1;
if(seq.Count > 1)
    while(seq.Next())
        frames = frames + 1;
        RawData(:,:,frames) = double(seq.ThermalImage.ImageProcessing.GetPixelsArray);
    end
end
%%%%%%%Subtract Background
%RawData = RawData-RawData(:,:,1); %%%I'm not sure if this is kosher
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

