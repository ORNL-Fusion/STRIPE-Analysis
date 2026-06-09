disp(viewType)
disp('Loading .seq files...')
tic
nshots = length(shot);
for s = 1:nshots
    % Load the Atlats SDK
    atPath = getenv('FLIR_Atlas_MATLAB');
    atImage = strcat(atPath,'Flir.Atlas.Image.dll');
    asmInfo = NET.addAssembly(atImage);
    %open the IR-file'
    PATHNAME = [a{addressLoc(s)},'\'];
    FILENAME = ['Shot ',num2str(shot(s)),'.seq'];
    videoFileName=[PATHNAME FILENAME];
    file = Flir.Atlas.Image.ThermalImageFile(videoFileName);

    % Define the seq file and associated thermal parameters
    seq{s} = file.ThermalSequencePlayer();
    seq{s}.ThermalImage.ThermalParameters.ExternalOpticsTransmission = 0.7;
    seq{s}.ThermalImage.ThermalParameters.AtmosphericTemperature = 24;
    seq{s}.ThermalImage.ThermalParameters.Distance = 1;
    seq{s}.ThermalImage.ThermalParameters.ExternalOpticsTemperature = 24;
    seq{s}.ThermalImage.ThermalParameters.ReferenceTemperature = 24;
    seq{s}.ThermalImage.ThermalParameters.Transmission = 1;
    seq{s}.ThermalImage.ThermalParameters.RelativeHumidity = 0;
    seq{s}.ThermalImage.ThermalParameters.ReflectedTemperature = 24;

    %Get the pixels
    img = seq{s}.ThermalImage.ImageProcessing.GetPixelsArray;
    im = double(img);

    intensityRaw{s}(:,:,1) = im;
    fr = 1;
    if(seq{s}.Count > 1)
        while(seq{s}.Next())
            img = seq{s}.ThermalImage.ImageProcessing.GetPixelsArray;
            im = double(img);
            intensityRaw{s}(:,:,fr) = im(end:-1:1,:);         
            fr = fr + 1;
        end
    end
end
initialVars{end+1} = 'intensityRaw';
toc
disp('.seq files loaded!!'); 
msgbox('.seq files loaded!!'); 