% Testing code: read the McPherson data:
close all
clear all

FileName = 'D2_1468_20um_11.SPE';
image = readSPE(FileName);
% we could use the dirPath by connecting directly to the Server and
% therefore no need to store the data in my computer.

% Fiber, wavelength, time

figure; hold on
for s = 1:50
    plot(image(1,:,s))
end
set(gca,'Yscale','log')