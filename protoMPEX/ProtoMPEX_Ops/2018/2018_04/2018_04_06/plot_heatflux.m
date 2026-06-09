%% plot heatflux
% note this script assumes you are in the folder the file is. Modify
% appropiately
shot =21047;
M = csvread(strcat('heatflux_',num2str(shot), '.csv'), 9);
M = M(:,3:end);
x_size=sqrt(size(M,1));
y_size=sqrt(size(M,1));
t_size = size(M,2);
HeatFlux = reshape(M, x_size, y_size, t_size);
HeatFlux = permute(HeatFlux, [2,1,3]);
figure
imagesc(HeatFlux(:,:,t_size/2)/1e6);
c=colorbar; c.Label.String = 'Heat Flux (MW/m^2)'; caxis([0,1.5]);  
xlabel('Pixel'); ylabel('Pixel');
colormap(flipud(hot));