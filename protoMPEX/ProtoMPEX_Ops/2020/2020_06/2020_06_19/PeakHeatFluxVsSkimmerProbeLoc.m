% Peak heat flux as a fuction of DLP 6.5 position:

rg       = [2.0,5.0,5.5,6.0,6.5 ,7.0 ];
heatFlux = [19 ,15 ,15 ,14 ,12.5,10.0];

figure('color','w')
plot(rg,heatFlux,'ko')
ylabel('Peak heat flux [MWm^{-2}]')
ylim([0,20])
xlim([0,15])
grid on