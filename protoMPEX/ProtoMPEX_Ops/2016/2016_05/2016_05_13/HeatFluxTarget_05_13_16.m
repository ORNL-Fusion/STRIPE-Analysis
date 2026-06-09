clear all 
close all

load('IRData_2016_05_13.mat','TR2','shotlist','Data')
% This data is only saved for frames 30 to 60

%% 
for n = 1:length(shotlist)
    N{n} = length(Data{n}(1,1,:))-1;
    for s = 2:N{n}
       dTdt_a{n}(:,:,s) = Data{n}(:,:,s)  - Data{n}(:,:,s-1);
       dTdt_b{n}(:,:,s) = Data{n}(:,:,s+1)- Data{n}(:,:,s-1);
    end
end
%%
close all

n = find(TR2 >= 320,1);

rngy = [280:480];
rngx = [100:310];
figure;
for s = 1:N{n}
   surf(dTdt_b{n}(rngx,rngy,s),'LineStyle','none')
   set(gcf,'Position',[600 300  230  200])
   view([0,90])
   title(num2str(s))
   zlim([-20,150])
    xlim([0,210])
    ylim([0,215])
    axis('square')
   caxis([0,50])
   colormap(flipud(hot))
       set(gca,'XTickLabel',[])
           set(gca,'YTickLabel',[])
           box on
    pause(0.05)

end

%%
figure
for s = 1:N{n}
   surf(Data{n}(rngx,rngy,s)-Data{n}(rngx,rngy,20),'LineStyle','none')
   set(gcf,'Position',[641.0000  299.0000  279.3333  218.6667])
   view([0,90])
   title(num2str(s))
   zlim([-10,300])
   caxis([0,200])
   colormap(flipud(hot))
   pause(0.01)
end

%%
close all
figure
[~,b] = sort(TR2);
rng = 90:110;

for n = 1:length(TR2)
     HeatFlux{n}(:,:,s) = dTdt_b{n}(rngx,rngy,s);
end


for n = 1:length(b)
    subplot(5,5,n); hold on
    s = 10;
   surf(HeatFlux{b(n)}(:,:,s),'LineStyle','none')
   plot3(rng,rng,HeatFlux{b(n)}(rng,rng,s),'g.')
    PeakHeatFlux(b(n)) = mean(mean(HeatFlux{b(n)}(rng,rng,s),1));
   view([0,90])
   title(num2str(TR2(b(n))))
   zlim([-20,150])
    xlim([0,210])
    ylim([0,215])
    axis('square')
   caxis([0,80])
   colormap(flipud(hot))
   set(gca,'XTickLabel',[])
   set(gca,'YTickLabel',[])
   box on
end

figure; plot(TR2(b),PeakHeatFlux(b))

%%
mdsconnect('mpexserver'); 
[heat_flux8, heat_peak, heat_peak_x, heat_peak_y, heat_power,heat_radius, T, X, Y, comments] = get_IR(18162);
figure; for s = 1:70; imagesc(X*100,Y*100, heat_flux8(:,:,s)*1e-6);
    caxis([0,0.4])
    colormap(flipud(hot))
    axis('square')
    xlim([-3,3])
    ylim([-3,3])
    grid on; 
    title(num2str(T(s)))
    drawnow; pause(0.1)
end

