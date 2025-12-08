function plot4(dataRanges,coord2DTor,coord2DPol,LayerData,plotqtt,plotTitle)

figure('Position',[0,0,1000,500]);
axes('Position',[0.05 0.05 0.4 0.9]);
for i=1:size(dataRanges,1)
    delTri = delaunay(coord2DTor(dataRanges(i,1):dataRanges(i,2)),coord2DPol(dataRanges(i,1):dataRanges(i,2)));
    trisurf(delTri,LayerData(dataRanges(i,1):dataRanges(i,2),1),LayerData(dataRanges(i,1):dataRanges(i,2),2),LayerData(dataRanges(i,1):dataRanges(i,2),3), ...
        plotqtt(dataRanges(i,1):dataRanges(i,2)), ...
        'EdgeColor','none');
    hold on;
end
caxis([min(plotqtt),max(plotqtt)])
hold off;
colorbar('southoutside')
axis equal;
shading interp;
camorbit(-90,0)
title(plotTitle)
for i=1:size(dataRanges,1)
    delTri = delaunay(coord2DTor(dataRanges(i,1):dataRanges(i,2)),coord2DPol(dataRanges(i,1):dataRanges(i,2)));
    if i==1
        axes('Position',[0.76 0.575 0.225 0.4]);
        trisurf(delTri,coord2DTor(dataRanges(i,1):dataRanges(i,2))-mean(coord2DTor(dataRanges(i,1):dataRanges(i,2))), ...
               -coord2DPol(dataRanges(i,1):dataRanges(i,2))+mean(coord2DPol(dataRanges(i,1):dataRanges(i,2))), ...
               coord2DTor(dataRanges(i,1):dataRanges(i,2))*0,plotqtt(dataRanges(i,1):dataRanges(i,2)),'EdgeColor','black');
        caxis([min(plotqtt),max(plotqtt)])
        %xlabel('Distance along limiter (~toroidal)')
        %ylabel('Distance along limiter (~poloidal)')
    elseif i==2
        axes('Position',[0.76 0.1 0.225 0.4]);
        trisurf(delTri,-coord2DTor(dataRanges(i,1):dataRanges(i,2))+mean(coord2DTor(dataRanges(i,1):dataRanges(i,2))), ...
               coord2DPol(dataRanges(i,1):dataRanges(i,2))-mean(coord2DPol(dataRanges(i,1):dataRanges(i,2))), ...
               coord2DTor(dataRanges(i,1):dataRanges(i,2))*0,plotqtt(dataRanges(i,1):dataRanges(i,2)),'EdgeColor','black');
        caxis([min(plotqtt),max(plotqtt)])
        xlabel('Distance along limiter (~toroidal)')
        %ylabel('Distance along limiter (~poloidal)')
    elseif i==3
        axes('Position',[0.5 0.575 0.225 0.4]);
        trisurf(delTri,-coord2DTor(dataRanges(i,1):dataRanges(i,2))+mean(coord2DTor(dataRanges(i,1):dataRanges(i,2))), ...
               -coord2DPol(dataRanges(i,1):dataRanges(i,2))+mean(coord2DPol(dataRanges(i,1):dataRanges(i,2))), ...
               coord2DTor(dataRanges(i,1):dataRanges(i,2))*0,plotqtt(dataRanges(i,1):dataRanges(i,2)),'EdgeColor','black');
        caxis([min(plotqtt),max(plotqtt)])
        %xlabel('Distance along limiter (~toroidal)')
        ylabel('Distance along limiter (~poloidal)')
    else
        axes('Position',[0.5 0.1 0.225 0.4]);
        trisurf(delTri,coord2DTor(dataRanges(i,1):dataRanges(i,2))-mean(coord2DTor(dataRanges(i,1):dataRanges(i,2))), ...
               coord2DPol(dataRanges(i,1):dataRanges(i,2))-mean(coord2DPol(dataRanges(i,1):dataRanges(i,2))), ...
               coord2DTor(dataRanges(i,1):dataRanges(i,2))*0,plotqtt(dataRanges(i,1):dataRanges(i,2)),'EdgeColor','black');
        caxis([min(plotqtt),max(plotqtt)])
        xlabel('Distance along limiter (~toroidal)')
        ylabel('Distance along limiter (~poloidal)')
    end
    view(2)
    axis tight;
    shading interp;
end

end