% Fit Preview, created on 2018_06_20
if FitShow
    for k = 1:12
         figure;
    for c = (1 + (k-1)*25):(25 + (k-1)*25);
        subplot(5,5,c- (k-1)*25 ); hold on
        try
        plot(DLPData.Vsweep{s}{c},DLPData.Isweep{s}{c}*1e3,'k')
        if GoodFits{s}(c) == 1
            plot(DLPData.Vsweep{s}{c},DLPData.Ifit{s}{c}*1e3,'r')
        else
            plot(DLPData.Vsweep{s}{c},DLPData.Ifit{s}{c}*1e3,'g')
        end
        ht = title(['I_s: ',num2str(DLPData.Isat_m{s}(c),2),' ,t = ',num2str(DLPData.time{s}(c),5),' s']);
        set(ht,'FontSize',6); grid on
        set(gcf,'color','w')
        set(gca,'FontSize',5)
        box on
        catch
            warning('error')
            break
        end
    end
    end
end