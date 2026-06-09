 load("ICRH_57877_8s.mat");
close all;

t_min=5;
t_max=10;
idx=find( (data.WDP.S57877.reflec.t>t_min).*(data.WDP.S57877.reflec.t<t_max) );

semilogy( data.WDP.S57877.reflec.position.r(:,idx) , data.WDP.S57877.reflec.ne(:,idx) , 'Color','#888888')
xlabel('R [m]')
ylabel('n_e [m^-3]')

nData=data.WDP.S57877.reflec.ne(:,idx);
rData=data.WDP.S57877.reflec.position.r(:,idx);
Rmin=min(rData(:));
Rmax=max(rData(:));

Rspace=linspace(Rmin,Rmax,100);
nData_interpolated=zeros([length(Rspace),length(idx)]);
for i=1:length(idx)
    nData_interpolated(:,i)=interp1(rData(:,i),nData(:,i),Rspace,'linear',nan); %extrapolation becomes NaN
end

nData_avg = mean(nData_interpolated,2,"omitnan"); %NaN not included in average
hold on;
semilogy(Rspace,nData_avg,'Color','#008800','LineWidth',3)
hold off;