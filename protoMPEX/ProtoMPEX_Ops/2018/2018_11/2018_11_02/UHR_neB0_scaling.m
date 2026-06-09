% Calculating how the magnetic field and plasma density  are related in the
% UHR frequency

clear all
close all


ne  = linspace(0.01e19,1.5e20,1e3);
f = [28,53,105]*1e9;
k = 0*2*pi*f./c_light;
Te = 6; 
vthe = sqrt(2*e_c*Te/m_e);
wUH = (2*pi*f  - 3*(k.^2)*(vthe^2) );

c1 = ((m_e/e_c)*wUH).^2;
c2 = ne*(m_e/e_0);


figure;
hold on
for s = 1:length(wUH)
    B{s} = sqrt( c1(s) - c2 );
    h(s) = plot(ne,B{s},'LineWidth',2);
end
ylim([0,4])

set(gcf,'color','w')
box on
legend(h,'28 GHz','53 GHz','105 GHz')
ylabel('B_0 [T]')
xlabel('n_e [m^{-3}]')
grid on
% set(gca,'XTick',[0:1:15]*1e19)