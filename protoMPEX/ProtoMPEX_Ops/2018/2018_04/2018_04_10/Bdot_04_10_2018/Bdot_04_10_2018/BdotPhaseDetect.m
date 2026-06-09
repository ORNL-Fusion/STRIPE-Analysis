function [phmeas] = BdotPhaseDetect(v1,v2,calfile,SS,wiggle,UL,LL)
phmeas = zeros(1,length(v1));
sheet = sprintf('S%i',SS);
data(1,:,:) = xlsread(calfile,sheet);
phase(1,:)  = data(1,:,1);
Vph1(1,:)   = data(1,:,3);
Vph2(1,:)   = data(1,:,4);
 
%%%Get the quadrants%%%%%%%%%%%%%
XI = (sign(diff(Vph1(1,:)))==-1);
x1n = [XI,XI(end)];
XI = (sign(diff(Vph1(1,:)))==1);
x1p = [XI,XI(end)];

XI = (sign(diff(Vph2(1,:)))==-1);
x2n = [XI,XI(end)];
XI = (sign(diff(Vph2(1,:)))==1);
x2p = [XI,XI(end)];
%%%Get the quadrants%%%%%%%%%%%%%

%%%%%Interpolate values from each quadrant
ph1_p = interp_pp(Vph1(x1p),phase(x1p),v1);
ph1_n = interp_pp(Vph1(x1n),phase(x1n),v1);

ph2_p = interp_pp(Vph2(x2p),phase(x2p),v2);
ph2_n = interp_pp(Vph2(x2n),phase(x2n),v2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Figure out where the phase and assign a phase
for ii = 1:length(v1)
   if abs(ph1_p(ii)-ph2_p(ii))<wiggle
       phmeas(ii) = (ph1_p(ii)+ph2_p(ii))./2;
   end
   if abs(ph1_p(ii)-ph2_n(ii))<wiggle
       phmeas(ii) = (ph1_p(ii)+ph2_n(ii))./2;
   end
   if abs(ph1_n(ii)-ph2_p(ii))<wiggle
       phmeas(ii) = (ph1_n(ii)+ph2_p(ii))./2;
   end
   if abs(ph1_n(ii)-ph2_n(ii))<wiggle
       phmeas(ii) = (ph1_n(ii)+ph2_n(ii))./2;
   end
   if v1(ii)>UL
       phmeas(ii) = ph2_p(ii);
   end
   if v1(ii)<LL
       phmeas(ii) = ph2_n(ii);
   end
   if v2(ii)>UL
       phmeas(ii) = ph1_n(ii);
   end
   if v2(ii)<LL
       phmeas(ii) = ph1_p(ii);
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

