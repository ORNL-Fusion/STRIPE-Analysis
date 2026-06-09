
%%%%%%Set Up Coil Parameters
shots = [21080:1:21082,21084:1:21085,21087:1:21091];
r     = [-4:0.5:-2,-1:0.5:1];
Atten = [10*ones(3,1)',zeros(7,1)'];
Coil  = (1:1:12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%Set up Phase Detector Calibration File and Detection characteristics
filename = 'Phase Calibration.xlsx';
wiggle = 4.0;  %%wiggle room for phase detection in degrees (typically resolution is better)
UL     = 1.85; %%Upper limit of phase detectors
LL     = 0.15; %%Lower limit of phase detectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Loop through shot data here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:length(shots)
    for jj=1:length(Coil)
        address = ['\MPEX::TOP.MACHOPS1:',sprintf('COIL1_%0.0f',Coil(jj))];   
        [~,~]  = mdsopen('MPEX',shots(ii)); 
        [Data(ii,jj,:),~]  = mdsvalue(address);
        bla = ['DIM_OF(',address,')'];
        if jj==1 && ii==1
            [t,~]  = mdsvalue(bla);
            t = t(1:end-1);
        end
    end
    %%%Conversion to Voltage Ratio and Phase
    
    for kk=1:4
        tic
        Vmag(ii,kk,:) = 10^(Atten(ii)/20).*Vratio('VMAG.csv',Data(ii,kk*3,:));
        toc
        v1 = Data(ii,1+3*(kk-1),:);
        v2 = Data(ii,2+3*(kk-1),:);
        SS = kk;
        if kk==3
            SS = 1; %%%Calibration for S3 sucks
        end
        phmeas(ii,kk,:) = BdotPhaseDetect(v1,v2,filename,SS,wiggle,UL,LL);
        toc
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure
% plot(t,squeeze(phmeas(6,4,:)),'o')

%%%%Write the data to spreadsheets sheet 1 = amplitude ratio sheet 2 =
%%%%phase S1,S2,S3,S4 corresponds to the channel. S1=9.5-Bphi, S2=9.5Br,
%%%% S3 = broken, S4 = 10.5-Br
for ii=1:4
    title = sprintf('S%i.xls',ii);
    %%Write Amplitude Ratio
    xlswrite(title,squeeze(Vmag(:,1,:))',1,'B4')
    xlswrite(title,r,1,'B2')
    xlswrite(title,shots,1,'B1')
    xlswrite(title,{'radius [cm]'},1,'A2')
    xlswrite(title,{'shot number'},1,'A1')
    xlswrite(title,{'Amplitude ratio V1/V2'},1,'A3')
    %%Write Phase
    xlswrite(title,squeeze(phmeas(:,1,:))',2,'B4')
    xlswrite(title,r,2,'B2')
    xlswrite(title,shots,2,'B1')
    xlswrite(title,{'radius [cm]'},2,'A2')
    xlswrite(title,{'shot number'},2,'A1')
    xlswrite(title,{'Phase [degrees]'},2,'A3')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



