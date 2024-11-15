% Plot single RFs from SWUS project

networks= {'AE'};
stations= {'U15A'};
event= 'EQ_2014_192_19_22';

for intwk=1:length(networks);
    netwk= networks{intwk};
    for istn= 1:length(stations)
        stn= stations{istn};
	filenameSp= ['/Users/evagolos/Research/ReceiverFunctions/RFcodes_Junlin/Data/Projects/SWUS/',netwk,'/',stn,'/','PreCalculated_RFs_Sp_4_100s_05.10.20_1_2.mat'];
    load(filenameSp);
    RF_Sp= PreCalculated_RFs.Sp.Individual_RFs.(event);
    T_Sp= PreCalculated_RFs.Sp.Time;
    filenamePs= ['/Users/evagolos/Research/ReceiverFunctions/RFcodes_Junlin/Data/Projects/SWUS/',netwk,'/',stn,'/','PreCalculated_RFs_Ps_4_100s_05.10.20_1.mat'];
    load(filenamePs);
    T_Ps= PreCalculated_RFs.Ps.Time;
    RF_Ps= PreCalculated_RFs.Ps.Individual_RFs.(event);
    neg_Ps= RF_Ps; ineg= RF_Ps>-0.001;
    neg_Ps(ineg)= -0.001;
    pos_Ps= RF_Ps; ipos= RF_Ps<0.001;
    pos_Ps(ipos)= 0.001;
    end
end


figure; hold on;
for i=1:1
    %area(T_Ps,RF_Ps',0.001);
    %area(T_Ps,neg_Ps',-0.001);
    %fill([zeros(length(T_Ps),1),T_Ps'],[RF_Ps,T_Ps'],[0 0.4470 0.7410]);
    fill([0.001; pos_Ps],[0; T_Ps'],[0 0.4470 0.7410]);
    fill([-0.001; neg_Ps],[0; T_Ps'],[0.8500 0.3250 0.0980]);
    plot([0.001 0.001],[min(T_Ps') max(T_Ps)],'w','linewidth',3);
    plot([-0.001 -0.001],[min(T_Ps') max(T_Ps)],'w','linewidth',3);
    %area(neg_Ps,-0.001);
    plot(RF_Ps,T_Ps','Linewidth',2,'Color','k');
end
set(gca,'Ydir','reverse');
ylim([-2 32]);
set(gca,'FontSize',14);
xlabel('RF'); ylabel('Time (s)');