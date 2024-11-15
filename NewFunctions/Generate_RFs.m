function Generate_RFs(Project, basedir, NETs, Phases, tag, varargin)

% ********* Function Description *********
%
% Take prepped waveform data (i.e., picked
% and windowed Ps and Sp waveforms), bin
% it, calculate the receiver function in
% time, and then migrate to depth.
%
%
% ****************************************
% *                                      *
% *  Written by David L. Abt - May 2008	 *
% *            (edited E. Hopper, 2014)  *
% *  Email: David_Abt@brown.edu          *
% *                                      *
% ****************************************

clc
Proj_Dir        = [basedir 'Data/Projects/' Project '/'];
addpath([basedir '/Functions/']);

load([Proj_Dir 'Networks.mat']);
% Choose networks (if not already hardwired in)
if isempty(NETs)
    nets=fieldnames(Networks);
    %disp('Available networks: '); disp(nets);
    %net=input('Type in network to generate RFs for, or hit ''return'' for all:  ','s');
    %if isempty(net); NETs=nets; else NETs=cellstr(net); end
    NETs=nets;%flipud(nets);%
end

% Choose stations (if multiple networks, all stations in
% those networks will be prepped)
if(length(Project)>6 && strcmp(Project(1:6),'Synth_'))
    if strcmp(tag,'SingleStaRF'); sta=cellstr('all');
    elseif length(varargin)==6; sta=varargin{end-1};
    else sta=cellstr('all');
    end
    version=varargin{end};
else
    if(length(NETs)==1 && strcmp(NETs,'all') || length(NETs)>1);
        sta=cellstr('all');
        NETs = fieldnames(Networks);
    else stas=fieldnames(Networks.(NETs{1})); clc
        disp('Available stations:  '); disp(stas);
        
        sta=input(['Type in station to generate RFs for,' ...
            ' or hit ''return'' for all:  '],'s');
        if isempty(sta); sta=cellstr('all'); else sta=cellstr(sta); end
    end
    if ~isempty(varargin) && strcmp(varargin{end}(1:5),'names')
        version=varargin{end};
        if ~exist([Proj_Dir version '.mat'],'file');
            %system(['cp ' Proj_Dir 'names.mat ' Proj_Dir version '.mat']);
            FileNames = [];
            save([Proj_Dir version '.mat'],'FileNames');
        end
    else version='names';
    end
end


switch tag
    case 'checkStatus'
        clc
        lowT=varargin{1};
        highT = varargin{2};
        Check_Prepped(Proj_Dir, NETs, sta, Networks, 'disp',lowT,highT);
        Check_RFs(Proj_Dir, NETs, sta, Networks, Phases, 'disp',lowT,highT);
        
    case 'generateHK'
        for iph=1:length(Phases)
            Phase=Phases{iph};
            lowT=varargin{1};
            highT = varargin{2};
            SvSh=varargin{3};
            
            % Check what has already been prepped/calculated
            [unprepped]=Check_Prepped(Proj_Dir, NETs, sta, Networks, 'n',lowT,highT, version);
            %             [matchrfs] = Check_RFs(Proj_Dir, NETs, sta, Networks, Phases, 'matchrfs', ...
            %                 lowT, CrustModel, MantleVpModel, MantleVsModel);
            matchrfs=[];
            clc
            
            
            
            if ~isempty(unprepped)
                fprintf('The following stations have not been prepped:\n');
                npn=fieldnames(unprepped); n=0;
                for k=1:length(npn); nps=fieldnames(unprepped.(npn{k}));
                    for l=1:length(nps);
                        n=n+1; fprintf('  %s.%s  ',npn{k},nps{l});
                        if rem(n,10)==0; fprintf('\n'); end
                    end
                end
            end
            
            if ~isempty(matchrfs)
                fprintf('\n\nThe following stations have already had RFs calculated with these parameters:\n');
                nrf=fieldnames(matchrfs); n=0;
                for k=1:length(nrf); srf=fieldnames(matchrfs.(nrf{k}));
                    for l=1:length(srf);
                        n=n+1; fprintf('  %s.%s  ',nrf{k},srf{l});
                        if rem(n,10)==0; fprintf('\n'); end
                    end
                end
            end
            
            if(~isempty(unprepped) || ~isempty(matchrfs))
                fprintf('\n\nUnprepped stations will be skipped, and RFs will be overwritten.\n');
                Continue='y';%[Continue]=input('Continue anyway (y/n)?  ','s');    clc
                if strcmp(Continue,'n'); return; end
            end
            
            [Params]=Set_RF_params(lowT, highT, Phase,'HK');
            Params.unprepped=unprepped; Params.NETs=NETs; Params.sta=sta; Params.Phase=Phase;
            Params.SvSh=SvSh;
            
            Calculate_RFs(Params, Networks,basedir, Proj_Dir, Project, 'HK', version);
            
        end
        
    case 'generate'
        for iph=1:length(Phases)
            Phase=Phases{iph};
            lowT=varargin{1};
            highT = varargin{2};
            SvSh=varargin{3};
            CrustModel=varargin{4};
            MantleVpModel=varargin{5}; MantleVsModel=varargin{6};
            isRPdiff = varargin{7};
            
            % Check what has already been prepped/calculated
            [unprepped]=Check_Prepped(Proj_Dir, NETs, sta, Networks, 'n',lowT,highT, version);
            %             [matchrfs] = Check_RFs(Proj_Dir, NETs, sta, Networks, Phases, 'matchrfs', ...
            %                 lowT, CrustModel, MantleVpModel, MantleVsModel);
            matchrfs=[];
            clc
            
            
            
            if ~isempty(unprepped)
                fprintf('The following stations have not been prepped:\n');
                npn=fieldnames(unprepped); n=0;
                for k=1:length(npn); nps=fieldnames(unprepped.(npn{k}));
                    for l=1:length(nps);
                        n=n+1; fprintf('  %s.%s  ',npn{k},nps{l});
                        if rem(n,10)==0; fprintf('\n'); end
                    end
                end
            end
            
            if ~isempty(matchrfs)
                fprintf('\n\nThe following stations have already had RFs calculated with these parameters:\n');
                nrf=fieldnames(matchrfs); n=0;
                for k=1:length(nrf); srf=fieldnames(matchrfs.(nrf{k}));
                    for l=1:length(srf)
                        n=n+1; fprintf('  %s.%s  ',nrf{k},srf{l});
                        if rem(n,10)==0; fprintf('\n'); end
                    end
                end
            end
            
            if(~isempty(unprepped) || ~isempty(matchrfs))
                fprintf('\n\nUnprepped stations will be skipped, and RFs will be overwritten.\n');
                Continue='y';%[Continue]=input('Continue anyway (y/n)?  ','s');    clc
                if strcmp(Continue,'n'); return; end
            end
            
            
            
            [Params]=Set_RF_params(lowT, highT, Phase,'RF',CrustModel,MantleVpModel,MantleVsModel,isRPdiff,basedir);
            Params.unprepped=unprepped; Params.NETs=NETs; Params.sta=sta; Params.Phase=Phase;
            Params.SvSh=SvSh;
            
            
            Calculate_RFs(Params, Networks,basedir, Proj_Dir, Project, 'RF', version);
            
            
            iffiltdat = 'no'; % Change to "yes" to use mantle cutoff criteria.
            if strcmp(iffiltdat,'yes')
                if ~exist([Proj_Dir 'RFrange_450.mat'],'file')
                    [RFrange.RF_min, RFrange.RF_max,RFrange.RF_med, RFrange.depth] = CompRefRF( Phase,Params,Networks,Proj_Dir,basedir,MantleVpModel,MantleVsModel);
                    save([Proj_Dir 'RFrange_450.mat'],'RFrange');
                else
                    load([Proj_Dir 'RFrange_450.mat'])
                end
                
                Dist_tot = [];
                Dist_tot_2 = [];
                Dist_tot_med = [];
                Self_tot_size = [];
                Self_shallow_tot_size = [];
                Self_deep_tot_size = [];
                Self_zero_tot_size = [];
                Dist_shallow_tot = [];
                Dist_deep_tot = [];
                Dist_zero_tot = [];
                Dist_shallow_tot_med = [];
                Dist_deep_tot_med = [];
                Self_shallow_Moho_size = [];
                Self_shallow_pv_size = [];
                Num_deep_out = [];
                Num_deep_out_2 = [];
                Max_out_range_2 = [];
                Max_deep_out_range_2 = [];
                
                ind_test = find(RFrange.depth<=450);
                norm_factor = sqrt(RFrange.RF_min(ind_test')'*RFrange.RF_min(ind_test'))*sqrt(RFrange.RF_max(ind_test')'*RFrange.RF_max(ind_test'));
                ind_shallow = find(RFrange.depth<=60 & RFrange.depth>=15);
                ind_deep = find(RFrange.depth>60 & RFrange.depth<=450);
                ind_zero = find(RFrange.depth>=0 & RFrange.depth<15);
                norm_factor_shallow = sqrt(RFrange.RF_min(ind_shallow')'*RFrange.RF_min(ind_shallow'))*sqrt(RFrange.RF_max(ind_shallow')'*RFrange.RF_max(ind_shallow'));
                norm_factor_deep = sqrt(RFrange.RF_min(ind_deep')'*RFrange.RF_min(ind_deep'))*sqrt(RFrange.RF_max(ind_deep')'*RFrange.RF_max(ind_deep'));
                norm_factor_zero = sqrt(RFrange.RF_min(ind_zero')'*RFrange.RF_min(ind_zero'))*sqrt(RFrange.RF_max(ind_zero')'*RFrange.RF_max(ind_zero'));
                for in = 1:length(NETs)
                    STAs=fieldnames(Networks.(NETs{in}));
                    for is=1:length(STAs)
                        if isfield(Params.unprepped,NETs{in})
                            if isfield(Params.unprepped.(NETs{in}),STAs{is})
                                continue;
                            end
                        end
                        disp(['Working on comparison of ' NETs{in} '.' STAs{is} '...'])
                        junk=dir([Proj_Dir NETs{in} '/' STAs{is} '/PreCalculated_RFs_*']);
                        if isempty(junk)
                            continue;
                        end
                        lj=length(junk); jd=zeros(lj,1);
                        for ij=1:lj; jd(ij)=datenum(junk(ij).date); end
                        [~,mi]=max(jd); precalc_name=junk(mi).name;
                        load([Proj_Dir NETs{in} '/' STAs{is} '/' precalc_name]);
                        load([Proj_Dir NETs{in} '/' STAs{is} '/Ray_Path_Data_',num2str(1/Params.hf),'_',num2str(1/Params.lf),'.mat']);
                        load([Proj_Dir NETs{in} '/' STAs{is} '/Migration_Models_',num2str(1/Params.hf),'_',num2str(1/Params.lf),'.mat']);
                        Migration_Data=Migration_Models.(Phase)(end);
                        evs=fieldnames(PreCalculated_RFs.(Phase).Individual_RFs);
                        for k=length(evs):-1:1
                            if isempty(PreCalculated_RFs.(Phase).Individual_RFs.(evs{k}));
                                evs(k)=[];
                            end
                        end
                        
                        pathstr=[Phase(1) '_Path'];
                        Mask_Dp=zeros(1,length(evs));
                        time=PreCalculated_RFs.(Phase).Time;
                        dZs   	= Migration_Data.Reference_Depth;
                        RPs     = Migration_Data.Relative_RPs;
                        Delay   = Migration_Data.Relative_DTs';
                        [rps,dzs] = meshgrid(RPs,dZs);
                        dzs=double(dzs);
                        dZs=double(dZs);
                        RFs=zeros(length(evs),length(time)); RP=zeros(1,length(evs));
                        for k = 1:length(evs)
                            RFs(k,:)=PreCalculated_RFs.(Phase).Individual_RFs.(evs{k});
                            RP(k)=All_Ray_Path_Data.(evs{k}).(pathstr).ray_parameter;
                            mdp=dZs(find(Migration_Data.Critical_RPs<=RP(k),1));
                            if isempty(mdp); Mask_Dp(k)=dZs(end); else Mask_Dp(k)=mdp; end
                        end
                        rfs_atdepth=zeros(length(evs),length(RFrange.depth));
                        for k=1:length(evs)
                            StaDelay=interp2(rps,dzs,Delay',RP(k),dZs);
                            if lowT==0; filtrfs=RFs(k,:);
                            else filtrfs=bpfilt(RFs(k,:),diff(time(1:2)),0.01,1/lowT,'lp');
                            end
                            rfs_atdepth_tmp=interp1(time,filtrfs,StaDelay+Migration_Data.Reference_DT);
                            rfs_atdepth_tmp(dZs>Mask_Dp(k))=0;
                            rfs_atdepth(k,:) = interp1(dZs,rfs_atdepth_tmp,RFrange.depth);
                        end
                        rfs_atdepth(isnan(rfs_atdepth)) = 0;
                        rfcompmin = rfs_atdepth-repmat(RFrange.RF_min',[size(rfs_atdepth,1),1]);
                        rfcompmax = rfs_atdepth-repmat(RFrange.RF_max',[size(rfs_atdepth,1),1]);
                        rfcompmed = rfs_atdepth-repmat(RFrange.RF_med',[size(rfs_atdepth,1),1]);
                        ind = find(rfcompmin.*rfcompmax<0);
                        rfcompmin(ind) = 0;
                        rfcompmax(ind) = 0;
                        rfcompmin_dist = sum(rfcompmin(:,ind_test).^2,2)/norm_factor;
                        rfcompmax_dist = sum(rfcompmax(:,ind_test).^2,2)/norm_factor;
                        rfcompmed_dist = sum(rfcompmed(:,ind_test).^2,2)/sum(RFrange.RF_med.^2);
                        rfcompmin_shallow_dist = sum(rfcompmin(:,ind_shallow).^2,2)/norm_factor_shallow;
                        rfcompmax_shallow_dist = sum(rfcompmax(:,ind_shallow).^2,2)/norm_factor_shallow;
                        rfcompmed_shallow_dist = sum(rfcompmed(:,ind_shallow).^2,2)/sum(RFrange.RF_med(ind_shallow').^2);
                        rfcompmin_deep_dist = sum(rfcompmin(:,ind_deep).^2,2)/norm_factor_deep;
                        rfcompmax_deep_dist = sum(rfcompmax(:,ind_deep).^2,2)/norm_factor_deep;
                        rfcompmed_deep_dist = sum(rfcompmed(:,ind_deep).^2,2)/sum(RFrange.RF_med(ind_deep').^2);
                        rfcompmin_zero_dist = sum(rfcompmin(:,ind_zero).^2,2)/norm_factor_zero;
                        rfcompmax_zero_dist = sum(rfcompmax(:,ind_zero).^2,2)/norm_factor_zero;
                        rf_self_size = sum(rfs_atdepth(:,ind_test).^2,2);
                        rf_shallow_self_size = sum(rfs_atdepth(:,ind_shallow).^2,2);
                        rf_deep_self_size = sum(rfs_atdepth(:,ind_deep).^2,2);
                        rf_zero_self_size = sum(rfs_atdepth(:,ind_zero).^2,2);
                        rf_nv_atdepth = rfs_atdepth;
                        rf_pv_atdepth = rfs_atdepth;
                        rf_nv_atdepth(rfs_atdepth>0) = 0;
                        rf_pv_atdepth(rfs_atdepth<0) = 0;
                        rf_shallow_self_moho_size = sum(rf_nv_atdepth(:,ind_shallow).^2,2);
                        rf_shallow_self_pv_size = sum(rf_pv_atdepth(:,ind_shallow).^2,2);
                        dist_sta = min([rfcompmin_dist,rfcompmax_dist],[],2);
                        dist_shallow_sta = min([rfcompmin_shallow_dist,rfcompmax_shallow_dist],[],2);
                        dist_deep_sta = min([rfcompmin_deep_dist,rfcompmax_deep_dist],[],2);
                        dist_zero_sta = min([rfcompmin_zero_dist,rfcompmax_zero_dist],[],2);
                        RF_mean = (RFrange.RF_min+RFrange.RF_max)/2;
                        RF_diff = (RFrange.RF_max-RFrange.RF_min)/2;
                        %                     rfcompmin_2 = rfs_atdepth-repmat((RF_mean-2*RF_diff)',[size(rfs_atdepth,1),1]);    %Fich Model
                        %                                         rfcompmax_2 = rfs_atdepth-repmat((RF_mean+2*RF_diff)',[size(rfs_atdepth,1),1]);  %Fich Model
                        rfcompmin_2 = rfs_atdepth-repmat((RF_mean-0.8*RF_diff)',[size(rfs_atdepth,1),1]);   %Nienke Model
                        rfcompmax_2 = rfs_atdepth-repmat((RF_mean+0.8*RF_diff)',[size(rfs_atdepth,1),1]);   %Nienke Model
                        ind2 = find(rfcompmin_2.*rfcompmax_2<0);
                        rfcompmin_2(ind2) = 0;
                        rfcompmax_2(ind2) = 0;
                        rfcompmin_2_d = rfcompmin_2;
                        rfcompmin_2_d(rfcompmin_2_d>0) = 0;
                        rfcompmax_2_d = rfcompmax_2;
                        rfcompmax_2_d(rfcompmax_2_d<0) = 0;
                        %                     rfcompmin_dist_2 = sum(rfcompmin_2(:,ind_test).^2,2)/norm_factor;
                        %                     rfcompmax_dist_2 = sum(rfcompmax_2(:,ind_test).^2,2)/norm_factor;
                        rfcompmin_dist_2 = sum(rfcompmin_2_d(:,ind_test).^2,2)/norm_factor;
                        rfcompmax_dist_2 = sum(rfcompmax_2_d(:,ind_test).^2,2)/norm_factor;
                        dist_sta_2 = rfcompmin_dist_2+rfcompmax_dist_2;
                        %dist_sta_2 = min([rfcompmin_dist_2,rfcompmax_dist_2],[],2);
                        max_out_range_2 = min([max(abs(rfcompmin_2(:,ind_test)),[],2),max(abs(rfcompmax_2(:,ind_test)),[],2)],[],2);
                        max_deep_out_range_2 = min([max(abs(rfcompmin_2(:,ind_deep)),[],2),max(abs(rfcompmax_2(:,ind_deep)),[],2)],[],2);
                        N_deep_out_range = sum(abs(sign(rfcompmin(:,ind_deep))),2);
                        N_deep_out_range_2 = sum(abs(sign(rfcompmin_2(:,ind_deep))),2);
                        CompRFres.(NETs{in}).(STAs{is}).eve = evs;
                        CompRFres.(NETs{in}).(STAs{is}).dist = dist_sta;
                        CompRFres.(NETs{in}).(STAs{is}).dist_2 = dist_sta_2;
                        CompRFres.(NETs{in}).(STAs{is}).dist_shallow = dist_shallow_sta;
                        CompRFres.(NETs{in}).(STAs{is}).dist_deep = dist_deep_sta;
                        CompRFres.(NETs{in}).(STAs{is}).dist_zero = dist_zero_sta;
                        CompRFres.(NETs{in}).(STAs{is}).dist_med = rfcompmed_dist;
                        CompRFres.(NETs{in}).(STAs{is}).dist_shallow_med = rfcompmed_shallow_dist;
                        CompRFres.(NETs{in}).(STAs{is}).dist_deep_med = rfcompmed_deep_dist;
                        CompRFres.(NETs{in}).(STAs{is}).RFsize = rf_self_size;
                        CompRFres.(NETs{in}).(STAs{is}).RFsize_shallow = rf_shallow_self_size;
                        CompRFres.(NETs{in}).(STAs{is}).RFsize_deep = rf_deep_self_size;
                        CompRFres.(NETs{in}).(STAs{is}).RFsize_zero = rf_zero_self_size;
                        CompRFres.(NETs{in}).(STAs{is}).RFsize_shallow_Moho = rf_shallow_self_moho_size;
                        CompRFres.(NETs{in}).(STAs{is}).RFsize_shallow_pv = rf_shallow_self_pv_size;
                        CompRFres.(NETs{in}).(STAs{is}).N_deep_out_range = N_deep_out_range;
                        CompRFres.(NETs{in}).(STAs{is}).N_deep_out_range_2 = N_deep_out_range_2;
                        CompRFres.(NETs{in}).(STAs{is}).Max_out_range_2 = max_out_range_2;
                        CompRFres.(NETs{in}).(STAs{is}).Max_deep_out_range_2 = max_deep_out_range_2;
                        
                        Dist_tot = [Dist_tot;dist_sta];
                        Dist_tot_2 = [Dist_tot_2;dist_sta_2];
                        Dist_shallow_tot = [Dist_shallow_tot;dist_shallow_sta];
                        Dist_deep_tot = [Dist_deep_tot;dist_deep_sta];
                        Dist_zero_tot = [Dist_zero_tot;dist_zero_sta];
                        Dist_tot_med = [Dist_tot_med;rfcompmed_dist];
                        Dist_shallow_tot_med = [Dist_shallow_tot_med;rfcompmed_shallow_dist];
                        Dist_deep_tot_med = [Dist_deep_tot_med;rfcompmed_deep_dist];
                        Self_tot_size = [Self_tot_size;rf_self_size];
                        Self_shallow_tot_size = [Self_shallow_tot_size;rf_shallow_self_size];
                        Self_deep_tot_size = [Self_deep_tot_size;rf_deep_self_size];
                        Self_zero_tot_size = [Self_zero_tot_size;rf_zero_self_size];
                        Self_shallow_Moho_size = [Self_shallow_Moho_size;rf_shallow_self_moho_size];
                        Self_shallow_pv_size = [Self_shallow_Moho_size;rf_shallow_self_pv_size];
                        Num_deep_out = [Num_deep_out;N_deep_out_range];   %number of depth out of the range
                        Num_deep_out_2 = [Num_deep_out_2;N_deep_out_range_2];
                        Max_out_range_2 = [Max_out_range_2;max_out_range_2];
                        Max_deep_out_range_2 = [Max_deep_out_range_2;max_deep_out_range_2];
                    end
                end
                
                Total_info.Dist_tot = Dist_tot;
                Total_info.Dist_tot_2 = Dist_tot_2;
                Total_info.Dist_shallow_tot = Dist_shallow_tot;
                Total_info.Dist_deep_tot = Dist_deep_tot;
                Total_info.Dist_zero_tot = Dist_zero_tot;
                Total_info.Dist_tot_med = Dist_tot_med;
                Total_info.Dist_shallow_tot_med = Dist_shallow_tot_med;
                Total_info.Dist_deep_tot_med = Dist_deep_tot_med;
                Total_info.Self_tot_size = Self_tot_size;
                Total_info.Self_shallow_tot_size = Self_shallow_tot_size;
                Total_info.Self_deep_tot_size = Self_deep_tot_size;
                Total_info.Self_zero_tot_size = Self_zero_tot_size;
                Total_info.Self_shallow_Moho_size = Self_shallow_Moho_size;
                Total_info.Self_shallow_pv_size = Self_shallow_pv_size;
                Total_info.Num_deep_out = Num_deep_out;
                Total_info.Num_deep_out_2 = Num_deep_out_2;
                Total_info.Max_out_range_2  = Max_out_range_2;
                Total_info.Max_deep_out_range_2  = Max_deep_out_range_2;
                
                save([Proj_Dir 'Similarity_RF_450.mat'],'CompRFres','Total_info')
                
            end
            
        end
        
    case 'HK'
        
        lowT=varargin{1};
        highT = varargin{2};
        tagHV=varargin{3};
        
        %HVstacking(Proj_Dir, Networks, lowT, highT, tagHV)
        HKstacking(Proj_Dir, Networks, lowT, highT, tagHV)
        
    case 'SingleStaRF'
        lowT=varargin{1}; highT = varargin{2};
        ifplot=varargin{3}; ifbaz=varargin{4};
        isRPdiff = varargin{5};
        
        iffiltdat = 'yes';
        if strcmp(iffiltdat,'yes')
            load([Proj_Dir 'RFrange_450.mat']);
            load([Proj_Dir 'Similarity_RF_450.mat'])
            dist_sort = sort((Total_info.Dist_tot_2));
            cutoff = dist_sort(round(0.9*length(Total_info.Self_shallow_Moho_size)));
            %cutoff = mean(Total_info.Dist_tot);
            moho_size_sort = sort((Total_info.Self_shallow_Moho_size));
            cutoff_shallow_size = moho_size_sort(round(0.15*length(Total_info.Self_shallow_Moho_size)));
            %cutoff_shallow_size = 0.01*median(Total_info.Self_shallow_Moho_size);    %ori
            pv_size_sort = sort((Total_info.Self_shallow_pv_size));
            cutoff_shallow_pv = pv_size_sort(round(0.85*length(Total_info.Self_shallow_Moho_size)));
            %cutoff_shallow_pv = 10*median(Total_info.Self_shallow_pv_size);   %ori
            %cutoff_deep_size = 0.8*median(Total_info.Self_deep_tot_size);
            %         cutoff_deep_dist = median(Total_info.Dist_deep_tot);
            deep_num_sort = sort((Total_info.Num_deep_out_2));
            cutoff_deep_num = deep_num_sort(round(0.5*length(Total_info.Num_deep_out_2)));
            zero_dist_sort = sort((Total_info.Dist_zero_tot));
            cutoff_zero_dist = zero_dist_sort(round(0.85*length(Total_info.Dist_zero_tot)));
            %cutoff_deep_num = 5*median(Total_info.Num_deep_out_2);   %ori
            %cutoff_deep_max = 1.2*median(Total_info.Max_deep_out_range_2);
            
            Phase=char(Phases{1}); % Need to change this if we do both phases in series.
            if strcmp(Phase,'Sp')
                cutoff = 5*median(Total_info.Dist_tot_2);
                cutoff_shallow_size = 0.2*median(Total_info.Self_shallow_Moho_size);
                cutoff_shallow_pv = 3*median(Total_info.Self_shallow_pv_size);
                cutoff_deep_num = 1.2*median(Total_info.Num_deep_out_2);
                %cutoff_shallow_size= 1.5;
                %cutoff_shallow_pv= 0.25;
            end
        end
        
        if ~exist([Proj_Dir '/SingleRF/'],'dir')
            mkdir([Proj_Dir '/SingleRF/']);
            
        end
        
        
        method = 'median';
        
        
        
        for iph=1:length(Phases)
            Phase=char(Phases{iph});
            
            
            for in=1:length(NETs)
                
                
                if strcmp(char(sta),'all')
                    STAs=fieldnames(Networks.(NETs{in}));
                else
                    STAs=sta;
                end
                
                Unis = [];
                marks_sta = zeros(size(STAs));
                num_c = 0;
                for is = 1:length(STAs)
                    if marks_sta(is) ~= 0
                        continue;
                    end
                    num_c = num_c+1;
                    Unis(num_c).name = STAs(is);
                    Unis(num_c).del = 0;
                    STA = STAs{is};
                    lstaname = length(STA);
                    marks_sta(is) = num_c;
                    for iss = 1:length(STAs)
                        if lstaname>=length(STAs{iss})
                            continue;
                        elseif strcmp(STA,STAs{iss}(1:lstaname))
                            Unis(num_c).name = [Unis(num_c).name,STAs(iss)];
                            if marks_sta(iss) == 0
                                marks_sta(iss) = num_c;
                            else
                                Unis(marks_sta(iss)).del = 1;
                                marks_sta(iss) = num_c;
                            end
                        end
                    end
                end
                
                for iu = 1:length(Unis)
                    clearvars rfs_sta
                    
                    if Unis(iu).del == 1
                        continue;
                    end
                    clearvars rfs_sta
                    
                    STAs = Unis(iu).name;
                    rfs_sta = [];
                    STAref = Unis(iu).name{1};
                    ind_dat = 0;
                    
                    disp(['Working on ' NETs{in} '.' STAref '...'])
                    
                    for is=1:length(STAs)
                        Mig_M = [];
                        rfs = [];
                        PreCalculated_RFs = [];
                        All_Waveform_Data = [];
                        All_Ray_Path_Data = [];
                        
                        
                        clearvars Mig_M PreCalculated_RFs All_Waveform_Data All_Ray_Path_Data rfs
                        clearvars -global
                        
                        junk=dir([Proj_Dir NETs{in} '/' STAs{is} '/PreCalculated_RFs_',...
                            Phase,'_',num2str(lowT),'_',num2str(highT),'*']);
                        if isempty(junk)
                            continue;
                        end
                        lj=length(junk); jd=zeros(lj,1);
                        for ij=1:lj; jd(ij)=datenum(junk(ij).date); end
                        [nd,mi]=max(jd);
                        
                        precalc_name=junk(mi).name;
                        load([Proj_Dir NETs{in} '/' STAs{is} '/' precalc_name]);
                        
                        load([Proj_Dir NETs{in} '/' STAs{is} '/Ray_Path_Data_',...
                            num2str(lowT),'_',num2str(highT),'.mat']);
                        load([Proj_Dir NETs{in} '/' STAs{is} '/Waveform_Data_',...
                            num2str(lowT),'_',num2str(highT),'.mat']);
                        load([Proj_Dir NETs{in} '/', STAs{is}, '/Migration_Models_',...
                            num2str(lowT),'_',num2str(highT),'.mat']);
                        
                        Migration_Data1= Migration_Models.(Phase)(end);
                        
                        if ~strcmp(iffiltdat,'yes')
                            eves = fieldnames(PreCalculated_RFs.(Phase).Individual_RFs);
                            evs = [];
                            for ie = 1:length(eves)
                                if isempty(PreCalculated_RFs.(Phase).Individual_RFs.(eves{ie}))
                                    continue;
                                end
                                if All_Waveform_Data.(eves{ie}).([Phase,'_Waveforms']).S2N_ratio<3 ||...
                                        abs(All_Waveform_Data.(eves{ie}).([Phase,'_Waveforms']).ZR_Xcorr_Coeff)<0
                                    continue;
                                end
                                
                                evs = [evs;eves(ie)];
                            end
                        else
                            evs =  CompRFres.(NETs{in}).(STAs{is}).eve;
                           evs(CompRFres.(NETs{in}).(STAs{is}).dist_2>cutoff | ...
                               CompRFres.(NETs{in}).(STAs{is}).RFsize_shallow_Moho<cutoff_shallow_size | ...
                               CompRFres.(NETs{in}).(STAs{is}).N_deep_out_range_2>cutoff_deep_num | ...
                               CompRFres.(NETs{in}).(STAs{is}).RFsize_shallow_pv>cutoff_shallow_pv) = [];
                            for ie = length(evs):-1:1
                                if isempty(PreCalculated_RFs.(Phase).Individual_RFs.(evs{ie}))
                                    evs(ie) = [];
                                % Combine to make elseif - EG 10/2020.
                                % For Sp in SWUS, use S2N_ratio<3.
                                elseif All_Waveform_Data.(evs{ie}).([Phase,'_Waveforms']).S2N_ratio<3 ||...
                                        abs(All_Waveform_Data.(evs{ie}).([Phase,'_Waveforms']).ZR_Xcorr_Coeff)<0
                                    evs(ie) = [];
                                end
                            end
                        end
                        
                        if isempty(evs)
                            continue;
                        end
                        
                        if strcmp(ifbaz,'yes')
                            bazs=zeros(length(evs),1); bazstr='_baz';
                            for k=1:length(evs);
                                bazs(k)=All_Ray_Path_Data.(evs{k}).S_Path.bAzimuth;
                            end
                            bazbins=-180:10:180; [~,mind]=max(smooth(histc(bazs,bazbins),3));
                            if mind==1
                                evs=evs(bazs<=bazbins(mind+2) | bazs>=bazbins(length(bazbins)-2));
                            elseif mind==2
                                evs=evs(bazs<=bazbins(mind+2) | bazs>=bazbins(length(bazbins)-1));
                            elseif mind==(length(bazbins)-1)
                                evs=evs(bazs<=bazbins(2) | bazs>=bazbins(mind-2));
                            elseif mind==(length(bazbins))
                                evs=evs(bazs<=bazbins(3) | bazs>=bazbins(mind-2));
                            else
                                evs=evs(bazs<=bazbins(mind+2) & bazs>=bazbins(mind-2));
                            end
                            bazbins(mind)
                            rps=zeros(length(evs),1);
                            for k=1:length(evs);
                                rps(k)=All_Ray_Path_Data.(evs{k}).S_Path.ray_parameter;
                            end
                            rpsbins = linspace(min(rps),max(rps),10);
                            [~,rind]=max(histc(rps,rpsbins));
                            evs=evs(rps<=rpsbins(min(rind+1,10)) & rps>=rpsbins(max(rind-1,1)));
                            rpmin = min(rps(rps<=rpsbins(min(rind+1,10)) & rps>=rpsbins(max(rind-1,1))));
                            rpmax = max(rps(rps<=rpsbins(min(rind+1,10)) & rps>=rpsbins(max(rind-1,1))));
                        else
                            rps=zeros(length(evs),1);
                            for k=1:length(evs)
                                rps(k)=All_Ray_Path_Data.(evs{k}).S_Path.ray_parameter;
                            end
                            rps_sort = sort(rps);
                            if length(rps_sort)<15
                                continue;
                            end
                            rpsbins = linspace(min(rps),max(rps),10);
                            [~,rind]=max(histc(rps,rpsbins));
                            rpsbins(min(rind+1,10)) = rps_sort(round(0.4*length(rps)));%0.110;
                            rpsbins(max(rind-1,1)) = rps_sort(1);%0.100;
                            %evs=evs(rps<=rpsbins(min(rind+1,10)) & rps>=rpsbins(max(rind-1,1)));
                            rpmin = min(rps(rps<=rpsbins(min(rind+1,10)) & rps>=rpsbins(max(rind-1,1))));
                            rpmax = max(rps(rps<=rpsbins(min(rind+1,10)) & rps>=rpsbins(max(rind-1,1))));
                            rpmin = rps_sort(1);%0.100;
                            rpmax = rps_sort(round(0.4*length(rps)));%0.110;
                            bazstr=''; bazbins=666; mind = 1;
                        end
                        
                        %                     if length(evs)<30
                        %                         continue;
                        %                     end
                        %
                        if iscell(Phase); ph=Phase{:}; else ph=Phase; end
                        
                        [depth,rfs] = Migrate_single_sta_RFs(PreCalculated_RFs, ...
                            Migration_Data1, All_Ray_Path_Data, ph, evs,lowT,isRPdiff);
                        
                        
                        % Set criteria for rejecting traces from the stack.
                        % I played around with this a lot, right now it's
                        % optimized for Sp data.
                        %                         ind_dep_um = find(depth >= 130 & depth <= 300);
                        %                         ind_dep_lm = find(depth >= 700 & depth <= 800);
                        ind_dep_m = find(depth >= 200);
                        ind_410 = find(depth <= 470 & depth >= 350);
                        ind_660 = find(depth <= 720 & depth >= 600);
                        %ind_dep_moho = find(depth >= Migration_Data1.Original_Models.Crustal_Model.Moho-10 & ...
                        %    depth <= Migration_Data1.Original_Models.Crustal_Model.Moho+10);
                        ind_dep_moho = find(depth >= 25 & depth <= 50);
                        ind_sfc= find(depth <= 5);
                        ind_dep_um = find(depth >= 60 & depth <= 150);
                        
                                                 rf_um = rfs(:,ind_dep_um);
                        %                         rf_lm = rfs(:,ind_dep_lm);
                        rf_m = rfs(:,ind_dep_m);
                        rf_moho = rfs(:,ind_dep_moho);
                        rf_410 = rfs(:,ind_410);
                        rf_660 = rfs(:,ind_660);
                        rf_sfc= rfs(:,ind_sfc);
                        
                        max_moho = nanmax(rf_moho,[],2);
                        max_410 = nanmax(rf_410,[],2);
                        max_660 = nanmax(rf_660,[],2);
                        max_sfc= nanmax(abs(rf_sfc),[],2);
                                         n_nnan_um = nansum(~isnan(rf_um),2);
                        %                 n_nnan_lm = nansum(~isnan(rf_lm),2);
                        n_nnan_m = nansum(~isnan(rf_m),2);
                        
                                         rms_um = sqrt(nansum(rf_um.^2,2)./n_nnan_um);
                        %                 rms_lm = sqrt(nansum(rf_lm.^2,2)./n_nnan_lm);
                        rms_m = sqrt(nansum(rf_m.^2,2)./n_nnan_m);
                        
                        %                 n_um = nansum(abs(rf_um)>0.05,2)./n_nnan_um;
                        %                 n_lm = nansum(abs(rf_lm)>0.05,2)./n_nnan_lm;
                        %n_m = nansum(abs(rf_m)>0.07,2)./n_nnan_m;
                        n_m = nansum(abs(rf_m)>0.1,2)./n_nnan_m;
                        %n_m2 = nansum(abs(rf_m)<0.003,2)./n_nnan_m;
                        
                        %                 ind_u = find(max_moho>0.02 & rms_um<0.1 & rms_lm<0.1 &...
                        %                     n_um<0.6 & (n_lm<0.6 | isnan(n_lm)));
                        %ind_u = find(max_moho>0.03 & rms_m<0.1 & n_m<0.4);
                        
                        %                         n_m = nansum(abs(rf_m)>0.07,2)./n_nnan_m;
                        %n_m2 = nansum(abs(rf_m)<0.003,2)./n_nnan_m;
                        
                        %                 ind_u = find(max_moho>0.02 & rms_um<0.1 & rms_lm<0.1 &...
                        %                     n_um<0.6 & (n_lm<0.6 | isnan(n_lm)));
                        %ind_u = find(max_moho>0.03 & rms_m<0.1 & n_m<0.4);
                        if strcmp(Phase,'Ps')
                            %ind_u = find(max_moho>0.02 & max_410>0.02 & max_660>0.02 & rms_m<0.15 & n_m<0.4);
                            if max(depth)>400
                                ind_u = find(max_moho>0.01 & max_410>0.015 & max_660>0.015 & rms_m<0.1 & n_m<0.3);
                            else
                                ind_u = find(max_moho>0.01 & rms_m<0.1 & n_m<0.3);
                            end
                        else
                         ind_u = find(max_moho>0.01 & rms_m<0.1 & n_m<0.3 & max_sfc<0.1); % Trial A
                              %ind_u = find(max_moho>0.06 & max_moho<0.3 & rms_m<0.1 & n_m<0.2 & max_sfc<0.1); % Trial F
                             %ind_u = find(max_moho>0.06 & rms_m<0.1 & rms_um>0.01 & n_m<0.3 & max_sfc<0.1); % Trial G
                             %ind_u = find(max_moho>0.06 & rms_m<0.1 & rms_um>0.02 & n_m<0.3 & max_sfc<0.1); % Trial H
                        end 
                        
                        rfs_sta = [rfs_sta;rfs(ind_u,:)];
                        
                        
                        ind_dat = ind_dat+1;
                    end
                    
                    if ind_dat == 0
                        continue;
                    end
                    
                    if isempty(rfs_sta); continue; end
                    
                    if strcmp(method,'median')
                        meanrf=nanmedian(rfs_sta,1);
                        
                        % Bootstrap
                        rng shuffle; nbt=100;
                        bt_rfs=zeros(nbt,length(meanrf));
                        for k=1:nbt
                            indx=randi(size(rfs_sta,1),size(rfs_sta,1),1);
                            bt_rfs(k,:)=nanmedian(rfs_sta(indx,:),1);
                            % this line was a test bt_rfs(k,:)=nanmean(rfs_sta(indx,:),1);
                        end
                        
                        % this line was a test bootstrap.mean=median(bt_rfs,1); bootstrap.std=std(bt_rfs);
                        bootstrap.mean=mean(bt_rfs,1); bootstrap.std=std(bt_rfs);
                        bootstrap.max=max(bt_rfs); bootstrap.min=min(bt_rfs);
                        meanrf=bootstrap.mean;
                    else
                        meanrf=nanmean(rfs_sta,1);
                        nrf = zeros(length(depth),1);
                        for i = 1:length(depth)
                            nrf(i) = length(find(~isnan(rfs_sta(:,i))));
                        end
                        bootstrap.mean = meanrf;
                        bootstrap.std = nanstd(rfs_sta,1)./sqrt(nrf');
                        
                    end
                    
                    RF_Depth=[depth meanrf'];
                    
                    Vptmp = Migration_Data1.Velocity_Model.Vp;
                    Vstmp = Migration_Data1.Velocity_Model.Vs;
                    Deptmp = Migration_Data1.Velocity_Model.Depth;
                    for i = 2:length(Deptmp)
                        if  Deptmp(i) == Deptmp(i-1)
                            Deptmp(i) = Deptmp(i)+0.01;
                        end
                    end
                    Vs_ray = interp1(Deptmp,Vstmp,depth);
                    Vp_ray = interp1(Deptmp,Vptmp,depth);
                    Moho_ray =  Migration_Data1.Original_Models.Crustal_Model.Moho;
                    dV_ray = gradient(Vs_ray);
                    
                    if strcmp(Phase,'Ps')
                        rp_rep = 0.06;
                    else
                        rp_rep = 0.105;
                    end
                    Period = 2;
                    recordid='RFcomp';
                    vs = Vs_ray(1);
                    vp = Vp_ray(1);
                    rho = 2.8;
                    zbot = [];
                    ztop = 0;
                    ind_na = 0;
                    for idep = 1:length(depth)
                        if Vp_ray(idep)*rp_rep>0.97
                            if zbot < depth(idep-1)
                                zbot = [zbot;depth(idep-1)];
                                ind_na = 1;
                            end
                            break;
                        end
                        if abs(Vs_ray(idep)-vs(end))/vs(end)>0.005
                            vs = [vs;Vs_ray(idep)];
                            vp = [vp;Vp_ray(idep)];
                            if depth(idep)<Moho_ray
                                rho = [rho;2.8];
                            else
                                rho = [rho;3.3];
                            end
                            zbot = [zbot;depth(idep)];
                            ztop = [ztop;depth(idep)];
                        end
                    end
                    if ind_na == 0
                        if ztop(end) ~= depth(end)
                            zbot = [zbot;depth(end)];
                        else
                            zbot = [zbot;depth(end)+1];
                        end
                    end
                    %                         write_propmat_syn(ztop,zbot,vp,vs,rho,rp_rep,Period,[recordid,'.mat'],Phase);
                    %                         SynWave = load([basedir,'NewFunctions/PROPMAT/',recordid,'.mat']);
                    %                         [SynRF,SynTime] = getrf(Phase,SynWave.R,SynWave.Z,SynWave.dt,vp(1),vs(1),rp_rep);
                    %
                    %                         [SynDepth,Synindts,Synindte] = getdepth(Phase,SynTime,zbot,vp,vs,rp_rep);
                    %
                    %                         Synrf_dep = interp1(SynDepth,SynRF(Synindts:Synindte),depth);
                    
                    if strcmp(ifplot,'yes') && size(rfs_sta,1)>15
                        
                        if strcmp(Phase,'Ps')
                            m_d = 1000;
                        else
                            m_d = 350;
                        end
                        
                        figure; %set(gcf,'position',[680 300 470 790])
                        subplot(1,3,1)
                        plot([0 0],[0 m_d],'k--'); hold on; set(gca,'ydir','reverse')
                        plot(RF_Depth(:,2),RF_Depth(:,1),'linewidth',2);
                        %plot(Synrf_dep,RF_Depth(:,1),'k:','linewidth',2);
                        mprf=(RF_Depth(:,2)-bootstrap.std'); trf=RF_Depth(:,1);
                        trf2=0:0.1:m_d; mprf=interp1(trf,mprf,trf2);
                        plot(mprf,trf2,'k--'); inds=find(mprf>=0);
                        if ~isempty(inds)
                            dinds=[find(diff(inds)>1) length(inds)]; sind=1;
                            for k=1:length(dinds); pinds=inds(sind):inds(dinds(k));
                                if mprf(pinds(1))>1e-4; MPRF=[0 mprf(pinds)];
                                    TRF2=[0 trf2(pinds)];
                                else MPRF=mprf(pinds); TRF2=trf2(pinds);
                                end
                                if mprf(pinds(end))>1e-4; MPRF=[MPRF 0]; TRF2=[TRF2 TRF2(end)]; end
                                h=fill(MPRF, TRF2,'r','linestyle','none');
                                set(h,'facealpha',0.5);sind=1+dinds(k);
                            end
                        end
                        pprf=(RF_Depth(:,2)+bootstrap.std'); pprf=interp1(trf,pprf,trf2);
                        plot(pprf,trf2,'k--');inds=find(pprf<=0);
                        if ~isempty(inds)
                            dinds=[find(diff(inds)>1) length(inds)]; sind=1;
                            for k=1:length(dinds); pinds=inds(sind):inds(dinds(k));
                                if pprf(pinds(1))<-1e-4; PPRF=[0 pprf(pinds)];
                                    TRF2=[0 trf2(pinds)];
                                else PPRF=pprf(pinds); TRF2=trf2(pinds);
                                end
                                if pprf(pinds(end))<-1e-4; PPRF=[PPRF 0]; TRF2=[TRF2 TRF2(end)]; end
                                h=fill(PPRF, TRF2,'b','linestyle','none');
                                set(h,'facealpha',0.5); sind=1+dinds(k);
                            end
                        end
                        xlim([-0.1,0.2]);
                        title([STAs{is} ': ' num2str(round(Networks.(NETs{in}).(STAs{is}).Latitude*10)/10) ...
                            '\circN, ' num2str(round(Networks.(NETs{in}).(STAs{is}).Longitude*10)/10) ...
                            '\circE (' num2str(size(rfs_sta,1)) ' waveforms)'])
                        %                         title([STAs{is} ': ' num2str(size(rfs,1)) ' waveforms'])
                        ylabel('Depth (km)'); xlabel('RF Amplitude'); ylim([0 m_d])
                        %ylim([0 100])  %Junlin
                        
                        
                        subplot(1,3,2)
                        plot([0 0],[0 m_d],'k--'); hold on; set(gca,'ydir','reverse')
                        plot(dV_ray,RF_Depth(:,1),'linewidth',2);
                        xlim([-0.05,0.1]);
                        
                        ylabel('Depth (km)'); xlabel('Diff Velocity'); ylim([0 m_d])
                        %ylim([0 100])  %Junlin
                        
                        subplot(1,3,3)
                        plot([0 0],[0 m_d],'k--'); hold on; set(gca,'ydir','reverse')
                        hold on
                        int = round(size(rfs_sta,1)/50);
                        int = max(int,1);
                        for ii = 1:int:size(rfs_sta,1)
                            plot(rfs_sta(ii,:),RF_Depth(:,1),'linewidth',0.2,'Color',[0.7,0.7,0.7]);
                        end
                        plot(RF_Depth(:,2),RF_Depth(:,1),'r','linewidth',2);
                        %plot(Synrf_dep,RF_Depth(:,1),'b','linewidth',2);
                        %title(['Lat = ',num2str(slat),'   Lon = ',num2str(slon)]);
                        title(['rp ',num2str(rpmin),' ',num2str(rpmax),' rprep ',num2str(rp_rep)]);
                        %                         saxislim = max(abs(min(min(-Synrf_dep),min(RF_Depth(:,2)))),max(max(-Synrf_dep),max(RF_Depth(:,2))));
                        %                         xlim([-2*saxislim,2*saxislim]);
                        xlim([-0.3,0.3])
                        ylim([0 m_d])
                        %ylim([0 100])  %Junlin
                        print([Proj_Dir '/SingleRF/' NETs{in} '_' STAref '_' Phase '_' num2str(bazbins(mind)) '_',num2str(lowT),'_',num2str(highT),'.png'],'-dpng','-r300');
                        
                        %print([Proj_Dir '/SingleRF/' NETs{in} '_' STAs{is} '_' Phase '_' num2str(bazbins(mind)) '.jpg'],'-djpeg','-r300');
                        %pause;
                        % close all
                        
                        close(gcf)
                    end
                    nrfs = size(rfs_sta,1);
                    save([Proj_Dir '/SingleRF/' NETs{in} '_' STAref '_' Phase '_' num2str(bazbins(mind)) '_',num2str(lowT),'_',...
                        num2str(highT),'_recent.mat'],'depth','bootstrap','nrfs');
                    
                    
                end
            end
            
        end
        
        
end



end


function [unprepped] = Check_Prepped(Proj_Dir, NETs, sta, Networks, tag, lowT,highT,version)

unprepped=[];
junk=strsplit(Proj_Dir,'Projects/');
if length(junk{2})>6 && strcmp(junk{2}(1:6),'Synth_')
    load([Proj_Dir version '.mat']); synthstr=FileNames.Synth;
else synthstr='';
end
for in=1:length(NETs)
    NET=NETs{in};
    if strcmp(sta,'all'); stas=fieldnames(Networks.(NET)); else stas=sta; end
    
    for is=1:length(stas)
        STA=stas{is};
        if ~exist([Proj_Dir '/' NET '/' STA '/Waveform_Data_' num2str(lowT) '_' num2str(highT) synthstr '.mat'],'file')
            if strcmp(tag, 'disp');
                disp(['Station ' NET '.' STA ' has not been prepped!'])
            end
            unprepped.(NET).(STA)=[];
        end
    end
end

end

function [matchrfs] = Check_RFs(Proj_Dir, NETs, sta, Networks, Phases, tag, varargin)

matchrfs=[];

for in=1:length(NETs)
    NET=NETs{in};
    if strcmp(sta,'all'); stas=fieldnames(Networks.(NET)); else stas=sta; end
    
    for is=1:length(stas)
        STA=stas{is};
        fID=fopen([Proj_Dir NET '/' STA '/RFstatus.txt'],'w');
        fprintf(fID,'*** RF Status (%s) ***\n', datestr(date,'dd.mm.yy'));
        fprintf(fID,'         %s.%s\n', NET,STA);
        if strcmp(tag,'disp'); fprintf('\n%s.%s', NET,STA); end
        
        for iph=1:length(Phases)
            Phase=Phases{iph};
            clear RFs_in_*
            if exist([Proj_Dir '/' NET '/' STA '/' Phase '_in_Depth.mat'],'file')
                load([Proj_Dir '/' NET '/' STA '/' Phase '_in_Depth.mat'])
                fprintf(fID,'\n\n%s RF\n', Phase);
                if strcmp(tag,'disp'); fprintf('\n\n%s RF\n', Phase); end
                
                for iRF=1:length(RFs_in_Depth)
                    hf=RFs_in_Depth(iRF).RF_Parameters.hf;
                    lf=RFs_in_Depth(iRF).RF_Parameters.lf;
                    Crust=RFs_in_Depth(iRF).Migration_Parameters.Original_Models.Crust;
                    MantleVp=RFs_in_Depth(iRF).Migration_Parameters.Original_Models.Mantle_Vp;
                    MantleVs=RFs_in_Depth(iRF).Migration_Parameters.Original_Models.Mantle_Vs;
                    
                    fprintf(fID,['  %g-%gs: %s (crust); %s (mantle Vp);' ...
                        ' %s (mantle Vs)\n'], round(1/hf), round(1/lf), ...
                        Crust, MantleVp, MantleVs);
                    switch tag
                        case 'disp'
                            fprintf(['  %g-%gs: %s (crust); ' ...
                                '%s (mantle Vp); %s (mantle Vs)\n'], ...
                                round(1/hf), round(1/lf), Crust, MantleVp, MantleVs);
                            matchrfs=[];
                        case 'matchrfs'
                            lowT=varargin{1}; CrustModel=varargin{2};
                            MantleVpModel=varargin{3}; MantleVsModel=varargin{4};
                            
                            if(lowT==hf && strcmp(CrustModel,Crust) && ...
                                    strcmp(MantleVp,MantleVpModel) && ...
                                    strcmp(MantleVs,MantleVsModel))
                                matchrfs.(NET).(STA)=Phase;
                            end
                    end
                end
            end
        end
        
        D=dir([Proj_Dir '/' NET '/'  STA '/PreCalculated_RFs_*']);
        if(~isempty(D) && strcmp(tag,'disp')); fprintf('\n ETMTM:\n'); end
        for id=1:length(D)
            out=strsplit(D(id).name,'_'); Phase=char(out(3));
            sT=str2double(out{4}); bT=str2double(out{5}(1:end-1));
            Crust=char(out(6)); MantleVp=char(out(end-2));
            MantleVs=char(out(end-1));
            
            fprintf(fID,['  %g-%gs (%s): %s (crust); %s (mantle Vp);' ...
                ' %s (mantle Vs)\n'], sT, bT, Phase, Crust, MantleVp, MantleVs);
            switch tag
                case 'disp'
                    fprintf(['  %g-%gs (%s): %s (crust); ' ...
                        '%s (mantle Vp); %s (mantle Vs)\n'], ...
                        sT, bT, Phase, Crust, MantleVp, MantleVs);
                    matchrfs.(NET).(STA)=[];
                case 'matchrfs'
                    lowT=varargin{1}; CrustModel=varargin{2};
                    MantleVpModel=varargin{3}; MantleVsModel=varargin{4};
                    
                    if(lowT==sT && strcmp(CrustModel,Crust) && ...
                            strcmp(MantleVp,MantleVpModel) && ...
                            strcmp(MantleVs,MantleVsModel))
                        matchrfs.(NET).(STA)=Phase;
                    end
            end
        end
        
        
        
        
        
        if(strcmp(tag,'disp') && (~isfield(matchrfs,NET) || ~isfield(matchrfs.(NET),STA)))
            fprintf('     - no RFs calculated\n'); fprintf(fID,'     - no RFs calculated\n');
        end
        fclose(fID);
        
    end
end
end


function [Params] = Set_RF_params(lowT, highT, Phase,tag,varargin)

% ********* Function Description *********
%
% Set all parameters here and check them
%
% ****************************************


% Set culling parameters
switch Phase
    case 'Sp' % Currently set to minimum possible values
        S2N=2; % Signal to noise ratio must be >= 0
        ZRx=0; % ZR xcorr coefficient must be between 0 - 1
        TaS2N=15; % TauP S2N misfit must be >=0
    case 'Ps'
        S2N=2; ZRx=0; TaS2N=15;
end

if S2N<0; S2N=0; end; Params.S2N=S2N;
if ZRx<0; ZRx=0; elseif ZRx>1; ZRx=1; end; Params.ZR_Xcorr=ZRx;
if TaS2N<0; TaS2N=0; end; Params.TauP_S2N=TaS2N;


% Set deconvolution parameters
switch Phase
    case 'Sp'
        hf=1/lowT; lf=1/highT; % High freq corner must be > low freq corner Junlin
        Params.PreFilter={'Yes'}; % Prefilter waveforms (Yes/No)?
        Params.Taper={'Both'}; % Taper parent phase (Beginning/End/Both)?
        Params.Weight_By={'SV'}; % Weight by S2N of (P/SV)?
        Params.Max_Depth=6371;% 300;
    case 'Ps'
        hf=1/lowT; lf=1/highT;
        Params.PreFilter={'Yes'};
        Params.Taper={'Beginning'};
        Params.Weight_By={'P'};
        Params.Max_Depth=6371;
end


Params.Bin='Distance'; % Bin waveforms by (Distance/Back Azimuth/Both)

if lf>hf; lf=hf-0.01; end; Params.lf=lf; Params.hf=hf;

if lf>hf; lf=hf-0.01; end; Params.lf=lf; Params.hf=hf;


if ~strcmp(tag,'HK')
    CrustModel = varargin{1};
    MantleVpModel = varargin{2};
    MantleVsModel = varargin{3};
    isRPdiff = varargin{4};
    basedir = varargin{5};
    
    Params.isRPdiff = isRPdiff;
    % Check names of migration models
    % N.B. This code only set up to migrate with predefined models
    if exist([basedir 'Data/Velocity_Models/Crust/' CrustModel '.mat'],'file') || strcmp(CrustModel,'HK')
        Params.Crust=CrustModel;
    elseif (sum(strcmp(CrustModel,{'None','none'})))
        Params.Crust='None';
    else fprintf(['\nThe specified crustal velocity model does not exist...\n' ...
            'Choose one of these instead:  \n']);
        D=dir([basedir 'Data/Velocity_Models/Crust/*.mat']);
        for j=1:length(D); fprintf('%s\n',D(j).name(1:end-4)); end
        Params.Crust=input('\nCrustal Model:  ','s');
    end
    
    if exist([basedir 'Data/Velocity_Models/Mantle/Vp/' MantleVpModel '.mat'],'file')
        Params.MantleVp=MantleVpModel;
    else fprintf(['\nThe specified mantle Vp velocity model does not exist...\n' ...
            'Choose one of these instead:  \n']);
        D=dir([basedir 'Data/Velocity_Models/Mantle/Vp/*.mat']);
        for j=1:length(D); fprintf('%s\n',D(j).name(1:end-4)); end
        Params.MantleVp=input('\nMantle Vp Model:  ','s');
    end
    
    if exist([basedir 'Data/Velocity_Models/Mantle/Vs/' MantleVsModel '.mat'],'file')
        Params.MantleVs=MantleVsModel;
    else fprintf(['\nThe specified mantle Vs velocity model does not exist...\n' ...
            'Choose one of these instead:  \n']);
        D=dir([basedir 'Data/Velocity_Models/Mantle/Vs/*.mat']);
        for j=1:length(D); fprintf('%s\n',D(j).name(1:end-4)); end
        Params.MantleVs=input('\nMantle Vs Model:  ','s');
    end
end

end



function [savename]=Calculate_RFs(Params, Networks,basedir, Proj_Dir, Project,tagdeconv,version)

% **
% *
if exist([Proj_Dir version '.mat'],'file'); load([Proj_Dir version '.mat']); end

if(length(Project)>6 && strcmp(Project(1:6),'Synth_'))
    tag='synth';
    savename=[FileNames.Synth '_' Params.Phase '_' num2str(round(1/Params.hf)) '_' ...
        num2str(round(1/Params.lf)) 's_' datestr(date,'dd.mm.yy') '.mat'];
    if strcmp(Project(end-3:end),'Test'); tag2='test'; end
else
    if strcmp(tagdeconv,'HK')
        savename=['_HK_' Params.Phase '_' num2str(round(1/Params.hf)) '_' ...
            num2str(round(1/Params.lf)) 's_' datestr(date,'dd.mm.yy') '.mat'];
    else
        savename=['_' Params.Phase '_' num2str(1/Params.hf) '_' ...
            num2str(1/Params.lf) 's_' datestr(date,'dd.mm.yy') '.mat'];
    end
    tag='data';
end


if strcmp(tag,'data'); FileNames.Synth=''; end

if ~strcmp(tagdeconv,'HK')
    HKfile = ['HKstacking_',num2str(round(1/Params.hf)),'_',num2str(round(1/Params.lf)),'.mat'];
end

NETs=Params.NETs; Phase = Params.Phase;
if isfield(Params,'SvSh'); SvSh=Params.SvSh; else SvSh='SV'; end

% Run through selected networks
for in=1:length(NETs)
    NET     = NETs{in};
    if strcmpi(Params.sta,'all')
        STAs=fieldnames(Networks.(NET));
        if isfield(Params.unprepped,NET)
            STAs(isfield(Params.unprepped.(NET),STAs))=[];
        end
    else
        STAs = Params.sta;
    end
    
    marks_sta = zeros(size(STAs));
    Unis = [];
    num_c = 0;
    for is = 1:length(STAs)
        if marks_sta(is) ~= 0
            continue;
        end
        num_c = num_c+1;
        Unis(num_c).name = STAs(is);
        Unis(num_c).del = 0;
        STA = STAs{is};
        lstaname = length(STA);
        marks_sta(is) = num_c;
        for iss = 1:length(STAs)
            if lstaname>=length(STAs{iss})
                continue;
            elseif strcmp(STA,STAs{iss}(1:lstaname))
                Unis(num_c).name = [Unis(num_c).name,STAs(iss)];
                if marks_sta(iss) == 0
                    marks_sta(iss) = num_c;
                else
                    Unis(marks_sta(iss)).del = 1;
                    marks_sta(iss) = num_c;
                end
            end
        end
    end
    
    
    % if strcmp(NET,'TA'); STAs=STAs([42 155 288]); end
    
    % Run through selected stations
    for is=1:length(STAs)
        STA = STAs{is}; Data_Dir = [Proj_Dir NET '/' STA '/'];
        
        
        
        if(in==1 && is==1); n=1;
            while exist([Data_Dir 'PreCalculated_RFs' savename],'file')
                savename=[savename(1:end-4) '_' num2str(n) savename(end-3:end)];
                n=n+1;
            end
            
            % Save filename for PreCalcRFs
            FileNames.(Params.Phase).PreCalcRFs=savename; save([Proj_Dir version '.mat'], 'FileNames');
        end
        
        clc;
        
        if ~strcmp(tagdeconv,'HK')
            
            clear RFs_in_*;
            if strcmp(tag,'synth')
                if(exist('tag2','var') && strcmp(tag2,'test'))
                    MG_File = [Proj_Dir 'ak135_Migration_Models.mat'];
                else
                    MG_File = [basedir 'Data/Projects/' Project(7:end) '/' NET ...
                        '/' STA '/Migration_Models.mat'];
                end
            else
                MG_File = [Data_Dir 'Migration_Models_',num2str(1/Params.hf),'_',num2str(1/Params.lf),'.mat'];
                tag='data';
            end
            %                     if exist(MG_File,'file') && strcmp(Phase,'Ps')
            %                         system(['rm ',Data_Dir 'Migration_Models_',num2str(round(1/Params.hf)),'_',num2str(round(1/Params.lf)),'.mat']);
            %                     end
            
            if(exist(MG_File,'file'))
                [Calc_MG,i_MG,Migration_Models] = Check_MG_Models(Params,MG_File,Phase);
                if strcmp(Calc_MG,'No')
                    Migration_Data = Migration_Models.(Phase)(i_MG);
                end
            else
                Calc_MG	= 'Yes'; i_MG = 1;
            end
            
            if strcmp(Calc_MG,'Yes')
                % Extract the migration velocity profile
                %load([Data_Dir 'Station_Data_',num2str(round(1/Params.hf)),'_',num2str(round(1/Params.lf)),'.mat'])
                load([Data_Dir 'Station_Data_',num2str(1/Params.hf),'_',num2str(1/Params.lf),'.mat'])
                SLAT = Station_Data.Latitude; SLON = Station_Data.Longitude;
                fID=fopen([Proj_Dir 'MigrationProblems.txt'],'a');
                if strcmp(Params.Crust,'HK')
                    [Migration_Data] = Extract_Migration_Model(Params,SLAT,...
                        SLON,Phase,NET,STA, basedir,fID,Proj_Dir,HKfile);
                else
                    [Migration_Data] = Extract_Migration_Model(Params,SLAT,...
                        SLON,Phase,NET,STA, basedir,fID);
                end
                fclose(fID);
                Migration_Models.(Phase)(i_MG,1) = Migration_Data;
                %save([Data_Dir 'Migration_Models_',num2str(round(1/Params.hf)),'_',num2str(round(1/Params.lf)),'.mat'],'Migration_Models')
                save([Data_Dir 'Migration_Models_',num2str(1/Params.hf),'_',num2str(1/Params.lf),'.mat'],'Migration_Models')
                for iu = 1:length(Unis(marks_sta(is)).name)
                    if ~strcmp(Unis(marks_sta(is)).name{iu},STA)
                        disp([STA,'_',Unis(marks_sta(is)).name{iu}]);
                        Data_Dir_n = [Proj_Dir NET '/' Unis(marks_sta(is)).name{iu} '/'];
                        %system(['cp ',Data_Dir 'Migration_Models_',num2str(round(1/Params.hf)),'_',num2str(round(1/Params.lf)),'.mat ',Data_Dir_n]);
                        system(['cp ',Data_Dir 'Migration_Models_',num2str(1/Params.hf),'_',num2str(1/Params.lf),'.mat ',Data_Dir_n]);
                    end
                end
            end
            
            
            
            
        end
        
        % Load all waveform data and cull according to the
        % specified parameters from the GUI
        %load([Data_Dir 'Waveform_Data_' num2str(round(1/Params.hf)) '_' num2str(round(1/Params.lf)) FileNames.Synth '.mat'])
        %load([Data_Dir 'Ray_Path_Data_' num2str(round(1/Params.hf)) '_' num2str(round(1/Params.lf)) '.mat'])
        load([Data_Dir 'Waveform_Data_' num2str(1/Params.hf) '_' num2str(1/Params.lf) FileNames.Synth '.mat'])
        load([Data_Dir 'Ray_Path_Data_' num2str(1/Params.hf) '_' num2str(1/Params.lf) '.mat'])
        
        RP_RNG = [0 1];%[min(Migration_Data.Relative_RPs) max(Migration_Data.Relative_RPs)];
        [Event_IDs,dt,Culled] = Cull_Data(All_Waveform_Data,All_Ray_Path_Data,...
            RP_RNG,Params,NET,STA,Phase);
        
        
        if ~isempty(Event_IDs)
            
            % Set binning parameters
            [Bin_Data] = Set_Bin_Parameters(Params.Bin,Phase);
            
            % Bin and normalize data
            if strcmp(tagdeconv,'HK')
                [Waveform_Data,~,Bin_Data,Culled] = Bin_and_Normalize(...
                    Event_IDs,All_Waveform_Data,All_Ray_Path_Data,Bin_Data,...
                    Params,NET,STA,Phase,tag,SvSh,Culled,tagdeconv);
            else
                [Waveform_Data,~,Bin_Data,Culled] = Bin_and_Normalize(...
                    Event_IDs,All_Waveform_Data,All_Ray_Path_Data,Bin_Data,...
                    Params,NET,STA,Phase,tag,SvSh,Culled,tagdeconv,Migration_Data);
            end
            save([Data_Dir 'CulledRFs' savename], 'Culled')
            
            if ~isempty(Waveform_Data)
                
                
                % Deconvolve data to get single receiver functions (ETMTM)
                Deconvolve_Data(Waveform_Data,Bin_Data,Params,...
                    dt,NET,STA,Phase,Data_Dir,1/Params.hf,1/Params.lf, tagdeconv,savename);
                
            else
                disp(['No data for ' NET '.' STA '!!!']); pause(2);
            end
            
        else
            disp(['No data for ' NET '.' STA '!!!']); pause(2);
        end
        
    end
    
end



end


function HKstacking(Proj_Dir, Networks, lowT, highT, tag)


VpVs = 1.5:0.01:2.2;
Cr = 5:0.5:50;

[VPVS,CR] = meshgrid(VpVs,Cr);

if strcmp(tag,'both') || strcmp(tag,'one')
    F = zeros(size(VPVS));
    fPStot = F;
    fPPPStot = F;
    fPPSStot = F;
    fSPtot = F;
    fSSSPtot = F;
    fSSPPtot = F;
    fPStot2 = F;
    fPPPStot2 = F;
    fPPSStot2 = F;
    fSPtot2 = F;
    fSSSPtot2 = F;
    fSSPPtot2 = F;
    Vpsum = 0;
    Vssum = 0;
    Mohosum = 0;
    Latsum = 0;
    Lonsum = 0;
    Elesum = 0;
    nSp_tot = 0;
    nPs_tot = 0;
end

MinS2N=2.5;
MinZRx=0.6;
MaxTaupS2N=15;

% w1 = 0.25;
% w2 = 0.125;
% w3 = 0.125;
% w4 = 0.3;
% w5 = 0.05;
% w6 = 0.15;

w1 = 0.275;
w2 = 0.125;
w3 = 0.275;
% w2 = 0.2;
% w3 = 0.2;
w4 = 0.175;
w5 = 0.05;
w6 = 0.1;

% w1 = 0.3;
% w2 = 0.2;
% w3 = 0.2;
% w4 = 1;
% w5 = 0;
% w6 = 0.7;

NETs = fieldnames(Networks);
for in = 1:length(NETs)
    NET = NETs{in};
    STAs=fieldnames(Networks.(NET));
    marks_sta = zeros(size(STAs));
    num_c = 0;
    for is = 1:length(STAs)
        if marks_sta(is) ~= 0
            continue;
        end
        num_c = num_c+1;
        Unis(num_c).name = STAs(is);
        Unis(num_c).del = 0;
        STA = STAs{is};
        lstaname = length(STA);
        marks_sta(is) = num_c;
        for iss = 1:length(STAs)
            if lstaname>=length(STAs{iss})
                continue;
            elseif strcmp(STA,STAs{iss}(1:lstaname))
                Unis(num_c).name = [Unis(num_c).name,STAs(iss)];
                if marks_sta(iss) == 0
                    marks_sta(iss) = num_c;
                else
                    Unis(marks_sta(iss)).del = 1;
                end
            end
        end
    end
    for iu = length(Unis):-1:1
        if Unis(iu).del == 1
            Unis(iu) = [];
        end
    end
    
    for iu = 1:length(Unis)
        STAs = Unis(iu).name;
        STAref = STAs{1};
        disp(['Working on ' NET '.' STAref '...'])
        
        crust = getcrust(Networks.(NET).(STAref).Latitude,Networks.(NET).(STAref).Longitude,'CRUST1.0');
        if crust.thk(1)~=0
            thk = crust.thk(6:8);
            vp = crust.vp(6:8);
            vs = crust.vs(6:8);
        else
            thk = crust.thk(3:8);
            vp = crust.vp(3:8);
            vs = crust.vs(3:8);
        end
        Vp = sum(thk)/nansum(thk./vp);
        Vs = sum(thk)/nansum(thk./vs);
        Moho = sum(crust.thk)+Networks.(NET).(STAref).Elevation;
        
        VP = ones(size(VPVS))*Vp;
        VS = VP./VPVS;
        
        fPS = zeros(size(VP));
        fPPPS = fPS;
        fPPSS = fPS;
        fSP = fPS;
        fSSSP = fPS;
        fSSPP = fPS;
        fPS2 = fPS;
        fPPPS2 = fPS;
        fPPSS2 = fPS;
        fSP2 = fPS;
        fSSSP2 = fPS;
        fSSPP2 = fPS;
        nPs = 0;
        nSp = 0;
        
        for is = 1:length(STAs)
            STA = STAs{is};
            junk1=dir([Proj_Dir NET '/' STA '/PreCalculated_RFs_HK_Ps*']);
            junk2=dir([Proj_Dir NET '/' STA '/PreCalculated_RFs_HK_Sp*']);
            if isempty(junk1) && isempty(junk2)
                continue;
            end
            
            if ~isempty(junk1)
                lj=length(junk1); jd=zeros(lj,1);
                for ij=1:lj; jd(ij)=datenum(junk1(ij).date); end
                [~,mi]=max(jd); precalc_name_Ps=junk1(mi).name;
                PsRF = load([Proj_Dir NET '/' STA '/' precalc_name_Ps]);
                if isfield(PsRF.PreCalculated_RFs.Ps,'Time')
                    if ~isempty(PsRF.PreCalculated_RFs.Ps.Time)
                        timePs = PsRF.PreCalculated_RFs.Ps.Time;
                    end
                end
            end
            
            if ~isempty(junk2)
                lj=length(junk2); jd=zeros(lj,1);
                for ij=1:lj; jd(ij)=datenum(junk2(ij).date); end
                [~,mi]=max(jd); precalc_name_Sp=junk2(mi).name;
                SpRF = load([Proj_Dir NET '/' STA '/' precalc_name_Sp]);
                if isfield(SpRF.PreCalculated_RFs.Sp,'Time')
                    if ~isempty(SpRF.PreCalculated_RFs.Sp.Time)
                        timeSp = SpRF.PreCalculated_RFs.Sp.Time;
                    end
                end
            end
            
            load([Proj_Dir NET '/' STA '/Ray_Path_Data_',num2str(lowT),'_',num2str(highT),'.mat']);
            load([Proj_Dir NET '/' STA '/Waveform_Data_',num2str(lowT),'_',num2str(highT),'.mat']);
            
            EQs = fieldnames(All_Waveform_Data);
            
            for ie = 1:length(EQs)
                EQ = EQs{ie};
                if ~isempty(junk1)
                    if ~isempty(PsRF.PreCalculated_RFs.Ps.Individual_RFs.(EQ))
                        if All_Waveform_Data.(EQ).Ps_Waveforms.ZR_Xcorr_Coeff >= MinZRx &&...
                                All_Waveform_Data.(EQ).Ps_Waveforms.S2N_ratio >= MinS2N &&...
                                abs(All_Waveform_Data.(EQ).Ps_Waveforms.Window.Taup_S2N_Misfit)<=MaxTaupS2N
                            
                            p = All_Ray_Path_Data.(EQ).P_Path.ray_parameter;
                            rf = double(PsRF.PreCalculated_RFs.Ps.Individual_RFs.(EQ));
                            
                            tPS = CR.*(sqrt(VS.^-2-p^2)-sqrt(VP.^-2-p^2));
                            tPPPS = CR.*(sqrt(VS.^-2-p^2)+sqrt(VP.^-2-p^2));
                            tPPSS = 2*CR.*sqrt(VS.^-2-p^2);
                            
                            fPS = fPS+interp1(timePs,rf,tPS);
                            fPPPS = fPPPS+interp1(timePs,rf,tPPPS);
                            fPPSS = fPPSS+interp1(timePs,rf,tPPSS);
                            
                            fPS2 = fPS2+(interp1(timePs,rf,tPS)).^2;
                            fPPPS2 = fPPPS2+(interp1(timePs,rf,tPPPS)).^2;
                            fPPSS2 = fPPSS2+(interp1(timePs,rf,tPPSS)).^2;
                            
                            nPs = nPs+1;
                            
                        end
                    end
                end
                
                if ~isempty(junk2)
                    if ~isempty(SpRF.PreCalculated_RFs.Sp.Individual_RFs.(EQ))
                        if All_Waveform_Data.(EQ).Sp_Waveforms.ZR_Xcorr_Coeff >= MinZRx &&...
                                All_Waveform_Data.(EQ).Sp_Waveforms.S2N_ratio >= MinS2N &&...
                                abs(All_Waveform_Data.(EQ).Sp_Waveforms.Window.Taup_S2N_Misfit)<=MaxTaupS2N
                            
                            p = All_Ray_Path_Data.(EQ).S_Path.ray_parameter;
                            rf = double(SpRF.PreCalculated_RFs.Sp.Individual_RFs.(EQ));
                            
                            tSP = -CR.*(sqrt(VS.^-2-p^2)-sqrt(VP.^-2-p^2));
                            tSSSP = CR.*(sqrt(VS.^-2-p^2)+sqrt(VP.^-2-p^2));
                            tSSPP = 2*CR.*sqrt(VP.^-2-p^2);
                            
                            fSP = fSP+interp1(timeSp,rf,tSP);
                            fSSSP = fSSSP+interp1(timeSp,rf,tSSSP);
                            fSSPP = fSSPP+interp1(timeSp,rf,tSSPP);
                            
                            fSP2 = fSP2+(interp1(timeSp,rf,tSP)).^2;
                            fSSSP2 = fSSSP2+(interp1(timeSp,rf,tSSSP)).^2;
                            fSSPP2 = fSSPP2+(interp1(timeSp,rf,tSSPP)).^2;
                            
                            nSp = nSp+1;
                            
                        end
                    end
                end
                
            end
        end
        
        if strcmp(tag,'both') || strcmp(tag,'station')
            %             fPSmax = nanmax(nanmax(nanmax(fPS)));
            %             fPPPSmax = nanmax(nanmax(nanmax(fPPPS)));
            %             fPPSSmax = nanmax(nanmax(nanmax(-fPPSS)));
            %             fSPmax = nanmax(nanmax(nanmax(-fSP)));
            %             fSSPPmax = nanmax(nanmax(nanmax(fSSPP)));
            %             fSSSPmax = nanmax(nanmax(nanmax(-fSSSP)));
            %HK.(NET).(STA).Fsta = w1*fPS/fPSmax + w2*fPPPS/fPPPSmax-w3*fPPSS/fPPSSmax-w4*fSP/fSPmax+w5*fSSPP/fSSPPmax-w6*fSSSP/fSSSPmax;
            %HK.(NET).(STA).Fsta = w1*fPS/nPs + w2*fPPPS/nPs-w3*fPPSS/nPs-w4*fSP/nSp+w5*fSSPP/nSp-w6*fSSSP/nSp;
            HK.(NET).(STAref).Fsta = w1*fPS + w2*fPPPS-w3*fPPSS-w4*fSP+w5*fSSPP-w6*fSSSP;
            if nPs~=0 && nSp~=0
                HK.(NET).(STAref).Fsta_std = sqrt(w1*(fPS2-fPS.^2/nPs)+w2*(fPPPS2-fPPPS.^2/nPs)+...
                    w3*(fPPSS2-fPPSS.^2/nPs)+w4*(fSP2-fSP.^2/nSp)+w5*(fSSPP2-fSSPP.^2/nSp)+...
                    w6*(fSSSP2-fSSSP.^2/nSp));
            elseif nPs==0
                HK.(NET).(STAref).Fsta_std = sqrt(w4*(fSP2-fSP.^2/nSp)+w5*(fSSPP2-fSSPP.^2/nSp)+...
                    w6*(fSSSP2-fSSSP.^2/nSp));
            elseif nSp==0
                HK.(NET).(STAref).Fsta_std = sqrt(w1*(fPS2-fPS.^2/nPs)+w2*(fPPPS2-fPPPS.^2/nPs)+...
                    w3*(fPPSS2-fPPSS.^2/nPs));
            else
                HK.(NET).(STAref).Fsta_std = nan*ones(size(fPS));
            end
            HK.(NET).(STAref).nSp = nSp;
            HK.(NET).(STAref).nPs = nPs;
            if nSp+nPs < 30
                HK.(NET).(STAref).Vs = Vs;
                HK.(NET).(STAref).Vp = Vp;
                HK.(NET).(STAref).Cr = Moho;
                HK.(NET).(STAref).tag = 'All_CRUST1.0';
            else
                %E = -log(HK.(NET).(STA).Fsta);
                maxF = max(max(HK.(NET).(STAref).Fsta));
                %Eexp = -log(maxF);
                [ind_Cr,ind_VpVs] = find(HK.(NET).(STAref).Fsta==maxF);
                if Cr(ind_Cr)<=8
                    HK.(NET).(STAref).Vs = Vs;
                    HK.(NET).(STAref).Vp = Vp;
                    HK.(NET).(STAref).Cr = Moho;
                    HK.(NET).(STAref).tag = 'All_CRUST1.0';
                else
                    stdF = HK.(NET).(STAref).Fsta_std(ind_Cr,ind_VpVs);
                    ind_po = find(HK.(NET).(STAref).Fsta>=maxF-stdF);
                    VpVs_Range = max(VPVS(ind_po))-min(VPVS(ind_po));
                    Cr_Range = max(CR(ind_po))-min(CR(ind_po));
                    
                    if VpVs_Range>0.3 || Cr_Range>8
                        HK.(NET).(STAref).Vp = Vp;
                        HK.(NET).(STAref).Vs = Vs;
                        Fc = interp2(VPVS,CR,HK.(NET).(STAref).Fsta,Vp/Vs*ones(size(Cr)),Cr);
                        stdFc = interp2(VPVS,CR,HK.(NET).(STAref).Fsta_std,Vp/Vs*ones(size(Cr)),Cr);
                        [maxFc,ind_c] = max(Fc);
                        ind_pc = find(Fc>=maxFc-stdFc(ind_c));
                        Crc_Range = max(Cr(ind_pc))-min(Cr(ind_pc));
                        if Crc_Range>8
                            HK.(NET).(STAref).Cr = Moho;
                            HK.(NET).(STAref).tag = 'All_CRUST1.0';
                        else
                            HK.(NET).(STAref).Cr = Cr(ind_c);
                            HK.(NET).(STAref).CrRange = [min(Cr(ind_pc)),max(Cr(ind_pc))];
                            HK.(NET).(STAref).tag = 'Depth_HK';
                        end
                    else
                        HK.(NET).(STAref).Vp = Vp;
                        HK.(NET).(STAref).Vs = Vp/VpVs(ind_VpVs);
                        HK.(NET).(STAref).Cr = Cr(ind_Cr);
                        HK.(NET).(STAref).CrRange = [min(CR(ind_po)),max(CR(ind_po))];
                        HK.(NET).(STAref).VpVsRange = [min(VPVS(ind_po)),max(VPVS(ind_po))];
                        HK.(NET).(STAref).tag = 'All_HK';
                    end
                end
            end
        end
        
        for is = 1:length(STAs)
            STA = STAs{is};
            HK.(NET).(STA) = HK.(NET).(STAref);
        end
        
        if strcmp(tag,'both') || strcmp(tag,'one')
            Vpsum = Vpsum + (nPs+nSp)*Vp;
            Vssum = Vssum + (nPs+nSp)*Vs;
            Mohosum = Mohosum + (nPs+nSp)*Moho;
            Latsum = Latsum+(nPs+nSp)*Networks.(NET).(STAref).Latitude;
            Lonsum = Lonsum+(nPs+nSp)*Networks.(NET).(STAref).Longitude;
            Elesum = Elesum+(nPs+nSp)*Networks.(NET).(STAref).Elevation;
            F = F+w1*fPS + w2*fPPPS-w3*fPPSS-w4*fSP+w5*fSSPP-w6*fSSSP;
            fPStot = fPStot+fPS;
            fPPPStot = fPPPStot+fPPPS;
            fPPSStot = fPPSStot+fPPSS;
            fSPtot = fSPtot+fSP;
            fSSSPtot = fSSSPtot+fSSSP;
            fSSPPtot = fSSPPtot+fSSPP;
            fPStot2 = fPStot2+fPS2;
            fPPPStot2 = fPPPStot2+fPPPS2;
            fPPSStot2 = fPPSStot2+fPPSS2;
            fSPtot2 = fSPtot2+fSP2;
            fSSSPtot2 = fSSSPtot2+fSSSP2;
            fSSPPtot2 = fSSPPtot2+fSSPP2;
            nSp_tot = nSp_tot+nSp;
            nPs_tot = nPs_tot+nPs;
            
        end
        
        
    end
end

if strcmp(tag,'both') || strcmp(tag,'one')
    %     fPSmax = nanmax(nanmax(nanmax(fPStot)));
    %     fPPPSmax = nanmax(nanmax(nanmax(fPPPStot)));
    %     fPPSSmax = nanmax(nanmax(nanmax(-fPPSStot)));
    %     fSPmax = nanmax(nanmax(nanmax(-fSPtot)));
    %     fSSPPmax = nanmax(nanmax(nanmax(fSSPPtot)));
    %     fSSSPmax = nanmax(nanmax(nanmax(-fSSSPtot)));
    %F = w1*fPStot/fPSmax + w2*fPPPStot/fPPPSmax-w3*fPPSStot/fPPSSmax-w4*fSPtot/fSPmax+w5*fSSPPtot/fSSPPmax-w6*fSSSPtot/fSSSPmax;
    %F = w1*fPStot/nPs_tot + w2*fPPPStot/nPs_tot-w3*fPPSStot/nPs_tot-w4*fSPtot/nSp_tot+w5*fSSPPtot/nSp_tot-w6*fSSSPtot/nSp_tot;
    HK.F = F;
    if nPs_tot~=0 && nSp_tot~=0
        HK.F_std = sqrt(w1*(fPStot2-fPStot.^2/nPs_tot)+w2*(fPPPStot2-fPPPStot.^2/nPs_tot)+...
            w3*(fPPSStot2-fPPSStot.^2/nPs_tot)+w4*(fSPtot2-fSPtot.^2/nSp_tot)+w5*(fSSPPtot2-fSSPPtot.^2/nSp_tot)+...
            w6*(fSSSPtot2-fSSSPtot.^2/nSp_tot));
    elseif nPs_tot==0
        HK.F_std = sqrt(w4*(fSPtot2-fSPtot.^2/nSp_tot)+w5*(fSSPPtot2-fSSPPtot.^2/nSp_tot)+...
            w6*(fSSSPtot2-fSSSPtot.^2/nSp_tot));
    elseif nSp_tot==0
        HK.F_std = sqrt(w1*(fPStot2-fPStot.^2/nPs_tot)+w2*(fPPPStot2-fPPPStot.^2/nPs_tot)+...
            w3*(fPPSStot2-fPPSStot.^2/nPs_tot));
    else
        HK.F_std = nan*ones(size(fPStot));
    end
    HK.nSp_tot = nSp_tot;
    HK.nPs_tot = nPs_tot;
    HK.Latitude = Latsum/(nSp_tot+nPs_tot);
    HK.Longitude = Lonsum/(nSp_tot+nPs_tot);
    HK.Elevation = Elesum/(nSp_tot+nPs_tot);
    if nSp_tot+nPs_tot < 30
        HK.Vs = Vssum/(nSp_tot+nPs_tot);
        HK.Vp = Vpsum/(nSp_tot+nPs_tot);
        HK.Cr = Mohosum/(nSp_tot+nPs_tot);
        HK.tag = 'All_CRUST1.0';
    else
        maxF = max(max(F));
        [ind_Cr,ind_VpVs] = find(F==maxF);
        if Cr(ind_Cr)<=8
            HK.Vs = Vssum/(nSp_tot+nPs_tot);
            HK.Vp = Vpsum/(nSp_tot+nPs_tot);
            HK.Cr = Mohosum/(nSp_tot+nPs_tot);
            HK.tag = 'All_CRUST1.0';
        else
            stdF = HK.F_std(ind_Cr,ind_VpVs);
            ind_po = find(F>=maxF-stdF);
            VpVs_Range = max(VPVS(ind_po))-min(VPVS(ind_po));
            Cr_Range = max(CR(ind_po))-min(CR(ind_po));
            
            if VpVs_Range>0.3 || Cr_Range>8
                HK.Vp = Vpsum/(nSp_tot+nPs_tot);
                HK.Vs = Vssum/(nSp_tot+nPs_tot);
                Fc = interp2(VPVS,CR,F,HK.Vp/HK.Vs*ones(size(Cr)),Cr);
                stdFc = interp2(VPVS,CR,HK.F_std,HK.Vp/HK.Vs*ones(size(Cr)),Cr);
                [maxFc,ind_c] = max(Fc);
                ind_pc = find(Fc>=maxFc-stdFc(ind_c));
                Crc_Range = max(Cr(ind_pc))-min(Cr(ind_pc));
                if Crc_Range>8
                    HK.Cr = Mohosum/(nSp_tot+nPs_tot);
                    HK.tag = 'All_CRUST1.0';
                else
                    HK.Cr = Cr(ind_c);
                    HK.CrRange = [min(Cr(ind_pc)),max(Cr(ind_pc))];
                    HK.tag = 'Depth_HK';
                end
            else
                HK.Vp = Vpsum/(nSp_tot+nPs_tot);
                HK.Vs = HK.Vp/VpVs(ind_VpVs);
                HK.Cr = Cr(ind_Cr);
                HK.CrRange = [min(CR(ind_po)),max(CR(ind_po))];
                HK.VpVsRange = [min(VPVS(ind_po)),max(VPVS(ind_po))];
                HK.tag = 'All_HK';
            end
        end
        HK.VPVS = VPVS;
        HK.CR = CR;
    end
end

save([Proj_Dir 'HKstacking_',num2str(lowT),'_',num2str(highT),'.mat'],'HK');

end

function [Calc_MG,i_MG,Migration_Models] = Check_MG_Models(Params,MG_File,Phase)

tmp=load(MG_File); Migration_Models=tmp.Migration_Models;

if isfield(Migration_Models,Phase)
    Calc_MG = 'Yes';  i_MG = length(Migration_Models.(Phase))+1;
    
    % For previously calculated migration models, check if names match
    for imi=1:length(Migration_Models.(Phase))
        Crust_im    = 0; Vp_im       = 0; Vs_im       = 0;
        Orig_Model  = Migration_Models.(Phase)(imi).Original_Models;
        
        if strcmp(Params.Crust,Orig_Model.Crust); Crust_im=imi; end
        if strcmp(Params.MantleVp,Orig_Model.Mantle_Vp); Vp_im=imi; end
        if strcmp(Params.MantleVs,Orig_Model.Mantle_Vs); Vs_im=imi; end
        
        if(sum([Crust_im Vp_im Vs_im])~=0 && Crust_im==Vp_im && Vp_im==Vs_im)
            Calc_MG = 'No'; i_MG = Vp_im;
            break
        end
    end
    
else
    Calc_MG = 'Yes';  i_MG = 1;
    
end
%Calc_MG = 'Yes';
end

function [Migration_Data] = Extract_Migration_Model(Params,SLAT,SLON,Phase,NET,STA, basedir, fID,varargin)

% **
% *

disp(['Extracting migration model for ' STA ' (' NET ')...'])

% Load velocity models
if (sum(strcmp(Params.Crust,{'None','none'})))
    Velocity_Models.Crust = [];
elseif strcmp(Params.Crust,'HK')
    Proj_Dir = varargin{1};
    HKfile = varargin{2};
    tmp = load([Proj_Dir,HKfile]);
    HK = tmp.HK;
    clear tmp
    if strcmp(NET,'ALL')
        Velocity_Models.Crust.Moho = HK.Cr;
        Velocity_Models.Crust.Vp = HK.Vp;
        Velocity_Models.Crust.Vs = HK.Vs;
        Velocity_Models.Crust.Name = 'HK';
        Velocity_Models.Crust.Latitude = [];
    else
        indtmp = 0;
        if ~isfield(HK,NET)
            indtmp=1;
        else
            if ~isfield(HK.(NET),STA)
                indtmp=1;
            end
        end
        if indtmp == 1
            Velocity_Models.Crust.Moho = HK.Cr;
            Velocity_Models.Crust.Vp = HK.Vp;
            Velocity_Models.Crust.Vs = HK.Vs;
            Velocity_Models.Crust.Name = 'HK';
            Velocity_Models.Crust.Latitude = [];
        else
            Velocity_Models.Crust.Moho = HK.(NET).(STA).Cr;
            Velocity_Models.Crust.Vp = HK.(NET).(STA).Vp;
            Velocity_Models.Crust.Vs = HK.(NET).(STA).Vs;
            Velocity_Models.Crust.Name = 'HK';
            Velocity_Models.Crust.Latitude = [];
        end
    end
else
    Velocity_Models=load([basedir '/Data/Velocity_Models/Crust/' Params.Crust '.mat']);
end
junkvp=load([basedir '/Data/Velocity_Models/Mantle/Vp/' Params.MantleVp '.mat']);
junkvs=load([basedir '/Data/Velocity_Models/Mantle/Vs/' Params.MantleVs '.mat']);
Velocity_Models.Mantle.Vp=junkvp.Vp_Model;
Velocity_Models.Mantle.Vs=junkvs.Vs_Model;

Er = 6371;              % Earth's radius
% Keep longitudes between 0 and 360 for the mantle
if SLON<0; SLONm = SLON+360; else SLONm = SLON; end
isRPdiff = Params.isRPdiff;

% ***************************************************************
% * Gather Velocity Models and Extract a Profile at the Station *
% ***************************************************************
if strcmp(Params.Crust,Params.MantleVp) && ...
        strcmp(Params.Crust,Params.MantleVs)
    Depth   = Velocity_Models.Mantle.Vp.Depth;
    Vp      = Velocity_Models.Mantle.Vp.Vp;
    Vs      = Velocity_Models.Mantle.Vp.Vs;
    Crust   = [];
else
    % Crust
    Temp = Velocity_Models.Crust;
    
    if isempty(Temp)
        Crust.Vp=6; Crust.Vs=3.5; Crust.Moho=-1;
    elseif isempty(Temp.Latitude)
        Crust.Vp  	= Temp.Vp;
        Crust.Vs  	= Temp.Vs;
        Crust.Moho	= Temp.Moho;
    else
        % For SESAME, regional models do not extend far enough...
        clear SLOLA str minmax gl nSLAT nSLON
        for directs={'east','west','north','south'}
            direct=directs{1};
            switch direct
                case 'east'; SLOLA='SLON'; str='Longitude'; minmax='max';
                case 'west'; minmax='min';  % otherwise keep from East case
                case 'north'; SLOLA='SLAT'; str='Latitude'; minmax='max';
                case 'south'; minmax='min';  % otherwise keep from North case
            end
            switch minmax;
                case 'min'; gl='<='; pm='+';
                case 'max'; gl='>='; pm='-';
            end
            
            if(eval([SLOLA gl minmax '(Temp.' str ')']))
                disp(['Crustal velocity model does not extend far enough ' ...
                    direct '!  Using closest available ' lower(str) '.'])
                eval(['n' SLOLA ' = ' minmax '(Temp.' str ')' pm '0.1;']);
            end
        end
        repl=0;
        if ~exist('nSLAT','var'); nSLAT=SLAT; else repl=1; end
        if ~exist('nSLON','var'); nSLON=SLON; else repl=1; end
        if(repl~=0); fprintf(fID,['\n%s.%s (%gN, %gE) has crustal migration model' ...
                ' from (%gN, %gE)'], NET,STA,SLAT,SLON,nSLAT,nSLON);
        end
        
        if(sum(strcmp(Params.Crust,{'None','none'})))
            Crust.Vp=6; Crust.Vs=3.5; Crust.Moho=-1;
        else
            Crust.Vp 	= interp2(double(Temp.Longitude),double(Temp.Latitude), ...
                double(Temp.Vp),double(nSLON),double(nSLAT));
            Crust.Vs  	= interp2(double(Temp.Longitude),double(Temp.Latitude), ...
                double(Temp.Vs),double(nSLON),double(nSLAT));
            Crust.Moho	= round(interp2(double(Temp.Longitude),double(Temp.Latitude),...
                double(Temp.Moho),double(nSLON),double(nSLAT)));
        end
    end
    
    % Mantle
    Vs = {'Vp','Vs'};
    for iv=1:length(Vs)
        V = Vs{iv};
        Temp = Velocity_Models.Mantle.(V);
        Mantle.([V '_Depth']) = Temp.Depth;
        if isempty(Temp.Latitude)
            VI = Temp.(V);
        else
            
            % For SESAME, regional models do not extend far enough...
            clear SLOLA str minmax gl nSLAT nSLON
            for directs={'east','west','north','south'}
                direct=directs{1};
                switch direct
                    
                    case 'east'
                        if strcmp(Temp.Name,'GLADM15')
                            SLOLA='SLONm'; str='Longitude'; minmax='max'; gl='>';
                        else
                            SLOLA='SLON'; str='Longitude'; minmax='max'; gl='>';
                        end
                    case 'west'; minmax='min'; gl='<'; % otherwise keep from East case
                    case 'north'; SLOLA='SLAT'; str='Latitude'; minmax='max'; gl='>';
                    case 'south'; minmax='min'; gl='<'; % otherwise keep from North case
                end
                switch minmax;
                    case 'min'; gl='<='; pm='+';
                    case 'max'; gl='>='; pm='-';
                end
                
                if(eval([SLOLA gl minmax '(Temp.' str ')']))
                    disp(['Mantle ' V ' model does not extend far enough ' ...
                        direct '!  Using closest available ' lower(str) '.'])
                    eval(['n' SLOLA ' = ' minmax '(Temp.' str ')' pm '0.1;']);
                end
            end
            
            repl=0;
            if ~exist('nSLAT','var'); nSLAT=SLAT; else repl=repl+1; rep=1; end
            if ~exist('nSLON','var'); nSLON=SLON; else repl=repl+1; rep=1; end
            if(exist('rep','var')); fprintf(fID,['\n%s.%s (%gN, %gE) has mantle ' ...
                    V ' migration model from (%gN, %gE)'], NET,STA, ...
                    SLAT,SLON,nSLAT,nSLON);
            end
            
            [LonG,LatG,DepG] = meshgrid(Temp.Longitude,Temp.Latitude,Temp.Depth);
            LatI = ones(size(Temp.Depth))*nSLAT;
            if strcmp(Temp.Name,'UCB_NA_Vs_3.09') || strcmp(Temp.Name,'GLADM15')
                LonI = ones(size(Temp.Depth))*SLONm;
            else % if goes from -180:+180 longitude
                LonI = ones(size(Temp.Depth))*nSLON;
            end
            DepI = Temp.Depth;
            VI = interp3(double(LonG),double(LatG),double(DepG),...
                double(Temp.(V)),double(LonI),double(LatI),double(DepI));
            % Fix for East Coast project, where top few layers of
            % easternmost velocity profiles are NaNs - move inland (west)
            % until get a profile without NaNs
            if strcmp(Params.MantleVs,'SchmandtLin14')
                while (sum(isnan(VI)))
                    LonI=LonI-0.25;
                    if LatI(1)<27; LatI=LatI+0.25;
                    elseif LatI(1)>49; LatI=LatI-0.25;
                    end
                    VI = interp3(double(LonG),double(LatG),double(DepG),...
                        double(Temp.(V)),double(LonI),double(LatI),double(DepI));
                end
            end
        end
        if max(VI)>100; Mantle.(V) = VI/1000; else Mantle.(V) = VI; end
    end
    
    %if(repl~=0); fprintf(fID,'\n\n'); end
    
    
    if Mantle.Vp_Depth(1)>Mantle.Vs_Depth(1)
        junk=find(Mantle.Vs_Depth<Mantle.Vp_Depth(1));
        Mantle.Vp_Depth=[Mantle.Vs_Depth(junk); Mantle.Vp_Depth];
        if strcmp(Params.MantleVs,'WagnerSESAMEDec14')
            ga_crust=load([basedir '/Data/Velocity_Models/Crust/LaraDec14.mat']);
            k=ga_crust.Crust.Vp./ga_crust.Crust.Vs;  k(isnan(k))=mean(k(~isnan(k)));
            k=interp2(ga_crust.Crust.Longitude,ga_crust.Crust.Latitude,k,nSLON,SLAT);
        else
            k=1.77;
        end
        Mantle.Vp=[Mantle.Vs(junk)*k; Mantle.Vp];
    end
    
    if Mantle.Vp_Depth(end)<Mantle.Vs_Depth(end)
        Mantle_Depth    = Mantle.Vp_Depth;
        Mantle_Vp       = Mantle.Vp;
        Mantle_Vs       = interp1(double(Mantle.Vs_Depth),double(Mantle.Vs),...
            double(Mantle.Vp_Depth));
        for idi=length(Mantle_Depth):-1:1
            if isnan(Mantle_Vs(idi))
                Mantle_Vs(1:idi,1) = ones(idi,1)*Mantle_Vs(idi+1);
                break
            end
        end
    else
        for i = length(Mantle.Vs_Depth):-1:2
            if Mantle.Vs_Depth(i)==Mantle.Vs_Depth(i-1)
                Mantle.Vs_Depth(i) = Mantle.Vs_Depth(i)+0.001;
            end
        end
        for i = length(Mantle.Vp_Depth):-1:2
            if Mantle.Vp_Depth(i)==Mantle.Vp_Depth(i-1)
                Mantle.Vp_Depth(i) = Mantle.Vp_Depth(i)+0.001;
            end
        end
        Mantle_Depth    = Mantle.Vs_Depth;
        Mantle_Vp       = interp1(double(Mantle.Vp_Depth),double(Mantle.Vp),...
            double(Mantle.Vs_Depth));
        for idi=length(Mantle_Depth):-1:1
            if isnan(Mantle_Vp(idi))
                Mantle_Vp(1:idi,1) = ones(idi,1)*Mantle_Vp(idi+1);
                break
            end
        end
        Mantle_Vs       = Mantle.Vs;
    end
    
    % *************************************
    % * Combine Crustal and Mantle Models *
    % *************************************
    
    Grad_Man_Vs = gradient(Mantle_Vs)./gradient(Mantle_Depth);
    ind_Moho = find(Grad_Man_Vs(Mantle_Depth<70 & Mantle_Depth>4)>0.05,1,'last');
    if isempty(ind_Moho)
        [~,ind_Moho] = max(Grad_Man_Vs(Mantle_Depth<70 & Mantle_Depth>4));
    end
    ind_Moho = ind_Moho+1;
    Depth_cr = Mantle_Depth(Mantle_Depth<70 & Mantle_Depth>4);
    Moho_Mantle = Depth_cr(ind_Moho);
    ind_M = find(Mantle_Depth==Moho_Mantle);
    for idm=1:length(Mantle_Depth)
        if Moho_Mantle<=Crust.Moho
            if Mantle_Depth(idm)>Crust.Moho
                if Crust.Moho==-1; cr=[]; crvp=[]; crvs=[];
                else
                    cr=[0; Crust.Moho; Crust.Moho];
                    crvp=[Crust.Vp; Crust.Vp; Mantle_Vp(idm)];
                    crvs=[Crust.Vs; Crust.Vs; Mantle_Vs(idm)];
                end
                
                Depth   = [cr; Mantle_Depth(idm:end)];
                Vp      = [crvp; Mantle_Vp(idm:end)];
                Vs      = [crvs; Mantle_Vs(idm:end)];
                break
            end
        else
            if Mantle_Depth(idm)>Crust.Moho
                if Crust.Moho==-1; cr=[]; crvp=[]; crvs=[];
                else
                    cr=[0; Crust.Moho; Crust.Moho];
                    crvp=[Crust.Vp; Crust.Vp; Mantle_Vp(ind_M)];
                    crvs=[Crust.Vs; Crust.Vs; Mantle_Vs(ind_M)];
                end
                
                Depth   = [cr; Mantle_Depth(ind_M:end)];
                Vp      = [crvp; Mantle_Vp(ind_M:end)];
                Vs      = [crvs; Mantle_Vs(ind_M:end)];
                break
            end
        end
        
    end
end


% *************************************
% * Set a Maximum Depth for Migration *
% *************************************
if strcmp(Phase,'Ps')
    Test_Max_Dep  	= 1000; %VED CRUSTAL vs
else
    Test_Max_Dep  	= 550;
end
Max_Dep_i       = find(abs(Depth(:,1)-Test_Max_Dep)==min(abs(Depth(:,1)-Test_Max_Dep)),1);
Max_Dep         = ceil(Depth(Max_Dep_i));	% Maximum migration depth

if strcmp(isRPdiff,'yes')
    % *****************************************************
    % * Calculate Absolute and Relative Delay Time Tables *
    % *****************************************************
    
    if ~exist([basedir,'Data/Velocity_Models/TaupInfo.mat'],'file')
        diffRP(basedir,{'Ps','Sp'});
    end
    
    TMP = load([basedir,'Data/Velocity_Models/TaupInfo.mat']);
    TaupInfo = TMP.TaupInfo;
    clear TMP;
    
    [SRC,DIS,DEP] = meshgrid(TaupInfo.(Phase).srcdepths,TaupInfo.(Phase).distances,TaupInfo.(Phase).depths);
    [SRc,DIs] = meshgrid(TaupInfo.(Phase).srcdepths,TaupInfo.(Phase).distances);
    MAXDEP = Max_Dep*ones(size(SRc));
    for idi = 1:length(TaupInfo.(Phase).distances)
        for isc = 1:length(TaupInfo.(Phase).srcdepths)
            indtmp = find(~isnan(TaupInfo.(Phase).DIST(idi,isc,:)),1,'last');
            if TaupInfo.(Phase).depths(indtmp) <= Max_Dep+20
                tmpmaxdep = TaupInfo.(Phase).depths(indtmp)-20;
                Max_Dep_i_tmp       = find(Depth(:,1)-tmpmaxdep<=0,1,'last');
                MAXDEP(idi,isc) = Depth(Max_Dep_i_tmp);
            end
        end
    end
    
    BottomDIS = interp3(SRC,DIS,DEP,TaupInfo.(Phase).DIST,SRc,DIs,MAXDEP);
    
    
    % Reference Depth and Reference Ray Parameter Delay Times
    Reference_Depth(:,1) = (0:1:Max_Dep)';  %VED CRUSTAL VS MANTLE DZ
    depth = Depth;
    for idi=1:length(depth)-1
        if depth(idi)==depth(idi+1)
            depth(idi+1) = depth(idi)+1e-3;
        end
    end
    
    indxfinalzero = max(find(depth==0));
    if isempty(indxfinalzero)
        indxfinalzero = 1;
        depth = [0;depth];
        Depth = [0;Depth];
        Vp = [Vp(1);Vp];
        Vs = [Vs(1);Vs];
    end
    
    % Interpolate velocity model onto reference depth vector
    Reference_Vp = interp1(double(depth(indxfinalzero:end)), ...
        double(Vp(indxfinalzero:end)),double(Reference_Depth));
    Reference_Vs = interp1(double(depth(indxfinalzero:end)),...
        double(Vs(indxfinalzero:end)),double(Reference_Depth));
    
    
    RayPara_RNG = linspace(min(min(TaupInfo.(Phase).RayParam)),max(max(TaupInfo.(Phase).RayParam)),200); % Range of Recorded Ray Parameters at surface
    if strcmp(Phase,'Ps')
        [D_i,t_i,~]=jhua_shootray(Reference_Vp,Reference_Depth,RayPara_RNG/Er,Er);
        [D_c,t_c,~]=jhua_shootray(Reference_Vs,Reference_Depth,RayPara_RNG/Er,Er);
    else
        [D_i,t_i,~]=jhua_shootray(Reference_Vs,Reference_Depth,RayPara_RNG/Er,Er);
        [D_c,t_c,~]=jhua_shootray(Reference_Vp,Reference_Depth,RayPara_RNG/Er,Er);
    end
    D_i = [zeros(1,length(RayPara_RNG));D_i];
    t_i = [zeros(1,length(RayPara_RNG));t_i];
    D_c = [zeros(1,length(RayPara_RNG));D_c];
    t_c = [zeros(1,length(RayPara_RNG));t_c];
    
    
    TopDIS = BottomDIS+interp2(repmat(RayPara_RNG,[length(Reference_Depth),1]),...
        repmat(Reference_Depth,[1,length(RayPara_RNG)]),D_i,TaupInfo.(Phase).RayParam,MAXDEP);
    
    DIS_I = repmat(TopDIS,[1,1,length(Reference_Depth)])-interp2(repmat(RayPara_RNG,[length(Reference_Depth),1]),...
        repmat(Reference_Depth,[1,length(RayPara_RNG)]),D_i,repmat(TaupInfo.(Phase).RayParam,[1,1,length(Reference_Depth)])...
        ,repmat(permute(Reference_Depth,[3,2,1]),[size(TaupInfo.(Phase).RayParam,1),size(TaupInfo.(Phase).RayParam,2),1]));
    
    DIS_C = interp2(repmat(RayPara_RNG,[length(Reference_Depth),1]),...
        repmat(Reference_Depth,[1,length(RayPara_RNG)]),D_c,repmat(TaupInfo.(Phase).RayParam,[1,1,length(Reference_Depth)])...
        ,repmat(permute(Reference_Depth,[3,2,1]),[size(TaupInfo.(Phase).RayParam,1),size(TaupInfo.(Phase).RayParam,2),1]));
    
    for idi = 1:length(TaupInfo.(Phase).distances)
        for isc = 1:length(TaupInfo.(Phase).srcdepths)
            indtmp = find(Reference_Depth>MAXDEP(idi,isc));
            DIS_I(idi,isc,indtmp) = nan;
            DIS_C(idi,isc,indtmp) = nan;
        end
    end
    
    if strcmp(Phase,'Ps')
        Distances = 29:1.5:86;
    else
        Distances = 49:1.5:91;
    end
    
    [Src,Dis] = meshgrid(TaupInfo.(Phase).srcdepths,Distances);
    inds_N = find(~isnan(TaupInfo.(Phase).RayParam(:)) & ~isnan(TopDIS(:)));
    F = scatteredInterpolant(SRc(inds_N),TopDIS(inds_N),TaupInfo.(Phase).RayParam(inds_N),'linear','none');
    RayParam_model = F(Src,Dis);
    
    RP_depths = zeros(length(Distances),length(TaupInfo.(Phase).srcdepths),length(Reference_Depth));
    for isr = 1:length(TaupInfo.(Phase).srcdepths)
        for idep = 1:length(Reference_Depth)
            
            tmpdis = DIS_I(:,isr,idep)+DIS_C(:,isr,idep);
            ind_N_tmp = find(~isnan(tmpdis) & ~isnan(TaupInfo.(Phase).RayParam(:,isr)));
            if length(ind_N_tmp)>=2
                RP_depths(:,isr,idep) = interp1(tmpdis(ind_N_tmp),TaupInfo.(Phase).RayParam(ind_N_tmp,isr),Distances);
            else
                RP_depths(:,isr,idep) = nan;
            end
        end
    end
    
    DT = zeros(length(Distances),length(TaupInfo.(Phase).srcdepths),length(Reference_Depth));
    for isrc = 1:length(TaupInfo.(Phase).srcdepths)
        for idis = 1:length(Distances)
            RP_conv = permute(RP_depths(idis,isrc,:),[1,3,2]);
            if strcmp(Phase,'Ps')
                [D_i,t_i,~]=jhua_shootray(Reference_Vp,Reference_Depth,RayParam_model(idis,isrc)/Er,Er);
                [D_c,t_c,~]=jhua_shootray(Reference_Vs,Reference_Depth,RP_conv/Er,Er);
            else
                [D_i,t_i,~]=jhua_shootray(Reference_Vs,Reference_Depth,RayParam_model(idis,isrc)/Er,Er);
                [D_c,t_c,~]=jhua_shootray(Reference_Vp,Reference_Depth,RP_conv/Er,Er);
            end
            if isnan(D_i(1))
                D_i = [nan;D_i];
                t_i = [nan;t_i];
            else
                D_i = [0;D_i];
                t_i = [0;t_i];
            end
            D_c = [zeros(1,length(RP_conv));D_c];
            t_c = [zeros(1,length(RP_conv));t_c];
            D_c(1,isnan(D_c(2,:))) = nan;
            t_c(1,isnan(t_c(2,:))) = nan;
            t_c = diag(t_c);
            D_c = diag(D_c);
            Pdelta = -(RayParam_model(idis,isrc)+RP_conv)'.*(D_c-D_i)/180*pi/2;
            DT(idis,isrc,:) = t_c+Pdelta-t_i;
        end
    end
    
    % *****************************************************
    % * Calculate Critical Ray Parameters Down To Max_Dep *
    % *****************************************************
    if strcmp(Phase,'Ps')
        V1 = Reference_Vp;	V2 = Reference_Vs;
    elseif strcmp(Phase,'Sp')
        V1 = Reference_Vs; 	V2 = Reference_Vp;
    end
    for idp=1:length(Reference_Depth)-1
        Critical_Angle = asind(V1(idp+1)/V2(idp));
        if imag(Critical_Angle)==0
            Critical_RPs(idp,1) = (Er-Reference_Depth(idp))/Er*sind(Critical_Angle)/V1(idp+1);
        else
            Critical_RPs(idp,1) = 666;
        end
    end
    
    if ~sum(strcmp(Params.Crust,{'None','none'}))
        Migration_Data.Original_Models.Crust = Velocity_Models.Crust.Name;
        if ~isempty(Crust)
            Migration_Data.Original_Models.Crust            = Velocity_Models.Crust.Name;
            Migration_Data.Original_Models.Crustal_Model    = Crust;
        end
    else
        Migration_Data.Original_Models.Crust = 'None';
        Migration_Data.Original_Models.Crustal_Model = Crust;
    end
    Migration_Data.Original_Models.Mantle_Vp = Velocity_Models.Mantle.Vp.Name;
    Migration_Data.Original_Models.Mantle_Vs = Velocity_Models.Mantle.Vs.Name;
    Migration_Data.Velocity_Model.Vp   	= Vp;
    Migration_Data.Velocity_Model.Vs   	= Vs;
    Migration_Data.Er                   = Er;
    Migration_Data.Velocity_Model.Depth	= Depth;
    Migration_Data.Reference_Depth      = Reference_Depth;
    Migration_Data.Reference_SrcDepth   = TaupInfo.(Phase).srcdepths;
    Migration_Data.Reference_Distance   = Distances;
    Migration_Data.Reference_DT         = DT;
    Migration_Data.Reference_RP         = RP_depths;
    Migration_Data.Critical_RPs         = Critical_RPs;
    
else
    % *****************************************************
    % * Calculate Absolute and Relative Delay Time Tables *
    % *****************************************************
    nRPs = 20;
    if ceil(nRPs/2)~=nRPs/2
        nRPs = nRPs+1;
    end
    import edu.sc.seis.TauP.*	% Import the TauP package
    if strcmp(Phase,'Ps')
        Ref_Dist        = 60;       % Reference Ps distance
        DistRng         = [30 90];
        fctr            = 1;
    elseif strcmp(Phase,'Sp')
        Ref_Dist        = 65;       % Starting reference Sp distance
        DistRng         = [50 90];
        fctr            = -1;
    end
    % Use a nominal event depth of 100 km
    Temp_Ray_Path   = Matlab_TauP('Path','ak135',100,Phase(1),'sta',[0 0],'evt',[0 Ref_Dist]);
    Reference_RP    = Temp_Ray_Path.rayParam/Er;
    Temp_Ray_Path   = Matlab_TauP('Path','ak135',100,Phase(1),'sta',[0 0],'evt',[0 DistRng(1)]);
    Test_RP_1       = Temp_Ray_Path.rayParam/Er;
    Temp_Ray_Path   = Matlab_TauP('Path','ak135',100,Phase(1),'sta',[0 0],'evt',[0 DistRng(2)]);
    Test_RP_2       = Temp_Ray_Path.rayParam/Er;
    dRP             = (Test_RP_1-Test_RP_2)/nRPs;
    
    % Reference Depth and Reference Ray Parameter Delay Times
    Reference_Depth(:,1) = 0:1:Max_Dep;  %VED CRUSTAL VS MANTLE DZ
    depth = Depth;
    for idi=1:length(depth)-1
        if depth(idi)==depth(idi+1);
            depth(idi+1) = depth(idi)+1e-3;
        end
    end
    
    indxfinalzero = max(find(depth==0));
    if isempty(indxfinalzero)
        indxfinalzero = 1;
        depth = [0;depth];
        Depth = [0;Depth];
        Vp = [Vp(1);Vp];
        Vs = [Vs(1);Vs];
    end
    
    
    % Interpolate velocity model onto reference depth vector
    Reference_Vp = interp1(double(depth(indxfinalzero:end)), ...
        double(Vp(indxfinalzero:end)),double(Reference_Depth));
    Reference_Vs = interp1(double(depth(indxfinalzero:end)),...
        double(Vs(indxfinalzero:end)),double(Reference_Depth));
    
    
    % % Calculate Reference and Relative Ray Parameter(s) and Delay Times
    disp(['Calculating ' Phase ' delay times for ' STA ' (' NET ')...'])
    Reference_DT = fctr*real(Calc_dT(Reference_RP*Er,Reference_Depth,Vp,Vs,Depth)');
    %
    Relative_RPs = Reference_RP-nRPs/2*dRP:dRP:Reference_RP+nRPs/2*dRP;
    for irp=1:length(Relative_RPs)
        fprintf(' %g%% ',round(irp/(length(Relative_RPs)+1)*100))
        if rem(irp,15)==0; fprintf('\n'); end
        rp = Relative_RPs(irp)*Er;
        Relative_DT         = fctr*real(Calc_dT(rp,Reference_Depth,Vp,Vs,Depth)');
        Relative_DTs(:,irp)	= Relative_DT-Reference_DT;
    end
    disp('         ...done!')
    
    
    % *****************************************************
    % * Calculate Critical Ray Parameters Down To Max_Dep *
    % *****************************************************
    if strcmp(Phase,'Ps')
        V1 = Reference_Vp;	V2 = Reference_Vs;
    elseif strcmp(Phase,'Sp')
        V1 = Reference_Vs; 	V2 = Reference_Vp;
    end
    for idp=1:length(Reference_Depth)-1
        Critical_Angle = asind(V1(idp+1)/V2(idp));
        if imag(Critical_Angle)==0
            Critical_RPs(idp,1) = (Er-Reference_Depth(idp))/Er*sind(Critical_Angle)/V1(idp+1);
        else
            Critical_RPs(idp,1) = 666;
        end
    end
    
    if ~sum(strcmp(Params.Crust,{'None','none'}))
        Migration_Data.Original_Models.Crust = Velocity_Models.Crust.Name;
        if ~isempty(Crust)
            Migration_Data.Original_Models.Crust            = Velocity_Models.Crust.Name;
            Migration_Data.Original_Models.Crustal_Model    = Crust;
        end
    else
        Migration_Data.Original_Models.Crust = 'None';
        Migration_Data.Original_Models.Crustal_Model = Crust;
    end
    Migration_Data.Original_Models.Mantle_Vp = Velocity_Models.Mantle.Vp.Name;
    Migration_Data.Original_Models.Mantle_Vs = Velocity_Models.Mantle.Vs.Name;
    Migration_Data.Velocity_Model.Vp   	= Vp;
    Migration_Data.Velocity_Model.Vs   	= Vs;
    Migration_Data.Velocity_Model.Depth	= Depth;
    Migration_Data.Reference_Depth      = Reference_Depth;
    Migration_Data.Reference_RP         = Reference_RP;
    Migration_Data.Reference_DT         = Reference_DT;
    Migration_Data.Relative_RPs         = Relative_RPs;
    Migration_Data.Relative_DTs         = Relative_DTs;
    Migration_Data.Critical_RPs         = Critical_RPs;
    
    
end

    function [dT] = Calc_dT(p_srad,H,modelVp,modelVs,modelDepth)
        
        % Earth radius (km)
        r_e = 6371;
        
        % Radial (in the sense of spherical geometry) slowness
        U_r = @(r,p,v) (v.^-2-(p./r).^2).^(0.5);
        
        % Differential radial slowness for Ps
        dU_r = @(r,p) (U_r(r,p,Vs(r_e-r))-U_r(r,p,Vp(r_e-r)));
        
        % Find discontinuity depths in our model
        allDiscon = cat(min([2 size(modelDepth,2)]),0,...
            modelDepth((modelDepth(1:end-1)-modelDepth(2:end))==0));
        dT = zeros(length(p_srad),length(H));
        %       slow! re-evaulates the integral from the lowest bound every time.
        %        for ip = 1:length(p_srad)
        %            % Loop over conversion depths
        %            for iH=1:length(H)
        %                % Integrate dT = dUps_r dr to surface from depth H
        %                discon = allDiscon(allDiscon<H(iH));
        %                if isempty(discon)
        %                    dT(ip,iH) = 0;
        %                else
        %                    dT(ip,iH) = quadl(@(r)dU_r(r,p_srad(ip)),r_e-H(iH),r_e-discon(end));
        %                    for iDiscon=length(discon):-1:2
        %                        dT(ip,iH) = dT(ip,iH)+quadl(@(r)dU_r(r,p_srad(ip)),r_e-discon(iDiscon),r_e - discon(iDiscon - 1));
        %                    end
        %                end
        %            end
        %        end
        
        %       faster! integrates H cumulatively
        rad=r_e - H;
        for ip = 1:length(p_srad)
            % Loop over conversion depths
            for iH=2:length(H)
                dTnew  = quadl(@(r)dU_r(r,p_srad(ip)),rad(iH),rad(iH-1));
                dT(ip,iH) = dT(ip,iH-1)+dTnew;
            end
        end
        
        % Recovery of linearly-interpolated velocities
        function [a] = Vp(d)
            for ik=1:length(d)
                if d(ik)==modelDepth(1)
                    a(ik) = modelVp(1);
                else
                    id   = find(modelDepth>=d(ik),1);
                    a(ik) = modelVp(id-1)+(modelVp(id)-modelVp(id-1))* ...
                        (d(ik)-modelDepth(id-1))/(modelDepth(id)-modelDepth(id-1));
                end
            end
        end
        
        function [b] = Vs(d)
            for ik=1:length(d)
                if d(ik)==modelDepth(1)
                    b(ik) = modelVs(1);
                else
                    id   = find(modelDepth>=d(ik),1);
                    b(ik) = modelVs(id-1)+(modelVs(id)-modelVs(id-1)) * ...
                        (d(ik)-modelDepth(id-1))/(modelDepth(id)-modelDepth(id-1));
                end
            end
        end
        
    end
end

function [D,t,sn]=jhua_shootray(v,z,p,re)

% Modified by Junlin to return angular distance D and t for each layer.
% SHOOTRAY: similar to RAYFAN but with less error checking (faster)
%
% [x,t]=shootray(v,z,p)
%
% SHOOTRAY shoots a single ray or fan of rays through a stratified
% velocity model.
% It is assumed that v is a vector of interval velocities and that
% the ray goes from z(1) to z(length(z)). This routine does no error
% checking for efficiency. Be sure that v and z are both the same
% length and that p is a row vector. It returns inf if a critical
% refration occurs. This routine will get the same results as rayfan
% but is 25% faster
%
% v... COLUMN vector of interval velocities
% z... COLUMN vector of depths to the tops of the intervals
% NOTE: v and z must be at least length 2 or an abort will occur
% p... scalar or ROW vector ray parameters (Must be a row vector.)
% x ... scalar or vector of horizontal distances
% t ... scalar or vector of traveltimes
%
% G.F. Margrave, CREWES Project, July 1995
%
% NOTE: It is illegal for you to use this software for a purpose other
% than non-profit education or research UNLESS you are employed by a CREWES
% Project sponsor. By using this software, you are agreeing to the terms
% detailed in this software's Matlab source file.

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by
% its author (identified above) and the CREWES Project.  The CREWES
% project may be contacted via email at:  crewesinfo@crewes.org
%
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) Use of this SOFTWARE by any for-profit commercial organization is
%    expressly forbidden unless said organization is a CREWES Project
%    Sponsor.
%
% 2) A CREWES Project sponsor may use this SOFTWARE under the terms of the
%    CREWES Project Sponsorship agreement.
%
% 3) A student or employee of a non-profit educational institution may
%    use this SOFTWARE subject to the following terms and conditions:
%    - this SOFTWARE is for teaching or research purposes only.
%    - this SOFTWARE may be distributed to other students or researchers
%      provided that these license terms are included.
%    - reselling the SOFTWARE, or including it or any portion of it, in any
%      software that will be resold is expressly forbidden.
%    - transfering the SOFTWARE in any form to a commercial firm or any
%      other for-profit organization is expressly forbidden.
%
% END TERMS OF USE LICENSE

%check for critical angle
iprop=1:length(z)-1;
v_u = (v(iprop)+v(iprop+1))/2;
z_u = (z(iprop)+z(iprop+1))/2;

sn = v_u*re./(re-z_u)*p;

%
% sn is an n by m matrix where n is the length of iprop (the number of
%    layers propagated through) and m is the length of p (the number of
%    unique ray parameters to use). Each column of sn corresponds to a
%    single ray parameter and contains the sin of the vertical angle;
%

%ichk=find(sn>1);

%compute x and t
cs=sqrt(1-sn.*sn)+eps;
vprop=v_u*ones(1,length(p));
thk=abs(diff(z))*ones(1,length(p));
R = (re-z_u)*ones(1,length(p));
if(size(sn,1)>1)
    D=cumsum( (thk.*sn)./cs./R)/pi*180;
    t=cumsum(thk./(vprop.*cs));
else
    D=(thk.*sn)./cs./R/pi*180;
    t=thk./(vprop.*cs);
end
%assign infs
%if(~isempty(ichk))
ind = find(imag(cs)~=0 | imag(D)~=0 | imag(t)~=0 | sn>1);
D(ind)=NaN;
t(ind)=NaN;
sn(ind) = NaN;


%end


end


function [Event_IDs,dt, Culled] = Cull_Data(All_WV_Data,All_RP_Data,RP_RNG,Params,NET,STA,Phase)

% ********* Function Description *********
%
% Load in waveform data that meets the
% specified requirements.
%
%
% ****************************
% >                          <
% >  David L. Abt - 4.2009	 <
% >                          <
% ****************************

disp(['Culling ' STA ' (' NET ') ' Phase ' waveforms...'])
Event_IDs   = {};  dt          = [];
eq_ids      = fieldnames(All_WV_Data);
lpwd        = length(eq_ids); Culled=[];

for ipwd=1:lpwd % for each event
    
    eq_id  	= eq_ids{ipwd};
    S2N  	= All_WV_Data.(eq_id).([Phase '_Waveforms']).S2N_ratio;
    if S2N>=Params.S2N
        ZR_Xcorr = All_WV_Data.(eq_id).([Phase '_Waveforms']).ZR_Xcorr_Coeff;
        
        if isfield(All_WV_Data.(eq_id).([Phase '_Waveforms']).Window,'Taup_S2N_Misfit')
            TauP_S2N = abs(All_WV_Data.(eq_id).([Phase '_Waveforms']).Window.Taup_S2N_Misfit);
        else
            TauP_S2N = 0; display('Taup_S2N_Misfit not found in Waveforms.mat!');
        end
        
        Depth       = All_RP_Data.(eq_id).([Phase(1) '_Path']).srcDepth;
        RP          = All_RP_Data.(eq_id).([Phase(1) '_Path']).ray_parameter;
        
        % Check to see if this particular event/waveform passes the other filters
        %            disp(['TauP-S2N: ' num2str(TauP_S2N) '; Dpth: ' num2str(Depth) '; RP: ' num2str(RP)])
        if TauP_S2N>Params.TauP_S2N
            Culled.(eq_id).TauP_S2N = TauP_S2N;
        elseif Depth>Params.Max_Depth
            Culled.(eq_id).Depth = Depth;
        elseif RP<RP_RNG(1) || RP>RP_RNG(2)
            Culled.(eq_id).RP = RP;
        elseif abs(ZR_Xcorr)<Params.ZR_Xcorr
            Culled.(eq_id).Xcorr = ZR_Xcorr;
        else
            if isempty(dt)
                dt = All_WV_Data.(eq_id).([Phase '_Waveforms']).dt;
            end
            Event_IDs = [Event_IDs; eq_id];
        end
    end
    
end

disp('         ...done!')


end

function [Bin_Data] = Set_Bin_Parameters(Bin_By,Phase)

% **
% *

% This is where you can manually set the bin ranges
baz_min = 0;    baz_max = 360;  baz_inc = 60;
if strcmp(Phase,'Ps')
    dist_min = 35;  dist_max = 80;  dist_inc = 5;
elseif strcmp(Phase,'Sp')
    dist_min = 55;  dist_max = 85;  dist_inc = 5;        %distance limit
end

if strcmpi(Bin_By,'Distance')
    BAZs = [baz_min baz_max];
    DISTs = [dist_min:dist_inc:dist_max-dist_inc; ...
        dist_min+dist_inc:dist_inc:dist_max]';
elseif strcmpi(Bin_By,'bAz')
    BAZs = [baz_min:baz_inc:baz_max-baz_inc; ...
        baz_min+baz_inc:baz_inc:baz_max]';
    DISTs = [dist_min dist_max];
elseif strcmp(Bin_By,'both')
    BAZs = [baz_min:baz_inc:baz_max-baz_inc; ...
        baz_min+baz_inc:baz_inc:baz_max]';
    DISTs = [dist_min:dist_inc:dist_max-dist_inc; ...
        dist_min+dist_inc:dist_inc:dist_max]';
end

Bin_Data.Event_IDs = {};
for ibaz=1:length(BAZs(:,1))
    Bin_Data.BAZ_Bins(ibaz,1).BAZ_Range = BAZs(ibaz,:);
    Bin_Data.BAZ_Bins(ibaz,1).Event_IDs	= {};
    for idist=1:length(DISTs(:,1))
        Bin_Data.BAZ_Bins(ibaz,1).DIST_Bins(idist,1).DIST_Range = DISTs(idist,:);
        Bin_Data.BAZ_Bins(ibaz,1).DIST_Bins(idist,1).Event_IDs	= {};
    end
end

end

function [Waveform_Data,Time_Vector,Bin_Data,Culled] = Bin_and_Normalize(Event_IDs,All_WV_Data,...
    All_RP_Data,Bin_Data,Params,NET,STA,Phase,tag,SvSh,Culled,tagdeconv,varargin)

% ********* Function Description *********
%
% Bin data based on a specified parameter
% (e.g., distance, azimuth), the normalize
% and weight all binned waveforms
%
%
% ****************************************
% *                                      *
% *  Written by David L. Abt - May 2008	 *
% *                                      *
% *  Email: David_Abt@brown.edu          *
% *                                      *
% ****************************************

disp(['Binning and normalizing ' STA ' (' NET ') ' Phase '...'])

if ~isempty(Event_IDs)
    
    pre_time	= All_WV_Data.(Event_IDs{1}).([Phase '_Waveforms']).Window.pre_phase_time;
    post_time	= All_WV_Data.(Event_IDs{1}).([Phase '_Waveforms']).Window.post_phase_time;
    wl        	= pre_time+post_time;
    dt        	= All_WV_Data.(Event_IDs{1}).([Phase '_Waveforms']).dt;
    n2        	= 2^nextpow2(wl/dt);  	% nearest power of 2 for final window length
    tpr       	= 0.5;               	% Taper length (s)
    Time_Vector	= -pre_time:dt:-pre_time+(n2-1)*dt;
    iRFwv       = 0;
    
    lf = Params.lf;    hf = Params.hf;
    % Changed by Ved to deal with SCOOBA
    
    max_Amp_scl = 10;	% Maximum amplitude of the daughter component relative to the parent
    
    
    % Bandpass filter, normalize, and taper each binned waveform
    for ie=1:length(Event_IDs)
        eq_id = Event_IDs{ie}; binned='no';
        
        % N.B. SH component is the transverse signal divided by the
        % free surface reflection coefficient (2)
        if strcmpi(SvSh,'sh');
            sh = All_WV_Data.(eq_id).([Phase '_Waveforms']).T./2;
            All_WV_Data.(eq_id).([Phase '_Waveforms']).SV=sh;
        end
        
        % Determine the masking depth and time if they exist
        RP   	= All_RP_Data.(eq_id).([Phase(1) '_Path']).ray_parameter;
        
        if ~strcmp(tagdeconv,'HK')
            Migration_Data = varargin{1};
            mask_i	= find(RP>=Migration_Data.Critical_RPs);
            if isempty(mask_i)
                mask_i = length(Migration_Data.Reference_Depth);
            else
                mask_i = mask_i(1);
            end
            Mask_Depth	= Migration_Data.Reference_Depth(mask_i);
        end
        
        Time    = All_WV_Data.(eq_id).([Phase '_Waveforms']).Time;
        
        
        if(strcmpi(Params.PreFilter,'no'))
            % Don't Bandpass filter waveforms
            P   = All_WV_Data.(eq_id).([Phase '_Waveforms']).P;
            SV  = All_WV_Data.(eq_id).([Phase '_Waveforms']).SV;
        elseif strcmpi(Params.PreFilter,'yes')
            if strcmpi(tag,'synth')
                P=lpfilt(All_WV_Data.(eq_id).([Phase '_Waveforms']).P,dt,hf);
                SV=lpfilt(All_WV_Data.(eq_id).([Phase '_Waveforms']).SV,dt,hf);
            else
                % Bandpass filter waveforms
                P	= bpfilt(All_WV_Data.(eq_id).([Phase '_Waveforms']).P,dt,lf,hf,'bp');
                SV  = bpfilt(All_WV_Data.(eq_id).([Phase '_Waveforms']).SV,dt,lf,hf,'bp');
            end
        end
        %phase_st = All_WV_Data.(eq_id).([Phase '_Waveforms']).Window.phase_window(1)-2.02; %2.02 for SWUS
        phase_st = All_WV_Data.(eq_id).([Phase '_Waveforms']).Window.phase_window(1)-2.0; %1 s for MOOS in Alaska
        phase_et = All_WV_Data.(eq_id).([Phase '_Waveforms']).Window.phase_window(1)+20.40;
        
        % Taper the waveforms with the phase of interest and pad with zeros to
        % optimize the fourier transform (i.e., make the length=n2; a factor of 2)
        if strcmpi(Params.Taper,'beginning'); Taper = [phase_st -1];
        elseif strcmpi(Params.Taper,'end');   Taper = [-1 phase_et];
        elseif strcmpi(Params.Taper,'both');  Taper = [phase_st phase_et];
        end
        
        if strcmp(Phase,'Ps')
            P_post_tpr  = 5;
            P  = taper(P,Time,P_post_tpr,dt,phase_st,phase_et);
            max_Parent_Amp      = max(abs(P));
            max_Daughter_Amp    = max(abs(SV));
        elseif strcmp(Phase,'Sp')
            
            %             plot(SV);
            %             hold on
            
            %Taper(1) = Taper(1)-2;  %Junlin - don't use!
            SV	= taper(SV,Time,tpr,dt,Taper(1),Taper(2));
            
            %             plot(SV,'r');
            %             hold off
            %             drawnow
            
            max_Parent_Amp      = max(abs(SV));
            max_Daughter_Amp    = max(abs(P));
        end
        pad = n2-length(Time);
        switch tag
            case 'synth';
                switch Phase;
                    % Need to pad such that the parent phase arrives at
                    % about 25s from the start (for Ps) or end (for Sp) of
                    % the waveform as this is where the parent phase is
                    % windowed in the deconvolution!
                    case 'Sp'; P = [zeros(1,pad) P];  SV = [zeros(1,pad) SV];
                    case 'Ps'; P = [P zeros(1,pad)];  SV = [SV zeros(1,pad)];
                end
            case 'data';    P = [P zeros(1,pad)];  SV = [SV zeros(1,pad)];
        end
        
        if max_Daughter_Amp<=max_Amp_scl*max_Parent_Amp
            % Collect in appropriate bins
            BAZ     = All_RP_Data.(eq_id).([Phase(1) '_Path']).BAZ;
            DIST    = All_RP_Data.(eq_id).([Phase(1) '_Path']).distance;
            for iBb=1:length(Bin_Data.BAZ_Bins)
                if BAZ>=Bin_Data.BAZ_Bins(iBb).BAZ_Range(1) && ...
                        BAZ<=Bin_Data.BAZ_Bins(iBb).BAZ_Range(2)
                    for iDb=1:length(Bin_Data.BAZ_Bins(iBb).DIST_Bins);
                        if DIST>=Bin_Data.BAZ_Bins(iBb).DIST_Bins(iDb).DIST_Range(1) && ...
                                DIST<=Bin_Data.BAZ_Bins(iBb).DIST_Bins(iDb).DIST_Range(2)
                            if strcmp(Phase,'Ps') && ~strcmp(tagdeconv,'HK')
                                if (DIST<55 && DIST>51) || (DIST<64 && DIST>59) || DIST<37 %% avoid 660 interference with PP, and 410 & 660 interference with PcP
                                    continue;
                                end
                            end
                            
                            % Normalize waveforms by max amp of parent comp
                            P_norm	= P/max_Parent_Amp;
                            SV_norm	= SV/max_Parent_Amp;
                            % Weight by S2N of either parent or daughter
                            if strcmp(Params.Weight_By{1}(1),Phase(1))
                                S2N	= All_WV_Data.(eq_id).([Phase '_Waveforms']).S2N_ratio;
                            else
                                if strcmp(Phase,'Ps')
                                    S2N	= 1/std(SV_norm(1:450));
                                elseif strcmp(Phase,'Sp')
                                    S2N	= 1/std(P_norm(1:450));
                                end
                                if S2N>40
                                    S2N	= 40;
                                end
                            end
                            if iRFwv==0 || sum(strcmp(fieldnames(RF_Waveforms),eq_id))==0
                                RF_Waveforms.(eq_id).P          = P_norm*S2N;
                                RF_Waveforms.(eq_id).SV         = SV_norm*S2N;
                                RF_Waveforms.(eq_id).RP         = RP;
                                if ~strcmp(tagdeconv,'HK')
                                    RF_Waveforms.(eq_id).Mask_Depth = Mask_Depth;
                                    %RF_Waveforms.(eq_id).Mask_Time  = Mask_Time;
                                end
                                iRFwv = iRFwv+1;
                            end
                            % Add Event_ID to bin lists
                            Bin_Data.BAZ_Bins(iBb).Event_IDs = ...
                                [Bin_Data.BAZ_Bins(iBb).Event_IDs; eq_id];
                            Bin_Data.BAZ_Bins(iBb).DIST_Bins(iDb).Event_IDs = ...
                                [Bin_Data.BAZ_Bins(iBb).DIST_Bins(iDb).Event_IDs; eq_id];
                            
                            binned='yes';
                        end
                    end
                end
            end
        else
            fprintf(' Daughter component amplitude for %s is > %gx parent component amplitude!\n',...
                eq_id,max_Amp_scl)
        end
        if strcmp(binned,'no'); Culled.(eq_id).GCARC=9999; end
        %if strcmp(binned,'no'); Culled.(eq_id).GCARC=DIST; end  %origin one Junlin
    end
    disp('         ...done!')
    
    if iRFwv~=0
        % Get Event IDs of all waveforms that have been binned
        Unique_Event_IDs = unique(fieldnames(RF_Waveforms));
        n_Shared_RF_Waveforms = length(Unique_Event_IDs)-length(fieldnames(RF_Waveforms));
        if n_Shared_RF_Waveforms~=0
            fprintf(' %g events are found in more than one bin!!\n',n_Shared_RF_Waveforms)
            for iuid=1:length(Unique_Event_IDs)
                eq_id = Unique_Event_IDs{iuid};
                Unique_RF_Waveforms.(eq_id) = RF_Waveforms.(eq_id);
            end
        else
            Unique_RF_Waveforms = RF_Waveforms;
        end
        Bin_Data.Event_IDs  = Unique_Event_IDs;
        Waveform_Data       = Unique_RF_Waveforms;
    else
        Waveform_Data   = [];
        Time_Vector 	= [];
    end
else
    Waveform_Data   = [];
    Time_Vector 	= [];
end

    function [y] = taper(x,time,tpr,dt,t1,t2)
        
        % ********* Function Description *********
        %
        % TAPER  Taper a time series.
        %
        % TAPER(X,TIME,TPR,DT,T1,T2) takes time
        % series sampled at DT and tapers it with
        % a cosine taper TPR seconds long from
        % beginning point T1-TPR and with reverse
        % cosine taper from point T2 to point T2+
        % TPR. Points outside the range (T1-TPR,
        % T2+TPR) are zeroed. If T1/T2 is negative
        % then taper is not implemented at the
        % beginning/end. If X is an array of
        % seismograms, then the taper is applied
        % to each row of X.
        %
        %
        % ****************************************
        % *                                      *
        % *  Modified from Kate Rychert's        *
        % *  receiver function code - May 2008   *
        % *                                      *
        % *  Email: David_Abt@brown.edu          *
        % *                                      *
        % ****************************************
        
        nn      = length(x(1,:));
        nx      = length(x(:,1));
        taper   = ones(1,nn);
        it  	= [0:fix(tpr/dt)]*dt/tpr;
        ct      = 0.5-0.5*cos(pi*it);
        T1      = fix(time(1)/dt);             % Absolute sample point of first time step
        it1     = fix(t1/dt+1)-T1;
        it2     = fix(t2/dt+1)-T1;
        
        % Emily, 6th May 2013: problem with start time of phase being <100s, so
        % taper subscripts are negative.
        % Temporary workaround, set taper length to be shorter.
        % N.B. Only one event so far has had this issue!
        
        if t1>0
            if it1>fix(tpr/dt)
                taper(it1-fix(tpr/dt):it1)	= ct;
                taper(1:it1-fix(tpr/dt))    = zeros(size(1:it1-fix(tpr/dt)));
            else
                taper(1:it1) = ct(fix(tpr/dt)-it1+2:end);
                taper(1) = 0;
                disp('Bizarre taper!')
                
            end
        end
        if t2>0
            if t2>time(end)-tpr
                t2  = time(end)-tpr;
                it2	= fix(t2/dt)-T1;
            end
            taper(it2:it2+fix(tpr/dt))	= fliplr(ct);
            taper(it2+fix(tpr/dt):nn)	= zeros(size(taper(it2+fix(tpr/dt):nn)));
        end
        
        y = zeros(nx,nn);
        for ix=1:nx
            y(ix,:) = x(ix,:).*taper;
        end
        
    end

end


function Deconvolve_Data(WV_Data,Bin_Data,RFPs,dt,NET,...
    STA,Phase,Data_Dir,lowT,highT, tag,savename)

% ********* Function Description *********
%
% Deconvolve waveform(s)
%
%
% ****************************************
% *                                      *
% *  Written by David L. Abt - May 2008	 *
% *                                      *
% *  Email: David_Abt@brown.edu          *
% *                                      *
% ****************************************

% This calculates with ETMTM deconvolution ONLY, saved as
% PreCalculated_RFs_[params].mat

% *************
% * Single RF *
% *************
% Generate a single receiver function in time from all waveforms
All_Event_IDs = Bin_Data.Event_IDs;
disp(['Calculating single ' STA ' (' NET ') ' Phase ' RF...'])

clear P D RPs RFs
if ~isempty(All_Event_IDs)
    for ie=1:length(All_Event_IDs)
        eq_id           = All_Event_IDs{ie};
        RPs(ie)      	= WV_Data.(eq_id).RP;
        if strcmp(Phase,'Ps')
            P(ie,:)	= WV_Data.(eq_id).P;    % Parent
            D(ie,:)	= WV_Data.(eq_id).SV;   % Daughter
%             P(ie,:)	= WV_Data.(eq_id).Z;    % Parent
%             D(ie,:)	= WV_Data.(eq_id).R;   % Daughter
        elseif strcmp(Phase,'Sp')
            P(ie,:)	= WV_Data.(eq_id).SV;   % Parent
            D(ie,:)	= WV_Data.(eq_id).P;   	% Daughter
        end
    end
    
    Get_PreCalculated_RFs(Data_Dir,P,D,All_Event_IDs,...
        NET,STA,Phase,RFPs,dt,lowT,highT,tag, savename);
end

disp('         ...done')

end

function [PreCalculated_RFs] = Get_PreCalculated_RFs(Data_Dir,P,D, ...
    All_Event_IDs,NET,STA,Phase,RFPs,dt,lowT,highT,tagd,savename)

% Recalculate everything for purposes of CCP stack
load([Data_Dir 'PreCalculated_RFs_',num2str(lowT),'_',num2str(highT),'.mat'])
junk=dir([Data_Dir '/phasewidth_*mat']);
if ~isempty(junk); junk=strsplit(junk.name,'_'); pwl=str2double(junk{2}(1:end-5));
else pwl=[];
end
ie_ReCalc = 1:length(All_Event_IDs);

[PreCalculated_RFs] = Calculate_Individual_Time_Domain_RFs(P(ie_ReCalc,:),...
    D(ie_ReCalc,:),All_Event_IDs(ie_ReCalc),dt,NET,STA,Phase,PreCalculated_RFs,pwl,tagd);

save([Data_Dir 'PreCalculated_RFs' savename],'PreCalculated_RFs')
%save([Data_Dir 'PreCalculated_RFs_ZRT' savename],'PreCalculated_RFs')


    function [PreCalculated_RFs] = Calculate_Individual_Time_Domain_RFs(P,D, ...
            Event_IDs,dt,NET,STA,Phase,PreCalculated_RFs,pwl,tagd)
        
        % ********* Function Description *********
        %
        % Calculate time-domain receiver functions
        % for Parent(P) and Daughter(D) components
        %
        %
        % ****************************************
        % *                                      *
        % *  Written by David L. Abt - 5.08,4.09 *
        % *                                      *
        % *  Based on code written by Scott      *
        % *  French                              *
        % *                                      *
        % *  Email: David_Abt@brown.edu          *
        % *                                      *
        % ********* Function Description *********
        
        disp(['Pre-Calculating ' STA ' (' NET ') ' Phase ' Time RFs...'])
        % Compute our impulse function(I)
        I_HWt           = 1;                % Impulse half-width time (sec)
        I_W             = 2*round(I_HWt/dt)+1;
        i_IW            = (I_W-1)/2;
        i_I0            = round(I_HWt/dt)+1;
        I               = zeros(1,I_W);
        I(i_I0)         = 1;                % Delta function at t=0
        I_CF            = RFPs.hf;          % Impulse upper corner frequency (Hz)
        f_nyq           = 0.5/dt;
        [b_I,a_I]       = butter(2,I_CF/f_nyq,'low');
        I               = filtfilt(b_I,a_I,I);
        
        n2 = length(P(1,:));
        
        % Time-bandwidth product and number of tapers in the multitaper
        TB = 4; NT = 7;
        
        % Create a cosine taper for the edges of our records
        T_Wt            = 20;               % Taper width (sec)
        T_W         	= ceil(T_Wt/dt);
        T             	= ones(1,n2);
        T(1:T_W+1)    	= (1-cos((0:T_W)*pi/T_W))/2;
        T(end-T_W:n2)	= (1+cos((0:T_W)*pi/T_W))/2;
        if strcmp(Phase,'Ps')
            Mask        = (n2+1):(2*n2-1);
            t0          = -5;
            t1          = dt*(n2-1);
        elseif strcmp(Phase,'Sp')
            Mask        = 1:(n2-1);
            t0          = -dt*(n2-1);
            t1          = 25;
        end
        t_RF0           = -dt*(n2-1);
        t_RF1           = dt*(n2-1);
        time            = t_RF0:dt:t_RF1;
        
        i_RF0           = find(abs(time)==min(abs(time)));
        
        if str2double(Event_IDs{end}(4:7))>2050; tag='synth'; else tag='data';end
        PreCalculated_RFs.(Phase).Time = t0:dt:t1;
        for ie2=1:length(Event_IDs)
            EQ_ID       = Event_IDs{ie2};	% Event ID
            P_T         = P(ie2,:).*T;   	% Parent tapered
            D_T         = D(ie2,:).*T;    	% Daughter tapered
            
            %         % Added by Ved Aug 11, 2010
            if(isempty(find(isnan(P_T), 1)) && isempty(find(isnan(D_T), 1)))
                %[RF_Time] = ETMTM(P_T,D_T,time,TB,NT,t0,t1,Phase,pwl,tag); % VED AUG 11, 2012
                if strcmp(Phase,'Ps')
                    if strcmp(tagd,'HK')
                        [RF_Time,~] = IDRF(Phase,P_T,D_T,dt,t0,t0,0.2,0.01,0.005,50);
                    else
                        %[RF_Time,~] = IDRF(Phase,P_T,D_T,dt,t0,t0,1,0.01,0.005,50);
                        [RF_Time,~] = IDRF(Phase,P_T,D_T,dt,t0,t0,1,0.005,0.0025,100); % 2x longer
                    end
                else
                    if strcmp(tagd,'HK')
                        [RF_Time,~] = IDRF(Phase,P_T,D_T,dt,t1,t1,0.2,0.01,0.005,50);
                    else
                        %[RF_Time,~] = IDRF(Phase,P_T,D_T,dt,t1,t1,1,0.01,0.005,50); %orig
                        [RF_Time,~] = IDRF(Phase,P_T,D_T,dt,t1,t1,1,0.005,0.0025,100); % 2x longer
                    end
                end
            else
                RF_Time = zeros(1,length(find(time>=t0 & time<=t1)));
            end
            if ~isempty(RF_Time)
                PreCalculated_RFs.(Phase).Individual_RFs.(EQ_ID) = single(RF_Time);
            end
            
        end
        fprintf('         ...done!\n')
        
        % Iterative time-domain deconvolution sub-function
        function[RF_Time] = ETMTM(P,D,time,TB,NT,t0,t1,Faza,pwl,tag)
            
            % Findings: when P arrival is not in the center of the window, the
            % amplitudes are not unity at the beginning and decreasing from there on.
            % Instead they peak at the time shift which corresponds to the middle index
            % in the P time window.
            
            % As your TB
            % increases, the frequency smearing gets worse, which means that the RFs
            % degrate at shorter and shorter lag times. Therefore, as you increase TB,
            % you should also increase Poverlap.
            
            %TB = 4; NT = 7; %choise of TB = 4, NT = 3 is supposed to be optimal
            %t0 = -5; t1 = max(time);
            %        function [RF_Time] = MTMDecon_VedOptimized(P,D,TB,NT,t0,t1,Faza)
            % Ved wrote MTM for MATLAB, which has the added advantage of
            % finding the optimal damping parameter.
            % TB  = time bandwidth product (usually between 2 and 4)
            % NT  = number of tapers to use, has to be <= 2*TB-1
            
            % Flip time axis in case of Sp
            if(strcmp(Faza,'Sp'))
                D = fliplr(D); P = fliplr(P); tmp0 = t0; tmp1 = t1;
                t0 = -tmp1; t1 = -tmp0;
            else
                fakt = 1;
            end
            
            % Length of moving time window in seconds
            win_len = 50;
            Nwin = round(win_len/dt);
            
            % Fraction of overlap overlap between moving time windows. As your TB
            % increases, the frequency smearing gets worse, which means that the RFs
            % degrate at shorter and shorter lag times. Therefore, as you increase TB,
            % you should also increase Poverlap.
            Poverlap = 0.90; nso=(1-Poverlap)*Nwin; nso=round(nso);
            
            
            
            npad=zeros(1,nso*1); D=[npad D npad]; P=[P npad npad];
            time = 0:dt:dt*(length(P)-1);
            
            % Create moving time windows and daughter/parent snippets
            starts = 1:round((1-Poverlap)*Nwin):length(P)-Nwin+1; nd=0;
            for j = 1:length(starts)
                tmp_times(j,1:Nwin) = time(starts(j):starts(j)+Nwin-1);
                if(j==1) % ASSUME THAT PARENT PHASE IS CENTERED IN FIRST WINDOW!
                    Pwin = interp1(double(time),double(P),double(tmp_times(j,:)),'linear',0)';
                    
                end
                Dwin(1:Nwin,j) = interp1(double(time),double(D),double(tmp_times(j,:)),'linear',0);
                
                ltp=win_len/5; % taper before deconvolving (important for synthetics)
                Dwin(1:Nwin,j)=taper(Dwin(1:Nwin,j)',tmp_times(j,1:Nwin),ltp,dt,...
                    tmp_times(j,1+round((ltp+dt)/dt)),tmp_times(j,Nwin-round((ltp+dt)/dt)));
                nd1=length(find(Dwin(:,j)>1e-2)); if nd1>nd; nd=nd1; end
            end
            
            % Search range for optimal damping parameter alpha
            switch tag
                case 'data'
                    alphas = logspace(-2,2,20)*var(D(round(end/4):3*round(end/4)))*length(P);
                case 'synth'
                    alphas = logspace(-2,2,20)*var(D)*length(P);
            end
            % Figure out average times for each moving window
            t0_times = median(tmp_times,2);
            
            % Construct Slepians
            [E,~] = dpss(length(Pwin),TB);
            
            % Length of waveforms;
            nh = length(Pwin);
            %
            
            misfit = zeros(size(alphas)); magntd = zeros(size(alphas));
            % Now, calculate misfit and RF size for each alpha
            for kj = 1:length(alphas)
                
                % Now loop through different time windows for the daughter component
                for k = 1:size(tmp_times,1)
                    % Create multitaper estimates
                    for j = 1:NT
                        tmp1 = fft(E(:,j).*Pwin,nh);
                        tmp2 = fft(E(:,j).*Dwin(:,k),nh);
                        if j==1
                            NUM = conj(tmp1).*tmp2;
                            DEN = conj(tmp1).*tmp1;
                        else
                            NUM = NUM + conj(tmp1).*tmp2;
                            DEN = DEN + conj(tmp1).*tmp1;
                        end
                    end
                    
                    % Calculate optimal RF
                    tmp = real(ifft(NUM./(DEN + alphas(kj))));
                    
                    % Filter and Normalize optimal RF
                    nrm = max(real(ifft(DEN./(DEN + alphas(kj)))));
                    
                    tmp = reshape(fftshift(tmp./nrm),[1 nh]);
                    
                    % Time vector
                    vrijeme = dt*[-0.5*(nh-1):1:0.5*(nh-1)]+t0_times(k)-t0_times(1) ...
                        -dt.*length(npad);
                    
                    
                    RF_Time_win(:,k) = interp1(double(vrijeme),double(tmp),double(t0:dt:t1),'linear',NaN);
                    length_Dwin=time(Nwin)-time(1);
                end
                
                
                tmp = conv(nanmean(RF_Time_win,2),P); t0d=round(t0/dt);
                if size(tmp,1)<size(tmp,2); mfD=D; else mfD=D'; end
                misfit(kj) = nansum(abs(mfD - tmp(1-t0d:length(D)-t0d)));
                magntd(kj) = nansum(abs(tmp(1-t0d:length(D)-t0d)));
            end
            
            % Find optimal alpha
            [~,j2] = min((misfit./std(misfit)).^2+(magntd./std(magntd)).^2);
            
            % Now loop through different time windows for the daughter component
            for k = 1:size(tmp_times,1)
                % Create multitaper estimates
                for j = 1:NT
                    tmp1 = fft(E(:,j).*Pwin,nh);
                    tmp2 = fft(E(:,j).*Dwin(:,k),nh);
                    if j==1
                        NUM = conj(tmp1).*tmp2;
                        DEN = conj(tmp1).*tmp1;
                    else
                        NUM = NUM + conj(tmp1).*tmp2;
                        DEN = DEN + conj(tmp1).*tmp1;
                    end
                end
                
                % Find optimal alpha
                % [junk,j2] = min(abs(misfit./std(misfit))+abs(magntd./std(magntd)));
                
                % Calculate optimal RF
                tmp = real(ifft(NUM./(DEN + alphas(j2))));
                
                % Filter and Normalize optimal RF
                nrm = max(real(ifft(DEN./(DEN + alphas(j2)))));
                
                tmp = reshape(fftshift(tmp./nrm),[1 nh]);
                
                % Time vector
                vrijeme = dt*[-0.5*(nh-1):1:0.5*(nh-1)]+t0_times(k)-t0_times(1) ...
                    -dt.*length(npad);
                
                RF_Time_win(:,k) = interp1(double(vrijeme),double(tmp),double(t0:dt:t1),'linear',0);
            end
            
            RF_Time = nanmean(RF_Time_win,2);
            
            if(strcmp(Faza,'Sp')), RF_Time = flipud(RF_Time); end
            
        end
        
        function [RF, RF_Time] = IDRF(Phase,P,D,dt,t_for,t_trun,gauss_t,accept_mis1,accept_mis2,itmax)
            % Iterative Deconvolution and Receiver-Function Estimation in time domain
            
            
            misfit = 1;
            misfit_old = 9999999999999;
            misfit_ref = sqrt(sum(D.^2));
            
            %RF = zeros(length(P)*2-1,1);
            RF_tmp = zeros(length(P)*2-1,1);
            
            D_cur = D;
            
            itnum = 0;
            
            [~,t_corr] = xcorr(D_cur,P);
            RF_Time = t_corr*dt;
            
            while (misfit_old-misfit>accept_mis1*misfit || misfit>accept_mis2*misfit_ref) && itnum <= itmax
                [amp_corr,~] = xcorr(D_cur,P);
                auto_corr = xcorr(P);
                [~,ind] = max(abs(amp_corr));
                amp_rf = amp_corr(ind)/auto_corr((length(t_corr)+1)/2);
                %
                %                 if RF_Time(ind) > -9.5 && (amp_rf>0 || abs(amp_rf)<0.02)
                %                     RF_tmp(ind) = RF_tmp(ind)+amp_rf;
                %                 else
                %                     RF_tmp(ind) = RF_tmp(ind)+amp_rf;
                %                     RF(ind) = RF(ind)+amp_rf;
                %                 end
                RF_tmp(ind) = RF_tmp(ind)+amp_rf;
                D_sub = conv(P,RF_tmp,'same');
                D_cur = D - D_sub;
                %plot(D_cur)
                misfit_old = misfit;
                misfit = sqrt(sum(D_cur.^2))/misfit_ref;
                itnum = itnum+1;
            end
            
            RF = RF_tmp;
            
            if strcmp(Phase,'Sp')
                RF(RF_Time>t_trun)=0;
                RF = RF(RF_Time<=t_for);
                RF_Time = RF_Time(RF_Time<=t_for);
            else
                RF(RF_Time<t_trun)=0;
                RF = RF(RF_Time>=t_for);
                RF_Time = RF_Time(RF_Time>=t_for);
            end
            
            if gauss_t~=0
                gauss_sig = gauss_t/dt;
                x = linspace(-gauss_sig*4,gauss_sig*4,gauss_sig*8);
                Gauss_win = exp(-x.^2/(2*gauss_sig^2));
                RF = conv(RF,Gauss_win,'same');
                % gauss_len = length(RF);
                % Gauss_win = gausswin(gauss_len,gauss_sig*4);
                %                 wl_T = 1/centfrq('mexh');
                %                 times_T = gauss_t/wl_T;
                %                 dt_wl = dt/times_T;
                %                 bd_wl = round(50/dt_wl)*dt_wl;
                %                 N_wl = 2*round(50/dt_wl)+1;
                %                 WL_win = mexihat(-bd_wl/2,bd_wl/2,N_wl);
                %
                %                 n_fft_gau = 2^nextpow2(length(WL_win));
                %                 pad = ceil((n_fft_gau-length(WL_win))/2);
                %                 gau_tmp = [zeros(pad-1,1);WL_win';zeros(pad,1)];
                %                 gau_fft = fft(hilbert(gau_tmp));
                %                 gau_level = gau_fft(1);
                %                 gau_fft = fftshift(gau_fft);
                %                 gau_F_fft = 1/dt*((-n_fft_gau/2):(n_fft_gau/2-1))/n_fft_gau;
                %                 %gau_fft = gau_fft./((1i*2*pi*gau_F_fft').^(1));
                %                 gau_fft = -gau_fft.*((1i*2*pi*gau_F_fft').^(3/2)).*1i;
                %
                %                 gau_fft = ifftshift(gau_fft);
                %                 gau_fft(1) = gau_level;
                %                 %RF_fft(1) = 0;
                %                 gau_shift = ifft(gau_fft);
                %                 if length(gau_shift) ~= length(WL_win)
                %                     gau_shift(end-pad+1:end) = [];
                %                     if pad>1
                %                         gau_shift(1:pad-1) = [];
                %                     end
                %
                %                     %gau_shift(length(Gauss_win)+1:end) = [];
                %                 end
                %                 gau_shift = real(gau_shift');
                %                 %     gau_shift = interp1(x,gau_shift,x_dt);
                %                 %     Gauss_win = interp1(x,Gauss_win,x_dt);
                %
                %                 %     gauss_sig = gauss_t/dt;
                %                 %     x = linspace(-gauss_sig*4,gauss_sig*4,gauss_sig*8);
                %                 %     Gauss_win = exp(-x.^2/(2*gauss_sig^2));
                %
                %                 %RF = conv(RF,Gauss_win,'same');
                %
                %                 if strcmp(tag,'ori')
                %                     RF = conv(RF,WL_win,'same');
                %                 elseif strcmp(tag,'shift')
                %                     RF = conv(RF,gau_shift,'same');
                %                 end
            end
            %             if strcmp(Phase,'Sp')
            %                 figure(1)
            %                 plot(RF_Time,RF,'k','linewidth',0.3);
            %                 hold on
            %             end
            %RF = flipud(RF);
            %RF_Time = fliplr(RF_Time);
        end
        
        
    end

    function [y] = taper(x,time,tpr,dt,t1,t2)
        
        % ********* Function Description *********
        %
        % TAPER  Taper a time series.
        %
        % TAPER(X,TIME,TPR,DT,T1,T2) takes time
        % series sampled at DT and tapers it with
        % a cosine taper TPR seconds long from
        % beginning point T1-TPR and with reverse
        % cosine taper from point T2 to point T2+
        % TPR. Points outside the range (T1-TPR,
        % T2+TPR) are zeroed. If T1/T2 is negative
        % then taper is not implemented at the
        % beginning/end. If X is an array of
        % seismograms, then the taper is applied
        % to each row of X.
        %
        %
        % ****************************************
        % *                                      *
        % *  Modified from Kate Rychert's        *
        % *  receiver function code - May 2008   *
        % *                                      *
        % *  Email: David_Abt@brown.edu          *
        % *                                      *
        % ****************************************
        
        nn      = length(x(1,:));
        nx      = length(x(:,1));
        taper   = ones(1,nn);
        it  	= [0:fix(tpr/dt)]*dt/tpr;
        ct      = 0.5-0.5*cos(pi*it);
        T1      = fix(time(1)/dt);             % Absolute sample point of first time step
        it1     = fix(t1/dt+1)-T1;
        it2     = fix(t2/dt+1)-T1;
        
        % Emily, 6th May 2013: problem with start time of phase being <100s, so
        % taper subscripts are negative.
        % Temporary workaround, set taper length to be shorter.
        % N.B. Only one event so far has had this issue!
        
        if t1>0
            if it1>fix(tpr/dt)
                taper(it1-fix(tpr/dt):it1)	= ct;
                taper(1:it1-fix(tpr/dt))    = zeros(size(1:it1-fix(tpr/dt)));
            else
                taper(1:it1) = ct(fix(tpr/dt)-it1+2:end);
                taper(1) = 0;
                disp('Bizarre taper!')
                
            end
        end
        if t2>0
            if t2>time(end)-tpr
                t2  = time(end)-tpr;
                it2	= fix(t2/dt)-T1;
            end
            taper(it2:it2+fix(tpr/dt))	= fliplr(ct);
            taper(it2+fix(tpr/dt):nn)	= zeros(size(taper(it2+fix(tpr/dt):nn)));
        end
        
        y = zeros(nx,nn);
        for ix=1:nx
            y(ix,:) = x(ix,:).*taper;
        end
        
    end

end

function [y] = bpfilt(x,dt,lf,hf,tag)

% ********* Function Description *********
%
% Bandpass filter a time seriers.
%
% [Y] = bpfilt(X,DT,LF,HF)
%
% Take a time series, X, sampled at DT and
% filter it with a 2nd order, 2 pass
% butterworth filter between frequencies
% LF and HF. If X is a matrix, this will
% filter the individual rows of X.
%
% ****************************************
% *                                      *
% *  Written by David L. Abt - May 2008	 *
% *                                      *
% *  Taken from code written by Michael  *
% *  Bostock and Stephane Rondenay, and  *
% *  used by Kate Rychert.               *
% *                                      *
% *  Email: David_Abt@brown.edu          *
% *                                      *
% ****************************************

if strcmp(tag,'bp')
    
    nyq     = 0.5/dt;           % Nyquist Frequency
    wn      = [lf/nyq,hf/nyq];
    [b,a]   = butter(2,wn);
    for ix=1:length(x(:,1))
        y(ix,:) = filtfilt(b,a,double(x(ix,:)));  % Edited by Ved because waveforms
        % are single not double, but filtfilt
        % requires double
    end
elseif strcmp(tag,'lp')
    nyq     = 0.5/dt;           % Nyquist Frequency
    wn      = [hf/nyq];
    [b,a]   = butter(2,wn,'low');
    for ix=1:length(x(:,1))
        y(ix,:) = filtfilt(b,a,double(x(ix,:)));  % Edited by Ved because waveforms
        % are single not double, but filtfilt
        % requires double
    end
end

end

function [dZs,rfs_atdepth] = Migrate_single_sta_RFs(PreCalculated_RFs, ...
    Migration_Data, All_Ray_Path_Data, Phase, evs,lowT,isRFdiff)

%Delay gives the delay time between the recorded wave and a ray coming in
%from 65 deg away (for an Sp path): DT - reference DT
% For Sp, all of these delay times are negative
% So +ve: wave travels faster than reference DT at that depth (small RP)
%    -ve: wave travels slower than reference DT at that depth (large RP)


pathstr=[Phase(1) '_Path']; Mask_Dp=zeros(1,length(evs));
time=PreCalculated_RFs.(Phase).Time;

if strcmp(isRFdiff,'yes')
    % Reference
    dZs   	= Migration_Data.Reference_Depth;
    SrcZs   = Migration_Data.Reference_SrcDepth;
    Dists   = Migration_Data.Reference_Distance;
    Delay   = Migration_Data.Reference_DT;
    Er      = Migration_Data.Er;
    
    [SRC,DIST,DZ] = meshgrid(SrcZs,Dists,dZs);
    
    
    RFs=zeros(length(evs),length(time));
    dist=zeros(1,length(evs));
    src =zeros(1,length(evs));
    
    % Data
    for k=1:length(evs)
        RFs(k,:)=PreCalculated_RFs.(Phase).Individual_RFs.(evs{k});
        dist(k) = All_Ray_Path_Data.(evs{k}).(pathstr).distance;
        src(k) = All_Ray_Path_Data.(evs{k}).(pathstr).srcDepth;
        RP = interp3(SRC,DIST,DZ,Migration_Data.Reference_RP,src(k)*ones(size(dZs)),dist(k)*ones(size(dZs)),dZs);
        mdp=dZs(find(Migration_Data.Critical_RPs<=RP(1:end-1)/Er,1));
        if isempty(mdp); Mask_Dp(k)=dZs(end); else Mask_Dp(k)=mdp; end
    end
    %
else
    % Reference
    dZs   	= Migration_Data.Reference_Depth;
    RPs     = Migration_Data.Relative_RPs;
    Delay   = Migration_Data.Relative_DTs';
    
    RFs=zeros(length(evs),length(time)); RP=zeros(1,length(evs));
    % Data
    for k=1:length(evs)
        RFs(k,:)=PreCalculated_RFs.(Phase).Individual_RFs.(evs{k});
        RP(k)=All_Ray_Path_Data.(evs{k}).(pathstr).ray_parameter;
        mdp=dZs(find(Migration_Data.Critical_RPs<=RP(k),1));
        if isempty(mdp); Mask_Dp(k)=dZs(end); else Mask_Dp(k)=mdp; end
    end
    
end

%
rfs_atdepth=zeros(length(evs),length(dZs));

if strcmp(isRFdiff,'yes')
    for k=1:length(evs)
        StaDelay=interp3(SRC,DIST,DZ,Delay,src(k)*ones(size(dZs)),dist(k)*ones(size(dZs)),dZs);
        
        if lowT==0; filtrfs=RFs(k,:);
        else filtrfs=bpfilt(RFs(k,:),diff(time(1:2)),0.01,1/lowT,'lp');
        end
        %filtrfs=RFs(k,:); %Junlin
        
        rfs_atdepth(k,:)=interp1(time,filtrfs,StaDelay);
        rfs_atdepth(k,dZs>Mask_Dp(k))=nan;
    end
    
    if strcmp(Phase,'Sp'); rfs_atdepth=-rfs_atdepth;  end
else
    for k=1:length(evs)
        [rps,dzs] = meshgrid(RPs,dZs);
        dzs=double(dzs);
        dZs=double(dZs);
        StaDelay=interp2(rps,dzs,Delay',RP(k),dZs);
        if lowT==0; filtrfs=RFs(k,:);
        else filtrfs=bpfilt(RFs(k,:),diff(time(1:2)),0.01,1/lowT,'lp');
        end
        rfs_atdepth(k,:)=interp1(time,filtrfs,StaDelay+Migration_Data.Reference_DT);
        rfs_atdepth(k,dZs>Mask_Dp(k))=nan;
    end
    
    if strcmp(Phase,'Sp'); rfs_atdepth=-rfs_atdepth;  end
end




%meanrf=nanmedian(rfs_atdepth,1);



end

function [RF_min, RF_max, RF_med, depth] = CompRefRF( Phase, Params, Networks, Proj_Dir,basedir, VpModelname , VsModelname)
%Compute reference receiver functions based on velocity models, to help remove
%ones far different with reference ones

%First set the ray parameter range of the events
disp('Setting ray parameter range')
NETs = Params.NETs;
rp_max = 0;
rp_min = 1;
rp_sum = [];
rp_n = 0;
for in = 1:length(NETs)
    STAs=fieldnames(Networks.(NETs{in}));
    for is=1:length(STAs)
        if isfield(Params.unprepped,NETs{in})
            if isfield(Params.unprepped.(NETs{in}),STAs{is})
                continue;
            end
        end
        if exist([Proj_Dir NETs{in} '/' STAs{is} '/Ray_Path_Data_',num2str(1/Params.hf),'_',num2str(1/Params.lf),'.mat'],'file')
            load([Proj_Dir NETs{in} '/' STAs{is} '/Ray_Path_Data_',num2str(1/Params.hf),'_',num2str(1/Params.lf),'.mat']);
            eves = fieldnames(All_Ray_Path_Data);
            rps = zeros(length(eves),1);
            for ie = 1:length(eves)
                if strcmp(Phase,'Sp')
                    rps(ie) =  All_Ray_Path_Data.(eves{ie}).S_Path.ray_parameter;
                else
                    rps(ie) =  All_Ray_Path_Data.(eves{ie}).P_Path.ray_parameter;
                end
            end
            rp_sum = [rp_sum;rps];
            rp_n = rp_n + length(rps);
            minrp_tmp = min(rps);
            maxrp_tmp = max(rps);
            if minrp_tmp<rp_min
                rp_min = minrp_tmp;
            end
            if maxrp_tmp>rp_max
                rp_max = maxrp_tmp;
            end
        end
    end
end
rp_med = median(rp_sum);

rp_ref = [rp_min,rp_med,rp_max];

disp('Calculating reference synthetics')
%getting synthetics
tmp_vfile_s = dir([basedir 'Data/Velocity_Models/Mantle/Vs/' VsModelname '.mat']);
tmp_vfile_p = dir([basedir 'Data/Velocity_Models/Mantle/Vp/' VpModelname '.mat']);
orimodel_s = load([tmp_vfile_s.folder '/' tmp_vfile_s.name]);
orimodel_p = load([tmp_vfile_p.folder '/' tmp_vfile_p.name]);

depth = 0:450;

[LonG,LatG,DepG] = meshgrid(orimodel_s.Vs_Model.Longitude,orimodel_s.Vs_Model.Latitude,orimodel_s.Vs_Model.Depth);
latmin = orimodel_s.Vs_Model.Latitude(1);
latmax = orimodel_s.Vs_Model.Latitude(end);
lonmin = orimodel_s.Vs_Model.Longitude(1);
lonmax = orimodel_s.Vs_Model.Longitude(end);
lat_model = latmin:1:latmax;
lon_model = lonmin:1:lonmax;
[LonM,LatM,DepM] = meshgrid(lon_model,lat_model,depth');
Vs_M = interp3(double(LonG),double(LatG),double(DepG),...
    double(orimodel_s.Vs_Model.Vs),double(LonM),double(LatM),double(DepM));
Vp_M = interp3(double(LonG),double(LatG),double(DepG),...
    double(orimodel_p.Vp_Model.Vp),double(LonM),double(LatM),double(DepM));
max_V = nanmax(nanmax(nanmax(Vp_M)));
max_rp = sind(85)/max_V;
ind_ex = find(rp_ref>max_rp);
if ~isempty(ind_ex)
    rp_ref(ind_ex) = [];
    rp_ref = [rp_ref,max_rp];
end

Period = 13;
recordid='RFcomp';

rf_model = zeros(length(rp_ref),length(lat_model),length(lon_model),length(depth));

cd([basedir,'NewFunctions/PROPMAT/'])
for ilat = 1:length(lat_model)
    for ilon = 1:length(lon_model)
        Vp_ray = permute(Vp_M(ilat,ilon,:),[2,3,1]);
        Vs_ray = permute(Vs_M(ilat,ilon,:),[2,3,1]);
        dV_ray = gradient(Vs_ray)./gradient(depth);
        dV_ray(isnan(dV_ray)) = 0;
        indv1 = find(~isnan(Vs_ray),1);
        indv2 = find(~isnan(Vs_ray),1,'last');
        if isempty(indv1)
            Vs_ray(isnan(Vs_ray)) = 3.2;
            Vp_ray(isnan(Vp_ray)) = 5.6;
            Moho_ray = 100000;
        else
            Vs_ray(1:indv1) = Vs_ray(indv1);
            Vp_ray(1:indv1) = Vp_ray(indv1);
            if indv2>length(Vs_ray)
                Vs_ray(indv2:end) = Vs_ray(indv2);
                Vp_ray(indv2:end) = Vp_ray(indv2);
            end
            Moho_ray = depth(dV_ray==max(dV_ray(depth>30 & depth<60)));
            if length(Moho_ray)>1
                Moho_ray = Moho_ray(1);
            end
        end
        vs = Vs_ray(1);
        vp = Vp_ray(1);
        rho = 2.8;
        zbot = [];
        ztop = 0;
        for idep = 1:length(depth)
            if abs(Vs_ray(idep)-vs(end))/vs(end)>0.005
                vs = [vs;Vs_ray(idep)];
                vp = [vp;Vp_ray(idep)];
                if depth(idep)<Moho_ray
                    rho = [rho;2.8];
                else
                    rho = [rho;3.3];
                end
                zbot = [zbot;depth(idep)];
                ztop = [ztop;depth(idep)];
            end
        end
        if ztop(end) ~= depth(end)
            zbot = [zbot;depth(end)];
        else
            zbot = [zbot;depth(end)+1];
        end
        for irp = 1:length(rp_ref)
            write_propmat_syn(ztop,zbot,vp,vs,rho,rp_ref(irp),Period,[recordid,'.mat'],Phase);
            SynWave = load([basedir,'NewFunctions/PROPMAT/',recordid,'.mat']);
            if isnan(SynWave.R(100))
                rf_model(irp,ilat,ilon,:) = nan;
                continue;
            end
            [SynRF,SynTime] = getrf(SynWave.R,SynWave.Z,SynWave.dt,vp(1),vs(1),rp_ref(irp));
            [SynDepth,Synindts,Synindte] = getdepth(Phase,SynTime,zbot,vp,vs,rp_ref(irp));
            SynRF_dep = interp1(SynDepth,SynRF(Synindts:Synindte),depth);
            SynRF_dep(isnan(SynRF_dep)) = 0;
            rf_model(irp,ilat,ilon,:) = permute(SynRF_dep,[1,4,2,3]);
        end
    end
end
cd(basedir)
RF_med = permute(nanmedian(nanmedian(rf_model(2,:,:,:),2),3),[4,3,2,1]);
RF_min = permute(nanmin(nanmin(nanmin(rf_model,[],1),[],2),[],3),[4,3,2,1]);
RF_max = permute(nanmax(nanmax(nanmax(rf_model,[],1),[],2),[],3),[4,3,2,1]);

end

function diffRP(basedir,Phases)
% Find the corresponding ray parameter at different depths and distances
% for P and S

import edu.sc.seis.TauP.*	% Import the TauP package
for j = 1:length(Phases)
    if strcmp(Phases{j},'Ps')
        distances = 30:1:85;
    else
        distances = 50:1:90;
    end
    depths = 0:1000;
    srcdepths = 0:20:1000;
    RayParas = zeros(length(distances),length(srcdepths));
    DISTs = zeros(length(distances),length(srcdepths),length(depths));
    for i = 1:length(distances)
        for k = 1:length(srcdepths)
            Temp = Matlab_TauP('Path','ak135',srcdepths(k),Phases{j}(1),'sta',[0 0],'evt',[0 distances(i)]);
            [maxdep,imax] = max(Temp.path.depth);
            inds = find(depths<maxdep);
            RayParas(i,k) = Temp.rayParam;
            DISTs(i,k,inds) = interp1(Temp.path.depth(imax:end),Temp.path.distance(imax:end),depths(inds));
        end
    end
    DISTs(DISTs==0) = nan;
    %
    % DDISTs = diff(DISTs/180*pi);
    % DDISTs = [zeros(1,length(depths));DDISTs];
    % DDISTs(1,isnan(DDISTs(2,:))) = nan;
    % int_RayPara = [RayParas(1);RayParas(1:end-1)+diff(RayParas)/2];
    % int_pdelta = cumsum(DDISTs.*repmat(int_RayPara,[1,length(depths)]));
    
    % RayPara_comp = 420:-1:265;
    % int_pdelta_comp = interp1(RayParas,int_pdelta(:,701),RayPara_comp);
    % plot(RayPara_comp,int_pdelta_comp-int_pdelta_comp(1));
    % hold on
    TaupInfo.(Phases{j}).distances = distances;
    TaupInfo.(Phases{j}).depths = depths;
    TaupInfo.(Phases{j}).srcdepths = srcdepths;
    TaupInfo.(Phases{j}).RayParam = RayParas;
    TaupInfo.(Phases{j}).DIST = DISTs;
    
    
end

save([basedir,'Data/Velocity_Models/TaupInfo.mat'],'TaupInfo');

end

