function CCP_stack_slopeangle_RPdiff(Project, basedir, Phases, tag,vmodel,varargin)

% ********* Function Description *********
%
% Migrate RFs into 3D volume and stack based
% on common convergence point
%
% ****************************************

clc
Proj_Dir        = [basedir 'Data/Projects/' Project '/'];
if strcmp(Project(1:6),'Synth_'); data_tag='synth'; version=varargin{1};
else data_tag='normal'; 
    if ~isempty(varargin) && strcmp(varargin{1}(1:5),'names')
        version=varargin{1};
    else version='names';
    end
end


switch tag
    case 'calc_deplatlon'
        if(length(varargin)>1); strike=varargin{2}; else strike=[]; end
        calc_deplatlon(basedir, Project, Phases,data_tag,version,strike,vmodel)
    case 'stack'
        CCP_stack_from_precalc_RF(basedir, Proj_Dir,Phases,version,data_tag)
    case 'bootstrap'
        rng shuffle
        for iph=1:length(Phases)
            Phase=Phases(iph);
            load([Proj_Dir version '.mat']);
            numbtstrp=varargin{2};
            if isfield(FileNames,'DomT'); T=FileNames.DomT; 
            else  T=13.5;   %calcDomPeriod(Proj_Dir, Phase{:}, FileNames.PreCalcRFs,version);
                FileNames.DomT=T;
                eval(['save ' Proj_Dir version '.mat FileNames']);
            end
            CCP_stack_from_precalc_RF(basedir, Proj_Dir,Phase,version,'bootstrap',numbtstrp,T)
        end
    case 'meanbtstrap'
        for iph=1:length(Phases)
            Phase=Phases{iph};
            calcmeanbootstrap(basedir, Proj_Dir,version)
        end
end



end


function calc_deplatlon(basedir, project, Phases, tag, version,strike,vmodel)

% ********* Function Description *********
%
% Migrate RFs in time into depth and associate each depth with a lat,lon
% based on shooting at the parent ray parameter along the backazimuth.
% Write all the RFs to a data structure that can be read and used by
% CCP stacking programs. - Ved Aug 13, 2010
% 
% Also calculate the corresponding distance, converted wave rayparameter
% relationship
%
% ****************************************

dir0 = [basedir 'Data/Projects/' project '/'];
Er = 6371; km_per_deg  = pi*Er/180; Ellipsoid   = [Er 1]/km_per_deg;
load([dir0 'Networks.mat']);

iffiltdat = 'yes';
if strcmp(iffiltdat,'yes')
    load([dir0 'Similarity_RF_450.mat']);
    cutoff = 5*median(Total_info.Dist_tot_2);
    cutoff_shallow_size = 0.2*median(Total_info.Self_shallow_Moho_size);
    cutoff_shallow_pv = 3*median(Total_info.Self_shallow_pv_size);
    cutoff_deep_num = median(Total_info.Num_deep_out_2);
end

if exist([dir0 version '.mat'],'file'); load([dir0 version '.mat']);
else disp('No PreCalc RFs name saved!'); return;
end


if isempty(strike)
    minbaz=-180; maxbaz=180;
else
    minbaz=strike-20; maxbaz=strike+20;
end


for iph=1:length(Phases)
    phase=Phases{iph};
    name = FileNames.(phase).PreCalcRFs;
    eles = strsplit(name,'_');
    lowT = str2double(eles{3});
    highT = str2double(eles{4}(1:end-1));
    
    network_codes=fieldnames(Networks);
    numtot = 0;

    jk = 0;
    for j = 1:length(network_codes)
        stations = fieldnames(Networks.(network_codes{j}));
        for k = 1:length(stations)
            jk = jk + 1;
            folder{jk} = strcat(dir0,network_codes{j},'/',stations{k},'/');
        end
    end
    
   
    minXcorr = -1; minS2N = 2; % Minimum allowable Xcorr and s2n
    
    typewf = strcat(phase,'_Waveforms');
    switch(lower(phase)); case 'ps'; ime = 'P_Path'; case 'sp'; ime = 'S_Path'; end
    nj = 0; n = 0;
    
    % Initialise matrices
   if strcmp(phase,'Sp')
        maxdep=450; deps=0:0.5:maxdep;    %change bin size
    else
        maxdep=780; deps=0:0.5:maxdep;
   end
   CP_lat=zeros(76000,length(deps)); CP_lon=CP_lat; RF_depth=CP_lat; time_difference=CP_lat; RP_conv = CP_lat;
   SLAT=zeros(length(folder),1); SLON=SLAT; station_at_cp=zeros(1,length(CP_lat(:,1)));
   baz=station_at_cp; src=station_at_cp; dist=station_at_cp; RP=station_at_cp;
   slats=station_at_cp; slons=station_at_cp; s_inds=station_at_cp;
   

    
    
    % Calculations begin here
    for jkjk = 1:length(folder)
        
        folderName=textscan(folder{jkjk}, '%s', 'delimiter','/');
        net=folderName{1}{end-1};  sta=folderName{1}{end};
        
        
        display(['Working on station number ' num2str(jkjk) '...']);
        
        clear depth RPs M depths delay RPs_grid depths_grid quakes time Delay Moho
        clear RFs mask_i Mask_Depth Relative_DT Mask_Time dt T Z  RP_mig
        clear sort_by_dist sort_by_baz yd yb nRFs_by_dist nRFs_by_baz
        
       if strcmp(tag,'synth'); 
           junk=strsplit(folder{jkjk}, 'Synth_'); migr_folder=[junk{1} junk{2}];
       else migr_folder=folder{jkjk};
       end
        
        
        if(exist([folder{jkjk} '/PreCalculated_RFs' name],'file') && ...
                exist([migr_folder '/Migration_Models_',num2str(lowT),'_',num2str(highT),'.mat'],'file') && ...
                exist([folder{jkjk} '/Ray_Path_Data_',num2str(lowT),'_',num2str(highT),'.mat'],'file') && ...
                exist([folder{jkjk} '/Station_Data_',num2str(lowT),'_',num2str(highT),'.mat'],'file') && ...
                exist([folder{jkjk} '/Waveform_Data_',num2str(lowT),'_',num2str(highT) '.mat'],'file'))
            
            load([folder{jkjk} 'Station_Data_',num2str(lowT),'_',num2str(highT),'.mat']);
            SLAT(jkjk) = Station_Data.Latitude; SLON(jkjk) = Station_Data.Longitude;
            SELE = Station_Data.Elevation;
            
            load([migr_folder 'Migration_Models_',num2str(lowT),'_',num2str(highT),'.mat']);
            which_model=length(Migration_Models.(phase));
            
            
            for jm=length(Migration_Models.(phase)):-1:1
                if(strcmp(Migration_Models.(phase)(jm).Original_Models.Mantle_Vp, ...
                        vmodel.Vp) && ...
                        strcmp(Migration_Models.(phase)(jm).Original_Models.Crust, ...
                        vmodel.Cr) && ...
                        strcmp(Migration_Models.(phase)(jm).Original_Models.Mantle_Vs,...
                        vmodel.Vs))
                    which_model=jm; break
                end
            end
            
            if ~exist('which_model','var')
                disp([sta ' is migrated incorrectly']); continue
            end
            v=Migration_Models.(phase)(which_model).Velocity_Model;
            depth = v.Depth; %Reference_Depth;
            if max(depth)<maxdep; depth=[depth; maxdep];
                v.Depth=depth; v.Vp=[v.Vp; v.Vp(end)]; v.Vs=[v.Vs; v.Vs(end)];
            end
 
            dZs   	= Migration_Models.(phase)(which_model).Reference_Depth;
            SrcZs   = Migration_Models.(phase)(which_model).Reference_SrcDepth;
            Dists   = Migration_Models.(phase)(which_model).Reference_Distance;
            Delay   = Migration_Models.(phase)(which_model).Reference_DT;
            RP_mig  = Migration_Models.(phase)(which_model).Reference_RP;
            Moho = Migration_Models.(phase)(which_model).Original_Models.Crustal_Model.Moho;
            
            [SRC,DIST,DZ] = meshgrid(SrcZs,Dists,dZs);
            
            depths = min(depth):0.5:max(depth);    %change bin size
            % Define subregion in depth that we are actually interested in
            depth_indx = find(depths>=0 & depths<=maxdep); %250 200
            cp_depths = depths(depth_indx);
 
            ind_dep_m = find(cp_depths >= 200);
            ind_dep_moho = find(cp_depths >= Moho-SELE-10 & cp_depths <= Moho-SELE+10);
            if strcmp(phase,'Ps')
                ind_410 = find(depth <= 470 & depth >= 350);
                ind_660 = find(depth <= 720 & depth >= 600);
            end
            
            Velocity_Model(jkjk) = Migration_Models.(phase)(which_model).Velocity_Model;
            display('Using Standard Velocity Model');
            
            load([folder{jkjk} strcat('/PreCalculated_RFs', name)]);
            load([folder{jkjk} 'Ray_Path_Data_',num2str(lowT),'_',num2str(highT),'.mat']);
            load([folder{jkjk} 'Waveform_Data_',num2str(lowT),'_',num2str(highT),'.mat']);
            
            duz = length(cp_depths);
            
            if strcmp(iffiltdat,'yes')
                quakes =  CompRFres.(net).(sta).eve;
                
                quakes(CompRFres.(net).(sta).dist_2>cutoff | ...
                    CompRFres.(net).(sta).RFsize_shallow_Moho<cutoff_shallow_size | ...
                    CompRFres.(net).(sta).N_deep_out_range_2>cutoff_deep_num | ...
                    CompRFres.(net).(sta).RFsize_shallow_pv>cutoff_shallow_pv ) = [];
            else
                quakes = fieldnames(PreCalculated_RFs.(phase).Individual_RFs);
            end
            
            numtot = numtot+length(quakes);
            
            if isempty(quakes)
                continue;
            end
            
            indx_same = find(diff(Velocity_Model(jkjk).Depth)==0) ;
            Velocity_Model(jkjk).Depth(indx_same+1) = Velocity_Model(jkjk).Depth(indx_same+1) + 0.01;
            
            vs = interp1(Velocity_Model(jkjk).Depth,Velocity_Model(jkjk).Vs,cp_depths,'linear');
            vp = interp1(Velocity_Model(jkjk).Depth,Velocity_Model(jkjk).Vp,cp_depths,'linear');
      
            
            Inter_Velo(jkjk).vs = vs;
            Inter_Velo(jkjk).vp = vp;
            
            
            for j = 1:length(quakes)
                
                if(~isempty(PreCalculated_RFs.(phase).Individual_RFs.(quakes{j})) && ...
                        isfield(All_Ray_Path_Data,(quakes{j})))
                    
                    bAz=All_Ray_Path_Data.(quakes{j}).(ime).bAzimuth; 
                    if bAz>180; BAZ=bAz-360; else BAZ=bAz; end
                    if (minbaz<-180 && bAz>0); BAZ=bAz-360;
                    elseif (maxbaz>180 && bAz<0); BAZ=bAz+360;
                    elseif bAz>180; BAZ=bAZ-360;
                    end
                    if((BAZ>minbaz && BAZ<maxbaz))
                       
                        if isfield(All_Waveform_Data.(quakes{j}).(typewf),'Window')
                            if isfield(All_Waveform_Data.(quakes{j}).(typewf).Window,'Taup_S2N_Misfit')
                                TauP_S2N = All_Waveform_Data.(quakes{j}).(typewf).Window.Taup_S2N_Misfit;
                            else
                                TauP_S2N = 0;
                            end
                        else
                            TauP_S2N = 0;
                        end
                        
                        if(All_Waveform_Data.(quakes{j}).(typewf).ZR_Xcorr_Coeff > minXcorr ...
                                && All_Waveform_Data.(quakes{j}).(typewf).S2N_ratio > minS2N ...
                                && TauP_S2N < 100)
                            n = n + 1;
                            
                            station_at_cp(n) = jkjk;
                            RF_in_time = double(PreCalculated_RFs.(phase).Individual_RFs.(quakes{j}));
                            
%                             if(max(abs(RF_in_time))>0)   %Junlin
%                                 % Normalize
%                                 RF_in_time = RF_in_time./max(abs(RF_in_time));
%                             end
                            
                            
                            baz(n)  = All_Ray_Path_Data.(quakes{j}).(ime).BAZ;
                            src(n)  = All_Ray_Path_Data.(quakes{j}).(ime).srcDepth;
                            RP(n)   = All_Ray_Path_Data.(quakes{j}).(ime).ray_parameter;
                            dist(n) = All_Ray_Path_Data.(quakes{j}).(ime).distance;
                            s_inds(n) = jkjk;
                            
                            RP_conv(n,:) = interp3(SRC,DIST,DZ,RP_mig,src(n)*ones(size(cp_depths)),dist(n)*ones(size(cp_depths)),cp_depths);
                            
                            
                            CP_lat(n,1) = SLAT(jkjk); CP_lon(n,1) = SLON(jkjk); RF_depth(n,1) = 0;
                             
                            
                            clear D
                            D(1) = 0;
                            if strcmp(phase,'Sp')
                                [Dtmp,~,~] = jhua_shootray(vp',cp_depths',RP_conv(n,2:duz)/Er,Er);
                                D(2:duz) = diag(Dtmp);
                            else
                                [Dtmp,~,~] = jhua_shootray(vs',cp_depths',RP_conv(n,2:duz)/Er,Er);
                                D(2:duz) = diag(Dtmp);
                            end
                            
                            t = interp3(SRC,DIST,DZ-SELE,Delay,src(n)*ones(size(cp_depths)),dist(n)*ones(size(cp_depths)),cp_depths);
                            
                            can_use = ~isnan(t);
                            cant_use = isnan(t);
                            
                            rf_depth = interp1(PreCalculated_RFs.(phase).Time(1:length(RF_in_time)), RF_in_time, t,'pchip');
                            
                            rf_m = rf_depth(ind_dep_m);
                            rf_moho = rf_depth(ind_dep_moho);
                            
                            if strcmp(phase,'Ps')
                                rf_410 = rf_depth(ind_410);
                                rf_660 = rf_depth(ind_660);
                                max_410 = nanmax(rf_410);
                                max_660 = nanmax(rf_660);
                            end
                            
                            max_moho = nanmax(rf_moho);
                            n_nnan_m = nansum(~isnan(rf_m));
                            
                            rms_m = sqrt(nansum(rf_m.^2)/n_nnan_m);
                            %n_m = nansum(abs(rf_m)>0.07,2)/n_nnan_m;
                            n_m = nansum(abs(rf_m)>0.1,2)/n_nnan_m;
                            
                            if strcmp(phase,'Ps')
                                %if max_moho<=0.03 || rms_m>=0.1 || n_m>=0.4
                                %if max_moho<=0 || rms_m>=0.15 || n_m>=0.4
                                if max_moho<=0.02 || max_410<=0.02 || max_660<=0.02 || rms_m>=0.15 || n_m>=0.4
                                    RP_conv(n,:) = zeros(size(t));
                                    station_at_cp(n) = 0;
                                    baz(n)  = 0;
                                    src(n)  = 0;
                                    RP(n)   = 0;
                                    dist(n) = 0;
                                    n = n-1;
                                    continue;
                                    
                                end
                            else
                                if max_moho<=0.02 || rms_m>=0.1 || n_m>=0.3
                                    RP_conv(n,:) = zeros(size(t));
                                    station_at_cp(n) = 0;
                                    baz(n)  = 0;
                                    src(n)  = 0;
                                    RP(n)   = 0;
                                    dist(n) = 0;
                                    n = n-1;
                                    continue;
                                    
                                end
                            end
                                
                            time_difference(n,:)=t;
                            
                            % shootray checks for critical angle;
                            % will set t = inf if have passed it
%                             can_use = ~isinf(t);
%                             cant_use = isinf(t);
                            can_use = ~isnan(t);
                            cant_use = isnan(t);
                            
                            CP_lat(n,1) = SLAT(jkjk); CP_lon(n,1) = SLON(jkjk); RF_depth(n,1) = 0;
                            
                            [CP_lat(n,can_use), CP_lon(n,can_use)] = reckon('gc',SLAT(jkjk),SLON(jkjk),D(can_use),baz(n),'degrees');
                            RF_depth(n,can_use) = rf_depth(can_use);
                            
                            RF_depth(n,cant_use) = 0; CP_lat(n,cant_use) = -999;
                            CP_lon(n,cant_use) = -999;
                            
                        end
                    end
                end
            end
 
        end
    end
    
 
    name=strcat(project, name(1:end-4), '_Slope_', num2str(maxdep),'km_', num2str(strike), 'deg.mat');
    
    RP(n+1:end)=[];  baz(n+1:end)=[]; s_inds(n+1:end)=[]; slats(n+1:end)=[]; slons(n+1:end)=[];
    CP_lat(n+1:end,:)=[]; CP_lon(n+1:end,:)=[]; RF_depth(n+1:end,:)=[]; RP_conv(n+1:end,:)=[];
    time_difference(n+1:end,:)=[]; station_at_cp(n+1:end)=[]; src(n+1:end)=[]; dist(n+1:end)=[];
    
    
    if ~exist([basedir '/Data/Projects/',project,'/CCP'],'dir')
        mkdir([basedir '/Data/Projects/',project,'/CCP']); 
    end
    
%     eval(['save ' basedir '/Data/CCP/AllNtwrks_' phase '_' name ...
%         ' RF_depth CP_lat CP_lon RP baz SLAT SLON cp_depths depth_indx' ...
%         ' phase Velocity_Model s_inds slons slats folder station_at_cp time_difference name project Sca_info']) %depths  RP
    FileNames.(phase).CCP_Slope=name; save([dir0 version '.mat'],'FileNames');
    save([basedir '/Data/Projects/',project,'/CCP/AllNtwrks_' phase '_' name],'RF_depth','CP_lat','CP_lon','RP','baz','src','dist','SLAT','SLON','cp_depths','depth_indx',...
        'phase','Velocity_Model','s_inds','RP_conv','folder','station_at_cp','time_difference','name','project','Inter_Velo','numtot','-v7.3')
    
end

end

function CCP_stack_from_precalc_RF(basedir, Proj_Dir,Phases,version,tag,varargin)

load([Proj_Dir version '.mat']);

for iph=1:length(Phases)
    Phase=Phases{iph};
    if isfield(FileNames.(Phase),'CCP_Slope'); name=FileNames.(Phase).CCP_Slope;
    else disp('No dep/lat/lon matrix name saved!'); return
    end
    
    %load([Proj_Dir '/CCP/AllNtwrks_' Phase '_' name]);
    aaa=load([Proj_Dir '/CCP/AllNtwrks_' Phase '_' name]);
    cp_depths = aaa.cp_depths;
    depth_indx = aaa.depth_indx;
    CP_lat = aaa.CP_lat;
    CP_lon = aaa.CP_lon;
    SLAT = aaa.SLAT;
    SLON = aaa.SLON;
    phase = aaa.phase;
    s_inds = aaa.s_inds;
    %RP = aaa.RP;
    RP_conv = aaa.RP_conv;
    baz = aaa.baz;
    Inter_Velo = aaa.Inter_Velo;
    time_difference = aaa.time_difference;
    RF_depth = aaa.RF_depth;
    
    clear aaa
    
    switch tag
        case 'normal'; rfi=1; 
            if strcmp(Phase,'Sp')
                period = 10; %period=calcDomPeriod(Proj_Dir, Phase, FileNames.PreCalcRFs,version); %Junlin
            else
                period = 6;
            end
        case 'bootstrap'; rfi=varargin{1}; period=varargin{2};
        case 'synth'; rfi=1;
            switch Phase; 
                case 'Sp'; 
                    junk=strsplit(FileNames.PreCalcRFs,'_Sp_'); lowT=str2double(junk{2}(1));
                    if lowT==1; period=6; elseif lowT==2; period=7; 
                    elseif lowT==4; period=9; else period=9;
                    end
                case 'Ps'; period=5;
            end
    end

    FileNames.(Phase).stack_slope=['CCP_' Phase '_' num2str(period) '_' FileNames.(Phase).CCP_Slope(1:end-4) '_'];
    save([Proj_Dir version '.mat'],'FileNames');
    
    CP_lat = transpose(CP_lat); CP_lon = transpose(CP_lon); map_size = 4;
    LatLims = [min(SLAT(find(SLAT)))-map_size  max(SLAT(find(SLAT)))+map_size];
    LonLims = [min(SLON(find(SLON)))-map_size  max(SLON(find(SLON)))+map_size];
    LonLims_Temp = LonLims;
    LonLims_Temp(LonLims_Temp>180)  = LonLims(LonLims_Temp>180)-360;
    [tmplat,tmplon] = meshgrid(LatLims(1):0.1:LatLims(2), LonLims_Temp(1):0.1:LonLims_Temp(2));
    model_lats = tmplat(:); model_lons = tmplon(:);
    
    
%     for ijk = 1:length(Sca_info)
%         if ~isempty(Sca_info(ijk))
%             slat = SLAT(ijk);
%             slon = SLON(ijk);
%             Dis_model = 2*asind(sqrt(sind((model_lats-slat)/2).^2+...
%                 cosd(slat)*cosd(model_lats).*(sind((model_lons-slon)/2).^2)));
%             baz_model = atan2d((sind(model_lons-slon).*cosd(model_lats)),...
%                 (cosd(slat)*sind(model_lats)-sind(slat)*cosd(model_lats).*cosd(model_lons-slon)));
%             baz_model(baz_model<0) = baz_model(baz_model<0) + 360;
%             F = griddedInterpolant(Sca_info(2).Sp.depth,Sca_info(2).Sp.dist,Sca_info(2).Sp.inc_ang);
%             F_GS = griddedInterpolant(Sca_info(2).Sp.depth,Sca_info(2).Sp.dist,Sca_info(2).Sp.gs);
%             ind_c = find(Dis_model<=10);
%             [DEP,DIS] = ndgrid(cp_depths,Dis_model(ind_c)); 
%             Dis_sta = sqrt(DEP.^2+(DIS/180*pi.*(Er-DEP)).^2);
%             Incang(ijk).sca_ang = single(F(DEP,DIS));
%             Incang(ijk).gs = single(F_GS(DEP,DIS));
%             Incang(ijk).gs2 = single(1./sqrt(4*pi*Dis_sta.^2));
%             Incang(ijk).ind_dis = ind_c;
%             Incang(ijk).baz_model = baz_model;
%             
%         end
%     end
% 
%     clear Sca_info
    
 %   cd([basedir 'Data/Misc/']);
 %   [latTopo,lonTopo,Topo] = satbath(1,LatLims,LonLims_Temp);
    latTopo = [];
    lonTopo =[];
    Topo = [];
  %  cd(basedir)

    disp('Removing bad waveforms...')
    % Throw out RF_depth that are all zero (i.e. always beyond the critical
    % angle)
    for j = size(RF_depth,1):-1:1
        iscrit=find(imag(CP_lat(:,j))~=0,1);
        if sum(RF_depth(j,:))==0
            RF_depth(j,:) = [];   CP_lat(:,j) = [];
            CP_lon(:,j) = [];     station_at_cp(j) = [];
        elseif ~isempty(iscrit)
            % if the LVZ is low enough velocity, can get below the critical
            % angle again, but will have imaginary values - set these to
            % -999 to make sure nothing wrong gets in
            RF_depth(j,iscrit:end)=0; CP_lat(iscrit:end,j)=-999;
            CP_lon(iscrit:end,j)=-999;
        end
    end
    
    parnum = 1;
    %     parpool(parnum)
    %     parfor kkk=1:rfi
    for kkk=1
        parcalculate(kkk,parnum,Proj_Dir,FileNames,RF_depth,tag,s_inds,RP_conv,baz,Inter_Velo,model_lats,model_lons,CP_lat,CP_lon,cp_depths,time_difference,depth_indx,Phase,phase,tmplat,tmplon,latTopo,lonTopo,Topo,SLON,SLAT)

        %         D=dir([basedir 'Data/CCP/' FileNames.stack '*']);
        %         inds=zeros(length(D),1);
        %         for k=1:length(D);
        %             junk=strsplit(D(k).name,'_');
        %             if strcmp(junk{end}(end-2:end),'mat')
        %                 inds(k)=str2double(junk{end}(1:end-4));
        %             elseif strcmp(junk{end}(end-2:end),'kkk');
%                 inds(k)=str2double(junk{end}(1:end-7));
%             end
%         end
%         for num=1:100; if isempty(find(inds==num,1)); break; end; end
%         
%         savename=[basedir 'Data/CCP/' FileNames.stack num2str(num) '.mat']; 
%         
%         disp(savename);  %eval(['save ' savename ' kkk']); %placeholder
%         parsave(savename,'kkk',kkk);
%         
%         switch tag
%             case 'normal';        clean_indx = 1:size(RF_depth,1);
%             case 'bootstrap';
%                 [clean_indx] = randi(size(RF_depth,1),size(RF_depth,1),1);
%             case 'synth';        clean_indx = 1:size(RF_depth,1);
%         end
%         
%         % initiate variables
%         plane_RF = zeros(length(model_lats),length(cp_depths));
%         weighted_n_events=plane_RF; n_events=plane_RF;
%         
%         disp('Calculating average velocity profile...')
%         % Calculate mean velocity profile for all stations
%         allvels=zeros(size(Velocity_Model(1,1).Vp,1),length(Velocity_Model));
%         alldeps=zeros(size(allvels));
%         meanVelocity=zeros(size(Velocity_Model(1,1).Vp));
%         for i=1:length(Velocity_Model)
%             if length(Velocity_Model(1,i).Vp)==length(allvels(:,1))
%                 allvels(:,i)=Velocity_Model(1,i).Vp;
%                 alldeps(:,i)=Velocity_Model(1,i).Depth;
%             end
%         end
%         for i=length(allvels(1,:)):-1:1
%             if sum(allvels(:,i))==0
%                 allvels(:,i)=[];
%                 alldeps(:,i)=[];
%             end
%         end
%         meanVelocity=mean(allvels,2);
%         meanDepth=mean(alldeps,2);
%         
%         disp('Calculating Fresnel zones...');
%         avg_dists=zeros(length(depth_indx),1);
%         for k = 1:length(depth_indx)
%             %  We are going to make it proportional to the fresnel zone width
%             switch(lower(Phase))
%                 case 'ps'
%                     wavelength = 1 * interp1(Velocity_Model(1).Depth,Velocity_Model(1).Vs,cp_depths(k),'linear');
%                     avg_dists(k) = 0.5*sqrt((wavelength/4 + cp_depths(k))^2 - cp_depths(k)^2);
%                     cutoff=0;
%                 case 'sp'
%                     wavelength = period * interp1(meanDepth,meanVelocity,cp_depths(k),'linear');
%                     avg_dists(k) = 0.5*sqrt((wavelength/3 + cp_depths(k))^2 - cp_depths(k)^2);
%                     cutoff=35;
%             end
%         end
%         avg_dists = avg_dists./111.11;
%         
%         pause(2); clc
%         for j = 1:length(model_lats)
%             if(mod(j,25)==0 || j==1)
%                 display(['Working on ' num2str(j) 'th point of ' num2str(length(model_lats))]);
%             end
% 
%             for k = 1:length(depth_indx)
%                 % Correct distance would be calculated as:
%                 % dalj = distance(model_lats(j),model_lons(j),CP_lat(k,clean_indx),CP_lon(k,clean_indx));
%                 % but, it is much faster to calculate distance like this:
%                 dalj = sqrt((model_lats(j) - CP_lat(k,clean_indx)).^2 + ...
%                     (model_lons(j) - CP_lon(k,clean_indx)).^2);
%                 
%                 fakt = zeros(size(dalj));
%                 avg_dist = avg_dists(k);
%                 indx_infl = find(dalj<=2*avg_dist);
%                 n_events(j,k) = length(indx_infl);
%                 
%                
%                 if(length(indx_infl)>cutoff)  % CHANGE?! to 20 from 35
%                     fakt(indx_infl) = bspl(avg_dist,dalj(indx_infl));
%                     weighted_n_events(j,k)=sum(fakt);
%                     plane_RF(j,k) = fakt*RF_depth(clean_indx,depth_indx(k))./(sum(fakt));
%                 end
%             end
%             if mod(j,2000)==0
% %                 eval(['save ' savename ' phase plane_RF tmplat tmplon cp_depths '...
% %                     'latTopo lonTopo Topo SLON SLAT model_lats model_lons n_events weighted_n_events'])
%                 parsave(savename,'phase',phase,'plane_RF',plane_RF,'tmplat',tmplat,'tmplon',tmplon,'cp_depths',cp_depths,...
%                     'latTopo',latTopo,'lonTopo',lonTopo,'Topo',Topo,'SLON',SLON,'SLAT',SLAT,'model_lats',model_lats,...
%                     'model_lons',model_lons,'n_events',n_events,'weighted_n_events',weighted_n_events);
%             end
%         end
%         
% %         eval(['save ' savename ' phase plane_RF tmplat tmplon cp_depths '...
% %             'latTopo lonTopo Topo SLON SLAT model_lats model_lons n_events weighted_n_events'])
%         parsave(savename,'phase',phase,'plane_RF',plane_RF,'tmplat',tmplat,'tmplon',tmplon,'cp_depths',cp_depths,...
%             'latTopo',latTopo,'lonTopo',lonTopo,'Topo',Topo,'SLON',SLON,'SLAT',SLAT,'model_lats',model_lats,...
%             'model_lons',model_lons,'n_events',n_events,'weighted_n_events',weighted_n_events);
%         
%         
    end
end


end

function parcalculate(kkk,parnum,Proj_Dir,FileNames,RF_depth,tag,s_inds,RP_conv,baz,Inter_Velo,model_lats,model_lons,CP_lat,CP_lon,cp_depths,time_difference,depth_indx,Phase,phase,tmplat,tmplon,latTopo,lonTopo,Topo,SLON,SLAT)


D=dir([Proj_Dir '/CCP/' FileNames.(Phase).stack_slope '*.mat']);
inds=zeros(length(D),1);
for k=1:length(D)
    junk=strsplit(D(k).name,'_');
    if strcmp(junk{end}(end-2:end),'mat')
        inds(k)=str2double(junk{end}(1:end-4));
    elseif strcmp(junk{end}(end-2:end),'kkk');
        inds(k)=str2double(junk{end}(1:end-7));
    end
end
for num=1:100; if isempty(find(inds==num,1)); break; end; end

%num = 1;
savename=[Proj_Dir '/CCP/' FileNames.(Phase).stack_slope num2str(num) '.mat'];

disp(savename);  %eval(['save ' savename ' kkk']); %placeholder
parsave(savename,'kkk',kkk);

switch tag
    case 'normal';        clean_indx = (1:size(RF_depth,1))'; num_indx=ones(size(clean_indx));
    case 'bootstrap'
        [clean_indx_all] = randi(size(RF_depth,1),size(RF_depth,1),1);
        clean_indx = unique(clean_indx_all);
        num_indx = zeros(size(clean_indx));
        for itmp = 1:length(clean_indx)
            num_indx(itmp) = length(find(clean_indx_all==clean_indx(itmp)));
        end
    case 'synth';        clean_indx = 1:size(RF_depth,1);
end

%load([basedir 'Data/CCP/indtmps.mat']);

% initiate variables
plane_RF = zeros(length(model_lats),length(cp_depths));
weighted_n_events=plane_RF; n_events=plane_RF; RP_means=plane_RF;
w2sum=plane_RF; w2RF2sum = plane_RF; w2RFsum = plane_RF; 

nrpcs = 1000; %number of tested ray parameters for distance relationship
dist_model = 0:0.05:10;    %Range of distance for calculating conversion angle
Er = 6371;
if strcmp(Phase,'Sp')
    sigslo = 5;    %How much slope difference is considered significant
    gau_wid = 1;
else
    sigslo = 10;    %How much slope difference is considered significant
    gau_wid = 2;
end
maxGB = 4;

s_inds = s_inds(clean_indx);
RF_depth = RF_depth(clean_indx,:);
CP_lat = CP_lat(:,clean_indx);
CP_lon = CP_lon(:,clean_indx);
RP = RP_conv(clean_indx,:)/Er;
baz = baz(clean_indx);

for i = 1:length(Inter_Velo)
    
    if mod(i,20)==1 && i~=1
        parsave(savename,'phase',phase,'plane_RF',plane_RF,'RP_means',RP_means,'tmplat',tmplat,'tmplon',tmplon,'cp_depths',cp_depths,...
            'latTopo',latTopo,'lonTopo',lonTopo,'Topo',Topo,'SLON',SLON,'SLAT',SLAT,'model_lats',model_lats,'w2sum',w2sum,'w2RF2sum',w2RF2sum,'w2RFsum',w2RFsum,...
            'model_lons',model_lons,'n_events',n_events,'weighted_n_events',weighted_n_events,'cur_sta',i,'clean_indx',clean_indx,'num_indx',num_indx);
    end
    
    if isempty(Inter_Velo(i).vs)
        continue;
    end
    
    rf_ind = find(s_inds == i);
                                                        
    if isempty(rf_ind)
        continue;
    end
    

    disp(['Working on ' num2str(i) 'th station']);
    vs = Inter_Velo(i).vs;
    vp = Inter_Velo(i).vp;
    
    
    switch lower(Phase)
        case 'ps'
            maxrp = 1/vs(1);
            v_con = vs;
            v_inc = vp;
        case 'sp'
            maxrp = 1/vp(1);
            v_con = vp;
            v_inc = vs;
    end
    rp_tmp = linspace(0,maxrp,nrpcs+1);
    rp_cur = rp_tmp(1:end-1);
    [D_con,~,~] = jhua_shootray(v_con',cp_depths',rp_cur,Er);
    D_con = [zeros(size(rp_cur));D_con];
    [DIST,DEPTH] = meshgrid(dist_model,cp_depths);
    RAYPARA = zeros(size(DEPTH));
    rp_max_for_dep = (Er-cp_depths)./([v_con(1),v_con(1:end-1)]);
    rp_max = cummin(rp_max_for_dep)/Er;   %The actual maximum ray parameter for each layer, goes to inf
    RAYPARA(1,:)=nan;
    for idep = 2:length(cp_depths)
        Dists = D_con(idep,:);
        inds_u = find(~isnan(Dists) & rp_cur<rp_max(idep));
        rp_dep = rp_cur(inds_u);
        Dist_dep = Dists(inds_u);
        rp_dep = [rp_dep,rp_max(idep)];
        Dist_dep = [Dist_dep,10e6];
        RAYPARA(idep,:) = interp1(Dist_dep,rp_dep,dist_model);
    end
    SIN_ANG = RAYPARA*Er./(Er-DEPTH).*repmat([v_con(1),v_con(1:end-1)]',[1,length(dist_model)]);
    SIN_ANG(SIN_ANG>1) = 1;
    CON_ANG = asind(SIN_ANG);
%    SO_ANG = asind(RAYPARA*v_con(1));   % incident angle at surface
%    GRAD = gradient(RAYPARA)./gradient(deg2rad(DIST));
%     GS = RAYPARA/4/pi./((Er-DEPTH).^2)./sind(DIST)./cosd(CON_ANG)./cosd(SO_ANG).*abs(GRAD);  %geometrical spreading
%     GS(:,1) = 1/4/pi./((Er-DEPTH(:,1)).^2)./cosd(CON_ANG(:,1))./cosd(SO_ANG(:,1)).*((GRAD(:,1)).^2);
    
    slat = SLAT(i);
    slon = SLON(i);
    Dis_model = 2*asind(sqrt(sind((model_lats-slat)/2).^2+...
        cosd(slat)*cosd(model_lats).*(sind((model_lons-slon)/2).^2)));
    baz_model = atan2d((sind(model_lons-slon).*cosd(model_lats)),...
        (cosd(slat)*sind(model_lats)-sind(slat)*cosd(model_lats).*cosd(model_lons-slon)));
    baz_model(baz_model<0) = baz_model(baz_model<0) + 360;
    F = griddedInterpolant(DEPTH,DIST,CON_ANG);
    %F_GS = griddedInterpolant(DEPTH,DIST,GS);
    ind_c = find(Dis_model<=max(dist_model));
    baz_model = baz_model(ind_c);
    Dis_model = Dis_model(ind_c);
    [DEP,DIS] = ndgrid(cp_depths,Dis_model);
    Dis_sta = sqrt((DEP+(Er-DEP).*(1-cosd(DIS))).^2+(sind(DIS).*(Er-DEP)).^2); 
    Dis_weight = DEP./Dis_sta;
%     Dis_sta = sqrt(DEP.^2+(DIS/180*pi.*(Er-DEP)).^2);
    Sca_ang = F(DEP,DIS);
%     Gs = single(F_GS(DEP,DIS));
%     Gs2 = single(1./sqrt(4*pi*Dis_sta.^2));
    
    % Find the maximum array size (total<maxGB)
    
    size_bas = whos('DIS');
    n_rf = max(floor(maxGB/(size_bas.bytes/1e9*10*parnum)),1);
    for in = 1:ceil(length(rf_ind)/n_rf)
        sec_ind = rf_ind(((in-1)*n_rf+1):(min(in*n_rf,length(rf_ind))));
        
        num_sec = num_indx(sec_ind);  
        baz_sec = baz(sec_ind);
        time_diffsec = time_difference(sec_ind,:);
        RPse = (RP(sec_ind,:))';
        DEPse = (repmat(cp_depths,[length(sec_ind),1]))';
        Vinse = repmat(v_inc',[1,length(sec_ind)]);
        Incse = asind(Er*RPse./(Er-DEPse).*Vinse);
        time_grasec = (gradient(time_diffsec)./gradient(DEPse'))';
       
        Mr = repmat(permute(baz_model,[2,1,3]),[length(depth_indx),1,length(sec_ind)])-...
            repmat(permute(baz_sec,[1,3,2]),[length(depth_indx),length(ind_c),1]);    %model baz- event baz
        Mr(Mr<0) = Mr(Mr<0) + 360;
        
        Sca = repmat(Sca_ang,[1,1,length(sec_ind)]);
        Vcon = repmat(permute(v_con,[2,1,3]),[1,length(ind_c),length(sec_ind)]);
        Vinc = repmat(permute(v_inc,[2,1,3]),[1,length(ind_c),length(sec_ind)]);
        Inc = repmat(permute(Incse,[1,3,2]),[1,length(ind_c),1]);
        
        [LAT_Nsec,LON_Nsec] = reckon(0,0,repmat(Dis_model',[length(depth_indx),1,length(sec_ind)]),Mr);
        
        Slope = abs(atand(sqrt((Vinc.^2).*((sind(Sca)).^2)+(Vcon.^2).*((sind(Inc)).^2)-...
            2*Vinc.*Vcon.*sind(Inc).*cosd(Mr).*sind(Sca).*sign(LAT_Nsec))./(Vcon.*cosd(Inc)-Vinc.*cosd(Sca))));
        SlopeR = abs(atand((Vcon.*sind(Inc)-Vinc.*sind(Sca).*cosd(Mr).*sign(LAT_Nsec))./(Vinc.*cosd(Sca)-Vcon.*cosd(Inc))));
        SlopeT = abs(atand((Vinc.*sind(Sca).*sind(Mr))./(Vinc.*cosd(Sca)-Vcon.*cosd(Inc))));
        
        clear Sca Vcon Vinc Inc
           
        CP_latsec = repmat(permute(CP_lat(:,sec_ind),[1,3,2]),[1,length(ind_c),1]);
        CP_lonsec = repmat(permute(CP_lon(:,sec_ind),[1,3,2]),[1,length(ind_c),1]);
%         model_latsec = repmat(permute(model_lats(ind_c),[2,1,3]),[length(depth_indx),1,length(sec_ind)]);
%         model_lonsec = repmat(permute(model_lons(ind_c),[2,1,3]),[length(depth_indx),1,length(sec_ind)]);
        %         Dis_consec = 2*asind(sqrt(sind((model_latsec-CP_latsec)/2).^2+...
        %             cosd(CP_latsec).*cosd(model_latsec).*(sind((model_lonsec-CP_lonsec)/2).^2)));
        Lat_consec = 2*asind(sqrt(sind((slat-CP_latsec)/2).^2+...
            cosd(CP_latsec).*cosd(slat).*(sind((slon-CP_lonsec)/2).^2))); 
        
        clear CP_latsec CP_lonsec
       
       
        Dep_offsec = (abs(LAT_Nsec-Lat_consec).*tand(SlopeR)+abs(LON_Nsec.*cosd(LAT_Nsec)).*tand(SlopeT))...
            /180*pi.*(Er-repmat(cp_depths',[1,length(ind_c),length(sec_ind)]));
        
        clear Lat_consec LAT_Nsec LON_Nsec SlopeT SlopeR
        
        sigdep = gau_wid./repmat(permute(time_grasec,[1,3,2]),[1,length(ind_c),1]);
        Nevesec = Dep_offsec;
        Nevesec(Dep_offsec>sigdep*2) = 0;
        Nevesec(Dep_offsec<=sigdep*2) = 1;
        Nevesec(Slope>sigslo*2) = 0;
        Nevesec(Slope<=sigslo*2) = 1;
        Weightsec = exp(-(Dep_offsec.^2)/2./(sigdep.^2)).*exp(-(Slope.^2)/2/(sigslo^2)).*repmat(Dis_weight,[1,1,length(sec_ind)]);
        Weightsec(Weightsec<0.02) = 0;
        NSEC = repmat(permute(num_sec,[3,2,1]),[length(depth_indx),length(ind_c),1]);
        Weightsec = Weightsec./repmat(nansum(Weightsec,2),[1,length(ind_c),1]); 
        WRFsec = repmat(permute(RF_depth(sec_ind,:),[2,3,1]),[1,length(ind_c),1]).*Weightsec;
        RP_w = repmat(permute(RPse,[1,3,2]),[1,length(ind_c),1]).*Weightsec;
        plane_RF(ind_c,depth_indx) = plane_RF(ind_c,depth_indx)+(nansum(WRFsec.*NSEC,3))';
        weighted_n_events(ind_c,depth_indx) = weighted_n_events(ind_c,depth_indx)+(nansum(Weightsec.*NSEC,3))';
        w2sum(ind_c,depth_indx) = w2sum(ind_c,depth_indx)+(nansum((Weightsec.^2).*NSEC,3))';
        w2RFsum(ind_c,depth_indx) = w2RFsum(ind_c,depth_indx)+(nansum(WRFsec.*Weightsec.*NSEC,3))';
        w2RF2sum(ind_c,depth_indx) = w2RF2sum(ind_c,depth_indx)+(nansum(WRFsec.*WRFsec.*NSEC,3))';
        RP_means(ind_c,depth_indx) = RP_means(ind_c,depth_indx)+(nansum(RP_w.*NSEC,3))';
        n_events(ind_c,depth_indx) = n_events(ind_c,depth_indx)+(nansum(Nevesec.*NSEC,3))';
        
%         if ~isempty(weighted_n_events(isnan(weighted_n_events)))
%             aaa11=0;
%         end
%         if ~isempty(weighted_n_events(imag(weighted_n_events)~=0))
%             aaa11=0;
%         end
        
        clear sigdep Dep_offsec Nevesec Weightsec Slope WRFsec RP_w NSEC
        
    end 
    
    %pause(10)
    
end

wRFsum = plane_RF;
std_RF = sqrt(w2RF2sum.*(weighted_n_events.^2)+w2sum.*(wRFsum.^2)-2*weighted_n_events.*wRFsum.*w2RFsum)./(weighted_n_events.^2);
plane_RF = plane_RF./weighted_n_events;
RP_means = RP_means./weighted_n_events;
weighting_norm = weighted_n_events./repmat(nansum(weighted_n_events,1),size(weighted_n_events,1),1)*nansum(nansum(weighted_n_events,1))/size(weighted_n_events,2);

parsave(savename,'phase',phase,'plane_RF',plane_RF,'RP_means',RP_means,'tmplat',tmplat,'tmplon',tmplon,'cp_depths',cp_depths,'weighting_norm',weighting_norm,...
    'latTopo',latTopo,'lonTopo',lonTopo,'Topo',Topo,'SLON',SLON,'SLAT',SLAT,'model_lats',model_lats,'w2sum',w2sum,'w2RF2sum',w2RF2sum,'w2RFsum',w2RFsum,...
    'model_lons',model_lons,'n_events',n_events,'weighted_n_events',weighted_n_events,'wRFsum',wRFsum,'clean_indx',clean_indx,'num_indx',num_indx,'std_RF',std_RF);
end

function parsave(savename,varargin)
    tmp = ['save ',savename,];
    for i = 1:length(varargin)/2
        tmp_name = varargin{2*i-1};
        eval([tmp_name,'=varargin{2*i};']);
        tmp = [tmp,' ',tmp_name];
    end
    eval(tmp)
end

function [out]=calcmeanbootstrap(basedir, Proj_Dir, version)
load([Proj_Dir version '.mat']);
D=dir([basedir 'Data/CCP/' FileNames.stack '*']);
nbt=length(D);
if nbt~=50; disp(['You have ' num2str(nbt) ' bootstraps!!']); pause; end

load([basedir 'Data/CCP/' D(1).name])
all_plane_RF=zeros(length(model_lats),length(cp_depths),nbt);
all_weighted_n_events=all_plane_RF;


for ibt=1:nbt
    disp(['Working on bootstrap number ' num2str(ibt)]);
    if ibt==1;  rfs.plane_RF=plane_RF; rfs.weighted_n_events=weighted_n_events;
    else        rfs=load([basedir 'Data/CCP/' D(ibt).name]);       
    end
    
    all_plane_RF(:,:,ibt)=rfs.plane_RF;
    all_weighted_n_events(:,:,ibt)=rfs.weighted_n_events;
    clear rfs

end

disp('Calculating means...');
plane_RF=mean(all_plane_RF,3); weighted_n_events=mean(all_weighted_n_events,3);
std_RF=zeros(length(model_lats), length(cp_depths));
for id=1:length(cp_depths)
    std_RF(:,id,:)=std(squeeze(all_plane_RF(:,id,:))')';
end

out.plane_RF=plane_RF; out.std_RF=std_RF; out.tmplat=tmplat; out.tmplon=tmplon;
out.cp_depth=cp_depths; out.latTopo=latTopo; out.lonTopo=lonTopo; out.Topo=Topo;
out.SLON=SLON; out.SLAT=SLAT; out.model_lats=model_lats; out.model_lons=model_lons;
out.n_events=n_events; out.weighted_n_events=weighted_n_events;

FileNames.meanstack=['meanbtstrap_' FileNames.stack(1:end-1)];

%eval(['save' basedir 'Data/CCP/' FileNames.meanstack '.mat phase plane_RF std_RF tmplat tmplon cp_depths '...
%    'latTopo lonTopo Topo SLON SLAT model_lats model_lons n_events weighted_n_events']);
save([basedir,'Data/CCP/' FileNames.meanstack '.mat'],'phase', 'plane_RF', 'std_RF','tmplat','tmplon','cp_depths',...
    'latTopo','lonTopo','Topo','SLON','SLAT','model_lats','model_lons','n_events','weighted_n_events');

save([Proj_Dir version '.mat'],'FileNames');


end




function [x,t,sn]=ved2_shootray(v,z,p)
% Modified by Ved to return x and t for each layer.
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
sn = v(iprop)*p;
%
% sn is an n by m matrix where n is the length of iprop (the number of
%    layers propagated through) and m is the length of p (the number of
%    unique ray parameters to use). Each column of sn corresponds to a
%    single ray parameter and contains the sin of the vertical angle;
%

[ichk,pchk]=find(sn>1);

%compute x and t
cs=sqrt(1-sn.*sn)+eps;
vprop=v(iprop)*ones(1,length(p));
thk=abs(diff(z))*ones(1,length(p));
if(size(sn,1)>1)
    x=cumsum( (thk.*sn)./cs);
    t=cumsum(thk./(vprop.*cs));
else
    x=(thk.*sn)./cs;
    t=thk./(vprop.*cs);
end
%assign infs
if(~isempty(ichk))
    x(ichk)=inf*ones(size(pchk));
    t(ichk)=inf*ones(size(pchk));
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
sn = (v(iprop)*re./(re-z(iprop)))*p;
%
% sn is an n by m matrix where n is the length of iprop (the number of
%    layers propagated through) and m is the length of p (the number of
%    unique ray parameters to use). Each column of sn corresponds to a
%    single ray parameter and contains the sin of the vertical angle;
%

%ichk=find(sn>1);

%compute x and t
cs=sqrt(1-sn.*sn)+eps;
vprop=v(iprop)*ones(1,length(p));
thk=abs(diff(z))*ones(1,length(p));
R = (re-z(iprop))*ones(1,length(p));
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

function T=calcDomPeriod(dir0, phase, name,version)
disp('Calculating dominant period...')
ph=[phase '_Waveforms']; load([dir0 'Networks.mat']);
junk=strsplit(dir0,'Projects/'); 
if length(junk{2})>6 && strcmp(junk{2}(1:6),'Synth_')
    load([dir0 version '.mat']); synthstr=FileNames.Synth;
else synthstr='';
end


nets=fieldnames(Networks);
jk = 0;
for j = 1:length(nets)
    stations = fieldnames(Networks.(nets{j}));
    for k = 1:length(stations)
        jk = jk + 1;
        folder{jk} = strcat(dir0,nets{j},'/',stations{k},'/');
        sta{jk}=stations{k};
    end
end

Ts=zeros(1,40000); p=0;

for i = 1:length(folder)
    disp(sta{i})
    pause(0.02)
    clear F N RFs_in_Time dt fft_exp ind j k jk junk
    if exist([folder{i} 'PreCalculated_RFs' name], 'file')
        RFs = load([folder{i} 'PreCalculated_RFs' name]);
        load([folder{i} 'Waveform_Data' synthstr '.mat'])
        events = fieldnames(RFs.PreCalculated_RFs.(phase).Individual_RFs);

  
        for j = 1:length(events)
            S2N_cutoff=5; % 0; % 1; %
            
            if ~isempty(RFs.PreCalculated_RFs.(phase).Individual_RFs.(events{j})) ...
                    && All_Waveform_Data.(events{j}).(ph).S2N_ratio>S2N_cutoff
                p=p+1;
%                 F=dfft(RFs.PreCalculated_RFs.(phase).Individual_RFs.(events{j}));
%                 F(1:2)=F(1:2)/2;
%                 [~,ind] = max(F(2:end));

                
%                 Fs = 1/(time(2)-time(1));
%                 Nr = length(RFs.PreCalculated_RFs.(phase).Individual_RFs.(events{j}));
%                 F = fft(RFs.PreCalculated_RFs.(phase).Individual_RFs.(events{j}));
                staind = max(All_Waveform_Data.(events{j}).(ph).Window.indx_full_phase(1),...
                    All_Waveform_Data.(events{j}).(ph).Window.indx_full_analysis(1))-...
                    All_Waveform_Data.(events{j}).(ph).Window.indx_full_analysis(1)+1;
                
                endind = min(All_Waveform_Data.(events{j}).(ph).Window.indx_full_phase(2),...
                    All_Waveform_Data.(events{j}).(ph).Window.indx_full_analysis(2))-...
                    All_Waveform_Data.(events{j}).(ph).Window.indx_full_analysis(1)+1;
                
                if mod(endind-staind+1,2)~=0
                    endind = endind-1;
                end
                time = All_Waveform_Data.(events{j}).(ph).Time(staind:endind);
                Fs = 1/(time(2)-time(1));
                Nr = length(time);
                
                F = fft(All_Waveform_Data.(events{j}).(ph).SV(staind:endind));
                F=F'; 
                F = F(1:Nr/2+1);
                freq = Fs*(0:(Nr/2))/Nr;
                psdx = (1/(Fs*Nr))*abs(F).^2;
                psdx(2:end-1) = 2*psdx(2:end-1);
                psdx=psdx';
                [~,ind] = max(psdx(2:end-1));
                Ts(p) = 1/freq(ind+1);

            end
            
        end
    end
    
end

if length(Ts)>p; Ts(p+1:end)=[]; end

T=mean(Ts); clc

end

function f = dfft(y)

% Use standard Matlab routine to find Fourier transform of y.
% z contains the complex coefficients of the Fourier exponential series.
z=fft(y);

n = length(y);
half=n/2;

% Take the exponential series coefficients and derive the
% coefficients of the Fourier Sine and Cosine series.
% f(1) = zero frequency cosine
%          = average value of the function over period sampled
% f(2) = cosine coefficient at frequency n/2

for i = 1:2:n
    j=(i+1)/2;
    f(i)=real(z(j))/half;
    f(i+1)=-imag(z(j))/half;
end

j=n/2+1;
f(2)=real(z(j))/half;
end


function[b] = bspl(d0,d)   % Calculates value of b-spline
% Spline definition from Wang and Dahlen, 1995
for j = 1:length(d)
    if(d(j)<=d0) 
        b(j) = 0.75*(d(j)^3)/d0^3 - 1.5*(d(j)^2)./d0^2 + 1;
    else
        if(d(j)<=2*d0) 
            %As in WD1995, but we use a simpler expression
            %d1 = (d(j)-d0)/d0;
            %b(j) = -0.25*d1^3 + 0.75*d1^2 - 0.75*d1 + 0.25; 
            b(j) = 0.25*(2-d(j)/d0)^3;
        else
            b(j) = 0;
        end
    end
end
end
