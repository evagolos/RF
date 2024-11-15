function CCP_stack(Project, basedir, Phases, tag, varargin)

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
        calc_deplatlon(basedir, Project, Phases,data_tag,version,strike)
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


function calc_deplatlon(basedir, project, Phases, tag, version,strike)

% ********* Function Description *********
%
% Migrate RFs in time into depth and associate each depth with a lat,lon
% based on shooting at the parent ray parameter along the backazimuth.
% Write all the RFs to a data structure that can be read and used by
% CCP stacking programs. - Ved Aug 13, 2010
%
% ****************************************

dir0 = [basedir 'Data/Projects/' project '/'];
Er = 6371; km_per_deg  = pi*Er/180; Ellipsoid   = [Er 1]/km_per_deg;
load([dir0 'Networks.mat']);
load([dir0 'Similarity_RF_450_oldway.mat']);
dist_sort = sort((Total_info.Dist_tot_2));
cutoff = 5*median(Total_info.Dist_tot_2);
cutoff_shallow_size = 0.2*median(Total_info.Self_shallow_Moho_size);
cutoff_shallow_pv = 3*median(Total_info.Self_shallow_pv_size);
cutoff_deep_num = median(Total_info.Num_deep_out_2);

if exist([dir0 version '.mat'],'file'); load([dir0 version '.mat']);
else disp('No PreCalc RFs name saved!'); return;
end
name=FileNames.PreCalcRFs;
if isfield(FileNames,'Synth'); synthstr=FileNames.Synth; else synthstr=''; end 
mods=textscan(name,'%s','delimiter','_'); if length(mods{1}{end})==5; n=1; else n=0; end
fprintf(['Your migration models: \n  Crust      - %s\n  Mantle Vp  - ' ...
    '%s\n  Mantle Vs  - %s\n'],mods{1}{end-3-n},mods{1}{end-2-n},mods{1}{end-1-n});
fprintf('\nHigh frequency corner: %ss\n',mods{1}{3});

if isempty(strike); minbaz=-180; maxbaz=180;
else minbaz=strike-20; maxbaz=strike+20;
end

for iph=1:length(Phases)
    phase=Phases{iph};
    
    network_codes=fieldnames(Networks);
    
    jk = 0;
    for j = 1:length(network_codes)
        stations = fieldnames(Networks.(network_codes{j}));
        for k = 1:length(stations)
            jk = jk + 1;
            folder{jk} = strcat(dir0,network_codes{j},'/',stations{k},'/');
        end
    end
    
    how_sig = 2; % How many standard deviations away from 0 to be significant
    n_boot_iter = 1; % Number of bootstrapping iterations to carry out
    map_size = 7; % Size of map in degrees (on each side of station)
    minXcorr = -1; minS2N = 2; % Minimum allowable Xcorr and s2n
    
    typewf = strcat(phase,'_Waveforms');
    switch(lower(phase)); case 'ps'; ime = 'P_Path'; case 'sp'; ime = 'S_Path'; end
    nj = 0; n = 0;
    
    % Initialise matrices
    maxdep=450; deps=0:0.5:maxdep;    %change bin size
    CP_lat=zeros(76000,length(deps)); CP_lon=CP_lat; RF_depth=CP_lat; time_difference=CP_lat;
    SLAT=zeros(length(folder),1); SLON=SLAT; station_at_cp=zeros(1,length(CP_lat(:,1)));
    baz=station_at_cp; RP=station_at_cp;
    
    
    % Calculations begin here
    for jkjk = 1:length(folder)
        
        folderName=textscan(folder{jkjk}, '%s', 'delimiter','/');
        net=folderName{1}{end-1};  sta=folderName{1}{end};
        
        
        display(['Working on station number ' num2str(jkjk) '...']);
        
        clear depth RPs M depths delay RPs_grid depths_grid quakes time
        clear dist baz RP RFs mask_i Mask_Depth Relative_DT Mask_Time dt T Z
        clear sort_by_dist sort_by_baz yd yb nRFs_by_dist nRFs_by_baz
        
       if strcmp(tag,'synth'); 
           junk=strsplit(folder{jkjk}, 'Synth_'); migr_folder=[junk{1} junk{2}];
       else migr_folder=folder{jkjk};
       end
        
        if(exist([folder{jkjk} '/PreCalculated_RFs' name],'file') && ...
                exist([migr_folder '/Migration_Models.mat'],'file') && ...
                exist([folder{jkjk} '/Ray_Path_Data.mat'],'file') && ...
                exist([folder{jkjk} '/Station_Data.mat'],'file') && ...
                exist([folder{jkjk} '/Waveform_Data' synthstr '.mat'],'file'))
            
            load([folder{jkjk} 'Station_Data.mat']);
            SLAT(jkjk) = Station_Data.Latitude; SLON(jkjk) = Station_Data.Longitude;
            
            load([migr_folder 'Migration_Models.mat']);
            which_model=length(Migration_Models.(phase));
            
            
            for jm=length(Migration_Models.(phase)):-1:1
                if(strcmp(Migration_Models.(phase)(jm).Original_Models.Mantle_Vp, ...
                        mods{1}{end-2}) && ...
                        strcmp(Migration_Models.(phase)(jm).Original_Models.Crust, ...
                        mods{1}{end-3}) && ...
                        strcmp(Migration_Models.(phase)(jm).Original_Models.Mantle_Vs,...
                        mods{1}{end-1}))
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

            depths = min(depth):0.5:max(depth);    %change bin size
            % Define subregion in depth that we are actually interested in
            depth_indx = find(depths>=0 & depths<=maxdep); %250 200
            cp_depths = depths(depth_indx);
            
            Velocity_Model(jkjk) = Migration_Models.(phase)(which_model).Velocity_Model;
            display('Using Standard Velocity Model');
            
            load([folder{jkjk} strcat('/PreCalculated_RFs', name)]);
            load([folder{jkjk} 'Ray_Path_Data.mat']);
            load([folder{jkjk} 'Waveform_Data' synthstr '.mat']);
            
            duz = length(cp_depths);
            
            % Go through and select all the available pre-calculated RFs
            %quakes = fieldnames(PreCalculated_RFs.(phase).Individual_RFs);
            quakes =  CompRFres.(net).(sta).eve;
            quakes(CompRFres.(net).(sta).dist_2>cutoff | ...
                CompRFres.(net).(sta).RFsize_shallow_Moho<cutoff_shallow_size | ...
                CompRFres.(net).(sta).N_deep_out_range_2>cutoff_deep_num | ...
                CompRFres.(net).(sta).RFsize_shallow_pv>cutoff_shallow_pv ) = [];
            if isempty(quakes)
                continue;
            end
            
            vs = interp1(Velocity_Model(jkjk).Depth,Velocity_Model(jkjk).Vs,cp_depths,'linear');
            vp = interp1(Velocity_Model(jkjk).Depth,Velocity_Model(jkjk).Vp,cp_depths,'linear');
                            
            
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
                            
%                             if(max(abs(RF_in_time))>0)
%                                 % Normalize
%                                 RF_in_time = RF_in_time./max(abs(RF_in_time));
%                             end
                            
                            
                            baz(n)  = All_Ray_Path_Data.(quakes{j}).(ime).BAZ;
                            RP(n)   = All_Ray_Path_Data.(quakes{j}).(ime).ray_parameter;
                                                        
                            % First, we need to get the velocity structure at that
                            % point.
                            
                            indx_same = find(diff(Velocity_Model(jkjk).Depth)==0) ;
                            Velocity_Model(jkjk).Depth(indx_same+1) = Velocity_Model(jkjk).Depth(indx_same+1) + 0.01;
                            CP_lat(n,1) = SLAT(jkjk); CP_lon(n,1) = SLON(jkjk); RF_depth(n,1) = 0;
                            
                           
%                             clear xs xp ts tp x t;
%                             xs(1) = 0; xp(1) = 0; ts(1) = 0; tp(1) = 0;
%                             
%                             [xs(2:duz), ts(2:duz), ~] = ved2_shootray(vs',cp_depths',RP(n));
%                             [xp(2:duz), tp(2:duz), ~] = ved2_shootray(vp',cp_depths',RP(n));
                            
                            clear Ds Dp ts tp D t;
                            Ds(1) = 0; Dp(1) = 0; ts(1) = 0; tp(1) = 0;
                            [Ds(2:duz),ts(2:duz),~] = jhua_shootray(vs',cp_depths',RP(n),Er);
                            [Dp(2:duz),tp(2:duz),~] = jhua_shootray(vp',cp_depths',RP(n),Er);
                            
                            % Both the difference in travel time of P and S waves
                            % from a depth and the effect that the incident parent
                            % plane wave will impact the interface at different
                            % times (rp * L) have to be accounted for.
                            
%                             t = abs(ts-tp) + RP(n).*(abs(xp-xs));
%                             time_difference(n,:)=t;
%                             
%                             can_use = ~isinf(t);
%                             cant_use = isinf(t);
                            
                            t = ts-tp + RP(n)*Er.*(Dp-Ds)/180*pi;
                            time_difference(n,:)=t;
                            
                            can_use = ~isnan(t);
                            cant_use = isnan(t);
                            % shootray checks for critical angle;
                            % will set t = inf if have passed it
                            
%                             switch lower(phase)
%                                 case 'ps'
%                                     fakt = 1; %flipping or not time axis
%                                     x = xs;
%                                 case 'sp'
%                                     fakt = -1; %flipping or not time axis
%                                     x = xp;
%                             end 
%                             
%                             [CP_lat(n,can_use), CP_lon(n,can_use)] = reckon(SLAT(jkjk),SLON(jkjk),x(can_use)./km_per_deg,baz(n),'degrees');
%                             RF_depth(n,can_use) = interp1(PreCalculated_RFs.(phase).Time(1:length(RF_in_time)), RF_in_time, fakt*t(can_use),'pchip');
%                             
%                             RF_depth(n,cant_use) = 0; CP_lat(n,cant_use) = -999;
%                             CP_lon(n,cant_use) = -999;
%                             
                            
                            
                            switch lower(phase)
                                case 'ps'
                                    fakt = 1; %flipping or not time axis
                                    D = Ds;
                                case 'sp'
                                    fakt = -1; %flipping or not time axis
                                    D = Dp;
                            end
                            
                            [CP_lat(n,can_use), CP_lon(n,can_use)] = reckon('gc',SLAT(jkjk),SLON(jkjk),D(can_use),baz(n),'degrees');
                            RF_depth(n,can_use) = interp1(PreCalculated_RFs.(phase).Time(1:length(RF_in_time)), RF_in_time, fakt*t(can_use),'pchip');
                            
                            RF_depth(n,cant_use) = 0; CP_lat(n,cant_use) = -999;
                            CP_lon(n,cant_use) = -999;
                            
                        end
                    end
                end
            end
            
        end
    end
    
   % name=strcat(project, name(1:end-4), '_', num2str(maxdep),'km_Oldverision_', num2str(strike), 'deg.mat');
    name=strcat(project, name(1:end-4), '_', num2str(maxdep),'km_Allold_', num2str(strike), 'deg.mat');
    
    CP_lat(n+1:end,:)=[]; CP_lon(n+1:end,:)=[]; RF_depth(n+1:end,:)=[];
    time_difference(n+1:end,:)=[]; station_at_cp(n+1:end)=[];
    
    if ~exist([basedir '/Data/CCP'],'dir'); mkdir([basedir 'CCP']); end
    
    eval(['save ' basedir '/Data/CCP/AllNtwrks_' phase '_' name ...
        ' RF_depth CP_lat CP_lon SLAT SLON cp_depths depth_indx' ...
        ' phase Velocity_Model folder station_at_cp time_difference name project']) %depths  RP
    FileNames.CCP=name; save([dir0 version '.mat'],'FileNames');
    
    
end

end

function CCP_stack_from_precalc_RF(basedir, Proj_Dir,Phases,version,tag,varargin)

load([Proj_Dir version '.mat']);
if isfield(FileNames,'CCP'); name=FileNames.CCP;
else disp('No dep/lat/lon matrix name saved!'); return
end

for iph=1:length(Phases);
    Phase=Phases{iph};
    load([basedir '/Data/CCP/AllNtwrks_' Phase '_' name]);
    aaa=load([basedir '/Data/CCP/AllNtwrks_' Phase '_' name]);
    Velocity_Model = aaa.Velocity_Model;
    cp_depths = aaa.cp_depths;
    depth_indx = aaa.depth_indx;
    CP_lat = aaa.CP_lat;
    CP_lon = aaa.CP_lon;
    SLAT = aaa.SLAT;
    SLON = aaa.SLON;
    phase = aaa.phase;
    
    switch tag
        case 'normal'; rfi=1; period = 13; %period=calcDomPeriod(Proj_Dir, Phase, FileNames.PreCalcRFs,version); %Junlin
            
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

    FileNames.stack=['CCP_' Phase '_' num2str(period) '_' FileNames.CCP(1:end-4) '_'];
    save([Proj_Dir version '.mat'],'FileNames');
    
    CP_lat = transpose(CP_lat); CP_lon = transpose(CP_lon); map_size = 4;
    LatLims = [max(min(SLAT(find(SLAT))),34)-map_size  max(SLAT(find(SLAT)))+map_size];
    LonLims = [min(SLON(find(SLON)))-map_size  max(SLON(find(SLON)))+map_size];
    LonLims_Temp = LonLims;
    LonLims_Temp(LonLims_Temp>180)  = LonLims(LonLims_Temp>180)-360;
    [tmplat,tmplon] = meshgrid(LatLims(1):0.1:LatLims(2), LonLims_Temp(1):0.1:LonLims_Temp(2));
    model_lats = tmplat(:); model_lons = tmplon(:);
    
    cd([basedir 'Data/Misc/']);
    [latTopo,lonTopo,Topo] = satbath(1,LatLims,LonLims_Temp); 
    cd(basedir)

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
    
    
    %parpool(12)
    for kkk=1:rfi
%     for kkk=1:rfi
%     for kkk=1
        pause(kkk-1);
        parcalculate(kkk,basedir,FileNames,RF_depth,tag,model_lats,model_lons,CP_lat,CP_lon,cp_depths,Velocity_Model,period,depth_indx,Phase,phase,tmplat,tmplon,latTopo,lonTopo,Topo,SLON,SLAT)
        
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
%    delete(gcp)
end


end

function parcalculate(kkk,basedir,FileNames,RF_depth,tag,model_lats,model_lons,CP_lat,CP_lon,cp_depths,Velocity_Model,period,depth_indx,Phase,phase,tmplat,tmplon,latTopo,lonTopo,Topo,SLON,SLAT)
%for kkk=1


D=dir([basedir 'Data/CCP/' FileNames.stack '*']);
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

num=1;
savename=[basedir 'Data/CCP/' FileNames.stack num2str(num) '.mat'];

disp(savename);  %eval(['save ' savename ' kkk']); %placeholder
parsave(savename,'kkk',kkk);

switch tag
    case 'normal';        clean_indx = 1:size(RF_depth,1);
    case 'bootstrap';
        [clean_indx] = randi(size(RF_depth,1),size(RF_depth,1),1);
    case 'synth';        clean_indx = 1:size(RF_depth,1);
end

% initiate variables
plane_RF = zeros(length(model_lats),length(cp_depths));
weighted_n_events=plane_RF; n_events=plane_RF; std_RF = plane_RF;

disp('Calculating average velocity profile...')
% Calculate mean velocity profile for all stations
allvels=zeros(size(Velocity_Model(1,1).Vp,1),length(Velocity_Model));
alldeps=zeros(size(allvels));
meanVelocity=zeros(size(Velocity_Model(1,1).Vp));
for i=1:length(Velocity_Model)
    if length(Velocity_Model(1,i).Vp)==length(allvels(:,1))
        allvels(:,i)=Velocity_Model(1,i).Vp;
        alldeps(:,i)=Velocity_Model(1,i).Depth;
    end
end
for i=length(allvels(1,:)):-1:1
    if sum(allvels(:,i))==0
        allvels(:,i)=[];
        alldeps(:,i)=[];
    end
end
meanVelocity=mean(allvels,2);
meanDepth=mean(alldeps,2);

disp('Calculating Fresnel zones...');
avg_dists=zeros(length(depth_indx),1);
for k = 1:length(depth_indx)
    %  We are going to make it proportional to the fresnel zone width
    switch(lower(Phase))
        case 'ps'
            wavelength = 1 * interp1(Velocity_Model(1).Depth,Velocity_Model(1).Vs,cp_depths(k),'linear');
            avg_dists(k) = 0.5*sqrt((wavelength/4 + cp_depths(k))^2 - cp_depths(k)^2);
            cutoff=0;
        case 'sp'
            wavelength = period * interp1(meanDepth,meanVelocity,cp_depths(k),'linear');
            avg_dists(k) = 0.5*sqrt((wavelength/3 + cp_depths(k))^2 - cp_depths(k)^2);
            cutoff=35;
            %cutoff=20;
    end
end
avg_dists = avg_dists./111.11;
cur_sta = 1;
%load(savename)

pause(2); clc
for j = cur_sta:length(model_lats)
%for j = 12001:length(model_lats)    
    if(mod(j,25)==0 || j==1)
        display(['Working on ' num2str(j) 'th point of ' num2str(length(model_lats))]);
    end
    
    for k = 1:length(depth_indx)
        % Correct distance would be calculated as:
        % dalj = distance(model_lats(j),model_lons(j),CP_lat(k,clean_indx),CP_lon(k,clean_indx));
        % but, it is much faster to calculate distance like this:
        dalj = sqrt((model_lats(j) - CP_lat(k,clean_indx)).^2 + ...
            (model_lons(j) - CP_lon(k,clean_indx)).^2);
        
        fakt = zeros(size(dalj));
        avg_dist = avg_dists(k);
        indx_infl = find(dalj<=2*avg_dist);
        n_events(j,k) = length(indx_infl);
        
        
        if(length(indx_infl)>cutoff)  % CHANGE?! to 20 from 35
            fakt(indx_infl) = bspl(avg_dist,dalj(indx_infl));
            weighted_n_events(j,k)=sum(fakt);
%             weight = fakt(indx_infl);
%             RF_eff = RF_depth(indx_infl,depth_indx(k));
%             [RF_eff,ind_od] = sort(RF_eff);
%             weight = weight(ind_od');
%             weight_norm = weight/sum(weight);
%             weight_cum = cumsum(weight_norm);
%             ind_rf = find(abs(weight_cum-0.5)==min(abs(weight_cum-0.5)));
%             if length(ind_rf)==1
%                 plane_RF(j,k) = RF_eff(ind_rf);
%             else
%                 plane_RF(j,k) = mean(RF_eff(ind_rf));
%             end
            plane_RF(j,k) = fakt*RF_depth(clean_indx,depth_indx(k))./(sum(fakt));
            w2sum = sum(fakt.^2);
            w2RFsum = sum((fakt.^2)*RF_depth(clean_indx,depth_indx(k)));
            w2RF2sum = sum((fakt.^2)*((RF_depth(clean_indx,depth_indx(k))).^2));
            wRFsum = sum(fakt*RF_depth(clean_indx,depth_indx(k)));
            wsum = sum(fakt);
            std_RF(j,k) = sqrt(w2RF2sum*(wsum^2)+w2sum*(wRFsum^2)-2*wsum*wRFsum*w2RFsum)/(wsum^2);
        end
    end
    if mod(j,3000)==0
        %                 eval(['save ' savename ' phase plane_RF tmplat tmplon cp_depths '...
        %                     'latTopo lonTopo Topo SLON SLAT model_lats model_lons n_events weighted_n_events'])
        parsave(savename,'phase',phase,'plane_RF',plane_RF,'tmplat',tmplat,'tmplon',tmplon,'cp_depths',cp_depths,...
            'latTopo',latTopo,'lonTopo',lonTopo,'Topo',Topo,'SLON',SLON,'SLAT',SLAT,'model_lats',model_lats,...
            'model_lons',model_lons,'n_events',n_events,'weighted_n_events',weighted_n_events,'std_RF',std_RF,'clean_indx',clean_indx,'cur_sta',j);
    end
end

%         eval(['save ' savename ' phase plane_RF tmplat tmplon cp_depths '...
%             'latTopo lonTopo Topo SLON SLAT model_lats model_lons n_events weighted_n_events'])
parsave(savename,'phase',phase,'plane_RF',plane_RF,'tmplat',tmplat,'tmplon',tmplon,'cp_depths',cp_depths,...
    'latTopo',latTopo,'lonTopo',lonTopo,'Topo',Topo,'SLON',SLON,'SLAT',SLAT,'model_lats',model_lats,...
    'model_lons',model_lons,'n_events',n_events,'weighted_n_events',weighted_n_events,'std_RF',std_RF,'clean_indx',clean_indx);
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

[ichk,pchk]=find(sn>1,1);
if ~isempty(ichk)
    ichk = ichk:length(sn);
    pchk = ones(size(ichk));
end

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

Ts=zeros(1,4000); p=0;

i=1; RFs=load([folder{i} 'PreCalculated_RFs' name]);
time=RFs.PreCalculated_RFs.(phase).Time;

for i = 1:length(folder)
    disp(sta{i})
    pause(0.02)
    clear coeffs F N RFs_in_Time dt fft_exp ind j k jk junk
    if exist([folder{i} 'PreCalculated_RFs' name], 'file')
        RFs = load([folder{i} 'PreCalculated_RFs' name]);
        load([folder{i} 'Waveform_Data' synthstr '.mat'])
        events = fieldnames(RFs.PreCalculated_RFs.(phase).Individual_RFs);

        
        fft_exp = floor(log2(length(time)));
        N=2^fft_exp;
        coeffs(1:2)=[0,N];
        for j = 1:N/2-1; coeffs(2*j+1:2*j+2) = 2*j;   end
  
        for j = 1:length(events);
            S2N_cutoff=5; % 0; % 1; %
            
            if ~isempty(RFs.PreCalculated_RFs.(phase).Individual_RFs.(events{j})) ...
                    && All_Waveform_Data.(events{j}).(ph).S2N_ratio>S2N_cutoff
                p=p+1;
%                 F=dfft(RFs.PreCalculated_RFs.(phase).Individual_RFs.(events{j}));
%                 F(1:2)=F(1:2)/2;
%                 [~,ind] = max(F(2:end));

                
                Fs = 1/(time(2)-time(1));
                Nr = length(RFs.PreCalculated_RFs.(phase).Individual_RFs.(events{j}));
                F = fft(RFs.PreCalculated_RFs.(phase).Individual_RFs.(events{j}));
                
%                 time = All_Waveform_Data.(events{j}).(ph).Time;
%                 Fs = 1/(time(2)-time(1));
%                 Nr = length(time);
%                 F = fft(All_Waveform_Data.(events{j}).(ph).SV);
%                 F=F';
                if mod(Nr,2)==0    
                    F = F(1:Nr/2+1);
                else
                    F = F(1:(Nr+1)/2);
                end
                psdx = (1/(Fs*Nr))*abs(F).^2;
                psdx(2:end-1) = 2*psdx(2:end-1);
                psdx=psdx';
                %freq = 0:Fs/N:Fs/2;
                ind = sum(psdx(2:end-1).*[2:(length(psdx)-1)]/sum(psdx(2:end-1)));
                Ts(p) = 1/((ind-1)/(length(psdx)-1)*Fs/2);
                
                 %Ts(p) = 2*(max(time)-min(time))/coeffs(ind+1);
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
