function update_eu_Data(Proj_Dir, basedir)
% The station information is found on http://orfeus-eu.org/opencms/stationbook/search-stations/
% Events.mat is supposed to be updated to NEIC

load([basedir 'Data/Event_Lists/NEIC_Catalog.mat']);
%TMP=load([Proj_Dir 'Events.mat']); Events = TMP.Events; clear TMP
TMP=load([Proj_Dir 'Networks.mat']); oldNetworks = TMP.Networks; clear TMP
TMP=load([Proj_Dir 'Networks_before_iris_update.mat']); olderNetworks = TMP.Networks; clear TMP
system(['cp ' Proj_Dir 'Networks.mat ' Proj_Dir 'Networks_' datestr(date,'dd.mm.yy') '.mat']);
system(['cp ' Proj_Dir 'Events.mat ' Proj_Dir 'Events_' datestr(date,'dd.mm.yy') '.mat']);

junk=fieldnames(NEIC_Catalog);
%Remove Multiples in NEIC_Catalog
for iy = 1:length(junk)
    oritime_cur = double(NEIC_Catalog.(junk{iy}).Origin_Time);
    datenum_cur = datenum([oritime_cur(:,1:5),oritime_cur(:,6)]);
    [~,ind] = unique(datenum_cur);
    NEIC_Catalog.(junk{iy}).Origin_Time = NEIC_Catalog.(junk{iy}).Origin_Time(ind,:);
    NEIC_Catalog.(junk{iy}).Location = NEIC_Catalog.(junk{iy}).Location(ind,:);
    NEIC_Catalog.(junk{iy}).Magnitude = NEIC_Catalog.(junk{iy}).Magnitude(ind,:);
    oritime_cur = double(NEIC_Catalog.(junk{iy}).Origin_Time);
    datenum_cur = datenum([oritime_cur(:,1:5),oritime_cur(:,6)]);
    [~,ind] = sort(datenum_cur);
    NEIC_Catalog.(junk{iy}).Origin_Time = double(NEIC_Catalog.(junk{iy}).Origin_Time(ind,:));
    NEIC_Catalog.(junk{iy}).Location = double(NEIC_Catalog.(junk{iy}).Location(ind,:));
    NEIC_Catalog.(junk{iy}).Magnitude = double(NEIC_Catalog.(junk{iy}).Magnitude(ind,:));
end
save([basedir 'Data/Event_Lists/NEIC_Catalog.mat'],'NEIC_Catalog');


datenum_last = datenum(double(NEIC_Catalog.(junk{end}).Origin_Time(:,1:3)));
last_ind = find(datenum_last==max(datenum_last));
lastNEIC=NEIC_Catalog.(junk{end}).Origin_Time(last_ind(1),1:3);
datenum_first = datenum(double(NEIC_Catalog.(junk{1}).Origin_Time(:,1:3)));
first_ind = find(datenum_first==min(datenum_first));
firstNEIC=NEIC_Catalog.(junk{1}).Origin_Time(first_ind(1),1:3);



% if datenum(lastNEIC)-datenum(lastEQ)<=0
%     disp('Update Events.mat first')
%     return
% end



DepthRange=[0 1000];
MagRange=[5.8 10];
Dist_Min = 30; Dist_Max = 90; RL = 30;
paranum_max = 12;

clc;
filename=[basedir '/Data/Station_Lists/Node_for_net_orfeus.txt'];
fileID = fopen(filename,'r');
ind = 1;
nets = [];
while ind~=-1
    ind = fgets(fileID);
    if ind==-1
        break;
    end
    net_tmp = split(ind);
    nets = [nets;net_tmp];
    ind = fgets(fileID);
    ind = fgets(fileID);
    for in = 1:length(net_tmp)-1
        data_source.(net_tmp{in}) = ind(1:end-1);
    end
end
fclose(fileID);

filename=[basedir '/Data/Station_Lists/Orfeus_request.txt'];
fileID = fopen(filename,'r');
ind = fgets(fileID);
NETs = [];
STAs = [];
LATs = [];
LONs = [];
ELEs = [];
STARTs = [];
ENDs = [];
while ind~=-1
    ind = fgets(fileID);
    if ind==-1
        break;
    end
    info_tmp = split(ind);
    NETs = [NETs;info_tmp(1)];
    STAs = [STAs;info_tmp(2)];
    if strcmp(info_tmp{4},'N')
        LATs = [LATs;str2num(info_tmp{3})];
    else
        LATs = [LATs;-str2num(info_tmp{3})];
    end
    if strcmp(info_tmp{6},'E')
        LONs = [LONs;str2num(info_tmp{5})];
    else
        LONs = [LONs;-str2num(info_tmp{5})];
    end
    ELEs = [ELEs;str2num(info_tmp{7})];
    STARTs = [STARTs;str2num(info_tmp{end-3})];
    if strcmp(info_tmp{end-2},'-')
        ENDs = [ENDs;2599];
    else
        ENDs = [ENDs;str2num(info_tmp{end-2})];
    end
end
fclose(fileID);

Earliest_Start_Year = max(min(STARTs),firstNEIC(1));
Latest_End_Year = min(max(ENDs),lastNEIC(1));

oldNets = fieldnames(oldNetworks);
for in=1:length(oldNets)
    oldSTAs = fieldnames(oldNetworks.(oldNets{in}));
    for is = 1:length(oldSTAs)
        Earliest_Start_Year = min([Earliest_Start_Year,max([oldNetworks.(oldNets{in}).(oldSTAs{is}).Start_Year,firstNEIC(1)])]);
        Latest_End_Year = max(Latest_End_Year,min([oldNetworks.(oldNets{in}).(oldSTAs{is}).End_Year,lastNEIC(1)]));
    end
end

Start_YR    = Earliest_Start_Year; Start_MO = 1; Start_DY = 1;
End_YR      = Latest_End_Year; % In case Latest_End_Year is not 2599!!
if Latest_End_Year == lastNEIC(1)
    End_MO = lastNEIC(2); End_DY = lastNEIC(3);
else
    End_MO = 12; End_DY = 31;
end
Mag_Min     = MagRange(1);  Mag_Max = MagRange(2); Dep_Min = DepthRange(1); Dep_Max = DepthRange(2);
Er= 6371;

days_nly    = [31 28 31 30 31 30 31 31 30 31 30 31];    % days in non-leap year
days_ly     = [31 29 31 30 31 30 31 31 30 31 30 31];    % days in leap year

ieg = 0;
iyr = 0;
for yr=Start_YR:End_YR
    
    iyr = iyr+1;
    events = NEIC_Catalog.(['NEIC_' num2str(yr)]);
    
    % Find indices of events within prescribed magnitude and depth limits
    iMagDep = find(events.Magnitude>=Mag_Min-1e-4 & events.Magnitude<=Mag_Max+1e-4 & ...
        events.Location(:,3)>=Dep_Min & events.Location(:,3)<=Dep_Max);
    
    for iei=iMagDep'
        
        if (events.Origin_Time(iei,1)>Start_YR && events.Origin_Time(iei,1)<End_YR) || ...
                (events.Origin_Time(iei,1)==Start_YR && events.Origin_Time(iei,1)~=End_YR &&...
                ((events.Origin_Time(iei,2)==Start_MO && events.Origin_Time(iei,3)>=Start_DY) || ...
                events.Origin_Time(iei,2)>Start_MO)) || ...
                (events.Origin_Time(iei,1)==End_YR && events.Origin_Time(iei,1)~=Start_YR && ...
                ((events.Origin_Time(iei,2)==End_MO && events.Origin_Time(iei,3)<=End_DY) || ...
                events.Origin_Time(iei,2)<End_MO)) || ...
                (events.Origin_Time(iei,1)==Start_YR && events.Origin_Time(iei,1)==End_YR && ...
                ((Start_MO~=End_MO && ...
                ((events.Origin_Time(iei,2)>Start_MO  && events.Origin_Time(iei,2)<End_MO) || ...
                (events.Origin_Time(iei,2)==Start_MO && events.Origin_Time(iei,3)>=Start_DY) || ...
                (events.Origin_Time(iei,2)==End_MO && events.Origin_Time(iei,3)<=End_DY))) || ...
                (Start_MO==End_MO && events.Origin_Time(iei,2)==Start_MO && ...
                events.Origin_Time(iei,3)>=Start_DY && events.Origin_Time(iei,3)<=End_DY)))
            
            ieg = ieg+1;
            
            % Calculate Julian Day
            if events.Origin_Time(iei,2)>1 % For months other than Jan
                % leap years will be any year divisible by 4
                % as the only century within seismic records used is 2000,
                % which was a leap year despite being a century year
                if rem(events.Origin_Time(iei,1),4)==0; days = days_ly;
                else days = days_nly;
                end
                % Julian Day is the sum of the days in the months
                % preceeding the event + the date of the event
                JDY = sum(days(1:events.Origin_Time(iei,2)-1))+events.Origin_Time(iei,3);
            elseif events.Origin_Time(iei,2)==1
                JDY = events.Origin_Time(iei,3);
            end
            
            % Create strings for event_ID
            YR_str = num2str(events.Origin_Time(iei,1));
            if JDY<10;           JDY_str = ['00' num2str(JDY)];
            elseif JDY<100;      JDY_str = ['0' num2str(JDY)];
            else                 JDY_str = num2str(JDY);
            end
            
            if events.Origin_Time(iei,4)<10
                HR_str = ['0' num2str(events.Origin_Time(iei,4))];
            else
                HR_str = num2str(events.Origin_Time(iei,4));
            end
            if events.Origin_Time(iei,5)<10
                MN_str = ['0' num2str(events.Origin_Time(iei,5))];
            else
                MN_str = num2str(events.Origin_Time(iei,5));
            end
            event_ID = [YR_str '.' JDY_str '.' HR_str '.' MN_str];
            Event_ID = ['EQ_' YR_str '_' JDY_str '_' HR_str '_' MN_str];
            
            Events.(Event_ID).Event_ID     	= event_ID;
            Events.(Event_ID).Year          = events.Origin_Time(iei,1);
            Events.(Event_ID).Julian_Day	= JDY;
            Events.(Event_ID).Month     	= events.Origin_Time(iei,2);
            Events.(Event_ID).Day           = events.Origin_Time(iei,3);
            Events.(Event_ID).Hour      	= events.Origin_Time(iei,4);
            Events.(Event_ID).Minute       	= events.Origin_Time(iei,5);
            Events.(Event_ID).Second      	= events.Origin_Time(iei,6)+events.Origin_Time(iei,7)/1000;
            Events.(Event_ID).Latitude   	= double(events.Location(iei,1));
            Events.(Event_ID).Longitude   	= double(events.Location(iei,2));
            Events.(Event_ID).Depth        	= double(events.Location(iei,3));
            clear x y z
            [x,y,z] = sph2cart(pi/180*events.Location(iei,2),...
                pi/180*events.Location(iei,1),...
                1-events.Location(iei,3)/Er);
            Events.(Event_ID).Globe_xyz   	= [x y z];
            Events.(Event_ID).Magnitude     = events.Magnitude(iei);
            
        end
        
    end
    
end

Eventnames = fieldnames(Events);
datenum_min = datenum([2599,1,1]);
datenum_max = datenum([0000,1,1]);
for ie = 1:length(Eventnames)
    datenum_cur = datenum(double([Events.(Eventnames{ie}).Year,Events.(Eventnames{ie}).Month,Events.(Eventnames{ie}).Day]));
    if datenum_cur<datenum_min
        datenum_min = datenum_cur;
        firstEve = [Events.(Eventnames{ie}).Year,Events.(Eventnames{ie}).Month,Events.(Eventnames{ie}).Day];
    end
    if datenum_cur>datenum_max
        datenum_max = datenum_cur;
        lastEQ = [Events.(Eventnames{ie}).Year,Events.(Eventnames{ie}).Month,Events.(Eventnames{ie}).Day];
    end
end
firstEve = double(firstEve);
lastEQ = double(lastEQ);

Networks = oldNetworks;
for in = 1:length(NETs)
    net = NETs{in};
    sta = STAs{in};
    if ~isfield(Networks,net)
        Networks.(net)=[];
    end
    if ~isfield(Networks.(net),sta)
        Networks.(net).(sta).Latitude = LATs(in);
        Networks.(net).(sta).Longitude = LONs(in);
        Networks.(net).(sta).Elevation = ELEs(in)/1000;
        [x,y,z] = sph2cart(pi/180*LONs(in),...
            pi/180*LATs(in),...
            1+ELEs(in)/Er/1000);
        Networks.(net).(sta).Globe_xyz = [x y z];
        Networks.(net).(sta).Start_Year = STARTs(in);
        Networks.(net).(sta).Start_Month = 1;
        Networks.(net).(sta).Start_Day = 1;
        Networks.(net).(sta).End_Year = ENDs(in);
        Networks.(net).(sta).End_Month = 12;
        Networks.(net).(sta).End_Day = 31;
    elseif Networks.(net).(sta).End_Year < ENDs(in)
        Networks.(net).(sta).End_Year = ENDs(in);
        Networks.(net).(sta).End_Month = 12;
        Networks.(net).(sta).End_Day = 31;
    end
    
    lastSTA = [Networks.(net).(sta).End_Year,Networks.(net).(sta).End_Month,Networks.(net).(sta).End_Day];
    
    
%     if isfield(Networks.(net).(sta),'Events_Requested')      %skip waveforms already updated through iris
%         if ~isfield(olderNetworks,net)
%             continue;
%         elseif ~isfield(olderNetworks.(net),sta)
%             continue;
%         elseif ~isfield(olderNetworks.(net).(sta),'Events_Requested')
%             continue;
%         elseif length(olderNetworks.(net).(sta).Events_Requested) ~=...
%                 length(Networks.(net).(sta).Events_Requested)
%             continue;
%         end
%     end

    if ~isfield(Networks.(net).(sta),'Events_Requested')
        if Networks.(net).(sta).Start_Year<1980
            Syear = 1980;
            Smonth = 1;
            Sday = 1;
        else
            Syear = Networks.(net).(sta).Start_Year;
            Smonth = Networks.(net).(sta).Start_Month;
            Sday = Networks.(net).(sta).Start_Day;
        end
        if datenum([Syear,Smonth,Sday])<=datenum(firstEve)
            Syear = firstEve(1);
            Smonth = firstEve(2);
            Sday = firstEve(3);
        end
        ind_sta = 1;
    else
        eveLast = Networks.(net).(sta).Events_Requested{end};
        eveY = eveLast(1:4);
        evejd = eveLast(6:8);
        eveh = eveLast(10:11);
        evem = eveLast(13:14);
        EveF = ['EQ_',eveY,'_',evejd,'_',eveh,'_',evem];
        isF = strcmp(EveF,Eventnames);
        indF = find(isF==1);
        ind_sta = indF+1;
        if indF==length(Eventnames)
            continue;
        else
            Syear = Events.(Eventnames{indF+1}).Year;
            Smonth = Events.(Eventnames{indF+1}).Month;
            Sday = Events.(Eventnames{indF+1}).Day;
            firstEQ = [Syear,Smonth,Sday];
            if datenum(lastSTA)<=datenum(firstEQ)
                continue;
            end
        end
    end
    
    if datenum(lastEQ)<=datenum(lastSTA)
        Eyear = Events.(Eventnames{end}).Year;
        Emonth = Events.(Eventnames{end}).Month;
        Eday = Events.(Eventnames{end}).Day;
    else
        Eyear = Networks.(net).(sta).End_Year;
        Emonth =Networks.(net).(sta).End_Month;
        Eday = Networks.(net).(sta).End_Day;
    end
    
    lastReq = [Eyear,Emonth,Eday];
    firstReq = [Syear,Smonth,Sday];
    NewReq = [];
    Yreqs = [];
    Mreqs = [];
    JDreqs = [];
    Dreqs = [];
    Hreqs = [];
    MNreqs = [];
    Sreqs = [];
    
    for ind_cur = ind_sta:length(Eventnames)
        Curyear = Events.(Eventnames{ind_cur}).Year;
        Curmonth = Events.(Eventnames{ind_cur}).Month;
        Curday = Events.(Eventnames{ind_cur}).Day;
        CurDate = [Curyear,Curmonth,Curday];
        if datenum(CurDate)<datenum(firstReq)
            continue;
        end
        if datenum(CurDate)>datenum(lastReq)
            break;
        end
        if Events.(Eventnames{ind_cur}).Magnitude<MagRange(1) || ...
                Events.(Eventnames{ind_cur}).Magnitude>MagRange(2)
            continue;
        end
        sta_eve_dist = distance([Networks.(net).(sta).Latitude,...
            Networks.(net).(sta).Longitude],[Events.(Eventnames{ind_cur}).Latitude,...
            Events.(Eventnames{ind_cur}).Longitude]);
        if sta_eve_dist<Dist_Min || sta_eve_dist>Dist_Max
            continue;
        end
        NewReq = [NewReq;{Events.(Eventnames{ind_cur}).Event_ID}]; 
        Yreqs = [Yreqs;Events.(Eventnames{ind_cur}).Year];
        Mreqs = [Mreqs;Events.(Eventnames{ind_cur}).Month];
        JDreqs = [JDreqs;Events.(Eventnames{ind_cur}).Julian_Day];
        Dreqs = [Dreqs;Events.(Eventnames{ind_cur}).Day];
        Hreqs = [Hreqs;Events.(Eventnames{ind_cur}).Hour];
        MNreqs = [MNreqs;Events.(Eventnames{ind_cur}).Minute];
        Sreqs = [Sreqs;floor(Events.(Eventnames{ind_cur}).Second)];
    end
    
    paranum = min(paranum_max,round(length(NewReq)/3));
    
    system(['/usr/local/bin/obspyDMT --datapath ',Proj_Dir,net,'/',...
        sta,'/EU_update',' --data_source "',data_source.(net),'" --min_epi ',...
        num2str(Dist_Min),' --max_epi ',num2str(Dist_Max),' --min_date ',...
        num2str(Syear,'%.4d'),'-',num2str(Smonth,'%.2d'),'-',num2str(Sday,'%.2d'),...
        ' --max_date ',num2str(Eyear,'%.4d'),'-',num2str(Emonth,'%.2d'),'-',...
        num2str(Eday,'%.2d'),' --preset 0 --offset ',num2str(RL*60),...
        ' --waveform_format "sac" --net "',net,'" --sta "',sta,'" --cha "BH*,HH*" ',...
        '--event_catalog NEIC_USGS --min_depth ',num2str(DepthRange(1)),...
        ' --max_depth ',num2str(DepthRange(2)),' --min_mag ',num2str(MagRange(1)),...
        ' --max_mag ',num2str(MagRange(2)),' --req_parallel --req_np ',num2str(paranum),...
        ' --parallel_process --process_np ',num2str(paranum)]);
    
    Fail = [];
    for ir = 1:length(NewReq)
        name_dir = [Proj_Dir,net,'/',sta,'/EU_update/',num2str(Yreqs(ir),'%.4d'),num2str(Mreqs(ir),'%.2d'),...
            num2str(Dreqs(ir),'%.2d'),'_',num2str(Hreqs(ir),'%.2d'),...
            num2str(MNreqs(ir),'%.2d'),num2str(Sreqs(ir),'%.2d'),'.a/processed/'];
        if ~exist(name_dir,'dir')
            Fail = [Fail;NewReq(ir)];
            continue;
        end
        files = dir([name_dir,net,'*']);
        for icp = 1:length(files)
           system(['cp ',name_dir,files(icp).name,' ',Proj_Dir,net,'/',sta,...
               '/SAC_Files/',NewReq{ir},'.',num2str(Sreqs(ir),'%.2d'),'.0000.',...
               files(icp).name,'..SAC'])
        end
    end
    
    Networks.Events_Requested = [Networks.Events_Requested;NewReq];
    save([Proj_Dir,net,'/',sta,'/EUupdateFail.mat'],'Fail');
    save([Proj_Dir,'Networks.mat'],'Networks');
    
end
save([Proj_Dir,'Networks.mat'],'Networks');
save([Proj_Dir,'Events.mat'],'Events');
end
