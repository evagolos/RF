function Auto_Prep(Project, basedir, tag, varargin)

% ********* Function Description *********
%
% This is where the raw data is input and
% prepped (phase and analysis windows
% chosen) for the deconvolution and
% migration code.
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
if ~isdir(Proj_Dir);  system(['mkdir ' Proj_Dir]);  end

% Password for encrypted SEED files
% psswrd = input('Decryption pass phrase:   ');
psswrd = '';


switch tag
    case 'neic'
        Update_NEIC(basedir);
        
    case 'request'
        Request_SAC_Files(Project, basedir);
        
    case 'Orfeus_request'
        Orfeus_Request_SAC_Files(Project, basedir);
        
    case 'unpack'
        Unpack_SAC_Files(Proj_Dir,Project, psswrd);
        
    case 'multiples'
        CheckForMultiples(Proj_Dir);
        
    case 'updateData'
        Update_NEIC(basedir);
        updateData(Project, Proj_Dir, basedir);
        %Unpack_SAC_Files(Proj_Dir,psswrd,'update');
    case 'updateEU'
        update_eu_Data(Proj_Dir, basedir)
        
    case 'arrayPick'
        nets=varargin{1};
        lowT = varargin{2};
        highT = varargin{3};
        ArrayPicker(Proj_Dir, Project, nets,lowT,highT);
        
    case 'arrayPickReformat'
        lowT = varargin{1};
        highT = varargin{2};
        ArrayPickReformat(Proj_Dir, Project,lowT,highT);
        
    case 'prep'
        if length(varargin)<=7
            NETs=varargin{1};
            Vp=varargin{3}; Vs=varargin{4};
            lowT = varargin{5};
            highT = varargin{6};
            if length(varargin)==7
                Para = varargin{7};
            else
                Para = 'single';
            end
        else disp('Wrong number of inputs to Prep Function'); return;
        end
        
        load([Proj_Dir 'Networks.mat']);
        % Choose networks to prep (if not already hardwired in)
        
        if isempty(NETs)
            nets=fieldnames(Networks); clc
            disp('Available networks: '); disp(nets);
            net=input(['Type in network you wish to prep, or ' ...
                'hit ''return'' for all:  '],'s');
            if isempty(net); NETs=nets; else NETs=cellstr(net); end
        end
        clc
        
        % Choose stations to prep (if multiple networks, all stations in
        % those networks will be prepped)
        if(length(NETs)==1 && strcmp(NETs,'all') || length(NETs)>1)
            NETs = fieldnames(Networks);
            sta=cellstr('all');
        else stas=fieldnames(Networks.(NETs{1})); clc
            if strcmp(stas,'FAKE');
                sta='FAKE';
            else
                disp('Available stations:  '); disp(stas);
                sta=input(['Type in station you wish to prep,or ' ...
                    'hit ''return'' for all:  '],'s');
            end
            if isempty(sta); sta=cellstr('all'); else sta=cellstr(sta); end
        end
        clc
        
        % Free-Surface Correction Velocities
        
        if ~strcmp(Para,'para')
            Prep_Waveforms(Proj_Dir,NETs,sta,'data',lowT,highT)
        elseif length(NETs)==1
            disp('para only for more than one station');
            Prep_Waveforms(Proj_Dir,NETs,sta,'data',lowT,highT)
        else
            npcs = min(8,length(NETs));
            sumrec = zeros(length(NETs),1);
            for in=1:length(NETs)
                stas = fieldnames(Networks.(NETs{in}));
                for is = 1:length(stas)
                    sumrec(in) = sumrec(in) + length(Networks.(NETs{in}).(stas{is}).Events_Requested);
                end
            end
            
            neach = sum(sumrec)/npcs;
            nleft = npcs;
            nnow = 0;
            nsumnow = 0;
            marksst = zeros(npcs,1);
            marksen = zeros(npcs,1);
            marksst(1) = 1;
            marksen(end) = length(NETs);
            im = 1;
            for inpc = 1:length(NETs)
                nnow = nnow + sumrec(inpc);
                nsumnow = nsumnow + sumrec(inpc);
                if nnow>4/5*neach
                    im = im+1;
                    marksen(im-1) = inpc;
                    marksst(im) = inpc+1;
                    nleft = nleft-1;
                    neach = (sum(sumrec)-nsumnow)/nleft;
                    nnow = 0;
                    if nleft == 1
                        break;
                    end
                end
            end
            
            %
            %             npcs = 16;
            %             stas = fieldnames(Networks.KO);
            %             sumrec = zeros(length(stas),1);
            %             for is=1:length(stas)
            %                 sumrec(is) = sumrec(is) + length(Networks.KO.(stas{is}).Events_Requested);
            %             end
            %
            %             neach = sum(sumrec)/npcs;
            %             nleft = npcs;
            %             nnow = 0;
            %             nsumnow = 0;
            %             marksst = zeros(npcs,1);
            %             marksen = zeros(npcs,1);
            %             marksst(1) = 1;
            %             marksen(end) = length(NETs);
            %             im = 1;
            %             for inpc = 1:length(stas)
            %                 nnow = nnow + sumrec(inpc);
            %                 nsumnow = nsumnow + sumrec(inpc);
            %                 if nnow>neach
            %                     im = im+1;
            %                     marksen(im-1) = inpc;
            %                     marksst(im) = inpc+1;
            %                     nleft = nleft-1;
            %                     neach = (sum(sumrec)-nsumnow)/nleft;
            %                     nnow = 0;
            %                     if nleft == 1
            %                         break;
            %                     end
            %                 end
            %             end
            parpool('local',npcs)
            
            parfor ip = 1:npcs
                %for ip = 2:2
                javaaddpath([basedir 'Functions/taup/lib/TauP-1.1.7.jar']);
                javaaddpath([basedir 'Functions/taup/lib/log4j-1.2.8.jar']);
                javaaddpath([basedir 'Functions/taup/lib/seisFile-1.0.1.jar']);
                Prep_Waveforms(Proj_Dir,NETs,sta,'data',lowT,highT,'para',ip,[marksst(ip),marksen(ip)]);
            end
            
            %             spmd
            %                 Prep_Waveforms(Proj_Dir,NETs,sta,FS_VpVs,'data');
            %             end
            
            delete(gcp)
            
            %             for ip = 1:npcs
            %                 nettmp = load([Proj_Dir 'Networks_2_',num2str(ip),'.mat']);
            %                 for in = marksst(ip):marksen(ip)
            %                     Networks.KO.(stas{in}) = nettmp.Networks.KO.(stas{in});
            %                 end
            %             end
            %
            for ip = 1:npcs
                nettmp = load([Proj_Dir 'Networks_',num2str(ip),'.mat']);
                for in = marksst(ip):marksen(ip)
                    Networks.(NETs{in}) = nettmp.Networks.(NETs{in});
                end
            end
            
            save([Proj_Dir 'Networks.mat'],'Networks');
            
        end
        
        
    case 'prepsynth'
        NETs=varargin{1}; Phase=varargin{2}; version=varargin{3};
        load([basedir 'Data/Projects/Synth_' Project '/Networks.mat']);
        % Choose networks to prep (if not already hardwired in)
        if isempty(NETs)
            nets=fieldnames(Networks); clc
            disp('Available networks: '); disp(nets);
            net=input(['Type in network you wish to prep, or ' ...
                'hit ''return'' for all:  '],'s');
            if isempty(net); NETs=nets; else NETs=cellstr(net); end
            NETs=nets;
        end
        clc
        
        % Choose stations to prep (if multiple networks, all stations in
        % those networks will be prepped)
        if(length(NETs)==1 && strcmp(NETs,'all') || length(NETs)>1);
            sta=cellstr('all');
        else stas=fieldnames(Networks.(NETs{1})); clc
            disp('Available stations:  '); disp(stas);
            sta=input(['Type in station you wish to prep,or ' ...
                'hit ''return'' for all:  '],'s');
            if isempty(sta); sta=cellstr('all'); else sta=cellstr(sta); end
        end
        clc
        
        Prep_Waveforms(Proj_Dir,NETs,sta,[],'synth', Phase, version);
    case 'prepsynth_test'
        NETs=varargin{1}; Phase=varargin{3}; sta=varargin{2}; version=varargin{4};
        load([basedir 'Data/Projects/Synth_' Project '/Networks.mat']);
        
        Prep_Waveforms(Proj_Dir,NETs,sta,[],'synth', Phase, version,'test');
        
end

end

function Update_NEIC(basedir)
% Update event catalogue

load([basedir '/Data/Event_Lists/NEIC_Catalog.mat']);
junk=fieldnames(NEIC_Catalog); times=NEIC_Catalog.(junk{end}).Origin_Time;
last=times(end,1:3);
eval(['save ' basedir 'Data/Event_Lists/Accessory/NEIC_Catalog_old.mat NEIC_Catalog']);


today=datevec(date); today=today(1:3);

% Compare last event in NEIC Catalogue with today's date - if there is a
% gap, update the catalogue
if last(1)<today(1) || last(1)==today(1) && last(2)<today(2) || ...
        last(1)==today(1) && last(2)==today(2) && last(3)<today(3)
    disp('Go to http://earthquake.usgs.gov/earthquakes/search/')
    disp('to request data with the following parameters:')
    disp(['   Start date:  ' num2str(last)])
    disp(['   End date:    ' num2str(today)])
    fprintf('   Magnitude:   6.0-10 \n   Depth:       0-1000km\n')
    fprintf('\nOutput as CSV, oldest First\n')
    
    fprintf('Save as Data/Event_Lists/Accessory/new.txt, plain text\n'); junk='junk';
    while ~strcmp(junk, 'done');
        junk=input('and type ''done'' when you have:   ', 's');
    end
    
    % Load the text copied from the website - if their syntax changes, will
    % need to update this!
    filename=[basedir '/Data/Event_Lists/Accessory/new.txt'];
    % import dates and times
    delimiter={',','-','T',':','Z','.'}; start=2;
    formatSpec = '%s%s%s%s%s%s%s%[^\n\r]';
    fileID=fopen(filename,'r');
    dataArray=textscan(fileID, formatSpec,'Delimiter', delimiter, 'HeaderLines', start-1, 'ReturnOnError', false);
    fclose(fileID); dates = [dataArray{1:end-1}];
    
    % import lat, lon, depth, magnitude
    delimiter=','; formatSpec= '%*s%s%s%s%s%[^\n\r]'; fileID=fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines', start-1,  'ReturnOnError', false);
    fclose(fileID); locs = [dataArray{1:end-1}];
    
    events=[dates, locs]; events=cellfun(@(x)str2double(x), events);
    clearvars filename start delimiter formatSpec fileID dataArray dates locs ans junk;
    
    % Ensure no duplicate events are added to NEIC Catalogue
    k=1; while times(end,3)==events(k,3) && times(end,4)==events(k,4); k=k+1; end
    events=events(k:end,:);   yrs=unique(events(:,1));
    
    for iy=1:length(yrs)
        yrStr=strcat('NEIC_', num2str(yrs(iy)));  inds=find(events(:,1)==yrs(iy));
        if isfield(NEIC_Catalog, yrStr)
            start=length(NEIC_Catalog.(yrStr).Magnitude);
            NEIC_Catalog.(yrStr).Origin_Time(start+1:length(inds)+start,:)=events(inds,1:7);
            NEIC_Catalog.(yrStr).Location(start+1:length(inds)+start,:)=events(inds,8:10);
            NEIC_Catalog.(yrStr).Magnitude(start+1:length(inds)+start,:)=events(inds,11);
        else
            NEIC_Catalog.(yrStr).Origin_Time=events(inds,1:7);
            NEIC_Catalog.(yrStr).Location=events(inds,8:10);
            NEIC_Catalog.(yrStr).Magnitude=events(inds,11);
        end
    end
end

eval(['save ' basedir 'Data/Event_Lists/NEIC_Catalog.mat NEIC_Catalog'])

end

function Request_SAC_Files(Project, basedir)
% Generate request files, Networks.mat and Events.mat

% Load events
load([basedir 'Data/Event_Lists/NEIC_Catalog.mat']);
% Create station list from rawstas.txt

if ~exist([basedir '/Data/Station_Lists/rawstas_',Project,'.txt'],'file')
    system(['touch ' basedir '/Data/Station_Lists/rawstas_',Project,'.txt']);
end

clc
disp('Go to http://www.iris.edu/SeismiQuery/channel.htm to find available stations')
disp('You need to tick NETWORK, STATION, channel, START TIME, END TIME, LAT, LON, ELEVATION')
disp('Fill in the following:')
disp('Channel: BH%;  Start time: 1980; End time: today; desired lat/lon limits')
disp('Ctrl+a and copy into a txt file (make sure it is plain text doc BEFORE copying in!)')
disp('Ctrl+f to replace all the OPEN with 2599-12-31 23:59:59')

fprintf(['Save as Data/Station_Lists/rawstas_',Project,'.txt, ']); junk='junk';
while ~strcmp(junk, 'done'); junk=input('and type ''done'' when you have:   ', 's'); end

clc;
filename=[basedir '/Data/Station_Lists/rawstas_',Project,'.txt'];
% Import station and network names
delimiter={'\t',' '}; formatSpec='%s%s%s%[^\n\r]';fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, ...
    'HeaderLines',5, 'ReturnOnError', false,'MultipleDelimsAsOne',1);
fclose(fileID); stas = [dataArray{1:end-1}];
% Import start and end times
formatSpec = [repmat('%*s',1,3) repmat('%s',1,3) '%*s' repmat('%s',1,3) '%[^\n\r]'];
fileID = fopen(filename,'r'); delimiter={'\t',' ','-'};
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, ...
    'HeaderLines',5, 'ReturnOnError', false,'MultipleDelimsAsOne',1);
fclose(fileID); times = [dataArray{1:end-1}]; times=cellfun(@(x)str2double(x), times);
% Import lat, lon, elevation
formatSpec = [repmat('%*s',1,7) repmat('%s',1,3) '%[^\n\r]'];
fileID = fopen(filename,'r'); delimiter = {'\t','+',' '};
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines',5, ...
    'ReturnOnError', false,'MultipleDelimsAsOne',1);
fclose(fileID); locs = [dataArray{1:end-1}]; locs=cellfun(@(x)str2double(x), locs);
% Remove stations below sea level - sometimes elevation is funky though
% so more than OBSs are removed!
%stas(locs(:,3)<0,:)=[]; times(locs(:,3)<0,:)=[]; locs(locs(:,3)<0,:)=[];

% If a station is multiply requested in .txt file, merge into single line
for iS=length(times(:,1)):-1:2
    if strcmp(stas{iS,2}, stas{iS-1,2}) && strcmp(stas{iS,1}, stas{iS-1,1})
        locs(iS,:)=[]; stas(iS,:)=[];
        if datenum(times(iS-1,1:3)) > datenum(times(iS,1:3))
            times(iS-1,1:3) = times(iS,1:3);
        end
        if datenum(times(iS-1,4:6)) < datenum(times(iS,4:6))
            times(iS-1,4:6) = times(iS,4:6);
        end
        times(iS,:)=[];
    end
end

ind = 1;
HHinfo.HHStalist = {};
HHinfo.HHNetlist = {};
HHinfo.HHBHtime = [];
HHinfo.HHHHtime = [];

while ind<size(stas,1)
    rep = find(strcmp(stas{ind,2},stas(ind+1:end,2)));
    if ~isempty(rep)
        mark = [];
        for iS = 1:length(rep)
            if strcmp(stas{ind,1},stas{rep(iS)+ind,1})
                if strcmp(stas{ind,3}(1:2),'HH')
                    deind = ind;
                    prind = rep(iS)+ind;
                elseif strcmp(stas{rep(iS)+ind,3}(1:2),'HH')
                    deind = rep(iS)+ind;
                    prind = ind;
                else
                    continue;
                end
                HHinfo.HHStalist = [HHinfo.HHStalist,stas{ind,2}];
                HHinfo.HHNetlist = [HHinfo.HHNetlist,stas{ind,1}];
                HHinfo.HHBHtime = [HHinfo.HHBHtime;times(prind,:)];
                HHinfo.HHHHtime = [HHinfo.HHHHtime;times(deind,:)];
                if datenum(times(prind,1:3))>datenum(times(deind,1:3))
                    times(prind,1:3) = times(deind,1:3);
                end
                if datenum(times(prind,4:6))<datenum(times(deind,4:6))
                    times(prind,4:6) = times(deind,4:6);
                end
                mark = [mark,deind];
            end
        end
        if ~isempty(mark)
            times(mark,:) = [];
            stas(mark,:)=[];
            locs(mark,:)=[];
        end
    end
    ind = ind+1;
    continue;
end


for iS=1:length(times(:,1))
    stations(iS).Network=stas(iS,1);  stations(iS).Station = stas{iS,2};
    stations(iS).Latitude=locs(iS,1); stations(iS).Longitude = locs(iS,2);
    stations(iS).Elevation=locs(iS,3); stations(iS).Start_Year=times(iS,1);
    stations(iS).Start_Month=times(iS,2); stations(iS).Start_Day=times(iS,3);
    stations(iS).End_Year=times(iS,4); stations(iS).End_Month=times(iS,5);
    stations(iS).End_Day=times(iS,6); stations(iS).Channel = stas{iS,3}(1:2);
end
stations=stations';
clearvars iS locs stas ans dataArray delimiter fileID filename formatSpec iS
eval(['save ' basedir 'Data/Station_Lists/' Project datestr(date, 'dd.mm.yy') '.mat stations']);
system(['cp ' basedir 'Data/Station_Lists/rawstas_',Project,'.txt ' basedir 'Data/Station_Lists/' Project datestr(date, 'dd.mm.yy') '.txt']);

[Networks,Events]=MakeRequests(basedir,Project,stations,NEIC_Catalog,HHinfo,'normal');

Proj_Dir = [basedir 'Data/Projects/' Project '/'];
save([Proj_Dir '/Events.mat'],'Events')
save([Proj_Dir '/Networks.mat'],'Networks')

end

function Orfeus_Request_SAC_Files(Project, basedir)
% Generate request files, Networks.mat and Events.mat

% Load events
load([basedir 'Data/Event_Lists/NEIC_Catalog.mat']);
% Create station list from rawstas.txt
clc
fprintf(['Go to http://145.23.252.222/eida/webdc3/ to find available stations' ...
    '\nSelect an area, and then copy/paste into an excel spreadsheet.'])

fprintf('\nSave as Data/Station_Lists/orf_rawstas.csv, '); junk='junk';
while ~strcmp(junk, 'done'); junk=input('\nand type ''done'' when you have:   ', 's'); end
clc;

filename=[basedir '/Data/Station_Lists/orf_rawstas.csv'];
% Import station and network names
delimiter=','; formatSpec=['%*s' repmat('%s',1,4) '%*[^\n\r]'];fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID); data = [dataArray{:}]; data(strcmp(data(:,1),''),:)=[];
nets=unique(data(:,1)); n=0; locs=cellfun(@(x)str2double(x),data(:,3:4));
for iN=1:length(nets)
    fprintf(['\n\nWorking on network %s - take average network ' ...
        'location (v.s. station by station) (y/n)?'],nets{iN}); av=input('  ','s');
    switch av
        case 'n';  stas=unique(data(strcmp(data(:,1),nets{iN}),2));
        case 'y';  stas={'av'};
    end
    
    for iS=1:length(stas)
        n=n+1;
        switch av
            case 'y'
                str=[' the network ' nets{iN}];
                inds=find(strcmp(data(:,1),nets{iN}));
            case 'n'
                str=[' station ' stas{iS} ' (network ' nets{iN} ')'];
                inds=find(strcmp(data(:,1),nets{iN}),1)+iS-1;
        end
        
        fprintf('\nInput start/end times for%s:\n',str)
        starttime=strsplit(input('   Start time (d/m/y): ','s'),'/');
        starttime=cellfun(@(x)str2double(x),starttime);
        endtime=  strsplit(input('   End time (d/m/y): ','s'),'/');
        endtime=cellfun(@(x)str2double(x),endtime);
        
        rep={'starttime','endtime'};
        for j=rep; eval(['junk=' j{1} '(3);']);
            if(junk<1900 && junk>50); junk=junk+1900; % turn yy into yyyy
            elseif(junk<1900 && junk<50); junk=junk+2000;
            end
            eval([j{1} '(3)=junk;']);
        end
        
        stations(n).Network=nets{iN}; stations(n).Station=stas{iS};
        stations(n).Latitude=mean(locs(inds,1));
        stations(n).Longitude=mean(locs(inds,2)); stations(n).Elevation=0;
        stations(n).Start_Day=starttime(1);  stations(n).Start_Month=starttime(2);
        stations(n).Start_Year=starttime(3); stations(n).End_Day=endtime(1);
        stations(n).End_Month=endtime(2);    stations(n).End_Year=endtime(3);
        
    end
end
fprintf('\n\n');
stations=stations';
clearvars iS locs stas ans dataArray delimiter fileID filename formatSpec iS
eval(['save ' basedir 'Data/Station_Lists/' Project datestr(date, 'dd.mm.yy') '.mat stations']);
system(['cp ' basedir 'Data/Station_Lists/rawstas.txt ' basedir 'Data/Station_Lists/' Project datestr(date, 'dd.mm.yy') '.txt']);

[Networks,Events]=MakeRequests(basedir,Project,stations,NEIC_Catalog,'orf');

Proj_Dir = [basedir 'Data/Projects/' Project '/'];
save([Proj_Dir '/Events.mat'],'Events')
save([Proj_Dir '/Networks.mat'],'Networks')

end

function [Networks,Events]=MakeRequests(basedir,Project,stations,NEIC_Catalog,HHinfo,tag)

% Contact details
Contact.NAME='Eva Golos'; Contact.INST='University of Wisconsin - Madison';
Contact.STREET='Dept. of Geoscience'; Contact.ADDR='Madison, WI 53706';
Contact.EMAIL='golosdata@gmail.com'; Contact.PHONE=''; Contact.FAX='';
Contact.MEDIUM='0 Electronic (ftp-able)'; Contact.ALTMEDIUM='0 Electronic (ftp-able)';


Er=6371; nstations = length(stations);
for is=1:nstations
    NET = stations(is).Network;    STA = stations(is).Station;
    % bug fixed by Emily 11/11/13
    % Matlab is confused by strings that start with numbers, so
    % add prefix to any station names starting with numbers
    NET = char(NET); STA = char(STA);
    if(str2num(STA(1))+1); STA=['bodge' STA];   end
    if(str2num(NET(1))+1); NET=['bodge' NET]; end
    Networks.(NET).(STA).Latitude       = stations(is).Latitude;
    Networks.(NET).(STA).Longitude      = stations(is).Longitude;
    
    Networks.(NET).(STA).Elevation      = stations(is).Elevation/1000; % Convert elevation from m to km
    Networks.(NET).(STA).Channel        = stations(is).Channel;
    
    if ~isempty(find(strcmp(STA,HHinfo.HHStalist) & strcmp(NET,HHinfo.HHNetlist))) %#ok<EFIND>
        Networks.(NET).(STA).Channel = 'BHHH';
    end
    
    clear x y z
    [x,y,z] = sph2cart(pi/180*stations(is).Longitude,...
        pi/180*stations(is).Latitude,...
        1+Networks.(NET).(STA).Elevation/Er);
    Networks.(NET).(STA).Globe_xyz      = [x y z];
    Networks.(NET).(STA).Start_Year     = stations(is).Start_Year;
    Networks.(NET).(STA).Start_Month    = stations(is).Start_Month;
    Networks.(NET).(STA).Start_Day      = stations(is).Start_Day;
    if isfield(stations(is),'End_Year')
        Networks.(NET).(STA).End_Year	= stations(is).End_Year;
        Networks.(NET).(STA).End_Month  = stations(is).End_Month;
        Networks.(NET).(STA).End_Day    = stations(is).End_Day;
    else
        Networks.(NET).(STA).End_Year	= 2999; % Assume permanent station
        Networks.(NET).(STA).End_Month  = 12;
        Networks.(NET).(STA).End_Day    = 31;
    end
    Networks.(NET).(STA).Events_Requested = {};
    
    % Find range of years over which stations were active
    if is==1
        Earliest_Start_Year	= Networks.(NET).(STA).Start_Year;
        Latest_End_Year   	= Networks.(NET).(STA).End_Year;
    else
        Earliest_Start_Year	= min([Earliest_Start_Year Networks.(NET).(STA).Start_Year]);
        Latest_End_Year  	= max([Latest_End_Year Networks.(NET).(STA).End_Year]);
    end
end

if Earliest_Start_Year<1980;  Earliest_Start_Year = 1980; end

days_nly    = [31 28 31 30 31 30 31 31 30 31 30 31];    % days in non-leap year
days_ly     = [31 29 31 30 31 30 31 31 30 31 30 31];    % days in leap year

% Grab values from Request_Data
junk=fieldnames(NEIC_Catalog); junk=NEIC_Catalog.(junk{end});
Start_YR    = Earliest_Start_Year; Start_MO = 1; Start_DY = 1;
End_YR      = min([junk.Origin_Time(end,1) Latest_End_Year]); % In case Latest_End_Year is not 2599!!
End_MO      = junk.Origin_Time(end,2); End_DY = junk.Origin_Time(end,3);
Mag_Min     = 6.0;  Mag_Max = 10; Dep_Min = 0; Dep_Max = 1000;
n_yrs = End_YR-Start_YR+1;


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

% Filter events and stations
%  This function eliminates event-station combinations that are  outside of the
% specified distance range and writes out the text file that is the breqfast email
% request.  Copy the text from the file for each network and paste it into an
% email addressed to:  breq_fast@iris.washington.edu
Request_Data.DistanceRange=[30 90]; Request_Data.DepthRange=[0 1000];
Request_Data.MagRange=[6.0 10]; Request_Data.RecordLength=30;
Proj_Dir        = [basedir 'Data/Projects/' Project '/'];
EQs = fieldnames(Events); NETs = fieldnames(Networks);
Dist_Min = 30; Dist_Max = 90; RL = 30;

for in=1:length(NETs)
    NET     = NETs{in};
    if length(NET)>5 && strcmp(NET(1:5),'bodge'); net=NET(6:end); else net=NET; end
    disp(['Working on Network ' net])
    STAs    = fieldnames(Networks.(NET));
    
    if isdir([Proj_Dir NET])==0
        mkdir([Proj_Dir NET]); mkdir([Proj_Dir NET '/Breqfast_Requests'])
        mkdir([Proj_Dir NET '/Seed_Files'])
    end
    
    % Write out a data log text file
    if strcmp(tag,'normal')
        fid_1 = fopen([Proj_Dir NET '/Breqfast_Requests/' net '_Data_Request_Log.txt'],'w');
        request_email = [Proj_Dir NET '/Breqfast_Requests/' NET '_Request.txt'];
    elseif strcmp(tag,'update')
        fid_1 = fopen([Proj_Dir NET '/Breqfast_Requests/UPDATE_' net '_Data_Request_Log.txt'],'w');
        request_email = [Proj_Dir NET '/Breqfast_Requests/UPDATE_' NET '_Request.txt'];
    elseif strcmp(tag,'orf')
        fid_1 = fopen([Proj_Dir NET '/Breqfast_Requests/' net '_Data_Request_Log.txt'],'w');
    end
    fprintf(fid_1,'\n ***** %s Data Request Log *****\n\n',net);
    datetime = datestr(clock);
    for i=1:length(datetime);
        if strcmp(datetime(i),':'); indx = i-2; break; end
    end
    todaydate = datetime(1:indx-2); time = datetime(indx:end);
    fprintf(fid_1,' Date:           %s\n',todaydate);
    fprintf(fid_1,' Time:           %s\n',time);
    fprintf(fid_1,' Requester:      %s\n\n',Contact.NAME);
    fprintf(fid_1,' Distance Range:  %g-%g deg\n',Dist_Min,Dist_Max);
    fprintf(fid_1,' Depth Range:     %g-%g km\n',Dep_Min,Dep_Max);
    fprintf(fid_1,' Record Length:   %g min\n\n',RL);
    
    first_year  = Events.(EQs{1}).Year;
    last_year   = Events.(EQs{end}).Year;
    % Bug fixed for requesting from ORFEUS - Emily 18.10.13
    % (space removed between ' and . below)
    if ~strcmp(tag,'orf')
        fid_2 = fopen(request_email,'w');
        fprintf(fid_2,'.NAME %s\n',Contact.NAME);
        fprintf(fid_2,'.LABEL %s_%s_%g_%g\n',Project,net,first_year,last_year);
        fprintf(fid_2,'.INST %s\n',Contact.INST);
        fprintf(fid_2,'.STREET %s\n',Contact.STREET);
        fprintf(fid_2,'.ADDR %s\n',Contact.ADDR);
        fprintf(fid_2,'.EMAIL %s\n',Contact.EMAIL);
        fprintf(fid_2,'.PHONE %s\n',Contact.PHONE);
        fprintf(fid_2,'.FAX %s\n',Contact.FAX);
        fprintf(fid_2,'.MEDIUM %s\n',Contact.MEDIUM);
        fprintf(fid_2,'.ALTMEDIUM %s\n',Contact.ALTMEDIUM);
        fprintf(fid_2,'.END\n');
    end
    
    for is=1:length(STAs)
        STA     = STAs{is};
        if strcmp(tag,'orf');
            if is>1; fclose(fid_2); end
            request_email = [Proj_Dir NET '/Breqfast_Requests/' NET '_' STA '_Request.csv'];
            fid_2 = fopen(request_email,'w');
        end
        
        for ie=1:length(EQs)
            event   = Events.(EQs{ie});
            ELAT    = event.Latitude;
            ELON    = event.Longitude;
            % Only request data if the event is after the station's start time.
            if Networks.(NET).(STA).Start_Year<event.Year || ...
                    (Networks.(NET).(STA).Start_Year==event.Year && ...
                    Networks.(NET).(STA).Start_Month<=event.Month)
                % Bug fixed Emily 11/11/13:
                % (do not include if after end date)
                if Networks.(NET).(STA).End_Year>event.Year || ...
                        (Networks.(NET).(STA).End_Year==event.Year && ...
                        Networks.(NET).(STA).End_Month>=event.Month)
                    
                    SLAT = Networks.(NET).(STA).Latitude;
                    SLON = Networks.(NET).(STA).Longitude;
                    DIST = distance('gc',ELAT,ELON,SLAT,SLON);
                    if DIST>=Dist_Min && DIST<=Dist_Max
                        
                        if length(STA)>5 && strcmp(STA(1:5), 'bodge')
                            sta=STA(6:end);
                        elseif length(STA)==2; sta = [STA '  '];
                        elseif length(STA)==3; sta = [STA ' '];
                        else   sta = STA;
                        end
                        
                        % Get start time for requested waveform
                        YRa = int32(event.Year);
                        % Get end time for requested waveform
                        if int16(YRa/4)~=YRa/4
                            days = days_nly;
                        elseif int16(YRa/4)==YRa/4
                            days = days_ly;
                        end
                        MOa = int8(event.Month);
                        DYa = int8(event.Day);
                        HRa = int8(event.Hour);
                        MNa = int8(event.Minute);
                        SCa = event.Second;
                        
                        while SCa>=60; SCa = SCa-60; MNa = MNa+1; end
                        while MNa>=60; MNa = MNa-60; HRa = HRa+1; end
                        while HRa>=24; HRa = HRa-24; DYa = DYa+1; end
                        while DYa>days(MOa); DYa = DYa-days(MOa); MOa = MOa+1; end
                        while MOa>12; MOa = MOa-12; YRa = YRa+1; end
                        
                        % Get strings for start time
                        YR1 = num2str(YRa);
                        if MOa<10; MO1 = [' ' num2str(MOa)]; else MO1 = num2str(MOa); end
                        if DYa<10; DY1 = [' ' num2str(DYa)]; else DY1 = num2str(DYa); end
                        if HRa<10; HR1 = [' ' num2str(HRa)]; else HR1 = num2str(HRa); end
                        if MNa<10; MN1 = [' ' num2str(MNa)]; else MN1 = num2str(MNa); end
                        if SCa<10;
                            if int8(SCa)~=SCa; SEC1 = [' ' num2str(SCa)]; elseif int8(SCa)==SCa; SEC1 = [' ' num2str(SCa) '.0']; end
                        else
                            if int8(SCa)~=SCa; SEC1 = num2str(SCa); elseif int8(SCa)==SCa; SEC1 = [num2str(SCa) '.0']; end
                        end
                        
                        HRadd = 0;  DYadd = 0;	MOadd = 0;  YRadd = 0;
                        
                        % Check to see if the hour, day, month, or
                        % year of the end time need to be changed
                        SCb = SCa;
                        MNb = MNa+RL;
                        while MNb>=60; MNb = MNb-60; HRadd = HRadd+1; end
                        HRb = HRa+HRadd;
                        while HRb>=24; HRb = HRb-24; DYadd = DYadd+1; end
                        DYb = DYa+DYadd;
                        while DYb>days(MOa); DYb = DYb-days(MOa); MOadd = MOadd+1; end
                        MOb = MOa+MOadd;
                        while MOb>12; MOb = MOb-12; YRadd = YRadd+1; end
                        YRb = YRa+YRadd;
                        
                        % Get strings for end time
                        YR2 = num2str(YRb);
                        if MOb<10; MO2 = [' ' num2str(MOb)]; else MO2 = num2str(MOb); end
                        if DYb<10; DY2 = [' ' num2str(DYb)]; else DY2 = num2str(DYb); end
                        if HRb<10; HR2 = [' ' num2str(HRb)]; else HR2 = num2str(HRb); end
                        if MNb<10; MN2 = [' ' num2str(MNb)]; else MN2 = num2str(MNb); end
                        if SCa<10
                            if int8(SCb)~=SCb; SEC2 = [' ' num2str(SCb)]; elseif int8(SCb)==SCb; SEC2 = [' ' num2str(SCb) '.0']; end
                        else
                            if int8(SCb)~=SCb; SEC2 = num2str(SCb); elseif int8(SCa)==SCa; SEC2 = [num2str(SCb) '.0']; end
                        end
                        
                        
                        if strcmp(Networks.(NET).(STA).Channel,'BHHH')
                            indtmp = find(strcmp(STA,HHinfo.HHStalist) & strcmp(NET,HHinfo.HHNetlist));
                            date_eq = datenum(double([YRa,MOa,DYa]));
                            if date_eq>=datenum(HHinfo.HHBHtime(indtmp,1:3)) &&...
                                    date_eq<=datenum(HHinfo.HHBHtime(indtmp,4:6))
                                CH = 'BH';
                            elseif date_eq>=datenum(HHinfo.HHHHtime(indtmp,1:3)) &&...
                                    date_eq<=datenum(HHinfo.HHHHtime(indtmp,4:6))
                                CH = 'HH';
                            else
                                CH = 'BH';
                            end
                        else
                            CH = Networks.(NET).(STA).Channel;
                        end
                        
                        % Write out a line for each event-station pair
                        if strcmp(tag,'orf')
                            % want csv file of format time,lat,lon,depth
                            % e.g. 2011-03-11T05:46:23;38.23;142.53;15
                            fprintf(fid_2,'%g-%g-%gT%g:%g:%g;%g;%g;%g\n',...
                                YRa,MOa,DYa,HRa,MNa,SCa,ELAT,ELON,...
                                event.Depth);
                        else
                            fprintf(fid_2,'\n%s %s %s  %s %s %s %s %s  %s  %s %s %s %s %s  1    %s?',...
                                sta,net,YR1,MO1,DY1,HR1,MN1,SEC1,YR2,MO2,DY2,HR2,MN2,SEC2,CH);
                        end
                        
                        Networks.(NET).(STA).Events_Requested = [Networks.(NET).(STA).Events_Requested; {event.Event_ID}];
                        
                    end
                end
            end
        end
    end
    fprintf(fid_1,' Stations:\n');
    for is=1:length(STAs)
        STA = STAs{is};
        fprintf(fid_1,'  %s  (%g events)\n',STA,length(Networks.(NET).(STA).Events_Requested));
        % Make SAC File and RF_Data Directories
        if isdir([Proj_Dir  NET '/' STA])==0 ; mkdir([Proj_Dir  NET '/' STA]); end
        Networks.(NET).(STA).Events_Prepped = {};
    end
    fclose(fid_1);
    fclose(fid_2);
    
    if strcmp(tag,'update');
        if ~isdir([Proj_Dir 'RequestUpdates_' datestr(date,'dd.mm.yy')]);
            mkdir([Proj_Dir '/RequestUpdates_' datestr(date,'dd.mm.yy')]);
        end
        system(['cp ' request_email ' ' Proj_Dir 'RequestUpdates_' ...
            datestr(date,'dd.mm.yy') '/']);
    end
end


end

function Unpack_SAC_Files(Proj_Dir,Project, psswrd, varargin)
% Unpack SEED files and move to correct folders

if isempty(varargin);
    TMP=load([Proj_Dir 'Networks.mat']); Networks = TMP.Networks; clear TMP;
else TMP=load([Proj_Dir 'UpdateNetworks.mat']); Networks=TMP.Networks; clear TMP;
end
[~,funct]=system(['echo ' Proj_Dir ' | awk -F/Data/Project ''{print $1}''']);
funct=[funct(1:end-1) '/Functions/'];
NETs = fieldnames(Networks);
for in=1:length(NETs)
    NET = NETs{in};
    if length(NET)>5 && strcmp(NET(1:5),'bodge'); net=NET(6:end); else net=NET; end
    disp(['Unpacking seed file(s) for ' net ])
    STAs = fieldnames(Networks.(NET));
    
    cd([Proj_Dir NET '/Seed_Files'])
    %     delete('*.SAC')
    temp    = dir; isd     = 0;  seed_files = {};
    for it=1:length(temp)
        if ~strcmp(temp(it).name(1),'.') && ~temp(it).isdir && ...
                length(temp(it).name)>=6  && ~strcmp(temp(it).name(1:6),'rdseed') && ...
                strcmp(temp(it).name(1:length(Project)),Project)
            isd = isd+1;
            seed_files{isd,1}	= temp(it).name;
        end
    end
    
    if ~isempty(seed_files)
        if ~isempty(varargin) % i.e. want to update, not unpack everything
            disp('Your SEED files:'); disp(seed_files); seed_files={};
            seed_files{1}=input('\n\nCopy in the SEED file you want to unpack: ','s');
        end
        label = cell(1,length(seed_files));
        for isd=1:length(seed_files) % Dataless format seems different -- might need to revist this in future. For MOOS example, had to say isd=2.
            Seed_File = seed_files{isd};
            if ~strcmp(Seed_File(end-8:end),'.dataless')
            label{isd} = Seed_File(1:end-10);
            
            if ~strcmp(Seed_File(end-3:end),'nssl')
                system(['tar -vxf ',Seed_File,' 2> /dev/null 1> /dev/null']);
                datalesstar= ls('*.tar.dataless');
                datalesstar= datalesstar(1:end-1); % b/c the last character is a return symbol.
                system(['tar -vxf ',datalesstar,' 2> /dev/null 1> /dev/null']); % Dataless doesn't like to appear in Seed_File.
                tarlabel= datalesstar(1:end-13);
                datalessfile= ['IRISDMC-',Project,'_',net,'_1980_2023.dataless'];
                system(['cp ',tarlabel,'/',datalessfile,' ',label{isd}]);
                cd(['./',label{isd}]);
                for is=1:length(STAs)
                    STA = STAs{is};
                    if length(STA)>5 && strcmp(STA(1:5),'bodge'); sta=STA(6:end); else sta=STA; end
                    mseedfile = [sta,'.',net,'.mseed'];
                    %datalessfile= ['IRISDMC-',Project,'_',net,'_2006_2009.dataless']; %datalessfile = ['IRISDMC-',sta,'.',net,'.dataless'];
                    %metafile= ['IRISDMC-',sta,'.',net,'.meta'];
                    if exist(mseedfile,'file') && exist(datalessfile,'file')
                        %system([funct '/Users/evagolos/Research/ReceiverFunctions/lib/mseed2sac-main/mseed2sac ' mseedfile ' -m ' metafile]);
                        system([funct 'rdseedv5.3.1/rdseed.mac.x86_64 -d -z 3 -o 1 -f ' mseedfile ...
                            ' -g ',datalessfile,' 2> /dev/null 1> /dev/null']);
                        %%%                         system([funct 'rdseedv5.3.1/rdseed.mac.x86_64 -d -o 1 -f ' mseedfile ...
                        %%%                             ' -g ',datalessfile,' 2> /dev/null 1> /dev/null']);
                    end
                    cd('../');
                    if ~isdir(['../' STA]); mkdir(['../' STA]); end
                    if ~isempty(dir((['./',label{isd},'/*.' net '.' sta '.*.SAC'])))
                        if ~isdir(['../' STA '/SAC_Files'])
                            mkdir(['../' STA '/SAC_Files'])
                        end
                        cd(['../' STA '/SAC_Files'])
                        %delete('*.SAC')
                        [errAll,result] = system(['mv ../../Seed_Files/' label{isd} '/*.' net '.' sta '.*.SAC .']);
                        if errAll~=0
                            [errAlln,result] = system(['find ../../Seed_Files/' label{isd} '/ -name "*.' net '.' sta '.*.SAC" -exec mv ''{}'' ./ \;']);
                            if errAlln~=0
                                Current_Date = clock;
                                for yr=1980:Current_Date(1)
                                    [errYR,result] = system(['find ../../Seed_Files/' label{isd} '/ -name "', num2str(yr) ,'.*.' net '.' sta '.*.SAC" -exec mv ''{}'' ./ \;']);
                                    %[errYR,result] = system(['mv ../../Seed_Files/' label{isd} '/' num2str(yr) '.*.' net '.' sta '.*.SAC .']);
                                end
                            end
                        end
                        cd ../../Seed_Files
                    end
                    cd(['./',label{isd}]);
                end
                cd(['../']);
                % Changed by Junlin Apr 30 2020
            else
                system([funct '/openssl enc -d -des-cbc -salt -in ' Seed_File ...
                    ' -out ' Seed_File(1:end-8) ' -pass pass:' char(psswrd)]);
                system([funct '/rdseedv5.3.1/rdseed.mac.x86_64 -f ' Seed_File(1:end-8) ...
                    ' -d 2> /dev/null 1> /dev/null']);
            end
        end % End if clause for non-dataless files.
        end % End loop for all seed files in current network.
        
        disp(['Moving SAC files for ' net '...'])
        for is=1:length(STAs)
            STA = STAs{is};
            if length(STA)>5 && strcmp(STA(1:5),'bodge'); sta=STA(6:end); else sta=STA; end
            if ~isdir(['../' STA]); mkdir(['../' STA]); end
            for isd = 1:length(label)
                if ~isempty(dir((['./',label{isd},'/*.' net '.' sta '.*.SAC'])))
                    if ~isdir(['../' STA '/SAC_Files'])
                        mkdir(['../' STA '/SAC_Files'])
                    end
                    cd(['../' STA '/SAC_Files'])
                    %delete('*.SAC')
                    [errAll,result] = system(['find ../../Seed_Files/' label{isd} '/ -name "*.' net '.' sta '.*.SAC" -exec mv ''{}'' ./ \;']);
                    if errAll~=0
                        Current_Date = clock;
                        for yr=1980:Current_Date(1)
                            [errYR,result] = system(['find ../../Seed_Files/' label{isd} '/ -name "', num2str(yr) ,'.*.' net '.' sta '.*.SAC" -exec mv ''{}'' ./ \;']);
                            %[errYR,result] = system(['mv ../../Seed_Files/' label{isd} '/' num2str(yr) '.*.' net '.' sta '.*.SAC .']);
                        end
                    end
                    cd ../../Seed_Files 
                end
            end
        end % end stations.
        for isd = 1:length(label)
            system(['rm -r ./',label{isd}]);
        end
        %end % if clause for non-dataless.
    end % end if there are seed files.
    cd ../../../../../
    
    disp(['     ' net ' done!'])
end

disp('All unpacked!')

end

function CheckForMultiples(Proj_Dir)
% Check for multiples and missing files, and edit Networks.mat and delete
% folders accordingly

TMP=load([Proj_Dir 'Networks.mat']); Networks = TMP.Networks; clear TMP;
NETs = fieldnames(Networks);

save([Proj_Dir '/Networks_' datestr(date, 'dd.mm.yy') '.mat'], 'Networks');
filename=['Multiples_Log_' datestr(date,'dd.mm.yy') '.txt'];
fid=fopen([Proj_Dir filename],'a');


for in=1:length(NETs)
    NET = NETs{in};
    if length(NET)>5 && strcmp(NET(1:5),'bodge'); net=NET(6:end); else net=NET; end
    disp(['Checking for multiple sites for ' net ])
    fprintf(fid,'Checking for multiple sites for %s',net);
    STAs = fieldnames(Networks.(NET));
    for is=1:length(STAs)
        STA = STAs{is};
        if length(STA)>5 && strcmp(STA(1:5),'bodge'); sta=STA(6:end); else sta=STA; end
        
        if isdir([Proj_Dir '/' NET '/' STA '/SAC_Files'])
            cur_dir=pwd; cd([Proj_Dir '/' NET '/' STA '/SAC_Files'])
            
            % Remove if only have Z components
            %[~,numAll]=system('ls -l *SAC | wc -l'); [~,numZ]=system('ls -l *BHZ*SAC | wc -l');
            [~,numAll]=system('ls -l *SAC | wc -l'); [~,numZ]=system('ls -l | grep "BHZ\|HHZ" | grep SAC | wc -l');  %grep statements added by Nick
            
            if str2double(numAll)==str2double(numZ);
                disp([net '.' sta ' only has Z components...  ']);
                junk=input('Do you want to remove it (y/n)?  ','s');
                if strcmp(junk,'y')
                    system(['rm -r ' Proj_Dir NET '/' STA]);
                    Networks.(NET)=rmfield(Networks.(NET),STA);
                    disp(['    Station ' sta ' removed.']);
                    fprintf(fid,'\n    Station %s removed.',sta);
                end
            elseif isnan(str2double(numZ))
                disp([net '.' sta ' has no Z components...  ']);
                junk=input('Do you want to remove it (y/n)?  ','s');
                if strcmp(junk,'y')
                    system(['rm -r ' Proj_Dir NET '/' STA]);
                    Networks.(NET)=rmfield(Networks.(NET),STA);
                    disp(['    Station ' sta ' removed.']);
                    fprintf(fid,'\n    Station %s removed.',sta);
                end
            else
                %[~,out]=system('for m in `ls -l | grep BHZ | grep SAC`; do echo $m | awk -F. ''{print $9}''; done');
                [~,out]=system('ls -l | grep "BHZ\|HHZ" | grep SAC | awk -F. ''{print $9}'' '); %grep statements added by Nick
                out2=splitstring(out); sites=unique(out2);
                if length(sites)==1; continue
                else
                    for ic=2:length(sites)
                        disp([ STA sites{ic}])
                        mkdir(['../../' STA sites{ic}])
                        mkdir(['../../' STA sites{ic} '/SAC_Files'])
                        %system(['mv *' sta '.' sites{ic} '.*SAC ../../' STA sites{ic} '/SAC_Files/']);
                        system(['for m in `ls | grep .' sta '.' sites{ic} '. | grep SAC`; do mv $m ../../' STA sites{ic} '/SAC_Files/ ; done']);
                        disp(['    Station ' sta '.' sites{ic} ' created!']);
                        fprintf(fid,'\n    Station %s.%s created!',sta,sites{ic});
                        Networks.(NET).([STA sites{ic}])=Networks.(NET).(STA);
                    end
                end
            end
            cd(cur_dir)
        else disp([net '.' sta ' does not have any SAC files...  ']);
            junk=input('Do you want to remove it (y/n)?  ','s');
            if strcmp(junk,'y')
                Networks.(NET)=rmfield(Networks.(NET),STA);  % Remove if empty
                system(['rm -r ' Proj_Dir NET '/' STA]);
                disp(['    Station ' sta ' removed.']);
                fprintf(fid,'\n    Station %s removed.',sta);
            end
        end
    end
end

fprintf(fid,'\n\n\n');
for n=1:length(NETs)
    if isempty(fieldnames(Networks.(NETs{n})));
        disp([NETs{n} ' does not have any SAC files...  ']);
        junk=input('Do you want to remove it (y/n)?  ','s');
        if strcmp(junk,'y')
            Networks=rmfield(Networks,(NETs{n}));
            system(['rm -rf ' Proj_Dir '/' NETs{n}]);
            fprintf(fid,'\nNETWORK %s removed.',NETs{n});
        end
    end
end

save([Proj_Dir '/Networks.mat'], 'Networks')

fclose('all');

end

function updateData(Project, Proj_Dir, basedir)

% Load events, Events.mat and Networks.mat
load([basedir 'Data/Event_Lists/NEIC_Catalog.mat']);
TMP=load([Proj_Dir 'Events.mat']); oldEvents = TMP.Events; clear TMP
TMP=load([Proj_Dir 'Networks.mat']); oldNetworks = TMP.Networks; clear TMP
system(['cp ' Proj_Dir 'Networks.mat ' Proj_Dir 'Networks_' datestr(date,'dd.mm.yy') '.mat']);
system(['cp ' Proj_Dir 'Events.mat ' Proj_Dir 'Networks_' datestr(date,'dd.mm.yy') '.mat']);

% Check last events
junk=fieldnames(oldEvents); junk=oldEvents.(junk{end}); lastEQ=[junk.Year junk.Month junk.Day];
junk=fieldnames(NEIC_Catalog); lastNEIC=NEIC_Catalog.(junk{end}).Origin_Time(end,1:3);
if datenum(lastNEIC)-datenum(lastEQ)<=0; return; end % If no more events since last update

% Create station list from rawstas.txt
clc
disp('Go to http://www.iris.edu/SeismiQuery/channel.htm to find available stations')
disp('Select ALL available stations - overlap with previous requests will be removed')
disp('You need to tick NETWORK, STATION, START TIME, END TIME, LAT, LON, ELEVATION')
disp('Fill in the following:')
disp('Channel: BH%;  Start time: 1980; End time: today; desired lat/lon limits')
disp('Ctrl+a and copy into a txt file (make sure it is plain text doc BEFORE copying in!)')

fprintf('Save as Data/Station_Lists/rawstas.txt, '); junk='junk';
while ~strcmp(junk, 'done'); junk=input('and type ''done'' when you have:   ', 's'); end

clc;
filename=[basedir '/Data/Station_Lists/rawstas.txt'];
% Import station and network names
delimiter={'\t',' '}; formatSpec='%s%s%[^\n\r]';fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, ...
    'HeaderLines',4, 'ReturnOnError', false,'MultipleDelimsAsOne',1);
fclose(fileID); stas = [dataArray{1:end-1}];
% Import start and end times
delimiter={'\t',' ','-'}; formatSpec = '%*s%*s%s%s%s%*s%s%s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, ...
    'HeaderLines',4, 'ReturnOnError', false,'MultipleDelimsAsOne',1);
fclose(fileID); times = [dataArray{1:end-1}]; times=cellfun(@(x)str2double(x), times);
% Import lat, lon, elevation
delimiter = {'\t','+'}; formatSpec = '%*s%*s%*s%*s%*s%s%s%s%[^\n\r]'; fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines',4, 'ReturnOnError', false);
fclose(fileID); locs = [dataArray{1:end-1}]; locs=cellfun(@(x)str2double(x), locs);

% If a station is multiply requested in .txt file, merge into single line
for iS=length(times(:,1)):-1:2
    if strcmp(stas{iS,2}, stas{iS-1,2})
        locs(iS,:)=[]; stas(iS,:)=[];
        if datenum(times(iS-1,1:3)) > datenum(times(iS,1:3))
            times(iS-1,1:3) = times(iS,1:3);
        end
        if datenum(times(iS-1,4:6)) < datenum(times(iS,4:6))
            times(iS-1,4:6) = times(iS,4:6);
        end
        times(iS,:)=[];
    end
end

% Remove overlap with previously requested stations
for iS=length(times(:,1)):-1:1;
    staEnd=times(iS,4:6);
    if(str2num(stas{iS,2}(1))+1); sta=['bodge' stas{iS,2}]; else sta=stas{iS,2}; end
    if(str2num(stas{iS,1}(1))+1); net=['bodge' stas{iS,1}]; else net=stas{iS,1}; end
    if isfield(oldNetworks, net) && isfield(oldNetworks.(net), sta)
        if datenum(lastEQ)-datenum(staEnd)>=0
            stas(iS,:)=[]; locs(iS,:)=[]; times(iS,:)=[];
        else times(iS,1:3)=lastEQ;
        end
    end
end

clear stations
for iS=1:length(times(:,1))
    stations(iS).Network=stas(iS,1);  stations(iS).Station = stas{iS,2};
    stations(iS).Latitude=locs(iS,1); stations(iS).Longitude = locs(iS,2);
    stations(iS).Elevation=locs(iS,3); stations(iS).Start_Year=times(iS,1);
    stations(iS).Start_Month=times(iS,2); stations(iS).Start_Day=times(iS,3);
    stations(iS).End_Year=times(iS,4); stations(iS).End_Month=times(iS,5);
    stations(iS).End_Day=times(iS,6);
end
stations=stations';
clearvars iS locs stas ans dataArray delimiter fileID filename formatSpec iS
eval(['save ' basedir 'Data/Station_Lists/Update_' Project datestr(date, 'dd.mm.yy') '.mat stations']);
system(['cp ' basedir 'Data/Station_Lists/rawstas.txt ' basedir 'Data/Station_Lists/Update_' Project datestr(date, 'dd.mm.yy') '.txt']);

[Networks,Events]=MakeRequests(basedir,Project,stations,NEIC_Catalog,'update');
save([Proj_Dir '/UpdateEvents.mat'],'Events')
save([Proj_Dir '/UpdateNetworks.mat'],'Networks')

% Merge with previous Events.mat file
evs=fieldnames(Events);  evsO=fieldnames(oldEvents);
evs1=double([Events.(evs{1}).Year Events.(evs{1}).Month Events.(evs{1}).Day]);
evsL=double([Events.(evs{end}).Year Events.(evs{end}).Month Events.(evs{end}).Day]);
evsO1=double([oldEvents.(evsO{1}).Year oldEvents.(evsO{1}).Month oldEvents.(evsO{1}).Day]);
evsOL=[oldEvents.(evsO{end}).Year oldEvents.(evsO{end}).Month oldEvents.(evsO{end}).Day];
[~,ind1]=min([datenum(evs1),datenum(evsO1)]); [~,indL]=max([datenum(evsL),datenum(evsOL)]);
if ind1==indL; names={'Events'; 'oldEvents'};
    eval(['mergeEvents=' names{ind1} ';']);
elseif ind1==1;  mergeEvents=Events; j=find(strcmp(evsO,evs(end)));
    for k=j+1:length(evsO); mergeEvents.(evsO{k})=oldEvents.(evsO{k}); end
else mergeEvents=oldEvents; j=find(strcmp(evs,evsO(end)));
    for k=j+1:length(evs); mergeEvents.(evs{k})=Events.(evs{k}); end
end
Events=mergeEvents;

% Merge with previous Networks.mat file
nets=fieldnames(Networks); mergeNetworks=oldNetworks;
for j=1:length(nets);
    if ~isfield(mergeNetworks,nets{j}); mergeNetworks.(nets{j})=[]; end
    
    stas=fieldnames(Networks.(nets{j}));
    for k=1:length(stas)
        if ~isfield(mergeNetworks.(nets{j}),stas{k});
            mergeNetworks.(nets{j}).(stas{k})=Networks.(nets{j}).(stas{k});
        else mergeNetworks.(nets{j}).(stas{k}).Events_Requested = ...
                unique([mergeNetworks.(nets{j}).(stas{k}).Events_Requested; ...
                Networks.(nets{j}).(stas{k}).Events_Requested]);
        end
    end
end
Networks=mergeNetworks;

save([Proj_Dir '/Events.mat'],'Events')
save([Proj_Dir '/Networks.mat'],'Networks')


end

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
MagRange=[6.0 10];
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
        Networks.(net)=[e];
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



function ArrayPicker(Proj_Dir, Project, nets,lowT,highT)
% Automatically find corrections to TauP pick for windowing
% Code by Ved Lekic


TTdir=[Proj_Dir 'TTmats/'];
if ~exist(TTdir,'dir'); mkdir(TTdir); end
TMP = load([Proj_Dir 'SAC_Filenames.mat']);
SAC_Filenames = TMP.SAC_Filenames; clear TMP

%Alter search windows if you wish to
Psearchrange=3; Ssearchrange=6;

% SLAT, SLON specify location of the station being prepped. The reason they
% are needed is that we want to avoid calculating correlations between
% waveform observed at stations that are very distant due to radiation
% pattern effects.

% Maximum distance allowed between stations in degrees
Search_Dist_Rng = 2.5;

winP_hf = 150;    winS_hf = 150; TTinfo = [];

load([Proj_Dir 'Events.mat']);  load([Proj_Dir 'Networks.mat'])
NETs  = fieldnames(Networks);
if strcmp(nets, 'all');
    inet=1:length(NETs);   %Junlin
else inet=[];
    for j=1:length(nets); inet=[inet, find(strcmp(nets{j},NETs))]; end
    inet=unique(inet);
end

it='';%[];%'c'; %

for n1_id = inet
    savename=[TTdir 'TT_' Project '_' num2str(lowT) '_' num2str(highT) '_' NETs{n1_id} it ...
        '_' datestr(date,'dd.mm.yy')];
    
    STAs_mstr = fieldnames(Networks.(NETs{n1_id})); % station names
    % Cycle through all the stations
    len=length(STAs_mstr);
    
    staid = 1;
    for s1_id = staid:len  %Junlin
        disp(['Working on ' NETs{n1_id} ' station ' STAs_mstr{s1_id} '...']);
        
        % Initialize counter of successfully array-picked quakes
        % for particular station for P and S separately.
        st_indx_P = 0; st_indx_S = 0;
        
        % Get useful station identification information
        SLAT  = Networks.(NETs{n1_id}).(STAs_mstr{s1_id}).Latitude;
        SLON  = Networks.(NETs{n1_id}).(STAs_mstr{s1_id}).Longitude;
        
        station.Name = STAs_mstr{s1_id}; station.Network = NETs{n1_id};
        station.LAT  = SLAT; station.LON = SLON;
        
        eqs = Networks.(NETs{n1_id}).(STAs_mstr{s1_id}).Events_Requested;
        
        for e1 = 1:length(eqs)  %Junlin
            
            
            % Alter to match earthquake name in Events structure
            tmp = eqs(e1);
            for n = 1:length(tmp{1});
                if(strcmp(tmp{1}(n),'.')); tmp{1}(n) = '_'; end
            end
            eq_name = strcat('EQ_',tmp); Event = Events.(eq_name{1});
            
            if(Event.Magnitude > 0 && Event.Magnitude < 10.0)
                % Load the event data for the station of interest
                [P1,S1,nS1,nP1,~,~] = load_all_data_for_event(NETs(n1_id),...
                    STAs_mstr(s1_id),Networks,Event,Proj_Dir,SLAT,SLON,...
                    Search_Dist_Rng,winP_hf,winS_hf,lowT,highT,SAC_Filenames);
                
                % (vrijeme: time in Croatian)
                
                % If data for the station of interest for this
                % particular earthquake exists, then proceed to load
                % all available data for that particular earthquake
                
                nS=0; nP=0;
                if(nP1+nS1>0)
                    [P,S,nS,nP,vrijemeP,vrijemeS] = load_all_data_for_event(...
                        NETs,[],Networks,Event,Proj_Dir,SLAT,SLON,Search_Dist_Rng,...
                        winP_hf,winS_hf,lowT,highT,SAC_Filenames);
                end
                
                
                if(nP>=5 && nP1==1) % Require at least 5 qualifying stations to record that event
                    uzorP = trimmean(P,20);
                    
                    % Find offset between P in trimmean and observed at station p1_id
                    indcs = find(vrijemeP > -5 & vrijemeP < 5);
                    dt = vrijemeP(2) - vrijemeP(1);
                    [~, max_indxP] = max(xcorr(P1(indcs),uzorP(indcs),...
                        round(Psearchrange/dt)));
                    max_indxP = max_indxP - round(Psearchrange/dt) - 1;
                end
                if(nS>=5 && nS1==1)
                    uzorS = trimmean(S,20);
                    % Find offset between S in trimmean and observed at station s1_id
                    indcs = find(vrijemeS > -15 & vrijemeS < 15);
                    dt = vrijemeS(2) - vrijemeS(1);
                    [~, max_indxS] = max(xcorr(S1(indcs),uzorS(indcs),...
                        round(Ssearchrange/dt)));
                    max_indxS = max_indxS - round(Ssearchrange/dt) - 1;
                    % This gives the offset in terms of the index, bc
                    % the max index is at the best xcorr between the
                    % mean and the single signal, and the max offset is
                    % Ssearchrange (i.e. searches between -Ssearchrange
                    % to +Ssearchrange; so if max_indxS==1, then the
                    % best xcorr is at -Ssearchrange)
                end
                
                if((nP>=5 && nP1==1) || (nS>=5 && nS1==1))
                    % Now, re-load the variables, etc. for the station of interest.
                    % Get useful station identification information
                    
                    station.Name = STAs_mstr{s1_id}; station.Network = NETs{n1_id};
                    station.LAT  = SLAT; station.LON = SLON;
                    
                    SAC_Dir = [Proj_Dir NETs{n1_id} '/' STAs_mstr{s1_id} '/SAC_Files/'];
                    
                    % Load SAC files
                    [tmp,n_sac_files,Multiples_Exist,tmp_comp] = ...
                        Load_SAC_Files(Event,station,SAC_Dir,eqs,e1,SAC_Filenames.(NETs{n1_id}).(STAs_mstr{s1_id}));
                    
                    % If SAC files were loaded and there is no problem with
                    % multiple files
                    
                    if(n_sac_files==3 && Multiples_Exist==0)
                        
                        % Align components to origin time of quake
                        [tmp,origin_time,dt] = Align_Components(tmp_comp,tmp,Event,'data');
                        if isempty(tmp)
                            continue;
                        end
                        
                        
                        % Get ray path data including P and S travel times
                        [Ray_Paths] = Find_TauP_Paths(Event, ...
                            Networks.(NETs{n1_id}).(STAs_mstr{s1_id}),origin_time, 'ArrayPicker');
                        
                        % Rotate to GCP
                        [tmp,tmp_comp] = Rotate_NE_into_RT(tmp,tmp_comp);
                        
                        % Filter P and S separately
                        [bpf_amps] = Filter_Waveforms(dt,tmp,tmp_comp,STAs_mstr{s1_id},lowT,highT);
                        
                        
                        if nP>=5 && nP1==1
                            % Identify vertical component, to be used for picking P
                            c_id = 1; %find(strcmp('BHZ',tmp_comp));
                            
                            % Now, put in the trimmean P into the bpf_amps
                            indcs = find(tmp(c_id).sac_file(:,1) > min(Ray_Paths.P_Path.time + vrijemeP) & ...
                                tmp(c_id).sac_file(:,1) < max(Ray_Paths.P_Path.time + vrijemeP));
                            indcs = indcs(indcs>-max_indxP & indcs<(length(bpf_amps(c_id,1).bpf_amp) - max_indxP));
                            display(['Indices ' num2str(length(indcs)) ' eq = ' num2str(e1)]);
                            bpf_amps(c_id,1).bpf_amp(indcs+max_indxP) = interp1(Ray_Paths.P_Path.time + vrijemeP,...
                                uzorP,tmp(c_id).sac_file(indcs,1));
                        end
                        
                        if nS>=5 && nS1==1
                            % Identify transverse component, to be used for picking
                            % S. This is implicitly assuming that SH and SV arrive
                            % concurrently, because SV is really giving rise to
                            % the S-to-P conversions we are analyzing.
                            
                            c_id = 5; %find(strcmp('BHT',tmp_comp));
                            
                            % Now, put in the trimmean S into the bpf_amps
                            indcs = find(tmp(c_id).sac_file(:,1) > min(Ray_Paths.S_Path.time + vrijemeS) & ...
                                tmp(c_id).sac_file(:,1) < max(Ray_Paths.S_Path.time + vrijemeS));
                            
                            indcs = indcs(indcs>-max_indxS & indcs<(length(bpf_amps(c_id,2).bpf_amp) - max_indxS));
                            
                            bpf_amps(c_id,2).bpf_amp(indcs+max_indxS) = interp1(Ray_Paths.S_Path.time + vrijemeS,...
                                uzorS,tmp(c_id).sac_file(indcs,1));
                        end
                        
                        %Find the Ps and Sp analysis and phase windows from Z and R/T
                        [Waveform_Data] = Set_Anaylsis_Windows(dt,[],tmp,bpf_amps,'ZT',...
                            Ray_Paths,nS,nP,'ArrayPicker');
                        
                        if(Waveform_Data.Ps_Waveforms.S2N_ratio>-1)
                            % If signal to noise was high enough to allow a S2N-based pick
                            if(Waveform_Data.Ps_Waveforms.Window.Taup_S2N_Misfit~=0)
                                st_indx_P = st_indx_P + 1;
                                % Store the picked travel-time vs. TauP prediction
                                TTinfo.(NETs{n1_id}).(STAs_mstr{s1_id}).Parr(st_indx_P) = ...
                                    Waveform_Data.Ps_Waveforms.Window.Taup_S2N_Misfit;
                                % Store P-wave ray parameter and backazimuth
                                TTinfo.(NETs{n1_id}).(STAs_mstr{s1_id}).P_rp(st_indx_P) = ...
                                    Ray_Paths.P_Path.ray_parameter;
                                TTinfo.(NETs{n1_id}).(STAs_mstr{s1_id}).P_baz(st_indx_P) = ...
                                    Ray_Paths.P_Path.bAzimuth;
                                TTinfo.(NETs{n1_id}).(STAs_mstr{s1_id}).P_S2N(st_indx_P) = ...
                                    Waveform_Data.Ps_Waveforms.S2N_ratio;
                            end
                        end
                        
                        if(Waveform_Data.Sp_Waveforms.S2N_ratio>-1)
                            % If signal to noise was high enough to allow a S2N-based pick
                            if(Waveform_Data.Sp_Waveforms.Window.Taup_S2N_Misfit~=0)
                                st_indx_S = st_indx_S + 1;
                                % Store the picked travel-time vs. TauP prediction
                                TTinfo.(NETs{n1_id}).(STAs_mstr{s1_id}).Sarr(st_indx_S) = ...
                                    Waveform_Data.Sp_Waveforms.Window.Taup_S2N_Misfit;
                                % Store S-wave ray parameter and backazimuth and S2Nratio
                                TTinfo.(NETs{n1_id}).(STAs_mstr{s1_id}).S_rp(st_indx_S) = ...
                                    Ray_Paths.S_Path.ray_parameter;
                                TTinfo.(NETs{n1_id}).(STAs_mstr{s1_id}).S_baz(st_indx_S) = ...
                                    Ray_Paths.S_Path.bAzimuth;
                                TTinfo.(NETs{n1_id}).(STAs_mstr{s1_id}).S_S2N(st_indx_S) = ...
                                    Waveform_Data.Sp_Waveforms.S2N_ratio;
                                %TTinfo.(NETs{n1_id}).(STAs_mstr{s1_id}).Env(st_indx_S)=...
                                %    Waveform_Data.Sp_Waveforms.Window.Envelope;
                            end
                        end
                        
                        
                    end
                end
            end
        end
        eval(['save ' savename '.mat TTinfo']);
    end
    %     profile viewer
    %     pause
end



end

function ArrayPickReformat(Proj_Dir, Project,lowT,highT)
% Compile all TTarray files and reformat them for use in prepping

TTdir=[Proj_Dir 'TTmats/'];  TTstring=['TT_' Project '_' num2str(lowT) '_' num2str(highT) '_'];
cur_dir=pwd; cd(TTdir);

% Identify individual TT matrices
files=dir; names{length(files)}=[]; k=0;
for j=1:length(files)
    if(length(files(j).name)>length(TTstring) &&  ...
            strcmp(files(j).name(1:length(TTstring)),TTstring))
        k=k+1; names{k}=files(j).name;
    end
end
names=names(~cellfun('isempty',names));

% Compile all TT pick matrices into one variable, TTall
for j=1:length(names)
    load(names{j});
    if ~isempty(TTinfo)
        nets=fieldnames(TTinfo);
        
        for k=1:length(nets)
            if(j==1 || ~isfield(TTall,nets{k}));
                TTall.(nets{k})=TTinfo.(nets{k});
            else stas=fieldnames(TTinfo.(nets{k}));
                for m=1:length(stas)
                    TTall.(nets{k}).(stas{m})=TTinfo.(nets{k}).(stas{m});
                end
            end
        end
    end
end

% Check TTall against Networks to check for missing stations
clear TTinfo files j k m names nets stas; load([Proj_Dir 'Networks.mat']);
nets=fieldnames(Networks); n=0;
for j=1:length(nets)
    if ~isfield(TTall, nets{j})
        disp(['No TTinfo for ' nets{j} ' network']); continue;
    end
    
    stas=fieldnames(Networks.(nets{j}));
    for k=1:length(stas)
        if ~isfield(TTall.(nets{j}),stas{k});
            disp(['No TTinfo for station ' nets{j} '.' stas{k}]); n=n+1;
        end
    end
end
disp(['Stations without TT info: ' num2str(n)])

% Reprocess TTall and save
nets=fieldnames(TTall); baz_bins=-180:90:180;
for j = 1:length(nets)
    stas = fieldnames(TTall.(nets{j}));
    for k = 1:length(stas)
        disp(['Working on network ' nets{j} ' and station ' stas{k} '...']);
        
        if isfield(TTall.(nets{j}).(stas{k}),'Parr')
            x = TTall.(nets{j}).(stas{k}).P_baz;
            y = TTall.(nets{j}).(stas{k}).Parr;
            [~, bins] = histc(x,baz_bins);
            ty = sparse(1:length(x),bins,y);
            mu = full(sum(ty)./sum(ty~=0));
            if length(mu)<4; mu(end+1:4)=nan; end %mu(isnan(mu)) = 0;
            TTall.(nets{j}).(stas{k}).Pbins = baz_bins;
            TTall.(nets{j}).(stas{k}).Pavgs = mu;
        end
        
        if isfield(TTall.(nets{j}).(stas{k}),'Sarr')
            x = TTall.(nets{j}).(stas{k}).S_baz;
            y = TTall.(nets{j}).(stas{k}).Sarr;
            [~, bins] = histc(x,baz_bins);
            ty = sparse(1:length(x),bins,y);
            mu = full(sum(ty)./sum(ty~=0));
            if length(mu)<4; mu(end+1:4)=nan; end
            TTall.(nets{j}).(stas{k}).Sbins = baz_bins;
            TTall.(nets{j}).(stas{k}).Savgs = mu;
        end
    end
end

TTinfo = TTall;
cd(cur_dir);

if exist([Proj_Dir 'names_' num2str(lowT) '_' num2str(highT) '.mat'],'file')
    load([Proj_Dir 'names_' num2str(lowT) '_' num2str(highT) '.mat']);
end
FileNames.TT=['TT_' Project '_' num2str(lowT) '_' num2str(highT) '_reformatted_' datestr(date,'dd.mm.yy') '.mat'];
eval(['save ' Proj_Dir 'names_' num2str(lowT) '_' num2str(highT) '.mat FileNames']);

eval(['save ' Proj_Dir FileNames.TT ' TTinfo']);

end

function out=splitstring(str)

ind = find(isspace(str));

for n=1:length(ind)
    if n==1; k=n; out{k}=str(1:ind(n)-1);
    elseif ind(n)-ind(n-1)==1
        continue
    else k=k+1; out{k}=str(ind(n-1)+1:ind(n)-1);
    end
end

end

function Prep_Waveforms(Proj_Dir,NETs,sta, tag,lowT,highT, varargin)


switch tag
    case 'data'
        TMP=load([Proj_Dir '/Events.mat']); Events=TMP.Events; clear TMP
        TMP=load([Proj_Dir '/Networks.mat']); Networks = TMP.Networks; clear TMP
        TMP=load([Proj_Dir '/SAC_Filenames.mat']); SAC_Filenames = TMP.SAC_Filenames; clear TMP
        version=['names_' num2str(lowT) '_' num2str(highT)];
        
        
        if ~isempty(varargin)
            ip = varargin{2};
            marks = varargin{3};
        else
            ip = 0;
        end
        
    case 'synth'
        junk=strsplit(Proj_Dir,'Projects/'); dataProj_Dir=Proj_Dir;
        Proj_Dir=[junk{1} 'Projects/Synth_' junk{2}]; % synthetics directory
        if strcmp(varargin{end},'test'); dataProj_Dir=Proj_Dir; end
        TMP=load([dataProj_Dir '/Events.mat']); Events=TMP.Events; clear TMP
        TMP=load([Proj_Dir '/Networks.mat']); Networks = TMP.Networks; clear TMP
        Phase=varargin{1}; version=varargin{2};
end


% Load Array picker info
if exist([Proj_Dir  version '.mat'],'file'); load([Proj_Dir version '.mat']);
    if isfield(FileNames,'TT');
        TMP=load([Proj_Dir FileNames.TT]); TTinfo=TMP.TTinfo;
        disp(FileNames.TT);
    else TTinfo=[];
    end
else TTinfo=[];
end

tpr = 5;

if ip==0
    stan = 1;
    ston = length(NETs);
else
    stan = marks(1);
    ston = marks(2);
end

for in=stan:ston
    NET = NETs{in};
    if strcmpi(sta,'all'); STAs = fieldnames(Networks.(NET)); else STAs=sta; end
    
    
    for is=1:length(STAs) %Junlin
        
        st_BestVpVs = cputime;
        
        clear *_Waveform_Data *_Ray_Path_Data PreCalculated_RFs BadData Temp_wave Temp_header
        STA = STAs{is};
        clc;
        disp(['Working on station ' NET '.' STA '...']);
        
        Station             = Networks.(NET).(STA);
        station.Name        = STA;
        station.Network     = NET;
        station.LAT         = Station.Latitude;
        station.LON         = Station.Longitude;
        SAC_Dir             = [Proj_Dir NET '/' STA '/SAC_Files/'];
        Events_to_Prep      = Station.Events_Requested; % Events that have been requested
        Events_Prepped      = {};                       % Events for which this station has data
        
        %         if length(Events_to_Prep)>50
        %             continue;
        %         end
        
        % Look to see whether we have array-picker info for this
        % station.
        clear TTarray
        if(isfield(TTinfo,NET) && isfield(TTinfo.(NET),STA))
            display('Using TTarray');
            if isfield(TTinfo.(NET).(STA),'S_rp')
                TTarray.S = scatteredInterpolant(TTinfo.(NET).(STA).S_rp',...
                    TTinfo.(NET).(STA).S_baz',...
                    TTinfo.(NET).(STA).Sarr','natural');
            else
                TTarray.S = [];
            end
            if isfield(TTinfo.(NET).(STA),'P_rp')
                TTarray.P = scatteredInterpolant(TTinfo.(NET).(STA).P_rp',...
                    TTinfo.(NET).(STA).P_baz',...
                    TTinfo.(NET).(STA).Parr','natural');
            else
                TTarray.P = [];
            end
        else
            TTarray.P = []; TTarray.S = [];
        end
        
        st_sta = cputime;
        
        disp(['Prepping ' STA ' (' NET ') waveform data...'])
        fid1 = fopen([Proj_Dir '/' NET '/' STA '/Data_Prep_Log_Temp.txt'],'w');
        fprintf('\n\n _____ Initial _____\n\n')
        fprintf(' Network:  %s\n Station:  %s\n Events:   %g\n Time:     ~%g min\n\n',...
            NET,STA,length(Events_to_Prep),ceil(length(Events_to_Prep)*1.75/60))
        fprintf(fid1,'\n\n _____ Initial _____\n\n');
        fprintf(fid1,' Network:  %s\n Station:  %s\n Events:   %g\n Time:     ~%g min\n\n',...
            NET,STA,length(Events_to_Prep),ceil(length(Events_to_Prep)*1.75/60));
        
        
        if strcmpi(tag,'synth')
            [traces,tvec,dt,shift]=read_synth_sac_files(Proj_Dir,NET,STA,Phase,version);
            nev=size(traces,3);
        else
            nev = length(Events_to_Prep);
            
        end
        % Counter for all events with multiples
        Events_with_Multiples   = {};
        
        % Run through all the events for which data was requested for this station
        for ie=1:nev   %Junlin

            
            EQ_ID = Events_to_Prep{ie};
            if strcmpi(tag,'synth'); eq_id=EQ_ID; else
                eq_id = ['EQ_' EQ_ID(1:4) '_' EQ_ID(6:8) '_' EQ_ID(10:11) '_' EQ_ID(13:14)];
            end
            
            fprintf(' Prepping event %g from station %s.%s\t(%s)...',ie,NET,STA,EQ_ID)
            fprintf(fid1,' Prepping event %g from station %s.%s\t(%s)...',ie,NET,STA,EQ_ID);
            
            st_evt = cputime;
            
            % Try to load all of the requested SAC files and check for multiples
            Event = Events.(eq_id);
            
            if strcmpi(tag,'synth')
                [sac_files,n_sac_files,Multiples_Exist,Components] = ...
                    synth_for_GUI(traces,tvec,ie);
                junk=dir([Proj_Dir NET '/' STA '/phasewidth_*mat']);
                junk=strsplit(junk.name,'_'); pwl=str2double(junk{2}(1:end-5));
            elseif strcmpi(tag,'data')
                [sac_files,n_sac_files,Multiples_Exist,Components] = ...
                    Load_SAC_Files(Event,station,SAC_Dir,Events_to_Prep,ie,SAC_Filenames.(NET).(STA));
            end
            
            if n_sac_files==3
                if strcmpi(tag, 'data')
                    % Make all components the same length and align to origin time
                    [sac_files,origin_time,dt] = Align_Components(Components,sac_files,Event,'data');
                    if isempty(sac_files)
                        BadData.(eq_id).Ps='Time sequence wrong'; BadData.(eq_id).Sp=BadData.(eq_id).Ps;
                        continue;
                    end
                    
                    
                    
                    % Rotate the N and E components into R and T
                    [sac_files,Components] = Rotate_NE_into_RT(sac_files,Components);
                    
                    % Filter the waveforms
                    [bpf_amps] = Filter_Waveforms(dt,sac_files,Components,STA,lowT,highT);
                    
                    
                    
                elseif strcmpi(tag,'synth')% origin_time=0;
                    % Make all components the same length and align to or  igin time
                    [sac_files,origin_time,dt] = Align_Components(Components,...
                        sac_files,Event,'synth',Phase);
                end
                
                % Get ray paths and ray parameters
                [Ray_Paths] = Find_TauP_Paths(Event,Station,origin_time);
                
                if isempty(Ray_Paths.P_Path) || isempty(Ray_Paths.S_Path)
                    continue;
                end
                
                if strcmpi(tag,'data')
                    % Find the Ps and Sp analysis and phase windows from Z and R/T
                    [Waveform_Data, BadEvent] = Set_Anaylsis_Windows(dt,[],sac_files,bpf_amps,'ZR',...
                        Ray_Paths,TTarray,'prep');
                    BadData.(eq_id)=BadEvent;
                    
                    if Waveform_Data.Ps_Waveforms.S2N_ratio~=-666 || ...
                            Waveform_Data.Sp_Waveforms.S2N_ratio~=-666
                        % Rotate the Z and R components into P and SV
                        
                        Events_Prepped = [Events_Prepped; {EQ_ID}];
                        
                        Temp_Waveform_Data.(eq_id)	= Waveform_Data;
                        Temp_Ray_Path_Data.(eq_id)	= Ray_Paths;
                        
                    else
                        fprintf(fid1,'!! No waveform data during the phase arrival times !!');
                    end
                    
                elseif strcmpi(tag,'synth')
                    [Waveform_Data] = Set_Anaylsis_Windows(dt,[],sac_files,[],[],...
                        Ray_Paths,Phase,shift,pwl,'synth');
                    Events_Prepped = [Events_Prepped; {EQ_ID}];
                    
                    % Instead of check polarity reversal, set xcorr
                    Waveform_Data.([Phase '_Waveforms']).ZR_Xcorr_Coeff = 1;
                    
                    [Waveform_Data,Ray_Paths]	= Reduce_File_Size(Waveform_Data,Ray_Paths,Phase,'synth');
                    All_Waveform_Data.(eq_id)  	= Waveform_Data;
                    All_Ray_Path_Data.(eq_id)  	= Ray_Paths;
                    PreCalculated_RFs.Ps.Impulse_Upper_Corner_Frequency	= [];
                    PreCalculated_RFs.Sp.Impulse_Upper_Corner_Frequency	= [];
                    PreCalculated_RFs.Ps.Individual_RFs.(eq_id)       	= [];
                    PreCalculated_RFs.Sp.Individual_RFs.(eq_id)        	= [];
                end
                
            elseif n_sac_files>3
                BadData.(eq_id).Ps='Too many SAC files'; BadData.(eq_id).Sp=BadData.(eq_id).Ps;
                fprintf(fid1,'!#* Somehow %g sac files exist for this event *#!',n_sac_files);
            elseif Multiples_Exist==1
                BadData.(eq_id).Ps='Multiples exist'; BadData.(eq_id).Sp=BadData.(eq_id).Ps;
                fprintf(fid1,'** Multiple components exist for this event **');
                Events_with_Multiples = [Events_with_Multiples; {EQ_ID}];
            elseif n_sac_files>0 && n_sac_files<3 && Multiples_Exist==0
                BadData.(eq_id).Ps='Too few SAC files'; BadData.(eq_id).Sp=BadData.(eq_id).Ps;
                fprintf(fid1,'## Only %g sac files exist for this event ##',n_sac_files);
            elseif n_sac_files==0
                BadData.(eq_id).Ps='No SAC files'; BadData.(eq_id).Sp=BadData.(eq_id).Ps;
                fprintf(fid1,'!! No data exists for this event !!');
            elseif n_sac_files<-660
                BadData.(eq_id).Ps='Pegged Component'; BadData.(eq_id).Sp=BadData.(eq_id).Ps;
                fprintf(fid1,'!! Pegged component !!');
            end
            et_evt = cputime;
            fprintf('done (%g sec).\n',et_evt-st_evt);
            fprintf(fid1,' done (%g sec).\n',et_evt-st_evt);
        end
        
        if strcmpi(tag,'data')
            et_BestVpVs = cputime;
        end
        
        
        disp(['Prepping ' STA ' (' NET ') waveform data...done'])
        
        if ~isempty(Events_Prepped)
            
            if strcmpi(tag,'data')
                
                st_BestVpVs = cputime;
                
                %start to search for the free surface velocity
                
                disp('Calculating Free Surface Velocity');
                [FS_VpVs,Vs_info,Vp_info] = Find_VpVs(Temp_Waveform_Data,Temp_Ray_Path_Data,lowT,highT);
                
                Events_Prepped_1 = Events_Prepped;
                Events_Prepped = {};
                for ie1=1:length(Events_Prepped_1)
                    
                    EQ_ID = Events_Prepped_1{ie1};
                    eq_id = ['EQ_' EQ_ID(1:4) '_' EQ_ID(6:8) '_' EQ_ID(10:11) '_' EQ_ID(13:14)];
                    if ie1<10
                        fprintf(' Prepping event %g\t\t(%s)...',ie1,EQ_ID)
                    else
                        fprintf(' Prepping event %g\t(%s)...',ie1,EQ_ID)
                    end
                    fprintf(fid1,' Prepping event %g\t(%s)...',ie1,EQ_ID);
                    
                    st_evt = cputime;
                    
                    Waveform_Data = Temp_Waveform_Data.(eq_id);
                    Ray_Path_Data = Temp_Ray_Path_Data.(eq_id);
                    
                    % Load in the sac files again
                    Event = Events.(eq_id);
                    [sac_files,n_sac_files,~,Components] = ...
                        Load_SAC_Files(Event,station,SAC_Dir,Events_Prepped_1,ie1,SAC_Filenames.(NET).(STA));
                    
                    if n_sac_files==3
                        [sac_files,~,dt]          = Align_Components(Components,sac_files,Event,'data');
                        if isempty(sac_files)
                            contimue;
                        end
                        [sac_files,Components]              = Rotate_NE_into_RT(sac_files,Components);
                        
                        % Rotate Z-R into P-SV using the mean best Vp and Vs
                        [sac_files,Components]      = Rotate_ZR_into_PSV(...
                            Waveform_Data,sac_files,Components,Ray_Path_Data,FS_VpVs);
                        Events_Prepped                      = [Events_Prepped; {EQ_ID}];
                        [bpf_amps]                          = Filter_Waveforms(dt,sac_files,Components,STA,lowT,highT);
                        
                        % Find the Ps and Sp analysis and phase windows from P and SV
                        [Waveform_Data]                     = Set_Anaylsis_Windows(dt,Waveform_Data,sac_files,bpf_amps,'PSV',...
                            Ray_Path_Data,TTarray,'prep');
                        [Waveform_Data]                     = Check_Polarity_Reversal(Waveform_Data);
                        [Waveform_Data,Ray_Path_Data]       = Reduce_File_Size(Waveform_Data,Ray_Path_Data);
                        All_Waveform_Data.(eq_id)           = Waveform_Data;
                        All_Ray_Path_Data.(eq_id)           = Ray_Path_Data;
                        PreCalculated_RFs.Ps.Impulse_Upper_Corner_Frequency	= [];
                        PreCalculated_RFs.Sp.Impulse_Upper_Corner_Frequency	= [];
                        PreCalculated_RFs.Ps.Individual_RFs.(eq_id)       	= [];
                        PreCalculated_RFs.Sp.Individual_RFs.(eq_id)        	= [];
                        
                    elseif n_sac_files~=3
                        fprintf(' !#* There is a problem with this event *#!')
                        fprintf(fid1,' !#* There is a problem with this event *#!');
                    end
                    et_evt = cputime;
                    fprintf('done (%g sec).\n',et_evt-st_evt)
                    fprintf(fid1,'done (%g sec).\n',et_evt-st_evt);
                    
                end
                
                clear Events_Prepped_1
                et_BestVpVs = cputime;
                fprintf('\n Calculation Time Using Best Vp,Vs:  %g sec\n',et_BestVpVs-st_BestVpVs)
                fprintf(fid1,'\n Calculation Using Best Mean Vp,Vs:  %g sec\n',et_BestVpVs-st_BestVpVs);
                
            end
            
            disp(['Prepping ' STA ' (' NET ') waveform data...done'])
            disp(['Saving prepped ' STA ' (' NET ') waveform data'])
            
            
            Station_Data                        = Networks.(NET).(STA);
            Station_Data.Events_Prepped     	= Events_Prepped;
            Station_Data.Events_with_Multiples	= Events_with_Multiples;
            Station_Data.FS_VpVs                = FS_VpVs;
            Station_Data.Vp_info                = Vp_info;
            Station_Data.Vs_info                = Vs_info;
            Networks.(NET).(STA).Events_Prepped = Events_Prepped;
            
            if ie==length(Events_to_Prep)
                save([Proj_Dir NET '/' STA '/Station_Data_' num2str(lowT) '_' num2str(highT) '.mat'],'Station_Data')
                if strcmpi(tag,'data'); FileNames.Synth=[]; end
                save([Proj_Dir NET '/' STA '/Waveform_Data_' num2str(lowT) '_' num2str(highT) FileNames.Synth '.mat'],'All_Waveform_Data')
                save([Proj_Dir NET '/' STA '/Ray_Path_Data_' num2str(lowT) '_' num2str(highT) '.mat'],'All_Ray_Path_Data')
                save([Proj_Dir NET '/' STA '/PreCalculated_RFs_' num2str(lowT) '_' num2str(highT) '.mat'],'PreCalculated_RFs')
                if strcmpi(tag,'data')
                    save([Proj_Dir NET '/' STA '/BadData_' num2str(lowT) '_' num2str(highT) '.mat'],'BadData')
                end
                
                if ip==0
                    save([Proj_Dir 'Networks.mat'],'Networks');
                else
                    save([Proj_Dir 'Networks_',num2str(ip),'.mat'],'Networks');
                end
                
                fprintf('\n\n Events with data:\t%g\n',length(Events_Prepped))
                fprintf(' Events with multiples:\t%g\n\n',length(Events_with_Multiples))
                fprintf(' \n The data for %s should now be ready for the migration code.\n\n',STA)
                fprintf(fid1,'\n\n Events with data:\t%g\n',length(Events_Prepped));
                fprintf(fid1,' Events with multiples:\t%g\n\n',length(Events_with_Multiples));
                fprintf(fid1,' \n The data for %s should now be ready for the migration code.\n\n',STA);
                
            elseif ie<length(Events_to_Prep)
                fprintf(fid1,'\n\n !! Some of the data for station %s has not been prepped !!\n\n',STA);
            end
        end
        et_sta = cputime;
        fprintf('\n Total Calculation Time For %s:\t%g sec\n',STA,et_sta-st_sta)
        fprintf(fid1,'\n Total Calculation Time For %s:\t%g sec\n',STA,et_sta-st_sta);
        fclose(fid1);
    end
end

end

function [P,S,nS,nP,vrijemeP,vrijemeS] = load_all_data_for_event(NETs,STAs,Networks,Event,Proj_Dir,SLAT,SLON,Search_Dist_Rng,winP_hf,winS_hf,lowT,highT,SAC_Filenames)
% If STAs is blank, then we just get data for all stations for that
% network, otherwise, we get the data for just the station(s) specified


% Initialize counter and clear P and S structures just in case
nP = 0; nS=0; clear P S;
P=zeros(50,6001); S=zeros(20,6001);


if(isempty(STAs)), useallSTAsforNET = 1; else useallSTAsforNET = 0; end;

% Now, cycle through all the Networks
for n_id = 1:length(NETs)
    
    if(useallSTAsforNET)
        % Get station names for the network NTs{n_id}
        STAs = fieldnames(Networks.(NETs{n_id}));
    end
    
    % Now, cycle through all the stations
    for s_id = 1:length(STAs)
        %display(['Working on station # ' num2str(s_id) '...']);
        
        % List of events requested for this Network/Station
        Events_Requested = Networks.(NETs{n_id}).(STAs{s_id}).Events_Requested;
        % Find current quake in list
        ie = find(strcmp(Event.Event_ID,Events_Requested));
        
        if(~isempty(ie))
            % Get useful station identification information
            station.Name = STAs{s_id}; station.Network = NETs{n_id};
            station.LAT  = Networks.(NETs{n_id}).(STAs{s_id}).Latitude;
            station.LON  = Networks.(NETs{n_id}).(STAs{s_id}).Longitude;
            
            % Check whether the station falls within the search range
            [dist_km] = distance(SLAT,SLON,station.LAT,station.LON);
            
            if(dist_km<Search_Dist_Rng)
                SAC_Dir = [Proj_Dir NETs{n_id} '/' STAs{s_id} '/SAC_Files/'];
                
                % Load SAC files
                %                 disp(STAs{s_id});
                %                 if strcmp(STAs{s_id},'UMT')
                %                     indaaa=1;
                %                 end
                
                [tmp,n_sac_files,Multiples_Exist,tmp_comp] = ...
                    Load_SAC_Files(Event,station,SAC_Dir,Events_Requested,ie,SAC_Filenames.(NETs{n_id}).(STAs{s_id}));
                
                % If SAC files were loaded and there is no problem with
                % multiple files
                
                
                if(n_sac_files==3 && Multiples_Exist==0)
                    
                    % Align components to origin time of quake
                    
                    %                     [tmpYR,tmpJDY,tmpHR,tmpMN,tmpSEC,tmpMSEC,tmpB,tmpdt,tmpnpts] = lh(tmp(1).sac_file,...
                    %                         'NZYEAR','NZJDAY','NZHOUR','NZMIN','NZSEC','NZMSEC','B','DELTA','NPTS');
                    tmpYR = tmp(1).sac_file(71,3);
                    tmpJDY = tmp(1).sac_file(72,3);
                    tmpHR = tmp(1).sac_file(73,3);
                    tmpMN = tmp(1).sac_file(74,3);
                    tmpSEC = tmp(1).sac_file(75,3);
                    tmpMSEC = tmp(1).sac_file(76,3);
                    tmpB = tmp(1).sac_file(6,3);
                    tmpdt = tmp(1).sac_file(1,3);
                    tmpnpts = tmp(1).sac_file(80,3);
                    evetime = Find_Absolute_Time(tmpdt,Event.Year,Event.Julian_Day,Event.Hour,Event.Minute,Event.Second,0);
                    headtime = Find_Absolute_Time(tmpdt,tmpYR,tmpJDY,tmpHR,tmpMN,tmpSEC,tmpMSEC)+double(tmpB);
                    durationt = tmpnpts*tmpdt;     %Junlin
                    
                    if headtime+durationt<evetime
                        continue;
                    else
                        %                         disp(STAs{s_id});
                    end
                    
                    [tmp,origin_time,dt] = Align_Components(tmp_comp,tmp,Event,'data');
                    
                    % Get ray path data including P and S travel times
                    [Ray_Paths] = Find_TauP_Paths(Event,...
                        Networks.(NETs{n_id}).(STAs{s_id}),origin_time,'ArrayPicker');
                    
                    % Rotate to GCP
                    [tmp,tmp_comp] = Rotate_NE_into_RT(tmp,tmp_comp);
                    
                    % Filter P and S separately
                    [bpf_amps] = Filter_Waveforms(dt,tmp,tmp_comp,STAs{s_id},lowT,highT);
                    
                    % Time vectors for windows around P and  S arrivals
                    vrijemeP = -winP_hf:dt:winP_hf;
                    vrijemeS = -winS_hf:dt:winS_hf;
                    
                    % Identify vertical component, to be used for picking P
                    c_id = 1; %find(strcmp('BHZ',tmp_comp));
                    
                    if ~isempty(Ray_Paths.P_Path)
                        p = interp1(tmp(c_id).sac_file(:,1),...
                            bpf_amps(c_id,1).bpf_amp,Ray_Paths.P_Path.time + vrijemeP);
                    else
                        p = nan;
                    end
                    
                    % Identify transverse component, to be used for picking
                    % S. This is implicitly assuming that SH and SV arrive
                    % concurrently, because SV is really giving rise to
                    % the S-to-P conversions we are analyzing.
                    
                    c_id = 5; %find(strcmp('BHT',tmp_comp));
                    
                    if ~isempty(Ray_Paths.S_Path)
                        s = interp1(tmp(c_id).sac_file(:,1),...
                            bpf_amps(c_id,2).bpf_amp,Ray_Paths.S_Path.time + vrijemeS);
                    else
                        s = nan;
                    end
                    
                    if(sum(isnan(s)==0)); nS=nS+1; S(nS,:)=s; end
                    if(sum(isnan(p)==0)); nP=nP+1; P(nP,:)=p; end
                else
                    
                end
            end
        end
    end
    P=P(1:nP,:); S=S(1:nS,:);
end
if(nP==0), P = []; vrijemeP = []; end
if(nS==0), S = []; vrijemeS = []; end
if(nS==0 && nP==0); dt = []; end

    function [ABS_TIME] = Find_Absolute_Time(dt,YR,JDY,HR,MN,SEC,MSEC)
        
        % ********* Function Description *********
        %
        % This will determine the absolute start
        % time for a sac file from calendar year
        % 1980 and day 0 in seconds (i.e., from
        % Julian day 0 of year 1980 at hour, min,
        % and sec 0)
        %
        %
        % ****************************************
        % *                                      *
        % *  Written by David L. Abt - May 2008	 *
        % *                                      *
        % *  Email: David_Abt@brown.edu          *
        % *                                      *
        % ****************************************
        
        sec_per_yr      = 365*24*3600;                  % Seconds in a non-leap year
        sec_per_lpyr    = 366*24*3600;                  % Seconds in a leap year
        sec_per_4yrs    = 3*sec_per_yr+sec_per_lpyr;    % Seconds in a four-year period
        
        dt = double(dt);
        YR = double(YR)-1980;
        JDY = double(JDY);
        HR = double(HR);
        MN = double(MN);
        SEC = double(SEC);
        MSEC = double(MSEC);
        
        SEC = SEC+(MSEC*1e-3);
        SEC = round(SEC/dt)*dt;
        
        lpyr_cycles = floor((YR-1)/4);
        if lpyr_cycles<0
            lpyr_cycles = 0;
        end
        extra_years = YR-1-lpyr_cycles*4;
        if extra_years<0
            extra_years = 0;
        end
        
        ABS_TIME =  (sec_per_4yrs*lpyr_cycles)+ ... % Seconds through the last leap year
            (sec_per_yr*extra_years)+ ...   % Seconds from last leap year to the end of the previous year
            ((JDY-1)*24*3600)+ ...          % Seconds from the first of this year up to "yesterday"
            (HR*3600+MN*60+SEC);            % Seconds "today"
        
    end

end

function [sac_files,n_sac_files,Multiples_Exist,Components] = Load_SAC_Files(EQ,STA,SAC_Dir,Events,ie,filelist)

% ********* Function Description *********
%
% Load the raw SAC files.  Check for
% multiple files of the same components,
% merge them, and pick which one to use.
% Also, some useful/necessary header
% information is added to the sac files.
%
%
% ****************************************
% *                                      *
% *  Written by David L. Abt - May 2008	 *
% *                                      *
% *  Email: David_Abt@brown.edu          *
% *                                      *
% ****************************************


%Junlin for components

n_sac_files     = 0;
Components      = {'BHZ';'BHN';'BHE'};
Event_ID        = EQ.Event_ID(1:11);
Multiples_Exist = 0;

Event_id=Events{ie}(1:14);
if ie<length(Events); Event_idp=Events{ie+1}(1:14); end
if ie>1; Event_idm=Events{ie-1}(1:14);end
nmin=str2double(Event_id(13:14));

for icmp=1:length(Components)
    COMP = char(Components(icmp));
    if strcmp(COMP(2:3),'HN');     NUM = '1';
    elseif strcmp(COMP(2:3),'HE'); NUM = '2';
    else NUM='666';
    end
    
    % Check for "normal" file names (YR,JDY,and HR of Event_ID)
    %     file_inds = [find(strcmp([Event_id(1:11),COMP(2:3)],filelist.mark));...
    %         find(strcmp([Event_id(1:11),'H',NUM],filelist.mark))];
    %     Filesn = filelist.name(file_inds);
    Filesn = filelist.name(strcmp([Event_id(1:11),COMP(2:3)],filelist.mark) | ...
        strcmp([Event_id(1:11),'H',NUM],filelist.mark));
    %     Filesn = [dir([SAC_Dir Event_id(1:11) '*' STA.Name '*' COMP(2:3) '*.SAC']); ...
    %         dir([SAC_Dir Event_id(1:11) '*' STA.Name '*H' NUM '*.SAC'])];
    
    
    for k=length(Filesn):-1:1
        namef = Filesn{k};
        keep=1; min=str2double(namef(13:14)); %min=str2double(Filesn(k).name(13:14));
        if abs(nmin-min)>20; keep=0; end
        % Check if preceding and subsequent events are in the same hr, in which
        % case check which is closest
        if ie>1 && strcmp(Event_idm(1:11),Event_id(1:11));
            mmin=str2double(Event_idm(13:14));
            if abs(mmin-min)<abs(nmin-min); keep=0; end
        end
        if ie<length(Events) && strcmp(Event_idp(1:11),Event_id(1:11));
            pmin=str2double(Event_idp(13:14));
            if abs(pmin-min)<abs(nmin-min); keep=0; end
        end
        if keep==0; Filesn(k)=[]; end
    end
    
    Filesm = [];    Filesp = [];
    
    % Set the EventID for the previous hour
    YRm = EQ.Year;  JDYm = EQ.Julian_Day; HRm = EQ.Hour-1;
    if HRm<0;   HRm = 23; JDYm = JDYm-1; end
    if JDYm<1 && floor((YRm-1)/4)~=(YRm-1)/4
        JDYm = 365;
        YRm = YRm-1;
    elseif JDYm<1 && floor((YRm-1)/4)==(YRm-1)/4
        JDYm = 366;
        YRm = YRm-1;
    end
    YRms = num2str(YRm);
    if JDYm<10;                   JDYms = ['00' num2str(JDYm)];
    elseif JDYm>=10 && JDYm<100;  JDYms = ['0' num2str(JDYm)];
    elseif JDYm>=100;             JDYms = num2str(JDYm);
    end
    if HRm<10;        HRms = ['0' num2str(HRm)];
    elseif HRm>=10;   HRms = num2str(HRm);
    end
    Event_IDm = [YRms '.' JDYms '.' HRms];
    
    % Set the EventID for the following hour
    YRp = EQ.Year; JDYp = EQ.Julian_Day; HRp = EQ.Hour+1;
    if HRp>23;    HRp = 0; JDYp = JDYp+1; end
    if JDYp>365 && floor(YRp/4)~=YRp/4
        JDYp = 1;
        YRp = YRp+1;
    elseif JDYp>366 && floor(YRp/4)==YRp/4
        JDYp = 1;
        YRp = YRp+1;
    end;
    YRps = num2str(YRp);
    if JDYp<10;                  JDYps = ['00' num2str(JDYp)];
    elseif JDYp>=10 && JDYp<100; JDYps = ['0' num2str(JDYp)];
    elseif JDYp>=100;            JDYps = num2str(JDYp);
    end
    if HRp<10;      HRps = ['0' num2str(HRp)];
    elseif HRp>=10; HRps = num2str(HRp);
    end
    Event_IDp = [YRps '.' JDYps '.' HRps];
    
    if isempty(Filesn)
        % Check for files from preceeding hour
        if length(Events)>1;
            if ie==1;     Filesm = [];
            elseif ie>1
                if strcmp(Event_idm(1:11),Event_ID);     Filesm = [];
                else
                    %                     file_inds_m = [find(strcmp([Event_IDm,COMP(2:3)],filelist.mark));...
                    %                         find(strcmp([Event_IDm,'H',NUM],filelist.mark))];
                    %                     Filesm = filelist.name(file_inds_m);
                    Filesm = filelist.name(strcmp([Event_IDm,COMP(2:3)],filelist.mark) | ...
                        strcmp([Event_IDm,'H',NUM],filelist.mark));
                    %                     Filesm = [dir([SAC_Dir Event_IDm '*' STA.Name '*' COMP(2:3) '.*.SAC']); ...
                    %                         dir([SAC_Dir Event_IDm '*' STA.Name '*H' NUM '*.SAC'])];
                    if strcmp(Event_idm(1:11),Event_IDm)
                        mmin=str2double(Event_idm(13:14));
                        for k=length(Filesm):-1:1
                            namefm = Filesm{k};
                            keep=1; min=str2double(namefm(13:14));
                            if abs(mmin-min)<abs((nmin+60)-min); keep=0; end
                            if keep==0; Filesm(k)=[]; end
                        end
                    end
                end
            end
        elseif length(Events)==1
            Filesm = filelist.name(strcmp([Event_IDm,COMP(2:3)],filelist.mark));
            %Filesm = dir([SAC_Dir Event_IDm '*' STA.Name '*' COMP(2:3) '.*.SAC']);
        end
        % Check for files from following hour
        if length(Events)>1
            if ie==length(Events);   Filesp = [];
            elseif ie<length(Events)
                if strcmp(Event_idp(1:11),Event_ID);  Filesp = [];
                else
                    %                     file_inds_p = [find(strcmp([Event_IDp,COMP(2:3)],filelist.mark));...
                    %                         find(strcmp([Event_IDp,'H',NUM],filelist.mark))];
                    %                     Filesp = filelist.name(file_inds_p);
                    Filesp = filelist.name(strcmp([Event_IDp,COMP(2:3)],filelist.mark) | ...
                        strcmp([Event_IDp,'H',NUM],filelist.mark));
                    %                     Filesp = [dir([SAC_Dir Event_IDp '*' STA.Name '*' COMP(2:3) '.*.SAC']); ...
                    %                         dir([SAC_Dir Event_IDp(1:11) '*' STA.Name '*H' NUM '*.SAC'])];
                    if strcmp(Event_idp(1:11),Event_IDp)
                        pmin=str2double(Event_idp(13:14));
                        for k=length(Filesp):-1:1
                            namefp = Filesp{k};
                            keep=1; min=str2double(namefp(13:14));
                            if abs(pmin-min)<abs((60-nmin)-min); keep=0; end
                            if keep==0; Filesp(k)=[]; end
                        end
                    end
                end
                
            end
        elseif length(Events)==1
            Filesp = filelist.name(strcmp([Event_IDp,COMP(2:3)],filelist.mark));
            %Filesp = dir([SAC_Dir Event_IDp '*' STA.Name '*' COMP(2:3) '.*.SAC']);
        end
    end
    
    Files = [Filesn; Filesm; Filesp];
    
    if length(Files)==1
        %fprintf( 'loading %s\n', Files.name )
        sac_file = rsac([SAC_Dir Files{1}]);  % If only one file exists, read it in
        
        %         if strcmp(STA.Name,'RPZ')
        %             [zzs,pps,ggs]=read_sac_pole_zero...
        %                 ('/Users/jhua/Documents/NZnewer/NewZealand/Scattered_Waves/Data/Projects/NewZealand/NZ/RPZ/2003-01-03T00:00:012015-10-14T23:00:00.txt');
        %             sac_file(:,2)=rm_SACPZ(sac_file(:,2),zzs,pps,ggs,100);
        % %             [zzs,pps,ggs]=read_sac_pole_zero...
        % %                 ('/Users/jhua/Documents/NZnewer/NewZealand/Scattered_Waves/Data/Projects/NewZealand/NZ/RPZ/2003-01-03T00:00:012015-10-14T23:00:00QRZ.txt');
        % %             sac_file(:,2)=ad_SACPZ(sac_file(:,2),zzs,pps,ggs,100);
        %         end
        %         if strcmp(STA.Name,'RPZ')
        %             sac_file(:,2)=-sac_file(:,2);
        %         end
        if strcmp(COMP(end-1:end),'HZ') && sac_file(59,3)==180
            sac_file(:,2)=-sac_file(:,2);
            
        end
        if ~isempty(findstr(Files{1},'4A.DAR')) && ~isempty(findstr(Files{1},'HHZ'))
            sac_file(:,2)=-sac_file(:,2);
        end
        if ~isempty(findstr(Files{1},'NZ.RPZ'))
            sac_file(:,2)=-sac_file(:,2);
        end
        
        if sac_file(307,3)==1
            sac_file=[];
        end
        
    elseif length(Files)>1
        Multiples_Exist = 1;
        sac_file = [];  % If multiple files exist, make sac_file empty
    elseif isempty(Files)
        sac_file = [];  % No files exist
    end
    
    if ~isempty(sac_file)
        sac_std = std(sac_file(:,2)-mean(sac_file(:,2)));
        % Check for pegged components i.e. adjacent values in the
        % SAC file are identical
        jn=20; inds=round(length(sac_file(:,2))*rand(jn,1));
        inds=inds(inds>0 & ~isnan(inds) & inds<length(sac_file(:,2))-1); jn=length(inds);
        diffs=diff(sac_file(reshape([inds inds+1]',1,2*jn),2));
        
        if sac_std==0 && ~strcmp(STA.Name,'SYN') && ~strcmp(STA.Name,'FAKE')
            clear sac_file;   sac_file = [];
        elseif any(isnan(sac_file(:,2))==1)
            clear sac_file;   sac_file = [];
        elseif(length(find(diffs(1:2:end)==0))==jn) && ~strcmp(STA.Name,'FAKE') % remove pegged components
            clear sac_file; sac_file = [];  %this also rejects synthetics
            n_sac_files=-666;   %this also rejects synthetics
        else
            n_sac_files = n_sac_files+1;
            ELAT = EQ.Latitude;
            ELON = EQ.Longitude;
            EDEP = EQ.Depth;
            SLAT = STA.LAT;
            SLON = STA.LON;
            
            % Determine the eq-sta distance and azimuth
            DIST	= distance('gc',ELAT,ELON,SLAT,SLON);
            BAZ    	= azimuth('gc',SLAT,SLON,ELAT,ELON);
            AZ      = BAZ+180;
            if AZ>360
                AZ  = AZ-360;
            end
            
            % Update the header information
            sac_file = ch(sac_file,'EVLA',ELAT,'EVLO',ELON,'EVDP',EDEP,...
                'DIST',DIST,'AZ',AZ,'BAZ',BAZ,...
                'KCMPNM',COMP,'KSTNM',STA.Name);
            
            sac_file(:,2) = sac_file(:,2)-mean(sac_file(:,2));  % remove the mean
        end
    end
    
    sac_files(icmp,1).sac_file = sac_file;
    
end


if n_sac_files==3
    if sac_files(2).sac_file(58,3)~=0 && sac_files(3).sac_file(58,3)~=90
        % If BHN/1 and/or BHE/2 components are not oriented N-S and E-W
        % then rotate into N-S and E-W
        [sac_files,~,~] = Align_Components(Components,sac_files,EQ,'data');
        
        AZ_1 = lh(sac_files(2).sac_file,'CMPAZ');
        
        BH1 = sac_files(2).sac_file(:,2);
        BH2 = sac_files(3).sac_file(:,2);
        
        rot = -AZ_1*pi/180;
        ROT = [cos(rot) sin(rot); -sin(rot) cos(rot)];
        
        X = ROT*[BH1'; BH2'];
        BHN = X(1,:)';
        BHE = X(2,:)';
        
        sac_files(2,1).sac_file(:,2) = BHN;
        sac_files(2,1).sac_file = ch(sac_files(2).sac_file,'CMPAZ',0);
        
        sac_files(3,1).sac_file(:,2) = BHE;
        sac_files(3,1).sac_file = ch(sac_files(3).sac_file,'CMPAZ',90);
    elseif (sac_files(2).sac_file(58,3)==0 && sac_files(3).sac_file(58,3)~=90) || ...
            (sac_files(2).sac_file(58,3)~=0 && sac_files(3).sac_file(58,3)==90)
        display(['ERROR: ' num2str(sac_files(2).sac_file(58,3)) ' and ' num2str(sac_files(3).sac_file(58,3))]);
        fprintf('  > !! There is a non-orthogonality problem with the horizontals !! \n')
        sac_files_temp  = sac_files(1);
        clear sac_files
        sac_files       = sac_files_temp;
        n_sac_files     = 1;
    end
end

end

function [traces,tvec, dt, shift] = read_synth_sac_files(Proj_Dir, NET, STA,Phase,version)

format compact
format short g
load([Proj_Dir version '.mat']); junk=strsplit(FileNames.Synth, '_synth');
name=[NET '_' STA junk{2}];
filename=[Proj_Dir NET '/' STA '/' Phase '_inputfiles/' name '.tr'];
fid=fopen(filename,'r');

% Read header:
line=fgetl(fid);
while line(1) == '#'; line=fgetl(fid); end
dum=str2num(line);
ntr=dum(1); nsamp=dum(2); dt=dum(3); align=dum(4); shift=dum(5);
% # traces, # samples or timesteps; # timestep; align or not; trace shift

% Read each trace
fprintf('Reading synthetic traces...\n'); int=0.05;
line=fgetl(fid);
for itr=1:ntr
    if itr/ntr>int; fprintf('  %s%%  ', num2str(int*100));
        if abs(int-0.5)<10E-5; fprintf('\n'); end; int=int+0.05;
    end
    while line(1) == '#'; line=fgetl(fid); end
    
    for isamp=1:nsamp
        traces(:,isamp,itr)=str2num(line);
        line=fgetl(fid);
    end
end
fclose(fid);
fprintf('\n\nNumber of traces: %d;  Number of timesteps: %d\n\n', size(traces,3), size(traces,2))

tvec=-shift:dt:(dt*(nsamp-1)-shift);

end

function [sac_files,n_sac_files,Multiples_Exist,Components] = synth_for_GUI(traces,tvec,ie)

n_sac_files=size(traces,1);
if n_sac_files==3;
    sac_files(1).sac_file=[tvec;traces(1,:,ie)]'; % P  component
    sac_files(2).sac_file=[tvec;traces(2,:,ie)]'; % SV component
    sac_files(3).sac_file=[tvec;traces(3,:,ie)]'; % SH component
else
    sac_files(1)=[]; sac_files(2)=[]; sac_files(3)=[];
end

Multiples_Exist=0; Components={'P';'SV';'SH'};


end

function [sac_files,origin_time,dt] = Align_Components(Components,sac_files,EQ,varargin)

% ********* Function Description *********
%
% Calculate the maximum time range of all
% the components, pad them with zeros if
% necessary and align them to the origin
% time of the earthquake
%
%
% ****************************************
% *                                      *
% *  Written by David L. Abt - May 2008	 *
% *                                      *
% *  Email: David_Abt@brown.edu          *
% *                                      *
% ****************************************

% Initialise arrays
start_times_abs=[0 0 0]; Npts=zeros(1,length(Components));
switch varargin{1}
    case 'data'
        for icmp=1:length(Components)
            clear sac_file
            sac_file = sac_files(icmp).sac_file;
            %             [YR,JDY,HR,MN,SEC,MSEC,B,dt,npts] = lh(sac_file,...
            %                 'NZYEAR','NZJDAY','NZHOUR','NZMIN','NZSEC','NZMSEC','B','DELTA','NPTS');
            YR = sac_file(71,3);
            JDY = sac_file(72,3);
            HR = sac_file(73,3);
            MN = sac_file(74,3);
            SEC = sac_file(75,3);
            MSEC = sac_file(76,3);
            B = sac_file(6,3);
            dt = sac_file(1,3);
            npts = sac_file(80,3);
            % VED ADDED FOR SYNTHETICS:
            %YR = EQ.Year; JDY = EQ.Julian_Day; HR = EQ.Hour; MN = EQ.Minute;
            %SEC = EQ.Second; MSEC = 0;
            
            [start_times_abs(icmp)] = Find_Absolute_Time(dt,YR,JDY,HR,MN,SEC,MSEC)+double(B);
            sps = 1/dt;
            dt = 1/double(sps);
            Npts(icmp) = double(npts);
        end
        
        [origin_time] = Find_Absolute_Time(dt,EQ.Year,EQ.Julian_Day,...
            EQ.Hour,EQ.Minute,EQ.Second,0);
        
        
        start_times     = floor((start_times_abs-origin_time)/dt)*dt;
        max_npts        = max(Npts);
        min_start_time  = round(min(start_times)/dt)*dt;
        START_TIME      = round(min_start_time/dt)*dt;
        ABS_TIME        = START_TIME:dt:START_TIME+(max_npts-1)*dt;
        if (ABS_TIME(end)-ABS_TIME(1))<140
            time_pre_add_2  = ABS_TIME(1)-140:dt:ABS_TIME(1)-dt;
            time_post_add_2 = ABS_TIME(end)+dt:dt:ABS_TIME(end)+140;
            ABS_TIME        = [time_pre_add_2 ABS_TIME time_post_add_2]';
        end
        NPTS = length(ABS_TIME);
        
        for icmp=1:length(Components)
            
            clear sac_file sac_file_temp abs_time
            sac_file_temp = sac_files(icmp).sac_file;
            sac_file_temp(:,2) = sac_file_temp(:,2)-mean(sac_file_temp(:,2));
            npts = length(sac_file_temp(:,1));
            abs_time = start_times(icmp):dt:start_times(icmp)+(npts-1)*dt;
            
            clear indx*
            indx1 = find(abs(ABS_TIME-abs_time(1))<(dt/1e3));
            indx2 = find(abs(ABS_TIME-abs_time(end))<(dt/1e3));
            if ~isempty(indx1)
                if length(indx1)>1; indx(1) = indx1(end); else indx(1) = indx1; end
            else indx(1) = 1;
            end
            if ~isempty(indx2)
                if length(indx2)>1; indx(2) = indx2(1); else indx(2) = indx2; end
            else indx(2) = NPTS;
            end
            if indx(1)~=1
                time_pre_add = ABS_TIME(1):dt:(abs_time(1)-dt);
                amp_pre_add = zeros(size(time_pre_add));
            else
                time_pre_add = [];
                amp_pre_add = [];
            end
            if indx(2)~=NPTS
                time_post_add = (abs_time(end)+dt):dt:ABS_TIME(end);
                amp_post_add = zeros(size(time_post_add));
            else
                time_post_add = [];
                amp_post_add = [];
            end
            clear sac_time sac_amp
            sac_time = [time_pre_add abs_time time_post_add]';
            sac_amp = [amp_pre_add sac_file_temp(:,2)' amp_post_add]';
            
            sac_file(:,1) = sac_time;
            sac_file(:,2) = sac_amp;
            samp_mis = length(sac_file(:,1))-NPTS;
            if samp_mis<0
                sac_file_temp(NPTS,3) = 0;
                sac_file_temp(1:length(sac_file(:,1)),1:2) = sac_file;
                sac_file_temp(NPTS+samp_mis+1:NPTS,1) = sac_file(end,1)+dt:dt:sac_file(end,1)-samp_mis*dt;
                clear sac_file
                sac_file = sac_file_temp;
            elseif samp_mis>0
                sac_file_temp(NPTS,3) = 0;
                sac_file_temp(1:NPTS,1:2) = sac_file(1:NPTS,:);
                clear sac_file
                sac_file = sac_file_temp;
            end
            sac_file(1:length(sac_file_temp),3) = sac_file_temp(:,3);
            
            B = sac_file(1,1);
            E = sac_file(end,1);
            OSECc = num2str(EQ.Second);
            if EQ.Second<10;              OSEC = str2double(OSECc(1));
                if length(OSECc)==1;      OMSEC = 0;
                elseif length(OSECc)==3;  OMSEC = str2double([OSECc(3) '00']);
                elseif length(OSECc)==4;  OMSEC = str2double([OSECc(3:4) '0']);
                elseif length(OSECc)>=5;  OMSEC = str2double(OSECc(3:5));
                end
            elseif EQ.Second>=10;         OSEC = str2double(OSECc(1:2));
                if length(OSECc)==2;      OMSEC = 0;
                elseif length(OSECc)==4;  OMSEC = str2double([OSECc(4) '00']);
                elseif length(OSECc)==5;  OMSEC = str2double([OSECc(4:5) '0']);
                elseif length(OSECc)>=5;  OMSEC = str2double(OSECc(4:6));
                end
            end
            %     sac_file = ch(sac_file,'NPTS',NPTS,'B',B,'E',E,'NZYEAR',EQ.Year,'NZJDAY',EQ.Julian_Day,...
            %         'NZHOUR',EQ.Hour,'NZMIN',EQ.Minute,'NZSEC',OSEC,'NZMSEC',OMSEC);
            sac_file(80,3)=NPTS;
            sac_file(6,3)=B;
            sac_file(7,3)=E;
            sac_file(71,3)=EQ.Year;
            sac_file(72,3)=EQ.Julian_Day;
            sac_file(73,3)=EQ.Hour;
            sac_file(74,3)=EQ.Minute;
            sac_file(75,3)=OSEC;
            sac_file(76,3)=OMSEC;
            
            sac_files(icmp,1).sac_file = sac_file;
            
        end
        
        %[B,E] = lh(sac_file,'B','E');
        B = sac_file(6,3);
        E = sac_file(7,3);
        B = double(B);  E = double(E);
        time = sac_file(:,1);
        
        % Decimate the data to a sampling rate of 20 sps (dt_dec=0.05s)
        dt_dec = 0.05;
        Bd = double(int32(B/dt_dec))*dt_dec;
        Ed = double(int32(E/dt_dec))*dt_dec;
        time_dec = Bd:dt_dec:Ed;
        
        dift = diff(time);
        if min(dift)<=0
            tseq = -1;
        else
            tseq = 1;
        end
        
        if tseq == 1
            for icmp=1:length(Components)
                sac_files_temp(icmp,1).sac_file(:,1)     = time_dec;
                sac_files_temp(icmp,1).sac_file(:,2)     = interp1(time,sac_files(icmp).sac_file(:,2),time_dec);
                sac_files_temp(icmp,1).sac_file(1,2)     = sac_files(icmp).sac_file(1,2);
                sac_files_temp(icmp,1).sac_file(end,2)   = sac_files(icmp).sac_file(end,2);
                sac_files_temp(icmp,1).sac_file(1,3)     = dt_dec;
                sac_files_temp(icmp,1).sac_file(2:306,3) = sac_files(icmp).sac_file(2:306,3);
            end
        else
            sac_files_temp = [];
        end
        clear sac_files
        sac_files = sac_files_temp;
        dt = dt_dec;
        
    case 'synth'
        tsynth=sac_files(2).sac_file(:,1); tsynthinds=find(tsynth>-2 & tsynth<10);
        dt=abs(diff(tsynth(1:2))); origin_time=0;
        switch varargin{2}
            case 'Sp'; [~,maxind]=max(abs(sac_files(2).sac_file(tsynthinds,2)));
            case 'Ps'; [~,maxind]=max(abs(sac_files(1).sac_file(tsynthinds,2)));
        end
        origin_time  	= tsynth(tsynthinds(maxind));
        ntinds=abs(round(origin_time/dt));
        if origin_time>0
            for k=1:3;
                sac_files(k).sac_file(:,2)= ...
                    [sac_files(k).sac_file(ntinds+1:end,2); zeros(ntinds,1)];
            end
        elseif origin_time<0
            for k=1:3;
                sac_files(k).sac_file(:,2)=[zeros(ntinds,1);...
                    sac_files(k).sac_file(1:end-ntinds,2)];
            end
        end
        
        
end

    function [ABS_TIME] = Find_Absolute_Time(dt,YR,JDY,HR,MN,SEC,MSEC)
        
        % ********* Function Description *********
        %
        % This will determine the absolute start
        % time for a sac file from calendar year
        % 1980 and day 0 in seconds (i.e., from
        % Julian day 0 of year 1980 at hour, min,
        % and sec 0)
        %
        %
        % ****************************************
        % *                                      *
        % *  Written by David L. Abt - May 2008	 *
        % *                                      *
        % *  Email: David_Abt@brown.edu          *
        % *                                      *
        % ****************************************
        
        sec_per_yr      = 365*24*3600;                  % Seconds in a non-leap year
        sec_per_lpyr    = 366*24*3600;                  % Seconds in a leap year
        sec_per_4yrs    = 3*sec_per_yr+sec_per_lpyr;    % Seconds in a four-year period
        
        dt = double(dt);
        YR = double(YR)-1980;
        JDY = double(JDY);
        HR = double(HR);
        MN = double(MN);
        SEC = double(SEC);
        MSEC = double(MSEC);
        
        SEC = SEC+(MSEC*1e-3);
        SEC = round(SEC/dt)*dt;
        
        lpyr_cycles = floor((YR-1)/4);
        if lpyr_cycles<0
            lpyr_cycles = 0;
        end
        extra_years = YR-1-lpyr_cycles*4;
        if extra_years<0
            extra_years = 0;
        end
        
        ABS_TIME =  (sec_per_4yrs*lpyr_cycles)+ ... % Seconds through the last leap year
            (sec_per_yr*extra_years)+ ...   % Seconds from last leap year to the end of the previous year
            ((JDY-1)*24*3600)+ ...          % Seconds from the first of this year up to "yesterday"
            (HR*3600+MN*60+SEC);            % Seconds "today"
        
    end

end

function [Ray_Paths] = Find_TauP_Paths(EQ,STA,origin_time,varargin)

% ********* Function Description *********
%
% Find the ray path information predicted
% by TauP.  This uses the matTaup function
% taupPath, and matTaup must be installed
% on the computer for this to work.
%
% Go here to download matTaup:
% http://www.ess.washington.edu/SEIS/FMI/matTaup.htm
%
% matTaup written by Qin Li, and based on
% the java program Taup by Philip Crotwell
% http://www.seis.sc.edu/TauP/
%
%
% ****************************************
% *                                      *
% *  Written by David L. Abt - May 2008	 *
% *                                      *
% *  Email: David_Abt@brown.edu          *
% *                                      *
% ****************************************

import edu.sc.seis.TauP.*   % Import the TauP package and load velocity model

REM             = 'ak135'; % 'prem'; % 'iasp91';
Ray_Paths.REM   = REM;

BAZ = azimuth('gc',STA.Latitude,STA.Longitude,EQ.Latitude,EQ.Longitude);
AZ = BAZ+180;
if AZ>360
    AZ = AZ-360;
end

Phases = {'P','S'};
for iph=1:length(Phases)
    Phase = Phases{iph};
    
    if(~isempty(varargin) && strcmp(varargin{1},'ArrayPicker'))
        Taup_Paths = Matlab_TauP('Time',REM,EQ.Depth,Phase,...
            'sta',[STA.Latitude STA.Longitude],...
            'evt',[EQ.Latitude EQ.Longitude]);
        if ~isempty(Taup_Paths) && length(Taup_Paths)>1
            Taup_Paths = Matlab_TauP('Path',REM,EQ.Depth,Phase,...
                'sta',[STA.Latitude STA.Longitude],...
                'evt',[EQ.Latitude EQ.Longitude]);
        end
    else
        Taup_Paths = Matlab_TauP('Path',REM,EQ.Depth,Phase,...
            'sta',[STA.Latitude STA.Longitude],...
            'evt',[EQ.Latitude EQ.Longitude]);
    end
    
    if ~isempty(Taup_Paths)
        if length(Taup_Paths)==1
            Taup_Path = Taup_Paths;
        elseif length(Taup_Paths)>1 % && ~strcmp(varargin{1},'ArrayPicker')
            for itpi=1:length(Taup_Paths)
                latmis = abs(Taup_Paths(itpi).path.latitude(end)-STA.Latitude);
                lonmis = abs(Taup_Paths(itpi).path.longitude(end)-STA.Longitude);
                Paths_mis(itpi,1) = sqrt(latmis^2+lonmis^2);
            end
            Taup_Path = Taup_Paths(Paths_mis==min(Paths_mis));
        end
        
        % Taup_Path.rayParam is in sec/rad, so to compare with output from
        % e.g. taup_time (sec/deg), need to multiply by 180/pi
        
        % We then calculate Taup_Path.ray_parameter
        % should be *km/rad = *rad/deg*deg/km = 180/pi*1/111.19 = 1/6371
        % where deg/km=1/(km/deg)=1/(circumference/360)=1/(2*pi*r/360) = 1/(r * pi/180)
        Taup_Path.ray_parameter        	= Taup_Path.rayParam*pi/180/(6371*pi/180);	% Ray parameter in s/km
        Taup_Path.TauP_abs_arrival_time	= Taup_Path.time+origin_time;
        Taup_Path.AZ                   	= AZ;
        Taup_Path.BAZ                  	= BAZ;
    elseif isempty(Taup_Paths)
        Taup_Path = [];
    end
    Ray_Paths.([Phase '_Path'])	= Taup_Path;
end

end

function [sac_files,Components] = Rotate_NE_into_RT(sac_files,Components)

% ********* Function Description *********
%
% Rotate the N and E components into R
% and T based on the azimuth of the wave.
%
%
% ****************************************
% *                                      *
% *  Written by David L. Abt - May 2008	 *
% *                                      *
% *  Email: David_Abt@brown.edu          *
% *                                      *
% ****************************************

for icmp=1:length(Components)
    sac_file = sac_files(icmp).sac_file;
    if strcmp(Components(icmp),'BHN')
        Ncomp = sac_file(:,2);
        %AZ = lh(sac_file,'AZ','BAZ');
        AZ = sac_file(52,3);
    elseif strcmp(Components(icmp),'BHE')
        Ecomp = sac_file(:,2);
    end
end
rot = AZ*pi/180; % AZ should now be the same as BAZ+180 (see Load_SAC_Files and Find_TauP_Paths)
ROT = [cos(rot) sin(rot); -sin(rot) cos(rot)];

X = ROT*[Ncomp'; Ecomp'];
Rcomp = X(1,:)';
Tcomp = X(2,:)';

sac_files(4,1).sac_file = sac_file;
sac_files(4,1).sac_file(:,2) = Rcomp;

sac_files(5,1).sac_file = sac_file;
sac_files(5,1).sac_file(:,2) = Tcomp;

Components{4,:} = 'BHR';
Components{5,:} = 'BHT';

end

function [bpf_amps] = Filter_Waveforms(dt,sac_files,Components,STA,lowT,highT)

% ********* Function Description *********
%
% Filter the waveforms to here helps
% better locate the desired analysis and
% phase windows.
%
%
% ****************************************
% *                                      *
% *  Written by David L. Abt - May 2008	 *
% *                                      *
% *  Email: David_Abt@brown.edu          *
% *                                      *
% ****************************************

Phases 	= {'P','S'};
% Set up filter
nyq   = 0.5/dt; % Nyquist frequency

Initial_Filter.P = [1/highT,1/lowT];
Initial_Filter.S = Initial_Filter.P;


for iph=1:length(Phases)
    Phase = Phases{iph};
    lf = Initial_Filter.(Phase)(1);
    hf = Initial_Filter.(Phase)(2);
    wn = [lf/nyq hf/nyq];
    [b,a] = butter(2,wn);
    % Filter waveforms
    for icmp=1:length(Components)
        
        
        clear sac_file; sac_file = sac_files(icmp).sac_file;
        if ~isempty(sac_file)
            % bpf_amps -> bandpass filtered amplitudes
            bpf_amps(icmp,iph).bpf_amp = filtfilt(b,a,sac_file(:,2));
            bpf_amps(icmp,iph).lf = lf;
            bpf_amps(icmp,iph).hf = hf;
        end
        
    end
end
end

function [Waveform_Data, Bad] = Set_Anaylsis_Windows(dt,Waveform_Data,sac_files,...
    bpf_amps,Comps,Ray_Paths,varargin)

% ********* Function Description *********
%
% Set the analysis and phase windows that
% will be used in the deconvolution to
% generate the receiver function.
%
% The phase window is the window that
% contains just the direct phase.  This is
% what will be deconvolved from the
% analysis window.
%
%
% ****************************************
% *                                      *
% *  Written by David L. Abt - May 2008	 *
% *                                      *
% *  Email: David_Abt@brown.edu          *
% *                                      *
% ****************************************

tag=varargin{end};
if strcmp(tag,'prep')
    TTarray=varargin{1}; Phases = ['Ps'; 'Sp'];
elseif strcmpi(tag,'synth')
    Phases=varargin{1}; shift = varargin{2}; pwl=varargin{3};
elseif strcmpi(tag,'arraypicker')
    nS=varargin{1}; nP=varargin{2};
    if(nS>=5 && nP>=5); Phases = ['Ps'; 'Sp'];
    elseif(nS>=5); Phases = 'Sp';
        Waveform_Data.Ps_Waveforms.S2N_ratio = -666;
    elseif(nP>=5); Phases = 'Ps';
        Waveform_Data.Sp_Waveforms.S2N_ratio = -666;
    end
elseif strcmp(tag,'prepSp')
    TTarray=varargin{1}; Phases = 'Sp';
end

time = sac_files(1).sac_file(:,1)';

for iph=1:size(Phases,1)
    wrong='fine';
    Phase = Phases(iph,:);
    clear sac_file bpf_amp_dec *S2N* good_indx
    
    if strcmpi(tag,'synth')
        if(strcmp(Phase(1),'S'))
            Continue = 1;
            Taup_phase_arrival = 0;
            pre_phase_t = shift;     post_phase_t = shift/2;
            Min_S2N = 2e6;
            phase_wl = 2*pwl;  % "signal" window length (sec)
        elseif(strcmp(Phase(1),'P'))
            Continue = 1;
            Taup_phase_arrival = 0;
            pre_phase_t = shift;     post_phase_t = shift*2;
            Min_S2N = 2e6;
            phase_wl = 3;  % "signal" window length (sec)
        end
        
        Taup_analysis_window = Taup_phase_arrival+[-pre_phase_t post_phase_t];
        Taup_phase_window = Taup_phase_arrival+phase_wl*[-1 1];
    else
        % for prepping or array picking
        if strcmp(Comps,'ZR');      CMPs = [1 4];
        elseif strcmp(Comps,'ZT');  CMPs = [1 5];
        elseif strcmp(Comps,'PSV');
            if iph==1;              CMPs = [6 7];
            elseif iph==2;          CMPs = [8 9];
            end
        end
        
        % Ps window parameters
        if(strcmp(Phase(1),'P') &&  ...
                (strcmp(Comps(1),'Z') || isfield(Waveform_Data.Ps_Waveforms,'Window')))
            Continue = 1;
            bpf_amp         = bpf_amps(CMPs(1),iph).bpf_amp;
            pre_phase_t     = 150;    post_phase_t    = 150;
            Taup_phase_wl   = 30;    	% Arbitrary phase window length for Taup estimates
            swl             = 5;        % "signal" window length (sec)
            nwl             = 20;   	% "noise" window length (sec) Junlin
            s2n_phase_wl    = 20;      	% Arbitrary phase window length for S2N
            
            if(strcmpi(varargin{end},'prep') && ~isempty(TTarray.P))
                % We set the Min_S2N so large that the search is never carried out
                % for optimal time window. Instead, the average Parr and Sarr times
                % for this particular station and the backazimuth to this
                % particular earthquake are used.
                
                search_range    = 15;%0.1;
                tmp_ved         = TTarray.P(Ray_Paths.P_Path.ray_parameter,Ray_Paths.P_Path.bAzimuth);
                if(~isempty(tmp_ved))
                    if(isnan(tmp_ved) || isinf(tmp_ved)), tmp_ved = -666; end;
                else  tmp_ved = -666;
                end
                
                Taup_phase_arrival = Ray_Paths.P_Path.time - tmp_ved;
                Taup_analysis_window    = Taup_phase_arrival+[-pre_phase_t post_phase_t];
                Min_S2N         = 2e6; % HUGE! so that it never searches
            end
            
            if(strcmpi(varargin{end},'arraypicker') || exist('tmp_ved','var') && tmp_ved==-666 || ...
                    (strcmpi(varargin{end},'prep') && isempty(TTarray.P)))
                search_range    = 15;
                Taup_analysis_window    = Ray_Paths.P_Path.time+[-pre_phase_t post_phase_t];
                Taup_phase_arrival  	= Ray_Paths.P_Path.time;
                Min_S2N         = 2;
            end
        elseif strcmp(Phase(1),'P'); Continue = 0;
        end
        
        % Sp window parameters
        if(strcmp(Phase(1),'S') &&  ...
                (strcmp(Comps(1),'Z') || isfield(Waveform_Data.Sp_Waveforms,'Window')))
            Continue = 1;
            % Use the R/SV component
            bpf_amp = bpf_amps(CMPs(2),iph).bpf_amp;
            pre_phase_t     = 150;         post_phase_t    = 150;
            Taup_phase_wl   = 30;     	% Arbitrary phase window length for Taup estimates
            swl             = 5;   	% "signal" window length (sec)
            nwl             = 20;   	% "noise" window length (sec) Junlin
            s2n_phase_wl    = 20;      	% Arbitrary phase window length for S2N
            
            if ((strcmpi(varargin{end},'prep') || strcmpi(varargin{end},'prepSp')) && ~isempty(TTarray.S))
                % We set the Min_S2N so large that the search is never carried out
                % for optimal time window. Instead, the average Parr and Sarr times
                % for this particular station and the backazimuth to this
                % particular earthquake are used.
                search_range    = 15;
                tmp_ved         = TTarray.S(Ray_Paths.S_Path.ray_parameter,Ray_Paths.S_Path.bAzimuth);
                if(~isempty(tmp_ved))
                    if(isnan(tmp_ved) || isinf(tmp_ved)), tmp_ved = -666; end
                else  tmp_ved = -666;
                end
                Taup_phase_arrival = Ray_Paths.S_Path.time - tmp_ved;
                Taup_analysis_window    = Taup_phase_arrival+[-pre_phase_t post_phase_t];
                Min_S2N         = 2e6; % HUGE! So it snever searches
            end
            
            if(strcmpi(varargin{end},'arraypicker') || exist('tmp_ved','var') && tmp_ved==-666 || ...
                    (strcmpi(varargin{end},'prep') && isempty(TTarray.S)))
                search_range    = 15; %changed from 6 to 12 by njm
                Taup_analysis_window    = Ray_Paths.S_Path.time+[-pre_phase_t post_phase_t];
                Taup_phase_arrival      = Ray_Paths.S_Path.time;
                Min_S2N         = 2;
            end
        elseif strcmp(Phase(1),'S');  Continue=0;
        end
    end
    
    if Continue==1
        clear Window
        Window.pre_phase_time       = pre_phase_t;
        Window.post_phase_time      = post_phase_t;
        Window.Taup_phase_arrival	= Taup_phase_arrival;
        Window.Min_S2N              = Min_S2N;
        
        if strcmpi(tag,'synth')
            indx_full_phase       	= find(time>=Taup_phase_window(1) & ...
                time<=Taup_phase_window(2));
            indx_full_analysis    	= find(time>=Taup_analysis_window(1) & ...
                time<=Taup_analysis_window(2));
            Window.phase_arrival  	= Taup_phase_arrival;
            Window.phase_window    	= Taup_phase_window;
            Window.analysis_window 	= Taup_analysis_window+Window.phase_arrival;
            Window.Taup_S2N_Misfit	= 0;
            Window.indx_full_phase    	= [indx_full_phase(1) indx_full_phase(end)];
            Window.indx_full_analysis  	= [indx_full_analysis(1) indx_full_analysis(end)];
            
            Waveforms.dt        = dt;
            if iph==1
                % Get the time and components in the analysis window
                Waveforms.Time 	= time(indx_full_analysis);
                Waveforms.P     = sac_files(1).sac_file(indx_full_analysis,2)';
                Waveforms.SV    = sac_files(2).sac_file(indx_full_analysis,2)';
                Waveforms.SH    = sac_files(3).sac_file(indx_full_analysis,2)';
            end
            Waveforms.Window   	= Window;
            Waveforms.S2N_ratio	= 5; % arbitrary S2N for synthetics
            
        else
            % Get the TauP predicted arrival time and analysis window
            Taup_phase_window           = Taup_phase_arrival+[0 Taup_phase_wl];
            
            % Generate an envelope function around the taup window and
            % calculate signal:noise for a pair of moving windows to
            % objectively find the arrival time.
            S2N_range   = search_range;
            indx_envl   = find(time>=(Taup_phase_arrival-S2N_range-pre_phase_t) & ...
                time<=(Taup_phase_arrival+S2N_range+post_phase_t));
            indx_big   = find(time>=(Taup_phase_arrival-pre_phase_t) & ...
                time<=(Taup_phase_arrival+post_phase_t));
            indx_S2N    = find(time>=(Taup_phase_arrival-S2N_range) & ...
                time<=(Taup_phase_arrival+S2N_range));
            indx_search = find(time>=(Taup_phase_arrival-search_range) & ...
                time<=(Taup_phase_arrival+search_range));
            nwlp        = nwl/dt;         % number of sample points in noise window
            swlp        = swl/dt;         % number of sample points in signal window
            
            if ~isempty(indx_search) && indx_S2N(1)>indx_envl(1)+double(nwlp) && ...
                    indx_S2N(end)<indx_envl(end)-swlp+1 && ~isempty(bpf_amp)
                
                Hilbert    	= hilbert(bpf_amp(indx_envl));
                Envelope  	= abs(Hilbert);
                time_envl  	= time(indx_envl);
                indx_big    = [max(find(indx_envl==indx_big(1)),1+nwlp) min(find(indx_envl==indx_big(end)),length(indx_envl)-swlp+1)];
                %time_big    = time(indx_big(1):indx_big(2));
                %indx_S2N   	= [find(indx_envl==indx_S2N(1)) find(indx_envl==indx_S2N(end))];
                indx_search	= [find(indx_envl==indx_search(1)) find(indx_envl==indx_search(end))];
                S2N         = zeros(size(Envelope));
                %for iw=indx_S2N(1):indx_S2N(2)
                for iw=indx_big(1):indx_big(2)
                    Noise   = mean(Envelope(iw-nwlp:iw-1));
                    Signal  = mean(Envelope(iw:iw+swlp-1));
                    S2N(iw) = Signal/Noise;
                end
                
                
                HWL             = swl;         	% Hanning window length
                HW              = hann(HWL/dt);	% Hanning window to be convolved with the S2N time series
                S2N_hann_temp   = conv(S2N,HW);	% Filtered S2N
                hann_buff       = (HWL/2)/dt;
                S2N_hann        = S2N_hann_temp(hann_buff:end-hann_buff);
                S2N_filt_norm   = S2N_hann/max(S2N_hann)*max(S2N);
                
                % Find the max S2N within the test window (some distance around the
                % Taup predicted arrival)
                S2N_max_indx = find(S2N_filt_norm==max(S2N_filt_norm(indx_search(1):indx_search(2))));
                if length(S2N_max_indx)>1
                    if length(S2N_max_indx)==2 && S2N_max_indx(1)==S2N_max_indx(2)+1
                        S2N_max_indx = S2N_max_indx(1);
                    else
                        S2N_max_indx    = 1;
                        wrong='More than one S2N peak';
                        S2N_filt_norm   = -666;
                    end
                elseif isempty(S2N_max_indx)
                    S2N_max_indx    = 1;
                    wrong='No S2N peak found';
                    S2N_filt_norm   = -666;
                end
                S2N_max = S2N_filt_norm(S2N_max_indx);
                
                if S2N_filt_norm ~= -666
                    S2N_max_time_envl = time_envl(S2N_filt_norm==max(S2N_filt_norm));
                    S2N_max_envl = max(S2N_filt_norm);
                    if Taup_phase_arrival-S2N_max_time_envl>S2N_range && S2N_max_envl>2 %&& S2N_max_time_envl>time_envl(1)+post_phase_t-50
                        %                         S2N_max_indx    = 1;
                        %                         wrong='No S2N peak found';
                        %                         S2N_filt_norm   = -666;
                        %                         S2N_max = S2N_filt_norm(S2N_max_indx);
                        S2N_max = 1e-6;
                    end
                end
                
                
                
                
                if S2N_max>=Min_S2N %&& strcmp(STA,'COEN')
                    % The SN2 maximum coincides with the best arrival pick.
                    % This will be our S2N_phase_arrival, and we pick the
                    % end time of the phase window by finding the width of
                    % the S2N peak at its inflection points and adding it
                    % to the phase arrival time.  The inflection points are
                    % on either side of the peak where the 2nd derivative
                    % of the S2N function becomes positive.
                    S2N_peak_width  = s2n_phase_wl/dt;
                    S2N_phase_arrival       = time_envl(S2N_max_indx-2/dt);
                    if S2N_max_indx+S2N_peak_width>length(time_envl)
                        S2N_phase_end       = time_envl(end);
                    else
                        S2N_phase_end       = time_envl(S2N_max_indx+S2N_peak_width);
                    end
                    S2N_phase_window      	= [S2N_phase_arrival S2N_phase_end];
                    S2N_analysis_window    	= S2N_phase_arrival+[-pre_phase_t post_phase_t];
                    
                    indx_full_phase       	= find(time>=S2N_phase_window(1) & ...
                        time<=S2N_phase_window(2));
                    indx_full_analysis    	= find(time>=S2N_analysis_window(1) & ...
                        time<=S2N_analysis_window(2));
                    Window.phase_arrival  	= S2N_phase_arrival;
                    Window.phase_window   	= S2N_phase_window;
                    Window.analysis_window 	= S2N_analysis_window;
                    Window.Taup_S2N_Misfit 	= Taup_phase_arrival-S2N_phase_arrival;
                else
                    indx_full_phase       	= find(time>=Taup_phase_window(1) & ...
                        time<=Taup_phase_window(2));
                    indx_full_analysis    	= find(time>=Taup_analysis_window(1) & ...
                        time<=Taup_analysis_window(2));
                    Window.phase_arrival  	= Taup_phase_arrival;
                    Window.phase_window    	= Taup_phase_window;
                    Window.analysis_window 	= Taup_analysis_window;
                    Window.Taup_S2N_Misfit	= 0;
                end
                S2N_ratio = S2N_max;
                
                Window.indx_full_phase    	= [indx_full_phase(1) indx_full_phase(end)];
                Window.indx_full_analysis  	= [indx_full_analysis(1) indx_full_analysis(end)];
                
                Waveforms.dt        = dt;
                if strcmp(Comps(1),'P')
                    % Get the time and components in the analysis window
                    Waveforms.Time 	= time(indx_full_analysis);
                    Waveforms.Z   	= sac_files(1).sac_file(indx_full_analysis,2)';
                    Waveforms.N   	= sac_files(2).sac_file(indx_full_analysis,2)';
                    Waveforms.E    	= sac_files(3).sac_file(indx_full_analysis,2)';
                    Waveforms.R    	= sac_files(4).sac_file(indx_full_analysis,2)';
                    Waveforms.T    	= sac_files(5).sac_file(indx_full_analysis,2)';
                    Waveforms.P     = sac_files(6+2*(iph-1)).sac_file(indx_full_analysis,2)';
                    Waveforms.SV    = sac_files(7+2*(iph-1)).sac_file(indx_full_analysis,2)';
                elseif strcmp(Comps(1),'Z')
                    Waveforms.Time 	= time(indx_full_analysis);
                    Waveforms.Z   	= sac_files(1).sac_file(indx_full_analysis,2)';
                    Waveforms.N   	= sac_files(2).sac_file(indx_full_analysis,2)';
                    Waveforms.E    	= sac_files(3).sac_file(indx_full_analysis,2)';
                    Waveforms.R    	= sac_files(4).sac_file(indx_full_analysis,2)';
                end
                
                Waveforms.Window   	= Window;
                Waveforms.S2N_ratio	= S2N_ratio;
                if 0 %strcmp(Phase,'Ps')
                    [tmax,thalf1,thalf2] = plotWindow(Window,Envelope,time_envl,Waveforms);
                    disp(num2str((thalf2-tmax)/2));
                end
            else
                S2N_ratio       	= -666;
                wrong='Arrival not included in waveform';
                Waveforms.S2N_ratio	= S2N_ratio;
            end
            
            
        end
        if strcmp(Phase(1),'P')
            Waveform_Data.Ps_Waveforms = Waveforms;
            clear Waveforms; Waveform_Data.Ps_Waveforms.wrong=wrong;
        elseif strcmp(Phase(1),'S')
            Waveform_Data.Sp_Waveforms = Waveforms;
            clear Waveforms; Waveform_Data.Sp_Waveforms.wrong=wrong;
        end
    end
    Bad.(Phase)=wrong;
end
end

function [Waveform_Data] = Check_Polarity_Reversal(Waveform_Data)
% **
% *

Phases = {'Ps','Sp'};
for iph=1:length(Phases)
    Phase = Phases{iph};
    Waveforms = Waveform_Data.([Phase '_Waveforms']);
    if isfield(Waveforms,'Z') && isfield(Waveforms,'R')
        
        Z       = Waveforms.Z-mean(Waveforms.Z);
        R       = Waveforms.R-mean(Waveforms.R);
        dt      = Waveforms.dt;
        Time	= Waveforms.Time;
        P_window    = Waveforms.Window.phase_window;
        P_arrival   = Waveforms.Window.phase_arrival;
        A_window	= Waveforms.Window.analysis_window;
        indx_P      = find(Time>=P_window(1) & Time<=P_window(2));
        indx_A      = find(Time>=A_window(1) & Time<=A_window(2));
        Time_grid   = Time-P_arrival;
        indx_p      = find(abs(Time_grid)==min(abs(Time_grid)));
        
        indx_X      = find(Time>=P_window(1) & Time<=(P_window(1)+10));
        nn          = length(indx_X);
        Tcor        = -dt*(nn-1):dt:dt*(nn-1);
        NN          = length(Tcor);
        
        xCorrZR     = xcorr(Z(indx_X),R(indx_X),'coeff');
        xCorr_0lag  = xCorrZR(nn);
        
        % First derivative of the cross-correlation function
        xCorrZR_d1(1)    = xCorrZR(1)-xCorrZR(NN);
        xCorrZR_d1(2:NN) = xCorrZR(2:NN)-xCorrZR(1:NN-1);
        xCorrZR_d1       = xCorrZR_d1/max(abs(xCorrZR_d1));
        
        trough_lag  = nn*dt;
        trough_val  = 0;
        peak_lag    = nn*dt;
        peak_val    = 0;
        
        if xCorrZR_d1(nn)==0 || ...
                (xCorrZR_d1(nn-1)<0 && xCorrZR_d1(nn+1)>0) || ...
                (xCorrZR_d1(nn-1)>0 && xCorrZR_d1(nn+1)<0)
            
            if xCorr_0lag>=0
                peak_val    = xCorrZR(nn);
                peak_lag    = 0;
                % Find trough to the left of zero-lag
                trough_indxs = [1 NN];
                for ixcr=1:nn-1
                    if xCorrZR_d1(nn-ixcr)<0
                        trough_indxs(1) = nn-ixcr+1;
                        trough_lags(1) = ixcr*dt;
                        break
                    end
                end
                trough_vals(1) = xCorrZR(trough_indxs(1));
                % Find trough to the right of zero-lag
                for ixcr=1:nn-1
                    if xCorrZR_d1(nn+ixcr)<0
                        trough_indxs(2) = nn+ixcr-1;
                        trough_lags(2) = ixcr*dt;
                        break
                    end
                end
                trough_vals(2)  = xCorrZR(trough_indxs(2));
                max_indx        = find(abs(trough_vals)==max(abs(trough_vals)));
                trough_val      = trough_vals(max_indx);
                trough_lag      = trough_lags(max_indx);
            elseif xCorr_0lag<0
                trough_val	= xCorrZR(nn);
                trough_lag  = 0;
                % Find peak to the left of zero-lag
                peak_indxs = [1 NN];
                for ixcr=1:nn-1
                    if xCorrZR_d1(nn-ixcr)<0
                        peak_indxs(1) = nn-ixcr+1;
                        peak_lags(1) = ixcr*dt;
                        break
                    end
                end
                peak_vals(1) = xCorrZR(peak_indxs(1));
                % Find peak to the right of zero-lag
                for ixcr=1:nn-1
                    if xCorrZR_d1(nn+ixcr)<0
                        peak_indxs(2) = nn+ixcr-1;
                        peak_lags(2) = ixcr*dt;
                        break
                    end
                end
                peak_vals(2)	= xCorrZR(peak_indxs(2));
                max_indx        = find(abs(peak_vals)==max(abs(peak_vals)));
                peak_val        = peak_vals(max_indx);
                peak_lag        = peak_lags(max_indx);
            end
            
        else
            
            if xCorrZR_d1(nn)>=0
                % Find trough to the left of zero-lag
                trough_indx = 1;
                for ixcr=1:nn-1
                    if xCorrZR_d1(nn-ixcr)<0
                        trough_indx = nn-ixcr+1;
                        trough_lag = ixcr*dt;
                        break
                    end
                end
                trough_val = xCorrZR(trough_indx);
                % Find peak to the right of zero-lag
                peak_indx = NN;
                for ixcr=1:nn-1
                    if xCorrZR_d1(nn+ixcr)<0
                        peak_indx = nn+ixcr-1;
                        peak_lag = ixcr*dt;
                        break
                    end
                end
                peak_val = xCorrZR(peak_indx);
            elseif xCorrZR_d1(nn)<0
                % Find peak to the left of zero-lag
                peak_indx = 1;
                for ixcr=1:nn-1
                    if xCorrZR_d1(nn-ixcr)>0
                        peak_indx = nn-ixcr+1;
                        peak_lag = ixcr*dt;
                        break
                    end
                end
                peak_val = xCorrZR(peak_indx);
                % Find trough to the right of zero-lag
                trough_indx = NN;
                for ixcr=1:nn-1
                    if xCorrZR_d1(nn+ixcr)>0
                        trough_indx = nn+ixcr-1;
                        trough_lag = ixcr*dt;
                        break
                    end
                end
                trough_val = xCorrZR(trough_indx);
            end
            
        end
        
        lag_tolerance = 1.5;
        if trough_lag<=lag_tolerance && peak_lag<=lag_tolerance
            if abs(trough_val)>peak_val
                Waveform_Data.([Phase '_Waveforms']).ZR_Xcorr_Coeff   = trough_val;
            elseif peak_val>=abs(trough_val)
                Waveform_Data.([Phase '_Waveforms']).ZR_Xcorr_Coeff   = peak_val;
            end
        elseif trough_lag<=lag_tolerance && peak_lag>lag_tolerance
            Waveform_Data.([Phase '_Waveforms']).ZR_Xcorr_Coeff       = trough_val;
        elseif trough_lag>lag_tolerance && peak_lag<=lag_tolerance
            Waveform_Data.([Phase '_Waveforms']).ZR_Xcorr_Coeff       = peak_val;
        elseif trough_lag>lag_tolerance && peak_lag>lag_tolerance
            Waveform_Data.([Phase '_Waveforms']).ZR_Xcorr_Coeff       = xCorrZR(nn);
        end
    else
        Waveform_Data.([Phase '_Waveforms']).ZR_Xcorr_Coeff = 0;
    end
    
end

end


function [FS_VpVs,Vs_info,Vp_info] =Find_VpVs(All_Waveform_Data,All_Ray_Path_Data,lowT,highT)
% find free surface velocity
% Junlin Hua 2019

idim = 181;
jdim = 181;
vps = linspace(2.7,8.1,idim);
vss = linspace(1.5,4.5,jdim);
[b,a] = meshgrid(vss,vps);


eqs = fieldnames(All_Waveform_Data);

phase_win = 7;

vsout = [];
vpout = [];
snrout_P = [];
snrout_S = [];
snrout_P_c = [];
snrout_S_c = [];
snrout_P_d = [];
snrout_S_d = [];
residual_P = [];
vsout_d = [];
vsout_c = [];
residual_P_d = [];
residual_P_c = [];
residual_S = [];
vpout_d = [];
residual_S_d = [];
vpout_c = [];
residual_S_c = [];
corr_P = [];
corr_S = [];
corr_P_c = [];
corr_S_c = [];
corr_P_d = [];
corr_S_d = [];
ind_P = [];
ind_P_c = [];
ind_P_d = [];
ind_S = [];
ind_S_c = [];
ind_S_d = [];

tpr = 1;

dt = 0.05;
for i = 1:length(eqs)
    if All_Waveform_Data.(eqs{i}).Ps_Waveforms.S2N_ratio~=-666
        dt = All_Waveform_Data.(eqs{i}).Ps_Waveforms.dt;
        break;
    end
end
[b_f,a_f] = butter(2,[1/highT/(0.5/dt),1/lowT/(0.5/dt)]);  %filter


%% For P arrival
for ie = 1:length(eqs)
    
    if isempty(All_Waveform_Data.(eqs{ie}).Ps_Waveforms)
        continue;
    end
    if All_Waveform_Data.(eqs{ie}).Ps_Waveforms.S2N_ratio < 5
        continue;
    end
    
    
    R = All_Waveform_Data.(eqs{ie}).Ps_Waveforms.R;
    Z = All_Waveform_Data.(eqs{ie}).Ps_Waveforms.Z;
    inds = All_Waveform_Data.(eqs{ie}).Ps_Waveforms.Window.indx_full_phase(1)-...
        All_Waveform_Data.(eqs{ie}).Ps_Waveforms.Window.indx_full_analysis(1)+1;
    inde = min(All_Waveform_Data.(eqs{ie}).Ps_Waveforms.Window.indx_full_analysis(2),...
        All_Waveform_Data.(eqs{ie}).Ps_Waveforms.Window.indx_full_phase(2))-...
        All_Waveform_Data.(eqs{ie}).Ps_Waveforms.Window.indx_full_analysis(1)+1;
    dt = All_Waveform_Data.(eqs{ie}).Ps_Waveforms.dt;
    RayParamP = All_Ray_Path_Data.(eqs{ie}).P_Path.ray_parameter;
    
    if RayParamP*vps(end)>0.98
        continue;
    end
    
    R_P = filtfilt(b_f,a_f,double(R));
    Z_P = filtfilt(b_f,a_f,double(Z));
    orilen = length(R);
    Cut_P(1) = inds;
    Cut_P(2) = inde;
    R_P = (taper(R_P,1:orilen,tpr/dt,1,Cut_P(1),Cut_P(2)));
    Z_P = (taper(Z_P,1:orilen,tpr/dt,1,Cut_P(1),Cut_P(2)));
    n_pad = round(tpr/dt);
    R_P = R_P(inds-n_pad:min(inde+n_pad,orilen));
    Z_P = Z_P(inds-n_pad:min(inde+n_pad,orilen));
    Time = All_Waveform_Data.(eqs{ie}).Ps_Waveforms.Time(inds-n_pad:min(inde+n_pad,orilen));
    
    Z_env = abs(hilbert(Z_P));
    Z_env(Time>Time(1)+12+tpr) = 0;
    T_peak = Time(Z_env==max(Z_env));
    
    if length(T_peak)>1
        T_peak = T_peak(1);
    end
    
    Taper_P(1) = T_peak-phase_win/2;
    Taper_P(2) = T_peak+phase_win/2;
    
    R_P = R_P(Time>=Taper_P(1) & Time<=Taper_P(2));
    Z_P = Z_P(Time>=Taper_P(1) & Time<=Taper_P(2));
    %Time_P = Time(Time>=Taper_P(1) & Time<=Taper_P(2));
    
    corr = corrcoef(R_P,Z_P);
    corr = abs(corr(1,2));
    if corr<0.65
        continue;
    end
    
    
    kdim = length(R_P);
    R_P = repmat(permute(R_P,[3,1,2]),[idim,jdim,1]);
    Z_P = repmat(permute(Z_P,[3,1,2]),[idim,jdim,1]);
    
    qa_P  = sqrt((a.^-2)-(RayParamP^2));
    qb_P  = sqrt((b.^-2)-(RayParamP^2));
    
    
    A_P   = RayParamP*(b.^2)./a;
    B_P   = ((b.^2)*(RayParamP^2)-0.5)./(a.*qa_P);
    C_P   = (0.5-(b.^2)*(RayParamP^2))./(b.*qb_P);
    D_P   = RayParamP*b;
    
    P_P = repmat(A_P,[1,1,kdim]).*R_P-repmat(B_P,[1,1,kdim]).*Z_P;
    S_P = repmat(C_P,[1,1,kdim]).*R_P-repmat(D_P,[1,1,kdim]).*Z_P;
    
    PS_P = sum(P_P.*S_P,3)./sum(R_P.*Z_P,3);
    PP_P = sum(P_P.*P_P,3)./sum(R_P.*Z_P,3);
    SS_P = sum(S_P.*S_P,3)./sum(R_P.*Z_P,3);
    
    b_test = vss;
    a_test = b_test*1.7;
    p = RayParamP;
    
    
    diffP = zeros(1,jdim);
    diffP1 = zeros(1,jdim);
    diffP2 = zeros(1,jdim);
    
    for i=1:jdim
        
        b0 = b_test(i);
        a0 = a_test(i);
        the = asind(p*a0); phi = asind(p*b0);
        A = (b0^(-2)-2*p^2)^2;
        B = 4*p^2*cosd(the)/a0*cosd(phi)/b0;
        C = 4*p*cosd(the)/b0*(b0^(-2)-2*p^2);
        D = 4*p*cosd(phi)/a0*(b0^(-2)-2*p^2);
        M1 = 1/(A+B)*[2*B*sind(the)+C*cosd(phi),2*A*cosd(phi)+D*sind(the);-2*A*cosd(the)-C*sind(phi),2*B*sind(phi)+D*cosd(the)];
        
        S11 = A_P*M1(1,1)+B_P*M1(2,1);
        S21 = C_P*M1(1,1)+D_P*M1(2,1);
        
        PSRZP_test = S21.*S11/M1(1,1)/M1(2,1);
        PPRZP_test = S11.*S11/M1(1,1)/M1(2,1);
        SSRZP_test = S21.*S21/M1(1,1)/M1(2,1);
        
        diffP(i) = sqrt(sum(sum((PSRZP_test+PS_P).^2+(PPRZP_test+PP_P).^2+...
            (SSRZP_test+SS_P).^2)));
        diffP1(i) = sqrt(sum(sum((PSRZP_test+PS_P).^2)));
        diffP2(i) = sqrt(sum(sum((PPRZP_test+PP_P).^2)));
        
    end
    
    if length(min(diffP))==1
        [~,peakpos] = findpeaks(-diffP);
        if length(peakpos)==1
            vsout = [vsout,vss(peakpos)];
            snrout_P = [snrout_P,All_Waveform_Data.(eqs{ie}).Ps_Waveforms.S2N_ratio];
            residual_P = [residual_P,min(diffP)];
            corr_P = [corr_P,corr];
            ind_P = [ind_P,ie];
        end
    end
    
    if length(min(diffP1))==1
        [~,peakpos_c] = findpeaks(-diffP1);
        if length(peakpos_c)==1
            vsout_c = [vsout_c,vss(peakpos_c)];
            residual_P_c = [residual_P_c,min(diffP1)];
            corr_P_c = [corr_P_c,corr];
            snrout_P_c = [snrout_P_c,All_Waveform_Data.(eqs{ie}).Ps_Waveforms.S2N_ratio];
            ind_P_c = [ind_P_c,ie];
        end
    end
    
    if length(min(diffP2))==1
        [~,peakpos_d] = findpeaks(-diffP2);
        if length(peakpos_d)==1
            vsout_d = [vsout_d,vss(peakpos_d)];
            residual_P_d = [residual_P_d,min(diffP2)];
            corr_P_d = [corr_P_d,corr];
            snrout_P_d = [snrout_P_d,All_Waveform_Data.(eqs{ie}).Ps_Waveforms.S2N_ratio];
            ind_P_d = [ind_P_d,ie];
        end
    end
    
    
end

if length(vsout) > 3
    weight_P = corr_P.*snrout_P;
    FS_VpVs.Vs = sum(weight_P.*vsout)/sum(weight_P);
    FS_VpVs.Vs_std = sqrt(sum((vsout-FS_VpVs.Vs).^2.*weight_P)/sum(weight_P));
else
    weight_P = -666;
    FS_VpVs.Vs = 2.8;
    FS_VpVs.Vp_std = -666;
end

Vs_info.vs = vsout;
Vs_info.snr = snrout_P;
Vs_info.residual = residual_P;
Vs_info.corr = corr_P;
Vs_info.ind = ind_P;
Vs_info.vs_c = vsout_c;
Vs_info.snr_c = snrout_P_c;
Vs_info.residual_c = residual_P_c;
Vs_info.corr_c = corr_P_c;
Vs_info.ind_c = ind_P_c;
Vs_info.vs_d = vsout_d;
Vs_info.snr_d = snrout_P_d;
Vs_info.residual_d = residual_P_d;
Vs_info.corr_d = corr_P_d;
Vs_info.ind_d = ind_P_d;
Vs_info.weight = weight_P;

%% For S arrival
for ie = 1:length(eqs)
    
    if isempty(All_Waveform_Data.(eqs{ie}).Sp_Waveforms)
        continue;
    end
    if All_Waveform_Data.(eqs{ie}).Sp_Waveforms.S2N_ratio < 5
        continue;
    end
    
    
    R = All_Waveform_Data.(eqs{ie}).Sp_Waveforms.R;
    Z = All_Waveform_Data.(eqs{ie}).Sp_Waveforms.Z;
    inds = All_Waveform_Data.(eqs{ie}).Sp_Waveforms.Window.indx_full_phase(1)-...
        All_Waveform_Data.(eqs{ie}).Sp_Waveforms.Window.indx_full_analysis(1)+1;
    inde = min(All_Waveform_Data.(eqs{ie}).Sp_Waveforms.Window.indx_full_analysis(2),...
        All_Waveform_Data.(eqs{ie}).Sp_Waveforms.Window.indx_full_phase(2))-...
        All_Waveform_Data.(eqs{ie}).Sp_Waveforms.Window.indx_full_analysis(1)+1;
    dt = All_Waveform_Data.(eqs{ie}).Sp_Waveforms.dt;
    RayParamS = All_Ray_Path_Data.(eqs{ie}).S_Path.ray_parameter;
    
    if RayParamS*vps(end)>0.98
        continue;
    end
    
    R_S = filtfilt(b_f,a_f,double(R));
    Z_S = filtfilt(b_f,a_f,double(Z));
    orilen = length(R);
    Cut_S(1) = inds;
    Cut_S(2) = inde;
    R_S = (taper(R_S,1:orilen,tpr/dt,1,Cut_S(1),Cut_S(2)));
    Z_S = (taper(Z_S,1:orilen,tpr/dt,1,Cut_S(1),Cut_S(2)));
    n_pad = round(tpr/dt);
    R_S = R_S(inds-n_pad:min(inde+n_pad,orilen));
    Z_S = Z_S(inds-n_pad:min(inde+n_pad,orilen));
    Time = All_Waveform_Data.(eqs{ie}).Sp_Waveforms.Time(inds-n_pad:min(inde+n_pad,orilen));
    
    R_env = abs(hilbert(R_S));
    R_env(Time>Time(1)+12+tpr) = 0;
    T_peak = Time(R_env==max(R_env));
    
    if length(T_peak)>1
        T_peak = T_peak(1);
    end
    
    Taper_S(1) = T_peak-phase_win/2;
    Taper_S(2) = T_peak+phase_win/2;
    
    R_S = R_S(Time>=Taper_S(1) & Time<=Taper_S(2));
    Z_S = Z_S(Time>=Taper_S(1) & Time<=Taper_S(2));
    %Time_S = Time(Time>=Taper_S(1) & Time<=Taper_S(2));
    
    corr = corrcoef(R_S,Z_S);
    corr = abs(corr(1,2));
    if corr<0.65
        continue;
    end
    
    kdim = length(R_S);
    R_S = repmat(permute(R_S,[3,1,2]),[idim,jdim,1]);
    Z_S = repmat(permute(Z_S,[3,1,2]),[idim,jdim,1]);
    
    qa_S  = sqrt((a.^-2)-(RayParamS^2));
    qb_S  = sqrt((b.^-2)-(RayParamS^2));
    
    
    A_S   = RayParamS*(b.^2)./a;
    B_S   = ((b.^2)*(RayParamS^2)-0.5)./(a.*qa_S);
    C_S   = (0.5-(b.^2)*(RayParamS^2))./(b.*qb_S);
    D_S   = RayParamS*b;
    
    P_S = repmat(A_S,[1,1,kdim]).*R_S-repmat(B_S,[1,1,kdim]).*Z_S;
    S_S = repmat(C_S,[1,1,kdim]).*R_S-repmat(D_S,[1,1,kdim]).*Z_S;
    
    PS_S = sum(P_S.*S_S,3)./sum(R_S.*Z_S,3);
    PP_S = sum(P_S.*P_S,3)./sum(R_S.*Z_S,3);
    SS_S = sum(S_S.*S_S,3)./sum(R_S.*Z_S,3);
    
    b0 = FS_VpVs.Vs;
    a_test = vps;
    p = RayParamS;
    
    
    diffS = zeros(1,idim);
    diffS1 = zeros(1,idim);
    diffS2 = zeros(1,idim);
    
    for i=1:idim
        
        
        a0 = a_test(i);
        the = asind(p*a0); phi = asind(p*b0);
        A = (b0^(-2)-2*p^2)^2;
        B = 4*p^2*cosd(the)/a0*cosd(phi)/b0;
        C = 4*p*cosd(the)/b0*(b0^(-2)-2*p^2);
        D = 4*p*cosd(phi)/a0*(b0^(-2)-2*p^2);
        M1 = 1/(A+B)*[2*B*sind(the)+C*cosd(phi),2*A*cosd(phi)+D*sind(the);-2*A*cosd(the)-C*sind(phi),2*B*sind(phi)+D*cosd(the)];
        
        S12 = A_S*M1(1,2)+B_S*M1(2,2);
        S22 = C_S*M1(1,2)+D_S*M1(2,2);
        
        PSRZS_test = S12.*S22/M1(1,2)/M1(2,2);
        PPRZS_test = S12.*S12/M1(1,2)/M1(2,2);
        SSRZS_test = S22.*S22/M1(1,2)/M1(2,2);
        
        diffS(i) = sqrt(sum(sum((PSRZS_test+PS_S).^2+(PPRZS_test+PP_S).^2+...
            (SSRZS_test+SS_S).^2)));
        diffS1(i) = sqrt(sum(sum((PSRZS_test+PS_S).^2)));
        diffS2(i) = sqrt(sum(sum((SSRZS_test+SS_S).^2)));
        
    end
    
    if length(min(diffS))==1
        [~,peakpos] = findpeaks(-diffS);
        if length(peakpos)==1
            vpout = [vpout,vps(peakpos)];
            snrout_S = [snrout_S,All_Waveform_Data.(eqs{ie}).Sp_Waveforms.S2N_ratio];
            residual_S = [residual_S,min(diffS)];
            corr_S = [corr_S,corr];
            ind_S = [ind_S,ie];
        end
    end
    
    if length(min(diffS1))==1
        [~,peakpos_c] = findpeaks(-diffS1);
        if length(peakpos_c)==1
            vpout_c = [vpout_c,vps(peakpos_c)];
            snrout_S_c = [snrout_S_c,All_Waveform_Data.(eqs{ie}).Sp_Waveforms.S2N_ratio];
            residual_S_c = [residual_S_c,min(diffS1)];
            corr_S_c = [corr_S_c,corr];
            ind_S_c = [ind_S_c,ie];
        end
    end
    
    if length(min(diffS2))==1
        [~,peakpos_d] = findpeaks(-diffS2);
        if length(peakpos_d)==1
            vpout_d = [vpout_d,vps(peakpos_d)];
            snrout_S_d = [snrout_S_d,All_Waveform_Data.(eqs{ie}).Sp_Waveforms.S2N_ratio];
            residual_S_d = [residual_S_d,min(diffS2)];
            corr_S_d = [corr_S_d,corr];
            ind_S_d = [ind_S_d,ie];
        end
    end
    
    
end

if length(vpout) > 3
    weight_S = corr_S.*snrout_S;
    FS_VpVs.Vp = sum(weight_S.*vpout)/sum(weight_S);
    FS_VpVs.Vp_std = sqrt(sum((vpout-FS_VpVs.Vp).^2.*weight_S)/sum(weight_S));
else
    weight_S = -666;
    FS_VpVs.Vp = FS_VpVs.Vs*1.8;
    FS_VpVs.Vp_std = -666;
end

Vp_info.vp = vpout;
Vp_info.snr = snrout_S;
Vp_info.residual = residual_S;
Vp_info.corr = corr_S;
Vp_info.ind = ind_S;
Vp_info.vp_c = vpout_c;
Vp_info.snr_c = snrout_S_c;
Vp_info.residual_c = residual_S_c;
Vp_info.corr_c = corr_S_c;
Vp_info.ind_c = ind_S_c;
Vp_info.vp_d = vpout_d;
Vp_info.snr_d = snrout_S_d;
Vp_info.residual_d = residual_S_d;
Vp_info.corr_d = corr_S_d;
Vp_info.ind_d = ind_S_d;
Vp_info.weight = weight_S;

end


function [sac_files,Components] = Rotate_ZR_into_PSV(Waveform_Data,sac_files,Components,Ray_Paths,FS_VpVs)

% ********* Function Description *********
%
% Rotate the Z and R components into P and
% SV components by finding the surface
% velocities the result in the minimum P
% amplitude on the SV component (in the
% predicted P-phase window) and the
% minimum SV amplitude on the P component
% (in the predicted S-phase window) after
% the free surface correction.
%
%
% ****************************************
% *                                      *
% *  Written by David L. Abt - May 2008	 *
% *                                      *
% *  Email: David_Abt@brown.edu          *
% *                                      *
% ****************************************

Z_comp_sac = -sac_files(1).sac_file(:,2);
R_comp_sac = sac_files(4).sac_file(:,2);
%     T_comp_sac = sac_files(5).sac_file(:,2);

Phases = {'Ps','Sp'};
Continue = zeros(length(Phases),1);
for iph=1:2
    Phase = Phases{iph};
    %RP = Ray_Paths.P_Path.ray_parameter;
    RP = Ray_Paths.([Phase(1) '_Path']).ray_parameter;
    
    if isfield(Waveform_Data.([Phase '_Waveforms']),'Window')
        Continue(iph) = 1;
    else
        Continue(iph) = 0;
    end
    
    if Continue(iph)==1
        
        % Decompose vertical/radial into P/SV
        % (i.e., Kennett [1991] and Bostock [1998])
        %
        % !! NOTE !! Z in this equation is positive down !!
        % We simply make Z_comp_sac (above) negative
        
        a   = FS_VpVs.Vp;
        b   = FS_VpVs.Vs;
        
        qa  = sqrt((a^-2)-(RP^2));
        qb  = sqrt((b^-2)-(RP^2));
        A   = RP*(b^2)/a;
        B   = ((b^2)*(RP^2)-0.5)/(a*qa);
        C   = (0.5-(b^2)*(RP^2))/(b*qb);
        D   = RP*b;
        
        P_comp_sac  = A*R_comp_sac+B*Z_comp_sac;
        SV_comp_sac = C*R_comp_sac+D*Z_comp_sac;
        %             SH_comp_sac = E*T_comp_sac;
        
        % The following stores the waveforms rotated to the parent
        % and daughter ray paramaters respectively
        %             if strcmp(Phase,'P') || strcmp(Phase,'Sp')
        %                 Waveform_Data.([Phases{iph1,2} '_Waveforms']).P_comp = P_comp;
        %                 sac_files(6+2*(iph1-1),1).sac_file      = sac_files(1).sac_file;
        %                 sac_files(6+2*(iph1-1),1).sac_file(:,2) = P_comp_sac;
        %             elseif strcmp(Phase,'S') || strcmp(Phase,'Ps')
        %                 Waveform_Data.([Phases{iph1,2} '_Waveforms']).SV_comp = SV_comp;
        %                 sac_files(7+2*(iph1-1),1).sac_file      = sac_files(4).sac_file;
        %                 sac_files(7+2*(iph1-1),1).sac_file(:,2) = SV_comp_sac;
        %             end
        
        % The following stores the waveforms rotated to the parent
        % ray paramater only
        sac_files(6+2*(iph-1),1).sac_file      = sac_files(1).sac_file;
        sac_files(6+2*(iph-1),1).sac_file(:,2) = P_comp_sac;
        sac_files(7+2*(iph-1),1).sac_file      = sac_files(4).sac_file;
        sac_files(7+2*(iph-1),1).sac_file(:,2) = SV_comp_sac;
        Components{6+2*(iph-1),:} = [Phase '_P'];
        Components{7+2*(iph-1),:} = [Phase '_SV'];
        
    end
    
end
end


function [Waveform_Data_Reduced,Ray_Path_Data_Reduced] = Reduce_File_Size(Waveform_Data_Full,Ray_Path_Data_Full,varargin)

% ********* Function Description *********
%
% Reduce the size of the data file that is
% saved by making the waveform data and
% time single precision.
%
%
% ****************************************
% *                                      *
% *  Written by David L. Abt - May 2008	 *
% *                                      *
% *  Email: David_Abt@brown.edu          *
% *                                      *
% ****************************************

clear Temp
if(~isempty(varargin) && strcmp(varargin{end},'synth'))
    Phases=cellstr(varargin{1}); tag='synth';
else Phases = {'Ps'; 'Sp'}; tag='data';
end

for iph=1:length(Phases)
    PW = [Phases{iph} '_Waveforms'];
    if Waveform_Data_Full.(PW).S2N_ratio~=-666
        Temp.(PW).Time              = single(Waveform_Data_Full.(PW).Time);
        Temp.(PW).dt                = Waveform_Data_Full.(PW).dt;
        if strcmpi(tag,'data')
            Temp.(PW).Z                 = single(Waveform_Data_Full.(PW).Z);
            Temp.(PW).N                 = single(Waveform_Data_Full.(PW).N);
            Temp.(PW).E                 = single(Waveform_Data_Full.(PW).E);
            Temp.(PW).R                 = single(Waveform_Data_Full.(PW).R);
            Temp.(PW).T                 = single(Waveform_Data_Full.(PW).T);
        end
        Temp.(PW).P                 = single(Waveform_Data_Full.(PW).P);
        Temp.(PW).SV                = single(Waveform_Data_Full.(PW).SV);
        if isfield(Waveform_Data_Full.(PW),'SH');
            Temp.(PW).SH                = single(Waveform_Data_Full.(PW).SH);
        end
        Temp.(PW).Window            = Waveform_Data_Full.(PW).Window;
        Temp.(PW).ZR_Xcorr_Coeff    = Waveform_Data_Full.(PW).ZR_Xcorr_Coeff;
    end
    Temp.(PW).S2N_ratio         = Waveform_Data_Full.(PW).S2N_ratio;
    Temp.(PW).wrong=Waveform_Data_Full.(PW).wrong;
end
Waveform_Data_Reduced = Temp;

clear Temp
Temp.REM = Ray_Path_Data_Full.REM;
Phases = {'P','S'};
for iph=1:length(Phases)
    Phase = Phases{iph};
    if ~isempty(Ray_Path_Data_Full.([Phase '_Path']))
        Temp.([Phase '_Path'])                  = Ray_Path_Data_Full.([Phase '_Path']);
        Temp.([Phase '_Path']).path             = rmfield(Temp.([Phase '_Path']).path,'p');
        Temp.([Phase '_Path']).path.time        = single(Temp.([Phase '_Path']).path.time);
        Temp.([Phase '_Path']).path.distance    = single(Temp.([Phase '_Path']).path.distance);
        Temp.([Phase '_Path']).path.depth       = single(Temp.([Phase '_Path']).path.depth);
        Temp.([Phase '_Path']).path.latitude    = single(Temp.([Phase '_Path']).path.latitude);
        Temp.([Phase '_Path']).path.longitude   = single(Temp.([Phase '_Path']).path.longitude);
    end
end
Ray_Path_Data_Reduced = Temp;

end

function [Vp_wmed,Vs_wmed] = get_weighted_median(filename)
%Written by NJM, Aug 2016
d=load(filename);

log10_SNR=log10(d(:,2));
Vp=d(:,3);
Vs=d(:,4);
w=log10_SNR./sum(log10_SNR);
Vp_wmed=weightedMedian(Vp,w);
Vs_wmed=weightedMedian(Vs,w);

end

function wMed = weightedMedian(D,W)
% Pulled function from MathWorks website. -NJM
% ----------------------------------------------------------------------
% Function for calculating the weighted median
% Sven Haase
%
% For n numbers x_1,...,x_n with positive weights w_1,...,w_n,
% (sum of all weights equal to one) the weighted median is defined as
% the element x_k, such that:
%           --                        --
%           )   w_i  <= 1/2   and     )   w_i <= 1/2
%           --                        --
%        x_i < x_k                 x_i > x_k
%
%
% Input:    D ... matrix of observed values
%           W ... matrix of weights, W = ( w_ij )
% Output:   wMed ... weighted median
% ----------------------------------------------------------------------


if nargin ~= 2
    error('weightedMedian:wrongNumberOfArguments', ...
        'Wrong number of arguments.');
end

if size(D) ~= size(W)
    error('weightedMedian:wrongMatrixDimension', ...
        'The dimensions of the input-matrices must match.');
end

% normalize the weights, such that: sum ( w_ij ) = 1
% (sum of all weights equal to one)

WSum = sum(W(:));
W = W / WSum;

% (line by line) transformation of the input-matrices to line-vectors
d = reshape(D',1,[]);
w = reshape(W',1,[]);

% sort the vectors
A = [d' w'];
ASort = sortrows(A,1);

dSort = ASort(:,1)';
wSort = ASort(:,2)';

sumVec = [];    % vector for cumulative sums of the weights
for i = 1:length(wSort)
    sumVec(i) = sum(wSort(1:i));
end

wMed = [];
j = 0;

while isempty(wMed)
    j = j + 1;
    if sumVec(j) >= 0.5
        wMed = dSort(j);    % value of the weighted median
    end
end

% final test to exclude errors in calculation
if ( sum(wSort(1:j-1)) > 0.5 ) & ( sum(wSort(j+1:length(wSort))) > 0.5 )
    error('weightedMedian:unknownError', ...
        'The weighted median could not be calculated.');
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
