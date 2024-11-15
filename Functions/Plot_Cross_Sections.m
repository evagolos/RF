clear all
% Modified from Junlin's script. Some of his functionality is commented
% out, and I haven't tested everything. Displays 4 subplots.
Data = '/Users/evagolos/Research/ReceiverFunctions/RFcodes_Junlin'; 
fol = 'Data/Projects/SWUS_big/CCP';
cd(Data)
addpath([Data '/Functions'])
Ps= 0;
if Ps==1
    file='CCP_Ps_13_SWUS_Ps_4_100s_05.10.20_1_780km_deg_1.mat'; % File created within CCP stacking code.
else
    %file='CCP_Sp_13_SWUS_Sp_4_100s_05.10.20_1_2_450km_deg_1.mat';
    %file='CCP_Sp_13_SWUS_big_Sp_4_100s_18.05.21_450km_deg_1.mat';
    file='CCP_Sp_13_SWUS_big_Sp_2_100s_22.04.22_450km_deg_1.mat';
end
    
% Load CCP stack data and reference velocity model    
load([Data '/' fol '/' file]);
junk=strsplit(file,[phase '_']); junk=strsplit(junk{end},'_'); T=junk{1};
string = ['pics_' T 's_'];
VsModelName='MITPS_10kmZ'; % or 'DNA13', 'FichtnerSV13', etc.
load([Data '/Data/Velocity_Models/Mantle/Vs/' VsModelName '.mat']);
[~,~,Vs_Model.Vs] = gradient(Vs_Model.Vs,1,1,Vs_Model.Depth(2)-Vs_Model.Depth(1));

plane_RF=real(plane_RF);  %CHANGE
%std_RF=real(std_RF);
%% Interpolate data
min_dep = 0; 
if Ps==1
    max_dep = 780; 
else
    max_dep= 450;
end
DepNpts = (max_dep-min_dep)/5;
rplane_rf=reshape(plane_RF,[size(tmplat), length(cp_depths)]);
rplane_std=reshape(std_RF,[size(tmplat), length(cp_depths)]);
weighted_n_events = weighting_norm;
rweighted_n_events=reshape(weighted_n_events,[size(tmplat), length(cp_depths)]);
rplane_rf(rweighted_n_events<0.05) = NaN;
rplane_rf(rplane_std>0.2 & 3*abs(rplane_std)>abs(rplane_rf)) = NaN;
%rweighted_n_events=reshape(n_events,[size(tmplat), length(cp_depths)]);
%rstd_RF=reshape(std_RF,[size(tmplat),length(cp_depths)]);
deps_tmp = linspace(min_dep,max_dep,DepNpts); deps=double(cp_depths)';
[glo,gla,gdep]=ndgrid(tmplon(:,1),tmplat(1,:)',deps);
[plo,pla]=ndgrid(tmplon(:,1),tmplat(1,:)');
RF=griddedInterpolant(glo,gla,gdep,rplane_rf);
%STD=griddedInterpolant(glo,gla,gdep,rplane_std);
COV=griddedInterpolant(glo,gla,gdep,rweighted_n_events);



[vellat,vellon]=meshgrid(Vs_Model.Latitude,Vs_Model.Longitude);
[vlo,vla,vdep]=ndgrid(vellon(:,1),vellat(1,:)',Vs_Model.Depth);
vsmo=Vs_Model.Vs;
P=[2,1,3];
vlo=permute(vlo, P);
vla=permute(vla, P);
vdep=permute(vdep, P);
velu=interp3(vlo,vla,vdep,vsmo,glo,gla,gdep);
VEL=griddedInterpolant(glo,gla,gdep,velu);




%% set up cross section lines

%close all;
% Specify orientation of profile line. I haven't tested the other options.
sects='randomW';% Other choices: 'NW-SE', 'NW-SW', 'SE-NW', 'W-E'; %
numdeg= 0;

% Here's where you specify the start and end points of your profile line.
latlons=[40 -122 40 -100];
        

if numdeg~=0; for k=2:numdeg % In Junlin's version, numdeg/0.1
        latlons(k,:)=[latlons(1,1)+k*lasc*0.1 latlons(1,2)+k*losc*0.1 ...
            latlons(1,3)+k*lasc*0.1 latlons(1,4)+k*losc*0.1];
    end; else k=2;
end

kj=round(k/2); sc = 0.04;%0.15;%0.06;

% Geographic coordinates of profile line.
Npts=ceil(111.16*distance(latlons(kj,1),latlons(kj,2), ...
    latlons(kj,3), latlons(kj,4))/10);
latsAll=zeros(size(latlons,1),Npts); lonsAll=latsAll;
for k=1:size(latlons,1);
    lonsAll(k,:)= linspace(latlons(k,2),latlons(k,4),Npts);
    latsAll(k,:)=linspace(latlons(k,1),latlons(k,3),Npts);
end

PLOTNOW=1:size(latlons,1);

%%
for kj =1:length(PLOTNOW)
    clear x y z;
    clear x y z; clear discmdl n_eventsINTERP weighted_eventsINTERP Overlay;
    lons = lonsAll(kj,:); lats = latsAll(kj,:);
    [gdep1,glo1]=ndgrid(deps,lons); [~,gla1]=ndgrid(deps,lats);

     latlima=max(latlons(kj,1),latlons(kj,3));
     latlimi=min(latlons(kj,1),latlons(kj,3));
     lonlima=max(latlons(kj,2),latlons(kj,4));
     lonlimi=min(latlons(kj,2),latlons(kj,4));
     coplat=99;
     coplon=[];

     

    %interpolate plane_RF onto coordinates of cross section
    discmdl=RF(glo1,gla1,gdep1);
    %discmdlstd=STD(glo1,gla1,gdep1);
    discmdlvel=VEL(glo1,gla1,gdep1);
    weighted_eventsINTERP=COV(glo1,gla1,gdep1);
    weighted_eventsINTERP(weighted_eventsINTERP<0)= NaN; % Otherwise you'll get errors.
    
    
    dalj = 111.16*distance(latlons(kj,1),latlons(kj,2),lats,lons);
     if coplat~=99
         dalcop= 111.16*distance(latlons(kj,1),latlons(kj,2),coplat,coplon);
     end
    
    Cutoff=50;
    weighted_or_not = weighted_eventsINTERP;
%     
    fontsz=20;
    
    % plot cross sections and coverage
    %close all
    savefig='no';
    min_dep = 0;   
    max_dep= 350;
    maxlat=max(Vs_Model.Latitude)-0.1; %sc=0.25;
    close all; clear x y z;
    
    % The way this is set up now, if you specify multiple lines it creates
    % a new figure for each.
    f=figure(kj); set(f,'position',[1 41 1600 783], 'PaperPositionMode','auto');
    Title=axes('position',[0 0 1 1],'visible','off');
    switch (lower(phase)); case 'sp'; cmap=fliplr(colormap(jet)); case 'ps'; cmap=colormap(jet); end
    text(0.35, 0.95,['\fontsize{14} \bf ' num2str(latlons(kj,1)) 'N, ' num2str(latlons(kj,2)) 'E   to   ' ...
        num2str(latlons(kj,3)) 'N, ' num2str(latlons(kj,4)) 'W'])
    
    if max_dep==75; fpos=[60 350 1480 750]; f2pos=[60 120 1480 300];
    else f2pos=[985 555 920 307]; fpos=[55 265 920 820];
    end
    
    % Set up points along profile in plotting space
    switch sects
        case {'NW-SE','SW-NE','W-E'}; dalj2=lons(1):0.5:lons(end);
        case {'N-S'}; dalj2=lats(1):-0.5:lats(end);
        otherwise
            if lats(1)>lats(end)
                dalj2=lats(1):-0.5:lats(end);
            else
                dalj2=lats(1):0.5:lats(end);
            end
    end
    halfdeg=abs((dalj(end)-dalj(1))/diff(dalj2([1 end]))/2);
    lims=[min(dalj) max(dalj)];
    
    
    % You can specify text files that outline the boundaries of geologic
    % provinces
%     BNR= importdata('/Users/evagolos/Research/Misc/GeoBounds/BNR.txt');
%     Col= importdata('/Users/evagolos/Research/Misc/GeoBounds/Colorado.txt');
%     RGR= importdata('/Users/evagolos/Research/Misc/GeoBounds/RioGrande.txt');
%     SNB= importdata('/Users/evagolos/Research/Misc/GeoBounds/Sierra.txt');
    
    % Find intersections with the geologic provinces.
%     P_BNR= InterX(BNR',[lons; lats]); nBNR= size(P_BNR,2);
%     PointPlot= [];
%     for ii=1:nBNR
%         dlatlon= sqrt(((lons-P_BNR(1,ii)).^2 - (lats-P_BNR(2,ii)).^2));
%         [~,minD]= min(dlatlon);
%         PointPlot= [PointPlot; minD];
%     end
%     P_Col= InterX(Col',[lons; lats]); nCol= size(P_Col,2);
%         for ii=1:nCol
%         dlatlon= sqrt(((lons-P_Col(1,ii)).^2 - (lats-P_Col(2,ii)).^2));
%         [~,minD]= min(dlatlon);
%         PointPlot= [PointPlot; minD];
%     end
%     P_RGR= InterX(RGR',[lons; lats]); nRGR= size(P_RGR,2);
%         for ii=1:nRGR
%         dlatlon= sqrt(((lons-P_RGR(1,ii)).^2 - (lats-P_RGR(2,ii)).^2));
%         [~,minD]= min(dlatlon);
%         PointPlot= [PointPlot; minD];
%     end
%     P_SNB= InterX(SNB',[lons; lats]); nSNB= size(P_SNB,2);
%         for ii=1:nSNB
%         dlatlon= sqrt(((lons-P_SNB(1,ii)).^2 - (lats-P_SNB(2,ii)).^2));
%         [~,minD]= min(dlatlon);
%         PointPlot= [PointPlot; minD];
%     end
    
    figure(f)
    
    % This part plots the CCP stack volume along the profile
    a1=axes('position',[0.03 0.49 0.65 0.4]); hold on
    h = imagesc(dalj,deps,discmdl); caxis(sc*[-1.3 1.3]);
    set(h,'alphadata',~isnan(discmdl));
    ylabel('Depth (km)','fontsize',10);
    %plot(dalj,5*topo_prof,'color',[0 0 0],'linewidth',2);
    plot(a1,dalj(1:10:end),0.*dalj(1:10:end)-20,'ko','markersize',6,'markerfacecolor','g');
     if coplat ~= 99
         plot(a1,dalcop,0,'ko','markersize',6,'markerfacecolor','r');
     end
     % Uncomment if you are plotting geologic province boundaries
%      for jj=1:length(PointPlot)
%         plot(a1,[dalj(PointPlot(jj)) dalj(PointPlot(jj))],[-10 max(deps)],'--k','linewidth',2);
%      end
    set(a1,'ydir','reverse','ytick',0:50:max_dep,'yticklabel',0:50:max_dep); daspect([1 1 1]);
    set(a1,'xtick',dalj(1):halfdeg:dalj(end),'xticklabel',dalj2)
    set(a1,'ygrid','on','xgrid','on','layer','top','fontsize',10,'box','on');
    axis([lims -25 max_dep]);   ylabel('Depth (km)','fontsize',10);
    title(['a) CCP Stack at ' T '-100s Filter Band']);
    xlabel('X Distance (degrees)');
    colormap(cmap); c=colorbar('location', 'southoutside');
    c.Label.String= 'Stack Amplitude'; c.Label.FontSize=12;
    set(a1,'FontSize',12);
    
    % Plot gradients of reference model.
    a2=axes('position',[0.03 0.04 0.65 0.4]); hold on
    contourf(dalj,deps,-discmdlvel,100,'linestyle','none'); shading flat; caxis([-0.015 0.015]);
    ylabel('Depth (km)','fontsize',10);
    %plot(dalj,5*topo_prof,'color',[0 0 0],'linewidth',2);
    plot(a2,dalj(1:10:end),0.*dalj(1:10:end)-20,'ko','markersize',6,'markerfacecolor','g');
     if coplat ~= 99
         plot(a2,dalcop,0,'ko','markersize',6,'markerfacecolor','r');
     end
    %errorbar(a2,dalj(~isnan(disndep)),disndep(~isnan(disndep)),disnbdh(~isnan(disndep))/2,'ko','linewidth',1,'MarkerFaceColor','blue','MarkerSize',4);
    %errorbar(a2,dalj(~isnan(dispdep)),dispdep(~isnan(dispdep)),dispbdh(~isnan(dispdep))/2,'ko','linewidth',1,'MarkerFaceColor','red','MarkerSize',4);
    set(a2,'ydir','reverse','ytick',0:50:max_dep,'yticklabel',0:50:max_dep); daspect([1 1 1]);
    set(a2,'xtick',dalj(1):halfdeg:dalj(end),'xticklabel',dalj2)
    set(a2,'ygrid','on','xgrid','on','layer','top','fontsize',10,'box','on');
    axis([lims -25 max_dep]);   ylabel('Depth (km)','fontsize',10);
    title(['Depth Gradient of Reference Model']);
    xlabel('X Distance (degrees)');
    colormap(cmap); c=colorbar('location', 'southoutside');
    c.Label.String= 'Gradient Amplitude'; c.Label.FontSize=12;
    set(a2,'FontSize',12);
    

    % Map view of profile locations, colored according to CCP weighting
    a3=axes('position',[0.72 0.49 0.25 0.4]); hold on
    %lon=[-120 -100]; lat=[28 42];
    lon= [min([latlons(2) latlons(4)])-2, max([latlons(2) latlons(4)+2])]; 
    lat= [min([latlons(1) latlons(3)])-2, max([latlons(1) latlons(3)+2])]; 
    %lat=[min(lats)-2 max(lats)+2]; lon=[min(lons)-2 max(lons)+2];%[26 40];
    if abs(diff(lat))<10; lat=mean(lat)+[-5 5]; end
    if abs(diff(lon))<6; lon=mean(lon)+[-5 5]; end
    daspect([111.16, 111.16*distance(mean(lat),0,mean(lat),1), 1]);
    box on;  axis([lon, lat]);%[-90 -77, lat]);
    contourf(tmplon(:,1),tmplat(1,:),log(rweighted_n_events(:,:,200))','linestyle','none');
    set(a3,'xgrid','on','ygrid','on','layer','top')
    % Plot state lines and geologic province boundaries
%         S=shaperead('usastatelo','UseGeoCoords',true);
%         for i=3:51; plot(S(i).Lon,S(i).Lat,'k'); end
%         plot(Col(:,1),Col(:,2),'k','linewidth',2);
%         text(248-360,37,'Colo','FontWeight','bold','FontSize',12);
%         plot(RGR(:,1),RGR(:,2),'k','linewidth',2);
%         text(256-360,34,'RG','Fontweight','bold');
%         plot([255 256]-360,[33 33.5],'k');
%         plot(BNR(:,1),BNR(:,2),'k','linewidth',2);
%         text(241-360,40.1,'BNR','FontWeight','bold');
%         plot(SNB(:,1),SNB(:,2),'k','linewidth',2);
%         text(-120.6,38.1,'SNB','FontWeight','bold');
    plot(lons,lats,'k-','linewidth',2); hold on;
    plot(lons(1:10:end),lats(1:10:end),'ko','markerfacecolor','c','markersize',8);
    title('Map View, Weighting');
    coast=load('coast');
    hold on
    plot(coast.long,coast.lat)
    xlabel('Longitude (degrees)','FontSize',12); ylabel('Latitude (degrees)','FontSize',12);
    c=colorbar('location', 'southoutside'); 
    set(a3,'FontSize',12);
    c.Label.String= 'log(weighting) at 200 km'; c.Label.FontSize=12;
    
    % Plot weighting in vertical profile view.
    a4=axes('position',[0.72 0.04 0.25 0.4]); hold on
    imagesc(dalj,deps,log(weighted_eventsINTERP));% caxis([0 2.5]);
    ylabel('Depth (km)','fontsize',10);
    %plot(dalj,5*topo_prof,'color',[0 0 0],'linewidth',2);
    plot(a4,dalj(1:10:end),0.*dalj(1:10:end)-20,'ko','markersize',6,'markerfacecolor','g');
    if coplat ~= 99
         plot(a2,dalcop,0,'ko','markersize',6,'markerfacecolor','r');
     end
    set(a4,'ydir','reverse','ytick',0:50:max_dep,'yticklabel',0:50:max_dep); daspect([1 1 1]);
    set(a4,'xtick',dalj(1):halfdeg:dalj(end),'xticklabel',dalj2,'fontsize',11)
    set(a4,'ygrid','on','xgrid','on','layer','top','fontsize',10,'box','on');
    axis([lims -25 max_dep]);   ylabel('Depth (km)','fontsize',10);
    title(['Depth View, Weighting']);
    xlabel('X Distance (degrees)');
    c=colorbar('location', 'southoutside');
    c.Label.String='log(weighting)'; c.Label.FontSize=12;

    
    
%     if sum(strcmp(sects(1:end-1),{'random','cir'}))
%         figname=[fignamestr num2str(length(dir([fignamestr '*']))+1)];
%     else figname=[fignamestr num2str(kj)];
%     end
    
    if strcmp(savefig, 'yes')
        
        pause(1); figure(f); fext='png';
        if strcmp(fext,'ps'); pstr='epsc'; elseif strcmp(fext,'jpg'); pstr='jpeg'; else pstr=fext; end
        eval(['print -d' pstr ' ' figname '.' fext]);
        if strcmp(fext,'ps'); system(['ps2pdf -dEPSCrop ' figname '.' fext]);
            system(['rm ' figname '.' fext]);   end
    end
end

