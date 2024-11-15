clear all
% So these are the data coverage matrices that I saved after the fact with
% Ved's script.  They're the same size as plane_RF.
Data = '/Users/evagolos/Research/ReceiverFunctions/RFcodes_Junlin'; 
fol = 'Data/Projects/SWUS_big/CCP';
cd(Data)
addpath('/Users/evagolos/Research/Misc/Colormaps/');
%addpath([Data '/Functions'])
plotGrads= 1;
Ps= 0;
if Ps==1
    file='CCP_Ps_13_SWUS_Ps_4_100s_05.10.20_1_780km_deg_1.mat';
    phase= 'Ps';
else
    %file='CCP_Sp_13_SWUS_Sp_4_100s_05.10.20_1_2_450km_deg_1.mat';
    %file='CCP_Sp_13_SWUS_big_Sp_4_100s_18.05.21_450km_deg_1.mat';
    file='CCP_Sp_13_SWUS_big_Sp_2_100s_22.04.22_450km_deg_2.mat';
    phase= 'Sp';
    %%SWUS_big
end

if plotGrads
    %load([Data '/' fol '/moho_peak_2s_0822_' file]);
    load([Data '/' fol '/moho_peak_2s_0822_use']);
    load([Data '/' fol '/nvg_peak_2s_0822_cutoff009' file]);
end

%if ~exist('cp_depths','var')
load([Data '/' fol '/' file]);
%end
junk=strsplit(file,[phase '_']); junk=strsplit(junk{end},'_'); T=junk{1};
string = ['pics_' T 's_'];
VsModelName='MITPS_10kmZ';
load([Data '/Data/Velocity_Models/Mantle/Vs/' VsModelName '.mat']);

[~,~,Vs_Model.Vs] = gradient(Vs_Model.Vs,1,1,Vs_Model.Depth(2)-Vs_Model.Depth(1));
% 
% 
plane_RF=real(plane_RF);  %CHANGE
%% Interpolate data
min_dep = 0; 
if Ps==1
    max_dep = 780; 
else
    max_dep= 450;
end
DepNpts = (max_dep-min_dep)/5;
rplane_rf=reshape(plane_RF,[size(tmplat), length(cp_depths)]);
rplane_std=reshape(std_RF,[size(tmplat), length(cp_depths)]);weighted_n_events = weighting_norm;
rweighted_n_events=reshape(weighted_n_events,[size(tmplat), length(cp_depths)]);
rplane_rf(rweighted_n_events<0.05) = NaN;
rplane_rf(rplane_std>0.2 & 3*abs(rplane_std)>abs(rplane_rf)) = NaN;
deps_tmp = linspace(min_dep,max_dep,DepNpts); deps=double(cp_depths)';
[glo,gla,gdep]=ndgrid(tmplon(:,1),tmplat(1,:)',deps);
[plo,pla]=ndgrid(tmplon(:,1),tmplat(1,:)');
RF=griddedInterpolant(glo,gla,gdep,rplane_rf);
%STD=griddedInterpolant(glo,gla,gdep,rplane_std);
COV=griddedInterpolant(glo,gla,gdep,rweighted_n_events);


% Put the Moho picks onto the same grid as tmplat and tmplon for
% interpolation later.
if plotGrads
nmolat= length(unique(tmplat)); nmolon= length(unique(tmplon));
mohogrid= zeros(size(tmplat)); labgrid1= zeros(size(tmplat));
labgrid2= zeros(size(tmplat));
for i=1:nmolon
    mohogrid(i,:)= pvgdeps(i:nmolon:end);
    labgrid1(i,:)= nvgdeps(i:nmolon:end,1);
    labgrid2(i,:)= nvgdeps(i:nmolon:end,2);
end
end


%% set up cross section lines
xslabels= {'A - A''', 'B - B''', 'C - C''', 'D - D''', 'E - E'''}; 
latlons=[42 -116 42 -104; 39 -116 39 -104; 37 -116 37 -104; 34 -116 34 -104]; % (E-W)
%latlons=[41 -118 41 -104; 39 -118 39 -104; 36 -118 36 -104; 34 -118 34 -104]; % (E-W)
%latlons=[34 -116 34 -104; 42 -105 31 -105]
%latlons=[41 -116 41 -104; 38 -116 38 -104; 36 -116 36 -104; 33 -116 33 -104]; % (E-W)

NS= 0; % For proper x-axis labels on cross sections. Set to 1 for N to S section.
k= size(latlons,1);
kj=round(k/2); sc = 0.08;

Npts=ceil(111.16*distance(latlons(kj,1),latlons(kj,2), ...
    latlons(kj,3), latlons(kj,4))/10);
latsAll=zeros(size(latlons,1),Npts); lonsAll=latsAll;
for k=1:size(latlons,1);
    lonsAll(k,:)= linspace(latlons(k,2),latlons(k,4),Npts);
    latsAll(k,:)=linspace(latlons(k,1),latlons(k,3),Npts);
end

PLOTNOW=1:size(latlons,1);

%% Make the figure


fig1= figure();  
ppos= get(fig1,'Position');
set(fig1,'Position',[ppos(1)-500 ppos(2) ppos(3)*1.5 ppos(4)*2.4]);
    
BNR= importdata('/Users/evagolos/Research/Misc/GeoBounds/BNR.txt');
Col= importdata('/Users/evagolos/Research/Misc/GeoBounds/Colorado.txt');
RGR= importdata('/Users/evagolos/Research/Misc/GeoBounds/RioGrande.txt');
SNB= importdata('/Users/evagolos/Research/Misc/GeoBounds/Sierra.txt');
load('roma_sat.mat');


for kj =1:length(PLOTNOW)
    if kj==1
        NS= 0;
    else
        NS= 1;
    end
%for kj=18
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

     

    %interpolate plane_RF onto grid of cross sections
     discmdl=RF(glo1,gla1,gdep1);
    weighted_eventsINTERP=COV(glo1,gla1,gdep1);
    weighted_eventsINTERP(weighted_eventsINTERP<0)= NaN; % Otherwise get errors.
    
    
    dalj = 111.16*distance(latlons(kj,1),latlons(kj,2),lats,lons);
     if coplat~=99
         dalcop= 111.16*distance(latlons(kj,1),latlons(kj,2),coplat,coplon);
     end
    
    Cutoff=50;
    weighted_or_not = weighted_eventsINTERP;
    
% Interpolate Moho at that depth
    if plotGrads
        moholine= interp2(tmplat,tmplon,mohogrid,lats,lons);
        labline1= interp2(tmplat,tmplon,labgrid1,lats,lons);
        labline2= interp2(tmplat,tmplon,labgrid2,lats,lons);
    end
    fontsz=20;
    
    % plot cross sections and coverage
    %close all
    savefig='no';
    min_dep = 0;   
    max_dep= 230;
    maxlat=max(Vs_Model.Latitude)-0.1; %sc=0.25;
    clear x y z %close all; clear x y z;
    
    %f=figure(kj);set(f,'position',[1 41 1600 783], 'PaperPositionMode','auto');
    Title=axes('position',[0 0 1 1],'visible','off');
    switch (lower(phase)); case 'sp'; cmap=fliplr(colormap(jet)); case 'ps'; cmap=colormap(jet); end
    
    
    if max_dep==75; fpos=[60 350 1480 750]; f2pos=[60 120 1480 300];
    else f2pos=[985 555 920 307]; fpos=[55 265 920 820];
    end
    
    if NS==1
        if lats(1)>lats(end)
            dalj2=lats(1):-0.5:lats(end);
        else
            dalj2=lats(1):0.5:lats(end);
        end
    else
        if lons(1)>lons(end)
            dalj2=lons(1):-0.5:lons(end);
        else
            dalj2=lons(1):0.5:lons(end);
        end
    end

    halfdeg=abs((dalj(end)-dalj(1))/diff(dalj2([1 end]))/2);
    
    lims=[min(dalj) max(dalj)];
    
    
% Find intersection points with surface geology to plot.
    P_BNR= InterX(BNR',[lons; lats]); nBNR= size(P_BNR,2);
    PointPlot= [];
    for ii=1:nBNR
        dlatlon= sqrt(((lons-P_BNR(1,ii)).^2 - (lats-P_BNR(2,ii)).^2));
        [~,minD]= min(dlatlon);
        PointPlot= [PointPlot; minD];
    end
    P_Col= InterX(Col',[lons; lats]); nCol= size(P_Col,2);
        for ii=1:nCol
        dlatlon= sqrt(((lons-P_Col(1,ii)).^2 - (lats-P_Col(2,ii)).^2));
        [~,minD]= min(dlatlon);
        PointPlot= [PointPlot; minD];
    end
    P_RGR= InterX(RGR',[lons; lats]); nRGR= size(P_RGR,2);
        for ii=1:nRGR
        dlatlon= sqrt(((lons-P_RGR(1,ii)).^2 - (lats-P_RGR(2,ii)).^2));
        [~,minD]= min(dlatlon);
        PointPlot= [PointPlot; minD];
    end
    P_SNB= InterX(SNB',[lons; lats]); nSNB= size(P_SNB,2);
        for ii=1:nSNB
        dlatlon= sqrt(((lons-P_SNB(1,ii)).^2 - (lats-P_SNB(2,ii)).^2));
        [~,minD]= min(dlatlon);
        PointPlot= [PointPlot; minD];
    end
    hold on;
    
    if kj==1
        text(0.04,0.98,'c)','FontSize',16);
    end
    
    % Cross-section label
    text(0.1, 0.99-(kj-1)*0.2125,['\fontsize{14} \bf ' xslabels{kj} ': ' num2str(latlons(kj,1)) 'N, ' num2str(latlons(kj,2)) 'E   to   ' ...
        num2str(latlons(kj,3)) 'N, ' num2str(latlons(kj,4)) 'W'])
    
    % Begin plotting section
    a1=axes('position',[0.09 0.79-(kj-1)*0.21 0.8 0.2]); hold on
    h = imagesc(dalj,deps,discmdl); caxis(sc*[-1.3 1.3]);
    set(h,'alphadata',~isnan(discmdl));
    ylabel('Depth (km)','fontsize',10);
    %plot(dalj,5*topo_prof,'color',[0 0 0],'linewidth',2);
    plot(a1,dalj(1:10:end),0.*dalj(1:10:end)-20,'ko','markersize',6,'markerfacecolor','g');
     if coplat ~= 99
         plot(a1,dalcop,0,'ko','markersize',6,'markerfacecolor','r');
     end
     for jj=1:length(PointPlot)
        plot(a1,[dalj(PointPlot(jj)) dalj(PointPlot(jj))],[-10 max(deps)],'--k','linewidth',2);
     end
     % Add Moho and LAB (or other NVG) markers
     if plotGrads
        plot(dalj(moholine>0),moholine(moholine>0),'ko');
        plot(dalj(labline1>0),labline1(labline1>0),'kx');
        plot(dalj(labline2>0),labline2(labline2>0),'kx');
     end
    set(a1,'ydir','reverse','ytick',0:50:max_dep,'yticklabel',0:50:max_dep); daspect([1 1 1]);
    %set(a1,'xtick',dalj(1):halfdeg:dalj(end),'xticklabel',dalj2)
    set(a1,'ygrid','on','xgrid','on','layer','top','fontsize',10,'box','on');
    axis([lims -25 max_dep]);   ylabel('Depth (km)','fontsize',10);
    %title(['a) CCP Stack at ' T '-100s Filter Band']);
    if kj==4
        if NS==1
            set(a1,'xtick',dalj(1):halfdeg:dalj(end),'xticklabel',dalj2);
        else
            set(a1,'xtick',dalj(1):4*halfdeg:dalj(end),'xticklabel',dalj2(1):2:dalj2(end));
        end
    else
        set(a1,'xtick',dalj(1):4*halfdeg:dalj(end),'xticklabel',{})
        %set(a1,'xtick',lats(1:halfdeg:end),'xticklabel',{})
    end
    set(a1,'FontSize',14);
    
end
    if NS==1;
        xlabel('Latitude (degrees)');
    else
        xlabel('Longitude (degrees)');
    end
    c=colorbar('location', 'southoutside');
    c.Position= c.Position - [0 0.055 0 0];
    c.Label.String= 'Stack Amplitude'; c.Label.FontSize=14;
    colormap(roma_sat);
    
%% Plot map
% Convert x from CCP plot (dalj) to degrees
dal2x= (latlons(1,4)-latlons(1,2));
xmarkdal= dalj(1:10:end);
latfromdal= linspace(lats(end),lats(1),length(dalj));
xmarklon= lons(1);

lon= [-124 -100]; lat= [28 45];
fig2= figure; hold on;
daspect([111.16, 111.16*distance(mean(lat),0,mean(lat),1), 1]);
box on;%[-90 -77, lat]);
contourf(tmplon(:,1),tmplat(1,:),log(rweighted_n_events(:,:,150))','linestyle','none');
grid on; %set(a3,'xgrid','on','ygrid','on','layer','top')
    S=shaperead('usastatelo','UseGeoCoords',true);
    for i=3:51; plot(S(i).Lon,S(i).Lat,'k'); end
    plot(Col(:,1),Col(:,2),'k','linewidth',2);
    text(248-360,36,'CP','FontWeight','bold','FontSize',12);

    plot(RGR(:,1),RGR(:,2),'k','linewidth',2);
    text(252.8-360,32.5,'RG','Fontweight','bold');
    plot([253.5 254.5]-360,[32.8 33.3],'k','linewidth',1.5);
    plot(BNR(:,1),BNR(:,2),'k','linewidth',2);
    text(241-360,40.1,'BNR','FontWeight','bold');
    plot(SNB(:,1),SNB(:,2),'k','linewidth',2);
    text(-119.6,37.1,'SN','FontWeight','bold');

% Plot profile lines
% for ii=1:size(latlons,1)
%     lons = lonsAll(ii,:); lats = latsAll(ii,:);
%     plot(lons,lats,'k-','linewidth',2); hold on;
%     plot(lons(1:10:end),lats(1:10:end),'ko','markerfacecolor','c','markersize',8);
%     text(lons(1)-1,lats(1),xslabels{ii}(1)); 
%     if NS==1
%         text(lons(end)-1,lats(end),xslabels{ii}(5:6));
%     else
%         text(lons(end)+0.5,lats(1),xslabels{ii}(5:6));
%     end
% end

% Plot stations used in single-station RF stack plots.
% scatter(-106.4572,34.9459,90,'r','^','filled','MarkerEdgeColor','k');
%     text(-105.2,34.95,'ANMO','FontSize',12);
% scatter(-110.7847,32.3098,90,'b','^','filled','MarkerEdgeColor','k');
%     text(-110.7,33,'TUC','FontSize',12);
% scatter(-112.2279,32.7006,90,'c','^','filled','MarkerEdgeColor','k');
%     text(-114,33.4,'115A','FontSize',12);
% scatter(-114.9151,41.4157,90,'k','^','filled','MarkerEdgeColor','k');
%     text(-115.8,40.9,'M12A','FontSize',12);
% scatter(-108.4925,37.1996,90,'m','^','filled','MarkerEdgeColor','k');
%     text(-110.9,37.8,'MVCO','FontSize',12);
% scatter(-110.5238,39.1108,90,'w','^','filled','MarkerEdgeColor','k');
%     text(-110.2,39.7,'SRU','FontSize',12);
  

%text(-122,44.3,'b)','FontSize',16);
%title('Map View, Weighting');
coast=load('coast');
plot(coast.long,coast.lat)
axis([lon, lat]);
xlabel('Longitude (degrees)','FontSize',12); ylabel('Latitude (degrees)','FontSize',12);
c=colorbar('location', 'southoutside'); 
set(gca,'FontSize',12);
c.Label.String= 'log(weighting) at 100 km'; c.Label.FontSize=12;

if strcmp(savefig, 'yes')

    pause(1); figure(f); fext='png';
    if strcmp(fext,'ps'); pstr='epsc'; elseif strcmp(fext,'jpg'); pstr='jpeg'; else pstr=fext; end
    eval(['print -d' pstr ' ' figname '.' fext]);
    if strcmp(fext,'ps'); system(['ps2pdf -dEPSCrop ' figname '.' fext]);
        system(['rm ' figname '.' fext]);   end
end
%end

