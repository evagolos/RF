load('Networks.mat');
ax=worldmap([-48,-40],[166,175]);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
for i=1:length(Alpine)
    plotm(Alpine(i).lat,Alpine(i).lon,'r','Linewidth',2.5);
end
NETs={'IU','XB','XC','XM','XU','Y3','YR','ZT','NZ','bodge4A'};
colors=[1,0,0;1,0,1;0,1,1;0.4,0.6,0.6;0,0,1;1,1,1;0.5,0.6,0.2;0.2,0.2,0.7;0.8,0.3,0.4;0.1,0.2,0.3];
for in=1:length(NETs)
    net = char(NETs{in});
    STAs = fieldnames(Networks.(net));
    disp([num2str(in)]);
    for is=1:length(STAs)
        sta=char(STAs{is});
        lat=Networks.(net).(sta).Latitude;
        lon=Networks.(net).(sta).Longitude;
        plotm(lat,lon,'color',[colors(in,:)],'Marker','^','Markerfacecolor',[colors(in,:)],'markeredgecolor','k','markersize',8);
    end
    plotm(-40.4-in*0.2,166.7,'color',[colors(in,:)],'Marker','^','Markerfacecolor',[colors(in,:)],'markeredgecolor','k','markersize',8);
    textm(-40.4-in*0.2,166.9,net((length(net)-1):length(net)));
end
