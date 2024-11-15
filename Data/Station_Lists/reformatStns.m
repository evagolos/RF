% Reformat station list

fid1= fopen('rawstas_WholeUS_orig.txt','r');

fid2= fopen('rawstas_WholeUS.txt','w');

%while ~feof(fid1)
AllLines= textscan(fid1,'%s %s %s %s %s %s %s %s',5736,'headerlines',3,'Delimiter','|');


for i=1:length(AllLines{1})
net= AllLines{1}{i};
stn= AllLines{2}{i};
lat= AllLines{3}{i};
lon= AllLines{4}{i};
el= AllLines{5}{i};
tS= AllLines{7}{i};
tS1= tS(1:10); tS2= tS(12:end);
tE= AllLines{8}{i};
tE1= tE(1:10); tE2= tE(12:end);
fprintf(fid2,'%s \t %s \t BHE \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n',net,stn,tS1,tS2,tE1,tE2,lat,lon,el);
fprintf(fid2,'%s \t %s \t BHN \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n',net,stn,tS1,tS2,tE1,tE2,lat,lon,el);
fprintf(fid2,'%s \t %s \t BHZ \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n',net,stn,tS1,tS2,tE1,tE2,lat,lon,el);






end





fclose(fid1);
fclose(fid2);