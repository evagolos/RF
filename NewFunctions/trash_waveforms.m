function trash_waveforms
clear all; clc
cd /Volumes/bigdisk1/hopper/Scattered_Waves
Project='Test';
TrashSAC(Project)

end


function TrashSAC(Project)
basedir=pwd;
cd(['Data/Projects/' Project]);
load Networks.mat
NET='Z9';
stas=fieldnames(Networks.(NET));

for is=1:length(stas)
    sta=stas{is};
    cd([NET '/' sta]);
    clear All_Waveform_Data; load('Waveform_Data.mat');
    eqs=fieldnames(All_Waveform_Data);
    
    for ie=1:length(eqs)
        if All_Waveform_Data.(eqs{ie}).Sp_Waveforms.S2N_ratio~=-666
            
            f=figure; evcode=strrep(strrep(eqs{ie},'EQ_',''),'_','.');
            subplot(3,1,1); plot(All_Waveform_Data.(eqs{ie}).Sp_Waveforms.N,'k-'); title([evcode ', ' sta ': BHN']);
            subplot(3,1,2); plot(All_Waveform_Data.(eqs{ie}).Sp_Waveforms.E,'k-'); title([sta ': BHE']);
            subplot(3,1,3); plot(All_Waveform_Data.(eqs{ie}).Sp_Waveforms.Z,'k-'); title([sta ': BHZ']);
            trash=input('Trash this waveform (y/[n])?  ','s');
             close(f);
        else trash='y';
        end
                if strcmp(trash,'y')
                    if ~isdir('./SAC_Files/trash'); mkdir('./SAC_Files/trash'); end
                    system(['mv ./SAC_Files/' evcode '*SAC ./SAC_Files/trash/']);
                    All_Waveform_Data=rmfield(All_Waveform_Data,eqs{ie});
                end
               

    end
    save Waveform_Data.mat All_Waveform_Data
    
    cd ../../    

end

cd(basedir);

end
