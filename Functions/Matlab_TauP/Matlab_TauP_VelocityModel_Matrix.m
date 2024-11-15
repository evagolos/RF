function [myVelocityModel] = Matlab_TauP_VelocityModel_Matrix(modelName,ND_Matrix,ND_Discon)

% function Matlab_TauP_VelocityModel
%
% Arguments
% ---------
% modelName             character array
% ND_Matrix             columns: Depth Vp Vs [Optional: Rho Qp Qs]
% ND_Discon.depth(:)    array of discontinuity depths
% ND_Discon.name(:)     cell array of discontinuity names
%
% Returns
% ------
% TauP VelocityModel object on success, nothing on failure
%
%  -- INHERENT ISSUES -- 
% The layers field in the VelocityModel class is protected and no public
% accessor method exists. Thus, we cannot simply create a VelocityModel
% object and populate the layers. Possible solutions:
% 1.    Read in some base model and use the replaceLayers() method to populate
%       it with our model.
% 2.    Make a ND type temporary file for our model and read it in.
% Though both are hacks and require file i/o, the second option is easiest.

% Import the TauP package
import edu.sc.seis.TauP.*

% Define our temporary directory
TMP_DIR = '/tmp';

% Define our temporary file
TMP_FILE_NAME = ['taup_tmp_' modelName '.nd'];

% We require Depth, Vp, and Vs [Rho, Qp, and Qs are optional]
validNumParameters = [3 4 6];
if ~isempty(find(validNumParameters == size(ND_Matrix,2),1))
    
    % Open the temporary file
    fTmp = fopen([TMP_DIR '/' TMP_FILE_NAME],'w');
    
    % Write the ND-style input file
    for layer=1:size(ND_Matrix)
        ndLayerString = num2str(ND_Matrix(layer,:),'%10.4f');
        fprintf(fTmp,'%s\n',ndLayerString);
        discon = find(ND_Discon.depth == ND_Matrix(layer,1),1);
        if ~isempty(discon) && ND_Matrix(layer,1) == ND_Matrix(layer+1,1)
            fprintf(fTmp,'%s\n',char(ND_Discon.name(discon)));
        end
    end
    
    % Close the temporary file
    fclose(fTmp);
    
    % Create a VelocityModel object
    myVelocityModel = VelocityModel();
    
    % Establish that the model is of type ND
    myVelocityModel.setFileType('nd');
    
    % Read in the file
    myVelocityModel.readVelocityFile([TMP_DIR '/' TMP_FILE_NAME]);
    
    % Set the model name
    myVelocityModel.setModelName(modelName);
    
    % Test the model
    if ~myVelocityModel.validate()
        clear myVelocityModel
    end
    
    % Remove the temporary file
    delete([TMP_DIR '/' TMP_FILE_NAME]);
end

end
