function [myVelocityModel] = Matlab_TauP_VelocityModel_File(modelName,modelFileType,modelFilePath)

% Import the TauP package
import edu.sc.seis.TauP.*

% Create a VelocityModel object
myVelocityModel = VelocityModel();

% Establish that the model is of type modelFileType
myVelocityModel.setFileType(modelFileType);

% Read in the file
myVelocityModel.readVelocityFile([modelFilePath '/' modelName '.' modelFileType]);

% Set the model name
myVelocityModel.setModelName(modelName);

% Test the model
if ~myVelocityModel.validate()
    clear myVelocityModel
end

end
