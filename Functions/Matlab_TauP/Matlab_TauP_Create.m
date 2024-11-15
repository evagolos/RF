function [myTauModel] = Matlab_TauP_Create(Mode,varargin)

% function Matlab_TauP_Create
%
% Arguments
% ---------
% Matlab_TauP_Create('FromObject','VelocityModel',VelocityModelObject)
% OR:
% Matlab_TauP_Create('FromFile', ...
%                    'modelFilePath',pathToModelFile, ...
%                    'modelName',nameOfModel, ...
%                    'modelFileType',typeOfModelFile)
%       where pathToModelFile, nameOfModel, and typeOfModelFile are
%       character arrays. typeOfModelFile will be either 'nd' or 'tvel'
%
% Returns
% -------
% TauP TauModel object on success, nothing on failure

% Import the TauP package
import edu.sc.seis.TauP.*

% Parse the supplied arguments
if strcmp(Mode,'FromObject')
    
    % Walk varargin for supplied VelocityModel
    nstdargs = 1;
    for argument=1:2:nargin-nstdargs
        switch varargin{argument}
            case {'VelocityModel'}
                myVelocityModel = varargin{argument+1};
            otherwise
                fprintf('\nWarning: unsupported argument %s ignored\n\n',varargin{argument});
        end
    end
    
    % Validate our arguments
    if ~ exist('myVelocityModel','var')
        fprintf('Error: no VelocityModel object supplied!\n');
        clear myVelocityModel
    else
        if ~ strcmp(class(myVelocityModel),'edu.sc.seis.TauP.VelocityModel')
            fprintf('Error: the supplied object is not a VelocityModel!\n');
            clear myVelocityModel
        else
            if ~ myVelocityModel.validate()
                fprintf('Error: the supplied VelocityModel is *not* valid!\n');
                clear myVelocityModel
            end
        end
    end
    
elseif strcmp(Mode,'FromFile')

    % Walk varargin for the location of the saved model
    nstdargs = 1;
    for argument=1:2:nargin-nstdargs
        switch varargin{argument}
            case {'modelFilePath'}
                 modelFilePath = varargin{argument+1};
            case {'modelName'}
                 modelName = varargin{argument+1};
            case {'modelFileType'}
                 modelFileType = varargin{argument+1};
            otherwise
                fprintf('\nWarning: unsupported argument %s ignored\n\n',varargin{argument});
        end
    end
    
    % Validate our arguments
    if ~ exist('modelFilePath','var')
        fprintf('Error: no path to velocity model file supplied!\n');
    else
        if ~ exist('modelName','var')
            fprintf('Error: no velocity model name supplied!\n');
        else
            if ~ exist('modelFileType','var')
                fprintf('Error: no velocity model type supplied!\n');
            else
                if ~ exist([modelFilePath '/' modelName '.' modelFileType],'file')
                    fprintf('Error: %s does not exist!\n',[modelFilePath '/' modelName '.' modelFileType]);
                else
                    % Create a VelocityModel object
                    myVelocityModel = VelocityModel();
                    
                    % Set the model file type
                    myVelocityModel.setFileType(modelFileType);
                    
                    % Read in the model
                    myVelocityModel.readVelocityFile([modelFilePath '/' modelName '.' modelFileType]);
                    if ~ myVelocityModel.validate()
                        fprintf('Error: the specified VelocityModel is *not* valid!\n');
                        clear myVelocityModel
                    end
                end
            end
        end
    end
    
else
    
    fprintf('\nError: unsupported mode %s\n\n',Mode);

end

myTauModel = [];

if exist('myVelocityModel','var')
    
    % Create a SphericalSModel object
    mySlownessModel = SphericalSModel();

    % Set slowness model parameters (lifted directly from TauP_Create
    % defaults)
    mySlownessModel.setMinDeltaP(0.1); % s / rad
    mySlownessModel.setMaxDeltaP(8.0); % s / rad
    mySlownessModel.setMaxDepthInterval(115.0); % km
    mySlownessModel.setMaxRangeInterval(1.75); % deg
    mySlownessModel.setMaxInterpError(0.05); % s
    mySlownessModel.setAllowInnerCoreS(1); % bool

    % Calculate the slowness model
    mySlownessModel.createSample(myVelocityModel);

    % Create a TauModel object
    myTauModel = TauModel();

    % Calculate the tau model from the give slowness model
    myTauModel.calcTauIncFrom(mySlownessModel);
end

end