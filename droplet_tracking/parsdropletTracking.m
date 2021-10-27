function P = parsdropletTracking

% ==========================================================================
% ANAYLISIS FLAGS 

% segmentationTest = 0; % NOTE: set to '0' for tracking
saveResults      = 1;
displayTracking  = 1; 
saveVideoFlag    = 1;

% ==========================================================================
% EXPERIMENT INFO 

% ### Exp name and folders
name          = 'petra21_009a';%'ma4063_075';%'mt9140_055';%
filename      = 'petra21_009aEnzoTest';%'ma4063_075';%'mt9140_055';%
folderXrayobj = 'Y:\4000_SharedDataAnalysis\Xray_Objs'; %'Y:\4000_SharedDataAnalysis\Xray_Objs';
folderData    = 'Y:\2021\I20210153\0009_a';

% ### Frame range
frameRange = 2200:3000; % frameRange must be firstFrame:step:lastFrame

% ### Data acquisition frame rate
frameRate     = 50; % in Hz (avg frame rate)
resolution    = 1.1; % micron/pixel

% ==========================================================================
% ANALYSIS PARAMETERS 

% ### Filter flags
% morphOp             = 0;                                    
% removeRecRegion     = 0;
% propertiesFilter    = 0;
% removeNonPeakRegion = 0;
% meanSubtract        = 1;
% spikesFilter        = 1;

% ### Image processing 
imadjust_inRange    = [0 0.5];
imadjust_outRange   = [1 0]; % [1 0] negative image
imadjust_gamma      = 5; % vals > 1 enhanche contrast of droplets

% ### Optical flow 
opflowAlgorithm          = 'HS';
opflowSmoothness         = 1;
opflowMaxIteration       = 10;
opflowVelocityDifference = 0;

% ### Segmentation 
% segThresh           = 0.5; % segmentation threshold
segFlowThresh       = 0.0055; % segmentation threshold

% strel               = [3 3];
% areawindow          = [5 1500]; % in pixel
% ecc                 = 0.75;
% thblob              = 0.1;
% sizeblob            = 5;

% ### Tracking parameters
TR.costOfNonAssignment  = 10;
TR.invisibleForTooLong  = 4;
TR.ageThreshold         = 8;
TR.MotionModel          = 'ConstantVelocity';
TR.InitialEstimateError = [200, 100];
TR.MotionNoise          = [150, 150];
TR.MeasurementNoise     = 25;
TR.minVisibleCount      = 2;

% ==========================================================================
% RESULTS SAVING INFO 

% ### Folder path
resultFolder = 'Y:\4000_SharedDataAnalysis\flowAnalysis'; %'G:\02_DataAnalysis\FluidFlow\Al-Pb_monotectic\Results_bulkLiquid\20190404';

% ### File name
expCondition = '';
annotation   = '';
resultFileName = [name,...
    '_frNo_',num2str(frameRange(1)),'_',num2str(frameRange(end)),...
    expCondition,annotation];
resultFilePath =  [resultFolder,filesep,resultFileName,'.mat'];

% ### Video parameters
videoFolder    = resultFolder;%'Z:\DataAnalysis\FluidFlow\Al-Pb_monotectic\Videos\Tracking - bulk liquid';
videoFilePath  = [videoFolder,filesep,resultFileName,'.avi'];
videoFrameRate = 6;

% ======================================================================= %
%                       END OF EDITABLE PARAMETERS                        %
% ======================================================================= %

% Load the Xray object of the sequence to analyse
if ~exist('obj','var') || ~strcmpi(obj.Experiment.name,name)
    eval('load([folderXrayobj,filesep,filename]);');
    eval(['obj = ',name,';'])
    eval(['clear ',name,'']) % change PETRA21_xxxa to obj
end
obj.changefolder(folderData); % Have to change the data path from Z:\ to Y:\

% Pack all parameters in a structure
P = v2struct;