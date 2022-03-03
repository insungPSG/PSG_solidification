T.id = []
T.bbox = []
T.kalmanFilter = []
T.age = []
T.totalVisibleCount = []
T.consecutiveInvisibleCount = []
T.active = []
T.intensity = []
T.resolution = []
T.frameNo = []
T.time = []
T.centroid = []
T.eqdiam = []
T.vel = []
T.acc = []
%% 
pars = {...
    'FNumPyramidLevels'       ,  1  ,...
    'FPyramidScale'           ,  0.99,...
    'FNumIterations'          ,  5  ,...
    'FNeighborhoodSize'       ,  5,...
    'FFilterSize'             ,  5 ,...
    'FDecimationFactor'       ,  7  ,...
    'FScaleFactor'            ,  10};


% Initialise optic flow object
% opticFlow = opticalFlowFarneback(...
%     'NumPyramidLevels',pars{2}, 'PyramidScale' ,pars{4},...
%     'NumIterations' ,pars{6}, 'NeighborhoodSize' ,pars{8},...
%     'FilterSize',pars{10});
opticFlow = opticalFlowHS;
opticFlow.Smoothness = 1;
opticFlow.MaxIteration = 10;
opticFlow.VelocityDifference = 0;

m =1;
for frn = P.frameRange
    % Read current frame
    D.img = imshowxray(P.obj,frn,'nofigure');
    
    % Remove spike pixels
    D.img(D.img>1) = 1;
    
    % Filter image - NOTE: Add more filters here if necessary
    
    % Apply non linear intensity transformation to enhance droplets contrast
    D.imgp = imadjust(D.img,...
        P.imadjust_inRange,P.imadjust_outRange,P.imadjust_gamma);
    flowtemp = estimateFlow(opticFlow,D.imgp);
    flow{m} = flowtemp;
%     imshow(imoverlay(D.img,flow{m}.Magnitude>0.005));
    imshow(D.img);hold on;
    plot(flowtemp,'DecimationFactor',[1 1],'ScaleFactor',150);
    pause(0.1);
    frn
    m = m+1;
end