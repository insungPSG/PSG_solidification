function opticalflow = initializeOpticalFlow(P)

if strcmpi(P.opflowAlgorithm,'HS')
    opticalflow = opticalFlowHS;
    opticalflow.Smoothness         = P.opflowSmoothness;
    opticalflow.MaxIteration       = P.opflowMaxIteration;
    opticalflow.VelocityDifference = P.opflowVelocityDifference;
end