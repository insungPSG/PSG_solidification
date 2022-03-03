function [D,P] = detectDroplets(frn,opticalflow,P,D)

% ==========================================================================
% PROCESSING 

% Read current frame 
D.img = imshowxray(P.obj,frn,'nofigure');

% Remove spike pixels
D.img(D.img>1) = 1;

% Filter image - NOTE: Add more filters here if necessary

% Apply non linear intensity transformation to enhance droplets contrast
D.imgp = imadjust(D.img,...
    P.imadjust_inRange,P.imadjust_outRange,P.imadjust_gamma); 

% ==========================================================================
% SEGMENTATION

% Use Optical flow to detect any moving object
D.flow = estimateFlow(opticalflow,D.imgp);
D.imgseg = D.flow.Magnitude > P.segFlowThresh;
% D.imgseg = D.imgp > P.segThresh;

% ==========================================================================
% ANALYSIS

%{ 
Probably not Needed when using Optical flow
% Morphological opertion
% D.imgseg = bwmorph(D.imgseg,'clean');
% D.imgseg = imopen(D.imgseg,P.strel);
% D.imgseg = imclose(D.imgseg,P.strel);
% D.imgseg = imfill(D.imgseg,'holes');

% Filter out noise and
% if P.removeNonPeakRegion
%     if ~bandPassFilterFlag
%         ibpass = (bpass(id,bpasslnoise,bpassSize,bpassTh));
%     else
%         ibpass = id;
%     end
%     pk     = pkfnd(ibpass,thblob,sizeblob);
%     ctr    = cntrd(ibpass,pk,15);
%     if ~isempty(ctr); ibbw = bwselect(ibbw,ctr(:,1),ctr(:,2),8); end
% end % NOT YET
% if P.propertiesFilter
%     D.imgseg = filterregionproperties(D.imgseg,...
%         {'Area', @gt, P.areawindow(1)},...
%         {'Area', @lt, P.areawindow(2)},...
%         {'Eccentricity',@lt,P.ecc}); % Remove small particles
% end % NOT YET
% if P.removeRecRegion
%     if ~isfield(D,'L')
%         D.L    = bwlabel(D.imgseg);
%     end
%     D.L = recurrentRegions(D.imgseg,D.L);
%     D.imgseg          = D.imgseg - D.L>0;
% end % NOT YET
%}

% Extract the needed information and store them in the data structure D
S = regionprops(D.imgseg,'Centroid','BoundingBox','EquivDiameter');

D.centroid  = cat(1,S.Centroid);
D.bbox      = cat(1,S.BoundingBox);
D.diameter  = cat(1,S.EquivDiameter);
D.meanInt   = zeros(size(D.diameter,1),1); % Not sure what this is for 

end