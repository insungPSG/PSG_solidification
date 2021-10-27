function [frame] = readFrame(frn,P)

 i = imshowxray(P.obj,frn,'noimage');

% Calculate image difference
[id,inew] = imdiff(P.obj,frn,P.step);
id(id<0) = 0; % Remove negative pixels
id = normalise(id);
imshow(id,quantile(id(:),[0.1,0.9]))
% Front filter, mean subtraction and bandPass filter
% id = (imbilatfilt(id,'spatialsigma',spatialsigma));
% if meanSubtract
%     id = (id - mean2(id) * meanAlpha);
%     id(id<0) = 0; % Remove negative pixels
%     id = normalise(id);
% end % filter id
% if bandPassFilterFlag % filter id image
%     id = (bpass(id,bpasslnoise,bpassSize,bpassTh));
% end

% Segment the image and apply filters
% ibbw = imbinarize(id,'global');
% if removeNonPeakRegion
%     if ~bandPassFilterFlag
%         ibpass = (bpass(id,bpasslnoise,bpassSize,bpassTh));
%     else
%         ibpass = id;
%     end
%     pk     = pkfnd(ibpass,thblob,sizeblob);
%     ctr    = cntrd(ibpass,pk,15);
%     if ~isempty(ctr); ibbw = bwselect(ibbw,ctr(:,1),ctr(:,2),8); end
% end % filter ibbw
% if frn==1; L    = bwlabel(ibbw);end
% if removeRecRegion
%     L = recurrentRegions(ibbw,L);
%     ibbw          = ibbw - L>0;
% end % filter ibbw
% if propertiesFilter
%     ibbw = filterregionproperties(ibbw, {'Area', @gt, areawindow(1)},...
%         {'Area', @lt, areawindow(2)},{'Eccentricity',@lt,ecc}); % Remove small particles
%     ibbw = imfill(ibbw,'holes');
% end % filter ibbw
% if frontFilter && frn >= firstFrameFront
%     %                 ifilt = normalise(entropyfilt(id,true(entropyFiltSize)));
%     %                 ibbw = imbinarize(id.*(ifilt));
%     %                 imask = imfill(ibbw,'holes');
%     %                 imask = imfillborder(imask);
%     %                 % ibw = bwselect(ibw,size(ibw,2),size(ibw,1));
%     %                 imask = bwareafilt(imask,1);
%     %                 ibbw = ibbw.* ~imask;
%     xx = round(size(id,1) - frontVel* (frn-firstFrameFront)*deltaT(m));
%     if ~isempty(xx);xx = xx(1);end
%     if xx < 2; xx = 1;end
%     if xx > size(id,1); xx = size(id,1); end
%     ibbw(xx:end,:) = 0;
% end

% Create display image
% frame = cat(3,id,id,id);

% imshow(im2uint8(frame));