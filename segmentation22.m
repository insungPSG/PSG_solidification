load('X:\4000_SharedDataAnalysis\Xray_Objs\petra21_009aEnzoTest.mat')

for frame = 2300:4400
    fprintf('Saving frame # %d ... \n',frame);
    a = imshowxray(petra21_009a, frame, 'nofigure'); 

    a(a>1) = 1;
    a(a<0) = 0;

    imadjust_inRange    = [0 0.5];
    imadjust_outRange   = [1 0]; % [1 0] negative image
    imadjust_gamma      = 5;

    [~,a1,~] = imreducehaze(a);

    a11 = medfilt2(a1);

    a2 = imadjust(a11,imadjust_inRange,imadjust_outRange,imadjust_gamma); 
    a22 = a2.*a2;

    a3 = a22 > 0.1;

    a4 = bwareaopen(a3, 6); % remove white clusters composed of fewer than 5 pixels

    BW = regionprops(a4, 'area', 'BoundingBox', 'PixelIdxList','MajorAxisLength', 'MinorAxisLength', 'Eccentricity');

    BW2 = zeros(size(a4,1), size(a4,2));
    for i = 1:size(BW,1)
        if BW(i).Area > 2000 || BW(i).MajorAxisLength / BW(i).MinorAxisLength > 2.0 || BW(i).Eccentricity > 0.85
            BW2(BW(i).PixelIdxList) = 0;
        else
            BW2(BW(i).PixelIdxList) = 1;
        end
    end
    A = [a, BW2];
    imwrite(mat2gray(A), sprintf('X:\\Insung_\\shrinkage\\petra009_segmentation\\%04d.tiff', frame), ...
       'compression', 'none');

end
%     AA = imfuse(a,BW2);
    % BW3 = regionprops(a4, 'area', 'BoundingBox', 'PixelIdxList','MajorAxisLength', 'MinorAxisLength');
%     imwrite(mat2gray(aa), sprintf('Y:\\Insung_\\shrinkage\\hottear009_1300_end_1_1200\\%04d.tiff', frame), ...
%            'compression', 'none');
% end
% figure; imshow(AA,[]);
% 
% figure; imshow(a,[]);
% figure; imshow(BW2,[]);


% a2 = medfilt2(a, [7,7]);
% BW2 = cat(1,BW.Area);
% BW3 = sort(BW2);

% figure; imshow(a,[]);
% figure; imshow(a4,[]);
% 
% 
% a3 = estimateFlow(opticalFlowLK,a22);
% 
% a_mag = mat2gray(a3.Magnitude);
% a_ori = mat2gray(a3.Orientation);
% a_vx = mat2gray(a3.Vx);
% a_vy = mat2gray(a3.Vy);
% 
% 
% 
% 
% 
% figure; imshowpair(a2, a22, 'montage');
% 
% 



% 
% 







% figure; imshowpair(a,a2,'montage');


% a3 = estimateFlow(opticalFlowLK,a2); %opticalFlowFarneback, opticalFlowLK, opticalFlowLKDoG
% 
% a_mag = mat2gray(a3.Magnitude);
% a_ori = mat2gray(a3.Orientation);
% a_vx = mat2gray(a3.Vx);
% a_vy = mat2gray(a3.Vy);

% a_ori2 = round(a_ori,3); % display up to w.xyz
% 
% a_ori2(a_ori2~=0.5) = 1;
% a_ori2(a_ori2~=1) = 0;
% 
% a_ori3 = a_ori2 .* a2 ; 
% a_ori3(a_ori3>0) = 1;
% a_ori4 = imfill(a_ori3,'holes');



% 
% a_ori3 = imfill(a_ori2,'holes');
% a_ori3(a_ori3~=1) = 0;
% a_ori3 = logical(a_ori3);
% 
% BW = regionprops(a_ori3, 'area', 'BoundingBox', 'PixelIdxList','MajorAxisLength', 'MinorAxisLength');
% 
% BW2 = zeros(size(a_ori3,1), size(a_ori3,2));
% 
% for i = 1:size(BW,1)
%     if BW(i).Area > 1500
%         BW2(BW(i).PixelIdxList) = 0;
%     else
%         BW2(BW(i).PixelIdxList) = 1;
%     end
% end
% 
% figure; imshow(BW2,[]);

% figure; imshowpair(a,BW2,'montage');

% A = imfuse(a,BW2); figure; imshow(A,[]);

% a_ori2(a_ori2~=0.5) = 1;
% a_ori3 = imfill(a_ori2,'holes');

% figure; imshowpair(a,a_ori3,'montage');

% binary = regionprops(a_ori3, 'BoundingBox');





% 
% s = regionprops( BW, 'BoundingBox');
% AR = s.BoundingBox(4) / s.BoundingBox(3);
% 
% figure; imshow(a3.Magnitude,[]);