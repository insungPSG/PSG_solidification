clear all
tic
i = 2300; % first frame here...

a = imread(sprintf('X:\\Insung_\\shrinkage\\petra009_segmentation2\\%04d.tiff', i)); 
a1 = regionprops(logical(a), 'area', 'BoundingBox', 'centroid'); % 'PixelIdxList'
initsize_id = [1:size(a1,1)]';
C = num2cell(initsize_id); 
[a1.id] = C{:};
a2  = cat(1,a1.Centroid);
activelist = ones(size(a1,1),1);
CC = num2cell(activelist);
[a1.active] = CC{:};
framelist = i*ones(size(a1,1),1);
CC = num2cell(framelist);
[a1.frame] = CC{:};
distancelist = zeros(size(a1,1),1);
CC = num2cell(distancelist);
[a1.distance] = CC{:};
for ii = 2301:4400
    fprintf('%d\n',ii);
    
    b = imread(sprintf('X:\\Insung_\\shrinkage\\petra009_segmentation2\\%04d.tiff', ii)); 
    b1 = regionprops(logical(b), 'area', 'BoundingBox', 'centroid'); % 'PixelIdxList'
    b2  = cat(1,b1.Centroid);

    closest = [];
    
    for i = 1:size(a1,1)
        dist_list = [];
        if a1(i).active > 0
            for j = 1:size(b2,1)
                dist = sqrt((a1(i).Centroid(end,1)-b2(j,1))^2 + (a1(i).Centroid(end,2)-b2(j,2))^2); %sqrt((a2(i,1)-b2(j,1))^2 + (a2(i,2)-b2(j,2))^2);
                dist_list = [dist_list; [dist,i,j]]; % here i is equivalent of id
            end
        
            A = dist_list(:,1);
            AA = find(A==min(A));
    %         closest = [closest; [a2(i,1), a2(i,2), b2(AA,1), b2(AA,2),AA,min(A)]];

            if a1(i).active > 0 && A(AA) < 5

                a1(i).Centroid(end+1,:) = b1(AA).Centroid; %[closest(t,3), closest(t,4)];
                a1(i).frame(end+1) = ii;
                a1(i).Area(end+1) = b1(AA).Area;
                a1(i).distance(end+1) = dist_list(AA,1);
                a1(i).BoundingBox(end+1,:) = b1(AA).BoundingBox; % x,y coordinates of the top-left corner of the box are first two (on the image). 
%                 a1(i).PixelIdxList(end+1,:) = aa1(AA).PixelIdxList;
                a1(i).active = 1;
            elseif a1(i).active > 0 && A(AA) > 5
    %             a1(i).Area(end+1) = a1(i).Area(end);
    %             a1(i).BoundingBox(end+1,:) = a1(i).BoundingBox(end,:);
    %             a1(i).Centroid(end+1,:) = a1(i).Centroid(end,:);
                a1(i).active = 0; %a1(i).active - 1;

                % only used one end+1 because repeated end+1 creates more rows...
                a1(end+1).Centroid(:) = b1(AA).Centroid;
                a1(end).Area = b1(AA).Area;
                a1(end).BoundingBox(:) = b1(AA).BoundingBox;
                a1(end).active = 1;
                a1(end).id = a1(end-1).id + 1;
%                 a1(end).PixelIdxList(:) = b1(AA).PixelIdxList;
                a1(end).frame = ii;
                a1(end).distance = 0;
            else
            end
                        
            if size(a1(end).Centroid(:)) == size(a1(end-1).Centroid(:))
                if ~any(a1(end).Centroid(:) - a1(end-1).Centroid(:)) & a1(end).Area == a1(end-1).Area
                    a1(end) = [];
                else
                end
            else
            end
            
%             if a1(end).Centroid(:) == a1(end-1).Centroid(:) & a1(end).Area == a1(end-1).Area
%                 a1(end) = [];
%                 
%             else
%             end
                
                
%         tmp = {a1(:).Centroid}';    
        else
        end
%         a2  = cat(1,a1.Centroid);

    end

end

% rectangle('Position', a1(3333).BoundingBox(1,:),'edgecolor', 'red');
% how to visuaalize bbox, id on image (with frame number)


toc
% tmp = {a1(:).Centroid}'; 
% for k = 1:10
%     fprintf('%d\n',k);
%     if k < 5   
%     else
%         return; % still remaining stuffs in the loop
%     end
%     fprintf('%d\n',k+5);
% end







% vector = [];
% for k = 1:size(closest,1) 
%     vector = [vector; [closest(k,3) - closest(k,1), closest(k,4) - closest(k,2)]];
% end



%  vizBlobs(subt5Cu11,Es_1,sas1); hold on;
% arrows = quiver(Es_1(:,1),Es_1(:,2),vector(:,1),vector(:,2), 'color', [1,0,0]);
%     set(arrows,'AutoScale','on', 'AutoScaleFactor',0.45); hold off;


% s(1).frameNo(:,1) = 1;

% s = struct;
% s.id = [];
% s.age = [];
% s.frameNo = [];
% s.centroid = [];
% s.size = [];
% s.diameter = [];







