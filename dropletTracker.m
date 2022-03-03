classdef dropletTracker < handle
    %DROPLETTRACKER class to track Pb liquid droplets in X-ray radiography
    %sequences of solidifying Al alloy

    % =====================================================================
    %% PROPERTIES

    properties(GetAccess = private, SetAccess = private) % Internal use
        settings_
        tracks_
        assignments_
        unassignedTracks_
        unassignedDetections_
        nextID_
    end
    properties(Dependent) % For user
        settings
        droplets
        frames
        ntracks
    end
    properties(Dependent, Hidden) % Internal use only
        diameter
        velocity
        acceleration
    end

    % =====================================================================
    %% METHODS

    % CONSTRUCTOR
    methods
        function obj = dropletTracker(S)
            %DROPLETTRACKER construct a dropletTracker class
            %object.

            arguments
                S struct = []
            end

            % Initialise tracks
            obj.initializeTracks;

            % Initialize counter
            obj.nextID_ = 1;

            % Set input parameters
            obj.settings_ = S;

        end
    end

    % PROPERTY ACCESS
    methods
        function s = get.settings(obj)
            s = obj.settings_;
        end
        function d = get.droplets(obj)
            ntrk = obj.ntracks; % Get number of tracks

            % Extract already stored data
            [d(1:ntrk).id]           = obj.tracks_.id;
            [d(1:ntrk).age]          = obj.tracks_.age;
            [d(1:ntrk).frameNo]      = obj.tracks_.frameNo;
            [d(1:ntrk).time]         = obj.tracks_.time;
            [d(1:ntrk).centroid]     = obj.tracks_.centroid;
            [d(1:ntrk).diameter]     = deal(obj.diameter{:});
            [d(1:ntrk).velocity]     = deal(obj.velocity{:});
            [d(1:ntrk).acceleration] = deal(obj.acceleration{:});
        end
        function f = get.frames(obj)
            f = unique(cat(1,obj.droplets.frameNo));
        end
        function n = get.ntracks(obj)
            n = length(obj.tracks_);
        end
        function d = get.diameter(obj)
            for n = 1:obj.ntracks
               d{n,1} = [obj.tracks_(n).eqdiam] * obj.settings_.resolution;
            end
        end
        function v = get.velocity(obj)
            % v = [vel, vx, vy] all vel are in micron/s
            
            for n = 1:obj.ntracks % For each droplet

                % Get the centres in microns
                centretmp = [obj.tracks_(n).centroid(:,1:2)]...
                    * obj.settings_.resolution; % in microns
                timeStep  = obj.tracks_(n).time(2:end) -...
                    obj.tracks_(n).time(1:end-1); % in seconds
                vx  = [0;...
                    (centretmp(2:end,1) - centretmp(1:end-1,1)) ./ timeStep];
                vy  = [0;...
                    (centretmp(2:end,2) - centretmp(1:end-1,2)) ./ timeStep];
                vel = sqrt(vx .^2 + vy .^2);
                v{n,1} = [vel,vx,vy];
            end
        end
        function a = get.acceleration(obj)
            % a = [acc, ax, ay] all acc in micron/s^2
            
            vel = obj.velocity;
            for n = 1:obj.ntracks % For each droplet
                timeStep  = obj.tracks_(n).time(2:end) -...
                    obj.tracks_(n).time(1:end-1); % in seconds
                ax = [0;...
                    (vel{n}(2:end,2) - vel{n}(1:end-1,2)) ./ timeStep];
                ay = [0;...
                    (vel{n}(2:end,3) - vel{n}(1:end-1,3)) ./ timeStep];
                acc = sqrt(ax .^2 + ay .^2);
                a{n,1} = [acc,ax,ay];
            end
        end
    end

    % TRACKING
    methods
        function step(obj,frn,P,D)
            %STEP

            obj.predictNewLocationsOfTracks(frn,P);
            obj.detectionToTrackAssignment(frn,D.centroid,D.flow);
            obj.updateAssignedTracks(D.centroid,D.bbox,D.diameter,D.meanInt);
            obj.updateUnassignedTracks();
            %deleteLostTracks();
            obj.deactivateLostTracks();
            obj.createNewTracks(frn,P,D);
        end
    end

    % FUNCTIONALITIES
    methods
        function vOut = velocities(obj,idxIn,varargin)
            % VELOCITY return the velocity of droplets.
            %
            % v = velocity(obj,frames)
            % v = velocity(obj,times,'time') NOT YET
            % v = velocity(obj,drtIDs,'ID') NOT YET
            % v = velocity(obj,...'parameters',...) NOT YET

            % Set the form of the output data. Could be 'frame','time' or
            % 'ID'
            if ~isempty(varargin) && isinput('time')
                outputType = 'time';
            elseif ~isempty(varargin) && isinput('ID')
                outputType = 'ID';
            else
                outputType = 'frame'; % Default
            end

            droplets = obj.droplets;
            vel      = cat(1,droplets.velocity);
            centres  = cat(1,droplets.centroid);
            switch outputType
                case 'frame' % all the velocity in the input frames
                    detectionFrames = cat(1,obj.droplets.frameNo);
                    for n = 1:numel(idxIn)
                        idxtmp = detectionFrames == idxIn(n);
                        vOut{n,2} = idxIn(n);
                        vOut{n,1} = [vel(idxtmp,:),centres(idxtmp,1:2)];
                    end


            end




        end
    end

    % VISUALIZATION
    methods
        function displaytracking(obj,frn,D,H)
            % Convert the frame and the mask to uint8 RGB.
            frame = im2uint8(D.img);
            frame = insertText(frame,[50 size(frame,1)*0.9],...
                ['FrameNo: ',num2str(frn)],'FontSize',80,'BoxColor','y',...
                'BoxOpacity',0.4);imshow(frame)

            if ~isempty(obj.tracks_)

                % Noisy detections tend to result in short-lived tracks.
                % Only display tracks that have been visible for more than
                % a minimum number of frames.
                reliableTrackInds = ...
                    [obj.tracks_(:).totalVisibleCount] > obj.settings.minVisibleCount;
                activeTrackInds = [obj.tracks_(:).active] == 1;
                reliableTracks = obj.tracks_(reliableTrackInds & activeTrackInds);

                % Display the objects. If an object has not been detected
                % in this frame, display its predicted bounding box.
                if ~isempty(reliableTracks)
                    % Get bounding boxes.
                    bboxes = cat(1, reliableTracks.bbox);

                    % Get ids.
                    ids = int32([reliableTracks(:).id]);

                    % Create labels for objects indicating the ones for
                    % which we display the predicted rather than the actual
                    % location.
                    labels = cellstr(int2str(ids'));
                    predictedTrackInds = ...
                        [reliableTracks(:).consecutiveInvisibleCount] > 0;
                    isPredicted = cell(size(labels));
                    isPredicted(predictedTrackInds) = {' predicted'};
                    labels = strcat(labels, isPredicted);

                    % Draw the objects on the frame.
                    frame = insertObjectAnnotation(frame, 'rectangle', ...
                        bboxes, labels,'FontSize',24);

                    %                 % Draw the objects on the mask.
                    %                 mask = insertObjectAnnotation(mask, 'rectangle', ...
                    %                     bboxes, labels);

                    % Display the mask and the frame.
                    imshow(frame)
                    pause(0.1);
                else
                    imshow(frame);
                    pause(0.1);
                end

            else
                imshow(frame);
                pause(0.1);
            end
        end
    end
    
    % INTERNAL USE ONLY
    methods(Access = private)
        function initializeTracks(obj)
            obj.tracks_ = newtrack;
        end
        function predictNewLocationsOfTracks(obj,frn,P)
            for i = 1:length(obj.tracks_)

                if ~obj.tracks_(i).active; continue; end

                % Predict the current location of the track.
                predictedCentroid = predict(obj.tracks_(i).kalmanFilter);

                % Shift the bounding box so that its center is at
                % the predicted location.
                bbox = obj.tracks_(i).bbox;
                predictedCentroid = int32(predictedCentroid) - int32(bbox(3:4)) / 2;
                obj.tracks_(i).bbox              = [predictedCentroid, bbox(3:4)];
                obj.tracks_(i).frameNo(end+1,:)  = frn;
                obj.tracks_(i).time(end+1,:)     = P.time;
                obj.tracks_(i).centroid(end+1,:) = [predictedCentroid 0];
                obj.tracks_(i).eqdiam(end+1,:)   = 0;
            end
        end
        function detectionToTrackAssignment(obj,frn,centroids,flow)

            % NOTE: modified to ignore inactive tracks
            if isempty(obj.tracks_)
                idxactive = [];
            else
                idxactive   = find(cat(1,obj.tracks_.active));
            end
            nTracks     = length(idxactive);
            nDetections = size(centroids, 1);          

            % Match
            if strcmpi(obj.settings.algorithm,'hungarian')
                
                % Compute the cost of assigning each detection to each
                % track.
                cost = zeros(nTracks, nDetections);
                for i = 1:nTracks
                    cost(i, :) = distance(obj.tracks_(idxactive(i)).kalmanFilter,...
                        centroids);
                end

                % Solve the assignment problem.

                [assignmentstmp, unassignedTrackstmp, obj.unassignedDetections_] = ...
                    assignDetectionsToTracksLocal(cost, obj.settings.costOfNonAssignment);
            
            elseif strcmpi(obj.settings.algorithm,'knn')

                % Update old centre location with optical flow
                ctrOldMoved = double.empty;
                maxImgSize = size(flow.Vx);
                for i = 1:nTracks

                    % NOTE: the centroids for the last frame are the end-1
                    % row as the Kalman filter already predicted the new
                    % location
                    xytmp = uint16(obj.tracks_(idxactive(i)).centroid(end-1,[1,2]));

                    % If xytmp is out of the image use the kalman filter
                    if xytmp(2)-2 < 1 ||...
                            xytmp(2)+2 > maxImgSize(1) ||...
                            xytmp(1)-2 < 1 ||...
                            xytmp(1)+2 > maxImgSize(2)

                        xdisp = obj.tracks_(idxactive(i)).centroid(end,1);
                        ydisp = obj.tracks_(idxactive(i)).centroid(end,2);
                    else
                        % Using optical flow to predict the locations in the
                        % current frame
                        ctr_dx = mean2(flow.Vx(...
                            xytmp(2)-2:xytmp(2)+2,xytmp(1)-2:xytmp(1)+2));
                        xdisp = xytmp(1) + ctr_dx;
                        ctr_dy(i) = mean2(flow.Vy(...
                            xytmp(2)-2:xytmp(2)+2,xytmp(1)-2:xytmp(1)+2));
                        ydisp = xytmp(2) + ctr_dy(i);
                    end
                    
                    % Stored the predicted location
                    ctrOldMoved(i,:) = [xdisp, ydisp];
                end

                % Temporary Hack to get the predicted location from Kalman
                % filter
%                 ctrOldMoved = cell2mat(arrayfun(@(x) x.centroid(end,[1,2]),...
%                    obj.tracks_(idxactive)','UniformOutput',false));

                [assignmentstmp, unassignedTrackstmp, obj.unassignedDetections_] =...
                    knnsearchInt(centroids,ctrOldMoved,obj.settings.maxDist);

            else
                error('Tracking algorithm not recognised.')
            end

            obj.assignments_      = [idxactive(assignmentstmp(:,1)),...
                assignmentstmp(:,2)];
            obj.unassignedTracks_ = idxactive(unassignedTrackstmp);
        end
        function updateAssignedTracks(obj,centroids,bboxes,diamframe,meanInt)
            numAssignedTracks = size(obj.assignments_, 1);
            for i = 1:numAssignedTracks
                trackIdx     = obj.assignments_(i, 1);
                detectionIdx = obj.assignments_(i, 2);
                centroid     = centroids(detectionIdx, :);
                bbox         = bboxes(detectionIdx, :);
                eqdiam       = diamframe(detectionIdx,:);
                intensity    = meanInt(i,:);

                % Correct the estimate of the object's location
                % using the new detection.
                centroidc = correct(obj.tracks_(trackIdx).kalmanFilter, centroid);

                % Replace predicted bounding box with detected
                % bounding box.
                obj.tracks_(trackIdx).bbox = bbox;

                % Update track's age.
                obj.tracks_(trackIdx).age = obj.tracks_(trackIdx).age + 1;

                % Update visibility.
                obj.tracks_(trackIdx).totalVisibleCount = ...
                    obj.tracks_(trackIdx).totalVisibleCount + 1;
                obj.tracks_(trackIdx).consecutiveInvisibleCount = 0;

                % Update intensity
                obj.tracks_(trackIdx).intensity = intensity;

                % Store current location and diameter
                %             tracks(trackIdx).frameNo(end,:)   = frn;
                %             tracks(trackIdx).time(end,:)      = time;
                obj.tracks_(trackIdx).centroid(end,:) = [centroidc 1];
                obj.tracks_(trackIdx).eqdiam(end,:)    = eqdiam;
            end
        end
        function updateUnassignedTracks(obj)
            for i = 1:length(obj.unassignedTracks_)
                ind = obj.unassignedTracks_(i);
                obj.tracks_(ind).age = obj.tracks_(ind).age + 1;
                obj.tracks_(ind).consecutiveInvisibleCount = ...
                    obj.tracks_(ind).consecutiveInvisibleCount + 1;
            end
        end
        function deactivateLostTracks(obj)
            if isempty(obj.tracks_)
                return;
            end

            % Compute the fraction of the track's age for which it was visible.
            ages = [obj.tracks_(:).age];
            totalVisibleCounts = [obj.tracks_(:).totalVisibleCount];
            visibility = totalVisibleCounts ./ ages;

            % Find the indices of 'lost' tracks.
            lostInds = (ages < obj.settings_.ageThreshold & visibility < 0.6) | ...
                [obj.tracks_(:).consecutiveInvisibleCount] >=...
                obj.settings_.invisibleForTooLong;

            % Deactivate lost tracks.
            if ~any(lostInds)
                return;
            else
                idxx = find(lostInds);
                for nn = 1:numel(idxx)
                    obj.tracks_(idxx(nn)).active = 0;
                end
            end
        end
        function createNewTracks(obj,frn,P,D)
            %centroids = centroids(obj.unassignedDetections_, :);
            %bboxes = bboxes(obj.unassignedDetections_, :);
            centroids   = D.centroid(obj.unassignedDetections_,:);
            bboxes      = D.bbox(obj.unassignedDetections_, :);
            eqdiams     = D.diameter(obj.unassignedDetections_,:);
            intensities = D.meanInt(obj.unassignedDetections_,:);
            for i = 1:size(centroids, 1)

                % Retrieve current obj info
                centroid  = centroids(i,:);
                bbox      = bboxes(i, :);
                eqdiam    = eqdiams(i,:);
                intensity = intensities(i,:);

                % Create a Kalman filter object.
                kalmanFilter = configureKalmanFilter('ConstantVelocity', ...
                    centroid,...
                    obj.settings.InitialEstimateError,...
                    obj.settings.MotionNoise,...
                    obj.settings.MeasurementNoise);

                % Create a new track.
                newTrack = newtrack(...
                    'id', obj.nextID_, ...
                    'bbox', bbox, ...
                    'kalmanFilter', kalmanFilter, ...
                    'age', 1, ...
                    'totalVisibleCount', 1, ...
                    'consecutiveInvisibleCount', 0, ...
                    'active', 1, ...
                    'intensity', intensity,...%                     'resolution', P.resolution, ...
                    'frameNo',frn, ...
                    'time', P.time, ...
                    'centroid', [centroid, 1], ...
                    'eqdiam', eqdiam, ...
                    'vel', 0, ...
                    'acc', 0);

                % Add it to the array of tracks.
                obj.tracks_(end + 1) = newTrack;

                % Increment the nextID counter.
                obj.nextID_ = obj.nextID_ + 1;
            end
        end
    end
end

% CLASS INTERNAL FUNCTIONS
function track = newtrack(T)
arguments
    T.id = {}
    T.bbox = {}
    T.kalmanFilter = {}
    T.age = {}
    T.totalVisibleCount = {}
    T.consecutiveInvisibleCount = {}
    T.active = {}
    T.intensity = {}
    %     T.resolution = {}
    T.frameNo = {}
    T.time = {}
    T.centroid = {}
    T.eqdiam = {}
    T.vel = {}
    T.acc = {}
end

track = struct(...
    'id', T.id, ...
    'bbox', T.bbox, ...
    'kalmanFilter', T.kalmanFilter, ...
    'age', T.age, ...
    'totalVisibleCount', T.totalVisibleCount, ...
    'consecutiveInvisibleCount', T.consecutiveInvisibleCount, ...
    'active', T.active,...
    'intensity', T.intensity,...%     'resolution', T.resolution,...
    'frameNo', T.frameNo,...
    'time', T.time,...
    'centroid', T.centroid, ...
    'eqdiam', T.eqdiam, ...
    'vel', T.vel, ...
    'acc', T.acc);
end
function TF = isinput(inputs,argument)

if any(strcmpi((inputs),argument))
    TF = 1;
else
    TF = 0;
end

end
function v = velCalculator()

end
function [assignments, unassignedTrks, unassignedDets] = knnsearchInt(newDetections,oldCentres,maxDist)

% Initialize output
assignments    = uint32.empty(0,2);
unassignedTrks = uint32.empty(0,1); % Unassigned tracks
unassignedDets = uint32.empty(0,1); % Unassigned detections

% NOTE: 1) In the case of no new detections (newDetections = []) --> no
% assignment, all tracks are unassigned and no unassigned detections (if
% oldCentres is empty it means there are no existing tracks so
% "unassignedTrks" will be empty). 2) In the case of no existing tracks
% (oldCentres = []) --> no assignments and no unassigned tracks.
% "unassignedDets" = "newDetections"  all new unassigned detections.
if isempty(newDetections)
    unassignedTrks = [1:size(oldCentres)]';
    return
end
if isempty(oldCentres)
    unassignedDets = [1:size(newDetections)]';
    return
end

[IdxNN,dist] = knnsearch(newDetections,oldCentres,'k',1);

% Sort out unassigned tracks and detections

% If maxDist is a value then filter out bad assignment
if ~isempty(maxDist)
    IdxNN(dist>maxDist) = NaN;
end

unassignedTrks = find(isnan(IdxNN));
% IdxNN(isnan(IdxNN)) = []; % Remove bad assignment
assignments    = [find(~isnan(IdxNN)),IdxNN(~isnan(IdxNN))];
unassignedDets = find(~ismember(1:size(newDetections,1),assignments))'; 

% listOfTrackCrt = [];
% for nn = 1:numel(oldCentres)
%     if isempty(assignments(nn)) || assignments(nn) <1; continue; end
%     if dist(nn) > 15
%         listOfTrackCrt(end+1) = assignments(nn); % Store index
%         continue;
%     end
%     ctemp = ctrframe(assignments(nn),:);
%     diam  = diamframe(assignments(nn));
%     vxtmp = vx(assignments(nn),1);
%     vytmp = vy(assignments(nn),1);
% 
%     % Structure of trktmp (table format)
%     % [ x | y | frNo | vx | vy | eqDiam | deltaT | res | trNo]
%     trtmp = array2table([ctemp, frn, vxtmp, vytmp, diam, deltaT, resolution,nn],...
%         'VariableNames',{'X','Y','FrameNo','vX','vY','Diameter',...
%         'deltaT','Resolution','TrackNo'});
%     tracks{nn} = [tracks{nn};trtmp];
%     %                 tracks{nn}(end+1,:) = [ctemp, frn, nn, diam];
%     listOfTrackCrt(end+1) = assignments(nn); % Store index
% end
% % Create new tracks
% idxNewCtr = find(~ismember(1:size(ctrframe,1),listOfTrackCrt));
% for nnn = 1:numel(idxNewCtr)
%     noTrk = numel(tracks);
%     trtmp = array2table([ctrframe(idxNewCtr(nnn),:),...
%         frn, vx(idxNewCtr(nnn),1), vy(idxNewCtr(nnn),1),...
%         diamframe(idxNewCtr(nnn)), deltaT, resolution,noTrk+1],...
%         'VariableNames',{'X','Y','FrameNo','vX','vY','Diameter',...
%         'deltaT','Resolution','TrackNo'});
%     tracks = [tracks,...
%         {trtmp}];
%     %                 tracks = [tracks,...
%     %                     {[ctrframe(idxNewCtr(nnn),:),...
%     %                     frn,noTrk+1,diamframe(idxNewCtr(nnn))]}];
% end

end

%% DEVELOPMENT NOTES:

% LAST UPDATE: 20211017 - EL

% =========================================================================
% ##### TODO

% =========================================================================
% ##### JOURNAL

%   @ 20211017 EL: Wrote the backbone of the class. The code is base on the
%   algorithm in the scrip s_20181207_AutomaticDropletTracking_BulkLiquid
%   developed to track Pb droplets in mt9140 and ma4035
%   beamtime videos. Need to be tested.

%   @ 20220209 EL: adding second algorithm based on KNN search instead of
%   hungarian algorithm. This algorithm should be faster. 

% =========================================================================
% ##### VERSION HISTORY

% - version 0.1 (date, programmer name)
%   + Feature1.
%   + Feature2...