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
    properties(Dependent, Hidden) % Internal use only
        
    end
    
    properties(Dependent) % For user
        settings
        ntracks
        droplets
        velocity
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
        function d = get.droplets(obj)
            d = obj.tracks_;
        end
        function s = get.settings(obj)
            s = obj.settings_;
        end
        function n = get.ntracks(obj)
            n = length(obj.tracks_);
        end
        function v = get.velocity(obj)
        
            
        end
    end

    % FUNCIONALITIES
    methods
        function step(obj,frn,P,D)
            %STEP
           
            obj.predictNewLocationsOfTracks(frn,P);
            obj.detectionToTrackAssignment(frn,D.centroid);
            obj.updateAssignedTracks(D.centroid,D.bbox,D.diameter,D.meanInt);
            obj.updateUnassignedTracks();
            %deleteLostTracks();
            obj.deactivateLostTracks();
            obj.createNewTracks(frn,P,D);
        end
    end

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
        function detectionToTrackAssignment(obj,frn,centroids)

            % NOTE: modified to ignore inactive tracks
            if isempty(obj.tracks_)
                idxactive = [];
            else
                idxactive   = find(cat(1,obj.tracks_.active));
            end
            nTracks     = length(idxactive);
            nDetections = size(centroids, 1);

            % Compute the cost of assigning each detection to each track.
            cost = zeros(nTracks, nDetections);
            for i = 1:nTracks
                cost(i, :) = distance(obj.tracks_(idxactive(i)).kalmanFilter,...
                    centroids);
            end

            % Solve the assignment problem.
            
            [assignmentstmp, unassignedTrackstmp, obj.unassignedDetections_] = ...
                assignDetectionsToTracks(cost, obj.settings.costOfNonAssignment);
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
                    'intensity', intensity,...
                    'resolution', P.resolution, ...
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
end

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
    T.resolution = {}
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
            'intensity', T.intensity,...
            'resolution', T.resolution,...
            'frameNo', T.frameNo,...
            'time', T.time,...
            'centroid', T.centroid, ...
            'eqdiam', T.eqdiam, ...
            'vel', T.vel, ...
            'acc', T.acc);
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

% =========================================================================
% ##### VERSION HISTORY

% - version 0.1 (date, programmer name)
%   + Feature1.
%   + Feature2...