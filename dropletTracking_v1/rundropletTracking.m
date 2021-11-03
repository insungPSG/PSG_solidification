function tracker = rundropletTracking
% RUNDROPLETTRACKING is the main function to run the droplet tracknig
% algorithm to meassure fluid flow in X-ray radiography video sequence
% collected during I20210153 beamtime at Petra III (2021).

% ==========================================================================
% LOAD PARAMETERS

P = parsdropletTracking;

% ==========================================================================
% INITIALIZATION

tracker     = dropletTracker(P.TR);
D           = initializeDataStructure;  % To speed up calculation
opticalflow = initializeOpticalFlow(P); % To detect moving droplet   
H           = initializeGraphic(P); 
m = 1;

% ==========================================================================
% TRACKING

for frn = P.frameRange
%     sprintf('step #%i',P.frameRange+m);
    P.time = P.obj.timeframeno(frn,'time'); % Convert frame number to time 

    [D,P] = detectDroplets(frn,opticalflow,P,D);
    % D : image information saved, P: parameters
    tracker.step(frn,P,D);
    
    if P.displayTracking || P.saveVideoFlag
        tracker.displaytracking(frn,D,H);
    end
    if P.saveResults
        save(P.resultFilePath,'tracker');
    end
    if P.saveVideoFlag
        writeVideo(H.vid,getframebg(H.fig));
%         writeVideo(H.vid,getframe(H.fig)); % want to figure out getframe or getframbg?
    end
    
    waitbar(m / numel(P.frameRange),H.w,...
        sprintf('Frame no. %i',frn));
    
    m = m + 1;
    
end

% ==========================================================================
% SAVE RESULTS - TO IMPLEMENT
% export as a .mat file?
% ==========================================================================

if saveVideoFlag
    close(H.vid);
    set(H.fig,'Visible','on');
end
close(H.w);

%% NOTES

% @ 2021 10 21 - EL

% Completed a working beta version. Tracking works and some basic data
% plotting capabilities

% =========================================================================
% ##### TODO
