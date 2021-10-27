function H = initializeGraphic(P)

H.fig = figure('position',[100 100 1400 800]);
H.ax = axes('Position',[0 0 1 1]);
axis ij equal off

H.w = waitbar(0,'Tracking droplets...');

if P.saveVideoFlag
    H.vid = VideoWriter(P.videoFilePath,'Motion JPEG AVI');
    H.vid.FrameRate = P.videoFrameRate;
    open(H.vid);
    set(H.fig,'Visible','off');
end
