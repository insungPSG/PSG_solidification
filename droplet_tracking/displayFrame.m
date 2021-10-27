function displayFrame(D,H)

% imagesc(H.ax,'CData',cat(3,D.imgseg,zeros(size(D.img)),D.img));
imshow(imoverlay(D.img,D.imgseg));

pause(0.1)
