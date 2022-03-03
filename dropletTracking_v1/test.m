thresh = multithresh(normalise(i),5);
seg_I = imquantize(i,thresh);
% RGB = label2rgb(seg_I); 	 
figure;
% imagesc(seg_I.*(seg_I<=3)),colormap(jet(4))
imshowpair(i,seg_I<=2);
axis off
title('RGB Segmented Image')

%% 

i(:,:,1) =imdiff(P.obj,frn,1);
i(:,:,2) =imdiff(P.obj,frn,3);
i(:,:,3) =imdiff(P.obj,frn,5);
imshow(i,[])