imgNum = 5;
imgFolder = './';
for i = 1:imgNum
    full_imgs(:,:,:,i) = imread(sprintf([imgFolder, 'undis%03d.jpg'], i-1));
    gray_imgs(:,:,i) = rgb2gray(full_imgs(:,:,:,i));
end
imgSize = size(gray_imgs(:,:,1));

for i = 1:imgNum
    feature_objs(i) = featureFind(gray_imgs(:,:,i));
end

img_pairs = bestOf2NearestMatch(feature_objs);

warning('off')
homographys = homoRANSAC(feature_objs, img_pairs);
warning('on')

if ~checkPano(homographys, imgNum)
    return 
end

cameras = cameraEstimate(homographys, imgNum, imgSize);

hat_cams = bundleAdjustment(cameras, homographys);
%hat_cams = bundleAdjustment_v2(cameras, homographys);

hat_cams_waved = waveCorrectHoriz(hat_cams);

[canvas, warped_imgs, warped_masks, corners] = cylindricalWarping(hat_cams_waved, full_imgs);

class(canvas)
imshow(uint8(canvas))