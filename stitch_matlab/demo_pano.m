function [cameras, rScale, aScale, warped_imgs, warped_masks, corners, canvas, canvas_corner] = demo_pano(scale, imgNum, imgFolder)
addpath(genpath('maxflow'));

if (nargin == 1)
    imgFolder = './input/';
    imgNum = 16
elseif nargin == 2
    imgFolder = './input/';
end

for i = 1:imgNum
    full_imgs(:,:,:,i) = imresize(imread(sprintf([imgFolder, '%03d_l.PNG'], i)), scale);
    %full_imgs(:,:,:,i) = imread(sprintf([imgFolder, '%03d_r.PNG'], i));
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

% hat_cams = bundleAdjustment(cameras, homographys);
hat_cams = bundleAdjustment_v2(cameras, homographys);

hat_cams_waved = waveCorrectHoriz(hat_cams);

cameras = hat_cams_waved;

aScale = 20;
[aScale, rScale, warped_imgs, warped_masks, corners] = cylindricalWarping(hat_cams_waved, full_imgs, aScale);
warped_masks = logical(warped_masks);

[warped_masks] = gc_graphcut(warped_imgs, warped_masks, corners, aScale);

[canvas, canvas_corner] = imageblending(warped_imgs, warped_masks, corners, aScale);

canvas_clip = canvas(:, 1:360 * aScale - 1, :);
class(canvas_clip)
figure; imshow(uint8(canvas_clip))

end