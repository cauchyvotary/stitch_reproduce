function feature_obj = featureFind( im )

if size(im, 3) ~= 1
    im = rgb2gray(im);
end

points = detectSURFFeatures(im);
[features, points] = extractFeatures(im, points);

feature_obj.features = features;
feature_obj.points = points;