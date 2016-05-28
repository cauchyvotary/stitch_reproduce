function [image_pairs] = bestOf2NearestMatch(feature_objs);

imageNum = length(feature_objs);
th1 = 0.1;

cnt = 1;
for i = 1:imageNum
    features1 = feature_objs(i).features;
    
    for j = i+1:imageNum
        features2 = feature_objs(j).features;
        [index_corrs, corrs_metrics] = matchFeatures(features1, features2);
        
        if size(index_corrs, 1) < min(size(features1, 1), size(features2, 1)) * th1
            continue;
        end
        
        image_pair.index_corrs = index_corrs;
        image_pair.corrs_metrics = corrs_metrics;
        image_pair.idx1 = i;
        image_pair.idx2 = j;
        image_pairs(cnt) = image_pair;
        cnt = cnt + 1;
    end
end
