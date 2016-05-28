% the trails of tears and blood

%%
clear
im1 = imread('001.jpg');
im2 = imread('002.jpg');

feature_objs(1) = featureFind(im1);
feature_objs(2) = featureFind(im2);

img_pairs = bestOf2NearestMatch(feature_objs);


%% test ransac
clear
H1 = eye(3,3);
H2 = reshape([rand(8,1) * 2; 1], 3,3);
H3 = randn(3,3); H3(3,3) = 1;
ps1 = rand(80,2) * 100;

ps2 = H1 * [ps1, ones(80,1)]';
ps2 = ps2 ./ (ones(3,1) * ps2(3, :));
ps2 = ps2(1:2, :)';
back_ps2 = ps2;
ps2 = ps2 + randn(80, 2);
ps2 = max(ps2, 0.1);

ps3 = H2 * [ps2, ones(80, 1)]';
ps3 = ps3 ./ (ones(3,1) * ps3(3,:));
ps3 = ps3(1:2, :)';
perm_idx = randperm(80)';
ps3 = ps3(perm_idx, :);

fake_feature_objs(1).points = SURFPoints(ps1);
fake_feature_objs(2).points = SURFPoints(ps2);
fake_feature_objs(3).points = SURFPoints(ps3);
img_pairs(1).idx1 = 2;
img_pairs(1).idx2 = 1;
img_pairs(1).index_corrs = repmat([1:80]', 1, 2);

img_pairs(2).idx1 = 3;
img_pairs(2).idx2 = 2;
img_pairs(2).index_corrs = [[1:80]', perm_idx];

homo_estimate = homoRANSAC(fake_feature_objs, img_pairs);

%% test homography
tf = maketform('projective', double(homographys(1).H'));
% transform image according to the homography information
imshow(imtransform(gx, tf))

%% test solve K, R with fsolve
params1 = [1 640 719/2 1 640 719/2 reshape(eye(3), 1, 9)];
[solu, fvalue] = fsolve(@homo2RK, params1)


%%
hat_R2 = reshape(solu(7:15), 3, 3)
hat_w2_norm = acos((trace(hat_R2) - 1)/2);
hat_w2 = hat_w2_norm / (2 * sin(hat_w2_norm)) * ...
    [hat_R2(6) - hat_R2(8); hat_R2(7) - hat_R2(3); hat_R2(2) - hat_R2(4)]
[solu_exp, fvalue_exp] = fsolve(@homo2RK_exp, [solu(1:6) hat_w2'])
    
%%
for i =1:length(homographys)
    H = homographys(i).H;
    src_idx = homographys(i).src_idx;
    dst_idx = homographys(i).dst_idx;
    dst_K = eye(3,3); dst_K(1,1) = cameras(dst_idx).focal; 
    dst_K(2,2) = cameras(dst_idx).focal * cameras(dst_idx).aspect;
    src_K = eye(3,3); src_K(1,1) = cameras(src_idx).focal; 
    src_K(2,2) = cameras(src_idx).focal * cameras(src_idx).aspect;
    hat_H = dst_K * cameras(dst_idx).R * inv(cameras(src_idx).R) * inv(src_K);
    hat_H = hat_H / hat_H(3, 3)
    homographys(i).H
    fprintf('fdsfasfa\n')
end

%%  reprojection error after cameraEstimate
F = 0;
for i = 1:length(homographys)
    src_idx = homographys(i).src_idx;
    dst_idx = homographys(i).dst_idx;
    src_ps = double(homographys(i).src_points');
    dst_ps = double(homographys(i).dst_points');
    %H = homographys(i).H;
    src_cam = cameras(src_idx);
    src_K = eye(3);
    src_K(1, 1) = src_cam.focal; src_K(2,2) = src_cam.focal * src_cam.aspect;
    src_K(1, 3) = src_cam.cx; src_K(2, 3) = src_cam.cy;
    dst_cam = cameras(dst_idx);
    dst_K = eye(3);
    dst_K(1, 1) = dst_cam.focal; dst_K(2,2) = dst_cam.focal * dst_cam.aspect;
    dst_K(1, 3) = dst_cam.cx; dst_K(2, 3) = dst_cam.cy;
    [U1 S1 V1] = svd(dst_cam.R);
    [U2 S2 V2] = svd(src_cam.R);
    H = dst_K * U1 * V1' * V2 * U2' * inv(src_K)
    
    tilde_src_ps = [src_ps; ones(1, size(src_ps, 2))];
    tilde_dst_ps = H * tilde_src_ps;
    hat_dst_ps = tilde_dst_ps(1:2, :) ./ repmat(tilde_dst_ps(3, :), [2, 1]);
    
    norm_r = sqrt(sum((dst_ps - hat_dst_ps).^2, 1));
    F = F + sum(norm_r(norm_r < 2) .^2);
    F = F + sum(norm_r(norm_r >=2) * 2*2 - 2^2);
end

%% test bundle adjustment
for i = 1:4
    w = rand(3,1)
    ws(:, i) = w
    W = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0]
    Ws(:,:, i) = W
    R = eye(3) + W * (sin(norm(w))/norm(w)) + W*W * ((1-cos(norm(w)))/norm(w)^2)
    Rs(:,:,i) = R;
    fake_cams(i).focal = 1000; fake_cams(i).aspect = 1; fake_cams(i).cx = 640;
    fake_cams(i).cy = 360; fake_cams(i).R = R;
end

k = 1;
for i = 1:4
    dst_cam = fake_cams(i);
    dst_K = eye(3);
    dst_K(1, 1) = dst_cam.focal; dst_K(2,2) = dst_cam.focal * dst_cam.aspect;
    dst_K(1, 3) = dst_cam.cx; dst_K(2, 3) = dst_cam.cy;
    for j = i+1:4
        src_cam = fake_cams(j);
        src_K = eye(3);
        src_K(1, 1) = src_cam.focal; src_K(2,2) = src_cam.focal * src_cam.aspect;
        src_K(1, 3) = src_cam.cx; src_K(2, 3) = src_cam.cy;
        H = dst_K * dst_cam.R * src_cam.R' * inv(src_K);
        H = H ./ H(3,3);
        fake_Hs(k).H = H;
        fake_Hs(k).src_idx = j;
        fake_Hs(k).dst_idx = i;
        src_points = 1000 * rand(20, 2);
        fake_Hs(k).src_points = src_points;
        src_points = [src_points, ones(20, 1)];
        dst_points = (H * src_points')';
        dst_points = dst_points(:, 1:2) ./ repmat(dst_points(:,3), [1,2]);
        fake_Hs(k).dst_points = dst_points;
        
        k = k+1;
    end
end

for i = 1:4
    w = ws(:, i)
    w = w + rand(3,1) * 0.01
    W = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
    R = eye(3) + W * (sin(norm(w))/norm(w)) + W*W * ((1-cos(norm(w)))/norm(w)^2);
    fake_cams(i).R = R;
end

bundleAdjustment(fake_cams, fake_Hs);