function [homographys] = homoRANSAC(feature_objs, img_pairs)

setNum = 4;
trialNum = 500;
error_tol = 2; % error tolerance of 2 pixels
homo_cnt = 0;
homographys = [];

for i = 1:length(img_pairs)
    img_idx1 = img_pairs(i).idx1;
    img_idx2 = img_pairs(i).idx2;
    idx_corrs = img_pairs(i).index_corrs;
    m_points1 = feature_objs(img_idx1).points(idx_corrs(:, 1)).Location;
    m_points2 = feature_objs(img_idx2).points(idx_corrs(:, 2)).Location; 
    mNum = length(m_points1);
    
    inliersNum = 0;
    inlier_ps1 = [];
    inlier_ps2 = [];
    H = zeros(3,3);
    suc_count = 0;
    for pha = 1:trialNum
        sampled_idx = randperm(mNum, setNum);
        sampled_ps1 = m_points1(sampled_idx, :);
        sampled_ps2 = m_points2(sampled_idx, :);
        [sampled_H, suc_flag] = homography4points(sampled_ps1, sampled_ps2);
        
        if suc_flag
            % test inliers
            suc_count  = suc_count + 1;
            tmp_ps2 = [m_points2, ones(mNum, 1)]';
            tmp_ps1 = [m_points1, ones(mNum, 1)]';
            hat_ps1 = sampled_H * tmp_ps2;
            hat_ps1 = hat_ps1 ./ (ones(3,1) * hat_ps1(3, :));
            tmp_ps1 = tmp_ps1(1:2, :);
            hat_ps1 = hat_ps1(1:2, :);
            
            %tmp_inliersNum = sum(sqrt(sum((tmp_ps1 - hat_ps1).^2, 1)) <= error_tol);
            tmp_inlier_idx = find(sqrt(sum((tmp_ps1 - hat_ps1).^2, 1)) <= error_tol);
            tmp_inliersNum = length(tmp_inlier_idx);
            tmp_inlier_ps1 = m_points1(tmp_inlier_idx, :);
            tmp_inlier_ps2 = m_points2(tmp_inlier_idx, :);
            
            %tmp_inliersNum / mNum
            if tmp_inliersNum > inliersNum
                H = sampled_H;
                inliersNum = tmp_inliersNum;
                inlier_ps1 = tmp_inlier_ps1;
                inlier_ps2 = tmp_inlier_ps2;
            end
        end
    end
    suc_count
    
    % probablistic check
    alpha = 8;
    beta = 0.3;
    if inliersNum > alpha + beta * mNum
        homo_cnt = homo_cnt + 1;
        homographys(homo_cnt).inliersNum = inliersNum;
        homographys(homo_cnt).H = H;
        homographys(homo_cnt).src_idx = img_idx2;
        homographys(homo_cnt).dst_idx = img_idx1;
        homographys(homo_cnt).src_points = inlier_ps2;
        homographys(homo_cnt).dst_points = inlier_ps1;
    end
end