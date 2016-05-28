function [cameras] = bundleAdjustment_v2(cameras, homographys)

sigma = 2;
imgNum = length(cameras);
x0 = zeros(7, imgNum);
for j = 1:imgNum
    x0(1, j) = cameras(j).focal;
    x0(2, j) = cameras(j).aspect;
    x0(3, j) = cameras(j).cx;
    x0(4, j) = cameras(j).cy;
    [U S V] = svd(cameras(j).R);
    newR = U*V';
    if det(newR) < 0
        newR = -newR;
    end
    
    if newR == eye(3)
        x0(5:7, j) = zeros(3,1);
    else
        tmp1 = acos((trace(newR) - 1)/ 2);
        x0(5:7, j) = tmp1 * (1/(2*sin(tmp1))) * ...
            [newR(3,2) - newR(2, 3); newR(1, 3) - newR(3,1); newR(2,1)-newR(1,2)];
    end
end

Cp_inv = eye(7*imgNum);

%for i = 1:imgNum
%    Cp_inv((i-1)*7 + 5, (i-1)*7 + 5) = (16/pi)^2;
%    Cp_inv((i-1)*7 + 6, (i-1)*7 + 6) = (16/pi)^2;
%    Cp_inv((i-1)*7 + 7, (i-1)*7 + 7) = (16/pi)^2;
%    Cp_inv((i-1)*7 + 1, (i-1)*7 + 1) = (10/ cameras(1).focal)^2;
%end

x = bundle_lm_solver(@calFun, Cp_inv,  x0(:));

x = reshape(x, [7, imgNum]);
for j = 1:imgNum
    cameras(j).focal = x(1, j);
    cameras(j).aspect = x(2, j);
    cameras(j).cx = x(3, j);
    cameras(j).cy = x(4, j);
    
    if x(5:7, j) == 0
        cameras(j).R = eye(3);
    else
        tmp1 = x(5:7, j);
        tmp2 = [0 -x(7, j) x(6, j); x(7, j) 0 -x(5, j); -x(6, j) x(5, j) 0];
        cameras(j).R = eye(3,3) + (tmp2./norm(tmp1)) * sin(norm(tmp1)) + ...
            (tmp2*tmp2)./dot(tmp1, tmp1) * (1-cos(norm(tmp1)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [F, r, J] = calFun(x)
        num_matchs = 0;
        for i = 1:length(homographys)
            num_matchs = num_matchs + size(homographys(i).src_points, 1);
        end
        
        r = calErr(homographys, x, imgNum, num_matchs);
        J = calJ(homographys, x, imgNum, num_matchs);
         
        % calculate scalar whole errors
        F = 0;
        norm_r = sqrt(sum(r.^2, 1));
        F = F + sum(norm_r(norm_r < sigma) .^ 2);
        F = F + sum(norm_r(norm_r >= sigma) * 2 * sigma - sigma^2);
        
        fprintf('quadratic err %f,  other %f\n', norm_r * norm_r', F)
        r = r(:);
    end

end

function [err] = calErr(homographys, x, imgNum, num_matchs)
    err = zeros(2, num_matchs);
    mat_x = reshape(x, [7, imgNum]);
    
    r_cnt = 0;
    for i = 1:length(homographys)
        src_idx = homographys(i).src_idx;
        dst_idx = homographys(i).dst_idx;
        src_ps = homographys(i).src_points';
        dst_ps = homographys(i).dst_points';
           
        tilde_src_ps = [src_ps; ones(1, length(src_ps))];
        [src_K, src_R] = genKR(mat_x(:, src_idx));
        [dst_K, dst_R] = genKR(mat_x(:, dst_idx));
        hat_tilde_dst_ps = dst_K * dst_R * src_R' * inv(src_K) * tilde_src_ps;
        hat_dst_ps = hat_tilde_dst_ps(1:2, :) ./ repmat(hat_tilde_dst_ps(3,:), [2,1]);
        r_i = (dst_ps - hat_dst_ps);
        err(:, r_cnt+1:r_cnt+size(r_i,2)) = r_i;
        r_cnt = r_cnt + size(r_i, 2);
    end
    
end

function [J] = calJ(homographys, x, imgNum, num_matchs)
    delta = 1e-4;
    mat_x = reshape(x, [7, imgNum]);
    J = zeros(num_matchs * 2, 7 * imgNum);
    
    
    for j = 1:imgNum
        % i:1~7 focal aspect cx cy w1 w2 w3
        for i = 1:7
            x_plus = mat_x;  x_plus(i, j) = x_plus(i, j) + delta;
            err_plus = calErr(homographys, x_plus, imgNum, num_matchs);
            x_minus = mat_x; x_minus(i, j) = x_minus(i, j) - delta;
            err_minus = calErr(homographys, x_minus, imgNum, num_matchs);
            err_plus = err_plus(:);
            err_minus = err_minus(:);
            J(:, 7*(j-1)+i) = (err_plus - err_minus) ./ (2*delta);
        end
    end
    
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K, R] = genKR(x)
K = eye(3,3);
K(1, 1) = x(1); K(2,2) = K(1,1) * x(2);
K(1, 3) = x(3); K(2, 3) = x(4);
w = x(5:7);
W = [0      -w(3) w(2); 
     w(3)   0     -w(1);
     -w(2)   w(1)  0];
if w == 0
    R = eye(3);
else
R = eye(3) + W ./ norm(w) * sin(norm(w)) + W*W ./ dot(w,w) * (1 - cos(norm(w)));
end
end