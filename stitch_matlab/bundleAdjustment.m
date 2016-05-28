function [cameras] = bundleAdjustment(cameras, homographys)

sigma = 2;
imgNum = length(cameras);
x0 = zeros(7, imgNum);
for j = 1:imgNum
    x0(1, j) = cameras(j).focal;
    x0(2, j) = cameras(j).aspect;
    x0(3, j) = cameras(j).cx;
    x0(4, j) = cameras(j).cy;
    [U S V] = svd(cameras(j).R);
    newR = U*V'
    if newR == eye(3)
        x0(5:7, j) = zeros(3,1);
    else
        tmp1 = acos((trace(newR) - 1)/ 2)
        x0(5:7, j) = tmp1 * (1/(2*sin(tmp1))) * ...
            [newR(3,2) - newR(2, 3); newR(1, 3) - newR(3,1); newR(2,1)-newR(1,2)]
    end
end

reproj(x0)

options.Algorithm = 'levenberg-marquardt';
options.Display = 'iter';
options.LargeScale = 'off';
[x, fval] = lsqnonlin(@reproj, x0, [], [], options)   
sqrt(fval)

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


%%%%%%%%%%%%%%%%%%%%
    function f = reproj(x)
        f = 0;
        for i = 1:length(homographys)
            src_idx = homographys(i).src_idx;
            dst_idx = homographys(i).dst_idx;
            src_points = double(homographys(i).src_points);
            dst_points = double(homographys(i).dst_points);
            [src_K, src_R] = genKR(x(:, src_idx));
            [dst_K, dst_R] = genKR(x(:, dst_idx));
            
            src_points = [src_points, ones(size(src_points, 1), 1)];
            hat_dst_points = dst_K * dst_R * src_R' * inv(src_K) * src_points';
            hat_dst_points = hat_dst_points ./ (ones(3,1)*hat_dst_points(3,:));
            hat_dst_points = hat_dst_points(1:2, :)';
            residual = dst_points - hat_dst_points;
            norm_residual = sqrt(sum(residual.^2, 2));
            f = f + sum((norm_residual(norm_residual < sigma)).^2) + ...
                sum(2*sigma*norm_residual(norm_residual > sigma) - sigma^2);
        end
    end
end

function [K, R] = genKR(para)
% para = [focal, aspect, cx cy w1 w2 w3];
K = eye(3,3);
K(1,1) = para(1);
K(2,2) = para(2) * para(1);
K(1,3) = para(3);
K(2,3) = para(4);

w = para(5:7);
W = [0 -para(7) para(6); para(7) 0 -para(5); -para(6) para(5) 0];
if w == 0
    R = eye(3);
else
    R = eye(3,3) + (W./norm(w)) * sin(norm(w)) + (W*W)./dot(w,w) * (1-cos(norm(w)));
end
end