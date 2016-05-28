function [H, suc_flag] = homography4points(points1, points2)
% calculate H from coord 2 to coord 1

A = [];
y = [];

for i = 1:4
    x2 = points2(i, 1); y2 = points2(i, 2);
    x1 = points1(i, 1); y1 = points1(i, 2);
    row1 = [x2, 0, -x2*x1, y2, 0, -y2*x1, 1, 0];
    row2 = [0, x2, -x2*y1, 0, y2, -y2*y1, 0, 1];
    A = [A;row1;row2];
    y = [y; x1; y1];
end

if (rcond(A) < 1e-9)
    H = [];
    suc_flag = false;
else
    H = reshape([A\y; 1], 3, 3); 
    suc_flag = true;
end

