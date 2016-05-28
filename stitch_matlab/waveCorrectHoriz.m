function [cameras] = waveCorrectHoriz(cameras)
% current cameras(i).R is rotation matrix from world coordinate to camera
% reference, the first step is extract all rotation matrix from camera
% reference to world coordinate
imgNum = length(cameras);
rmats = zeros(3, 3, imgNum);
for i = 1:imgNum
    rmats(:,:, i) = cameras(i).R';
end

X = squeeze(rmats(:, 1, :));
moment = X*X';
[U D] = eigs(moment);

% the eigen vector corresponding to the smallest eigen value
new_Y = U(:, 3);

Z = sum(squeeze(rmats(:, 3, :)), 2);
new_X = cross(new_Y, Z);
if norm(new_X) < 1e-4
    % no need to do correct
    return
end

new_X = new_X ./ norm(new_X);
new_Z = cross(new_X, new_Y);

conf = sum(new_X' * squeeze(rmats(:,1,:)));
if conf < 0
    new_X = -new_X;
    new_Y = -new_Y;
end

% [new_X, new_Y, new_Z] new reference in world coordinate,
R_old2new = [new_X, new_Y, new_Z]';
for i = 1:length(cameras)
    cameras(i).R = (R_old2new * rmats(:,:,i))';
end


end