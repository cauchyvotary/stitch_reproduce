function [aScale, rScale, warped_imgs, warped_masks, corners] = cylindricalWarping(cameras, imgs, aScale)
    imgNum = length(cameras);
    
    focals = zeros(imgNum, 1);
    for i = 1:imgNum
        focals(i) = cameras(i).focal;
    end
    rScale = median(focals);
    
    if (nargin == 2)
        aScale = 10;
    end
    
    width = size(imgs, 2);
    height = size(imgs, 1);
    
%     % add some magic effect
%     rule_xy = [[1:width]', ones(width, 1);
%                [1:width]', 2*ones(width, 1);
%                [1:width]', 3*ones(width, 1);
%                [1:width]', (height-2)*ones(width, 1);
%                [1:width]', (height-1)*ones(width, 1);
%                [1:width]', height*ones(width, 1);
%                1*ones(height,1), (1:height)';
%                2*ones(height,1), (1:height)';
%                3*ones(height,1), (1:height)';
%                (width - 2)*ones(height, 1), (1:height)';
%                (width - 1)*ones(height, 1), (1:height)';
%                width*ones(height, 1), (1:height)'];
%     edge_idx1 = [rule_xy* [height; 1] - height;
%                 rule_xy* [height; 1] - height + (height*width) ;
%                 rule_xy* [height; 1] - height + 2*(height*width)];
%     edge_idx2 = [rule_xy* [height; 1] - height + 2*(height*width)];
%     for i = 1:imgNum
%         tmp = imgs(:,:,:,i);
%         tmp(edge_idx1) = 0;
%         tmp(edge_idx2) = 255;
%         imgs(:,:,:,i) = tmp;
%     end
%     % magic effect added
    
    U = [[1:width]', ones(width, 1); 
         [1:width]', height*ones(width, 1);
        ones(height-2,1), (2:height-1)';
        width*ones(height-2, 1), (2:height-1)'];    
    XDatas = zeros(imgNum, 2);
    YDatas = zeros(imgNum, 2);
    tdata.rScale = rScale;
    tdata.aScale = aScale;
    for i = 1:imgNum
        tdata.K = genK(cameras(i));
        tdata.R = cameras(i).R;
        T.tdata = tdata;
        [X] = cylindricalWarpIMPLE(U, T);
        
        XDatas(i, 1) = floor(min(X(:, 1)));
        XDatas(i, 2) = ceil(max(X(:, 1)));
        YDatas(i, 1) = floor(min(X(:, 2)));
        YDatas(i, 2) = ceil(max(X(:, 2)));
    end
   
    corners = [XDatas(:, 1), YDatas(:,1)];
    warped_width = max(XDatas(:, 2) - XDatas(:, 1)) + 1;
    warped_height = max(YDatas(:, 2) - YDatas(:, 1)) + 1;
    warped_imgs = zeros(warped_height, warped_width, 3, imgNum);
    warped_masks = zeros(warped_height, warped_width, imgNum);
    
    tdata.rScale = rScale;
    tdata.aScale = aScale;
    mask = ones(size(imgs, 1), size(imgs,2));
    for i = 1:imgNum
        tdata.K = genK(cameras(i));
        tdata.R = cameras(i).R;
        tf = maketform('custom', 2, 2, @cylindricalWarpIMPLE, ...
                        @invCylindricalWarpIMPLE, tdata);
        [tmp_img, XData, YData] = imtransform(imgs(:,:,:,i), tf,... 
                                                'XData', XDatas(i,:),...
                                                'YData', YDatas(i,:));
        warped_imgs(1:size(tmp_img,1),1:size(tmp_img, 2),:,i) = tmp_img;
        tmp_img = imtransform(mask, tf, ...
            'XData', XDatas(i, :), 'YData', YDatas(i, :));
        warped_masks(1:size(tmp_img, 1),1:size(tmp_img, 2), i) = tmp_img;
    end
end

function [Uw] = cylindricalWarpIMPLE(Xc, T)
    rScale = T.tdata.rScale;
    aScale = T.tdata.aScale;
    
    K = T.tdata.K;
    R = T.tdata.R;
    
    tildeXc = [Xc, ones(size(Xc, 1), 1)];
    tildeXw = (R' * inv(K) * tildeXc' )' ;
    
    rads = acos(tildeXw(:, 1) ./ sqrt(tildeXw(:, 1).^2 + tildeXw(:, 3).^2));
    bt90_mask = ( tildeXw(:, 3) < 0 );
    rads(bt90_mask) = 2*pi - rads(bt90_mask);
    
    % reverse the direction (from anti-clock to clockwise)
    rads = 2 * pi - rads;
    
    % make mapping smooth, i.e. if the angles pass 360-0, then the angles
    % greater than 360 will not restart from 0
    if (max(rads) - min(rads) > pi)
        pass0_mask = (rads < pi);
        rads(pass0_mask) = rads(pass0_mask) + 2*pi;
    end
    
    % Uw(:, 1) -- clockwise degree, can be bigger than 360 * aScale;
    % Uw(:, 2) -- h;
    Uw(:, 1) = rads .* (180 / pi) .* aScale;
    Uw(:, 2) = tildeXw(:, 2) ./ sqrt(tildeXw(:, 1).^2 + tildeXw(:, 3).^2) .* rScale;
    
end

function [Xc] = invCylindricalWarpIMPLE(Uw, T)
    
    % Uw(:, 1) -- clockwise degree, can be bigger than 360 * aScale
    % Uw(:, 2) -- h
    rScale = T.tdata.rScale;
    aScale = T.tdata.aScale;
    K = T.tdata.K
    R = T.tdata.R;
    
    % the range of the rads dont matters
    anti_rads = - Uw(:, 1) .* pi ./ (aScale * 180) ;
    normal_h = Uw(:, 2) ./ rScale;
    tildeXw = [cos(anti_rads), normal_h, sin(anti_rads)];
    tildeXc = (K * R * tildeXw')' ;
    
    assert(isempty(find(tildeXc(:, 3) <= 0) ));
    
    Xc = tildeXc(:, 1:2) ./ repmat(tildeXc(:, 3), [1, 2]);
end

% function [X] = cylindricalWarpIMPLI(U, T)
%     rScale = T.tdata.rScale;
%     aScale = T.tdata.aScale;
%     camera = T.tdata.camera;
%     K = K(camera);
%     R = camera.R;
%     
%     tilde_U = [U, ones(size(U, 1), 1)];
%     tilde_W = (R' * inv(K) * tilde_U')';
%     
%     rads = acos(tilde_W(:,1) ./ sqrt(tilde_W(:, 1).^2 + tilde_W(:, 3).^2));
%     rads(tilde_W(:, 3) < 0) = 2*pi - rads(tilde_W(:, 3) < 0);
%     X(:, 1) = rads * (180 / pi) * aScale;
%     if (max(X(:, 1)) - min(X(:, 1)) > 90 * aScale)
%         tmp_idx = X(:, 1) < 180 * aScale;
%         X(tmp_idx, 1) = X(tmp_idx, 1) + 360 * aScale;
%     end
%     
%     X(:, 2) = tilde_W(:, 2) ./ sqrt(tilde_W(:, 1).^2 + tilde_W(:, 3).^2) * rScale;
% end
% 
% function [U] = invCylindricalWarpIMPLI(X, T)
%     rScale = T.tdata.rScale;
%     aScale = T.tdata.aScale;
%     camera = T.tdata.camera;
%     K = K(camera);
%     R = camera.R;
%     
%     rads = X(:, 1) / aScale * (pi/180);
%     tilde_W = zeros(length(rads), 3);
%     tilde_W(:, 1) = cos(rads);
%     tilde_W(:, 3) = sin(rads);
%     tilde_W(:, 2) = X(:, 2) ./ rScale;
%     tilde_U =  (K * R * tilde_W')';
%     
%     % assert(isempty(find(tilde_U(:,3) <= 0)));
%     U = tilde_U(:, 1:2) ./ repmat(tilde_U(:, 3), [1, 2]);
%     
% end

function [K] = genK(camera)
    K = eye(3);
    K(1,1) = camera.focal;
    K(2,2) = camera.focal * camera.aspect;
    K(1,3) = camera.cx;
    K(2,3) = camera.cy;
end