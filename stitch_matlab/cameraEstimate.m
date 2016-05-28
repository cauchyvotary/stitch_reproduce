function [cameras] = cameraEstimate(Hs, imgNum, imgSize)
% Before call this function, call CHECKPANO to check all the images
% can be stitched into single panorama.

% step1 : estimate focals
cnt = 0;
for i = 1:length(Hs)
    [f1, f1ok, f2, f2ok] = focalEstimate(Hs(i).H);
    if (f1ok & f2ok)
        cnt = cnt + 1;
        eFocals(cnt) = sqrt(f1*f2);
    end
end
init_focal = median(sort(eFocals, 'ascend'));
%prepare cameras
% \lmabda \tilde x = K * R * X_w
for i = 1:imgNum
    cameras(i).focal = init_focal;
    cameras(i).aspect = 1;
    cameras(i).cx = 0;
    cameras(i).cy = 0;
    cameras(i).R = eye(3);
    cameras(i).t = zeros(3, 1);
end

% step2 : generate max spaning tree
for i = 1:length(Hs)
    edge.src_idx = Hs(i).src_idx;
    edge.dst_idx = Hs(i).dst_idx;
    edge.weight = Hs(i).inliersNum;
    edges(i) = edge;
end

[tree, center] = maxSpanningTree(edges, imgNum);
queue1 = zeros(imgNum, 1);
queue1(1) = center;
hp = 1; ep = 2;
while (hp < ep)
    imgIdx = queue1(hp);
    hp = hp + 1;
    
    childs = tree(imgIdx).childs;
    for i = 1:length(childs)
        childIdx = childs(i);
        queue1(ep) = childIdx;
        ep = ep + 1;
        
        % find out homography form childIdx to curIdx
        c2p_H = child2parentH(childIdx, imgIdx, Hs);
        cameras(childIdx).R = inv(genK(cameras(childIdx))) * ...
            inv(c2p_H) * genK(cameras(imgIdx)) * cameras(imgIdx).R;
    end
end

for i = 1:length(cameras)
    cameras(i).cx = imgSize(2) / 2;
    cameras(i).cy = imgSize(1) / 2;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K] = genK(camera)
K = zeros(3,3);
K(1,1) = camera.focal;
K(2,2) = camera.focal * camera.aspect;
K(1,3) = camera.cx;
K(2,3) = camera.cy;
K(3,3) = 1;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c2p_H] = child2parentH(childIdx, parentIdx, Hs)
c2p_H = [];
for i = 1:length(Hs)
    if ((Hs(i).src_idx == childIdx) && (Hs(i).dst_idx == parentIdx) )
        c2p_H = Hs(i).H;
    else
        if ((Hs(i).src_idx == parentIdx) && (Hs(i).dst_idx == childIdx))
            c2p_H = inv(Hs(i).H);
            c2p_H = c2p_H ./ c2p_H(3,3);
        end
    end
end