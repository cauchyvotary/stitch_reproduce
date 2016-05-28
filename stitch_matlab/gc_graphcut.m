function [masks] = gc_graphcut(imgs, masks, corners, aScale)

assert(size(imgs, 3) == 3);
imgs = double(imgs);
masks = logical(masks);
imgNum = size(imgs, 4);

for first = 1:imgNum - 1
    for i = 1:imgNum - first
        second = imgNum - i + 1;
        [masks(:, :, first), masks(:, :, second)] = cut_pair(imgs(:,:,:, first),  imgs(:, :, :, second), ...
                                                            masks(:, :, first), masks(:, :, second), ...
                                                            corners(first, :),  corners(second, :), ...
                                                            aScale); 
    end
    fprintf('% 6dth image completed\n', first);
end

end

function [mask1, mask2] = cut_pair(img1, img2, mask1, mask2, corner1, corner2, aScale) 

BOUND = 360 * aScale;


if (corner1(1)+size(img1, 2) < BOUND)
    if (corner2(1)+size(img2, 2) < BOUND)
        TYPE = 1;
    else
        TYPE = 2;
    end
else
    if (corner2(1)+size(img2, 2) < BOUND)
        TYPE = 3;
    else
        TYPE = 4;
    end
end

if TYPE == 1
    [mask1, mask2] = aux_func(img1, img2, mask1, mask2, corner1, corner2);
elseif TYPE == 2
    tmp_img2 = img2(:, BOUND - (corner2(1)-1) : end, :);
    tmp_mask2 = mask2(:, BOUND - (corner2(1)-1) : end);
    tmp_corner2 = [0, corner2(2)];
    [mask1, tmp_mask21] = aux_func(img1, tmp_img2, mask1, tmp_mask2, corner1, tmp_corner2);
    
    tmp_img2 = img2(:, 1:BOUND - 1 - (corner2(1)-1), :);
    tmp_mask2 = mask2(:, 1:BOUND - 1 - (corner2(1)-1));
    [mask1, tmp_mask22] = aux_func(img1, tmp_img2, mask1, tmp_mask2, corner1, corner2);
    
    mask2 = cat(2, tmp_mask22, tmp_mask21);
    
elseif TYPE == 3
    tmp_img1 = img1(:, BOUND - (corner1(1)-1) : end, :);
    tmp_mask1 = mask1(:, BOUND - (corner1(1)-1) : end);
    tmp_corner1 = [0, corner1(2)];
    [tmp_mask11, mask2] = aux_func(tmp_img1, img2, tmp_mask1, mask2, tmp_corner1, corner2);
    
    tmp_img1 = img1(:, 1:BOUND - 1 - (corner1(1)-1), :);
    tmp_mask1 = mask1(:, 1:BOUND - 1 - (corner1(1)-1));
    [tmp_mask12, mask2] = aux_func(tmp_img1, img2, tmp_mask1, mask2, corner1, corner2);
    
    mask1 = cat(2, tmp_mask12, tmp_mask11);
    
elseif TYPE == 4
    [mask1, mask2] = aux_func(img1, img2, mask1, mask2, corner1, corner2);
end

end
     
function [mask1, mask2] = aux_func(img1, img2, mask1, mask2, corner1, corner2)
    GAP = 1;
    
    off = min([corner1; corner2], [], 1);
    subwidth = max(size(mask1, 2)+corner1(1)-off(1), size(mask2, 2)+corner2(1)-off(1)) + 2*GAP;
    subheight = max(size(mask1, 1)+corner1(2)-off(2), size(mask2, 1)+corner2(2) -off(2)) + 2*GAP;
    
    submask1 = false(subheight, subwidth);
    submask1(corner1(2)-off(2)+1+GAP : corner1(2)-off(2)+size(mask1, 1)+GAP, ...
              corner1(1)-off(1)+1+GAP : corner1(1)-off(1)+size(mask1, 2)+GAP) = mask1;
    submask2 = false(subheight, subwidth);
    submask2(corner2(2)-off(2)+1+GAP : corner2(2)-off(2)+size(mask2, 1)+GAP, ...
               corner2(1)-off(1)+1+GAP : corner2(1)-off(1)+size(mask2, 2)+GAP) = mask2;
           
    subintersect = submask1 & submask2;
    if (isempty(find(subintersect)))
        return;
    end    
    
    subimg1 = zeros(subheight, subwidth, 3);
    subimg1(corner1(2)-off(2)+1+GAP : corner1(2)-off(2)+size(img1, 1)+GAP, ...
             corner1(1)-off(1)+1+GAP : corner1(1)-off(1)+size(img1, 2)+GAP, :) = img1;
    subimg2 = zeros(subheight, subwidth, 3);
    subimg2(corner2(2)-off(2)+1+GAP : corner2(2)-off(2)+size(img2, 1)+GAP, ...
              corner2(1)-off(1)+1+GAP : corner2(1)-off(1)+size(img2, 2)+GAP, :) = img2;
    edgecost = sum((subimg1 - subimg2) .^ 2, 3) + 1;
    
    EI = []; EV = [];
    tmp_idx = find(imerode(subintersect, [0, 1, 1]));
    EI = [EI;
          tmp_idx, tmp_idx + subheight;
          tmp_idx+subheight, tmp_idx];
    EV = [EV;
           edgecost(tmp_idx) + edgecost(tmp_idx + subheight);
           edgecost(tmp_idx) + edgecost(tmp_idx + subheight)];
    tmp_idx = find(imerode(subintersect, [0; 1; 1]));
    EI = [EI;
          tmp_idx, tmp_idx+1;
          tmp_idx+1, tmp_idx];
    EV = [EV;
           edgecost(tmp_idx) + edgecost(tmp_idx + 1);
           edgecost(tmp_idx+1) + edgecost(tmp_idx)];
    A = sparse(EI(:, 1), EI(:, 2), EV, subwidth*subheight, subwidth*subheight, 4 * subwidth*subheight);
    
    suboutline1 = xor(submask1, imerode(submask1, [0 1 0; 1 1 1; 0 1 0]));
    suboutline2 = xor(submask2, imerode(submask2, [0 1 0; 1 1 1; 0 1 0]));
    suboutline = xor(subintersect, imerode(subintersect, [0 1 0; 1 1 1; 0 1 0]));
    TI = []; TV = [];
    tmp_idx = find((submask1 & ~subintersect) | (suboutline & suboutline2));
    TI = [TI;
          tmp_idx, ones(length(tmp_idx), 1)];
    TV = [TV;
           inf(length(tmp_idx), 1)];
     tmp_idx = find((submask2 & ~subintersect) | (suboutline & ~suboutline2));
     TI = [TI;
            tmp_idx, 2*ones(length(tmp_idx), 1)];
     TV = [TV;
            inf(length(tmp_idx), 1)];
     tmp_idx = find(~(submask1 | submask2));
     TI = [TI;
           tmp_idx, 1*ones(length(tmp_idx), 1)];
     TV = [TV;
            inf(length(tmp_idx), 1)];
     T = sparse(TI(:, 1), TI(:, 2), TV, subwidth*subheight, 2, 2*subwidth*subheight);
     
     [flow, labels] = maxflow(A, T);
     
     labels = double(reshape(labels, [subheight, subwidth]));
     labels(~(submask1 | submask2)) = 0.5;
     submask1(labels ~= 0) = false;
     submask2(labels ~= 1) = false;
     mask1 = submask1(corner1(2)-off(2)+1+GAP : corner1(2)-off(2)+size(mask1, 1)+GAP, ...
              corner1(1)-off(1)+1+GAP : corner1(1)-off(1)+size(mask1, 2)+GAP);
     mask2 = submask2(corner2(2)-off(2)+1+GAP : corner2(2)-off(2)+size(mask2, 1)+GAP, ...
               corner2(1)-off(1)+1+GAP : corner2(1)-off(1)+size(mask2, 2)+GAP);
     
end