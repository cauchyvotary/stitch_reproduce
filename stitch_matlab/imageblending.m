function [canvas, canvas_corner, mask_canvas] = imageblending(warped_imgs, warped_masks, corners, aScale)

canvas_corner = min(corners, [], 1);
canvas_corner(1) = 0;

warped_height = size(warped_imgs, 1);
warped_width = size(warped_imgs, 2);
imgNum = size(warped_imgs, 4);

canvas = zeros(max(corners(:, 2)) - min(corners(:, 2)) + warped_height, ... 
                 max(corners(:, 1)) - 0+ warped_width, ...
                 3);
             
% canvas_weight = zeros(size(canvas, 1), size(canvas, 2));
mask_canvas = zeros(size(canvas, 1), size(canvas, 2));

for i = 1:imgNum
    tl_x = corners(i, 1) - canvas_corner(1) + 1;
    tl_y = corners(i, 2) - canvas_corner(2) + 1;
    br_x = tl_x + warped_width - 1;
    br_y = tl_y + warped_height - 1;

%    canvas_weight(tl_y:br_y, tl_x:br_x) = canvas_weight(tl_y:br_y, tl_x:br_x) ...
%        + warped_masks(:,:, i);
    mask_canvas(tl_y:br_y, tl_x:br_x) = mask_canvas(tl_y:br_y, tl_x:br_x) + warped_masks(:,:,i) .* i;
    for j = 1:3
        canvas (tl_y:br_y, tl_x:br_x, j) = canvas(tl_y:br_y, tl_x:br_x, j) +  ...
                        warped_imgs(:,:,j,i).* warped_masks(:,:,i);
    end
end

% update 2016/3/2, round >360degree to the head
% tmp_img = canvas(:, 360*aScale - canvas_corner(1) + 1:end, :);
% canvas(:, 1:size(tmp_img, 2), :) = canvas(:, 1:size(tmp_img, 2), :) + tmp_img;

end