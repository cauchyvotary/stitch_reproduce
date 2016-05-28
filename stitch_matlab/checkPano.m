function [canPano] = checkPano(Hs, imgNum)

adjMat = zeros(imgNum, imgNum);
for i = 1:length(Hs)
    adjMat(Hs(i).src_idx, Hs(i).dst_idx) = 1;
    adjMat(Hs(i).dst_idx, Hs(i).src_idx) = 1;
end

% bfs to check all images can be stitched
init_idx = 1;
queue1 = zeros(imgNum, 1);
hp = 1;
ep = 2;
queue1(hp) = init_idx;

marked = zeros(imgNum, 1);
marked(init_idx) = 1;
while (hp < ep)
    cur_idx = queue1(hp);
    hp = hp + 1;
    
    childs = find(adjMat(cur_idx, :));
    for i = 1:length(childs)
        if marked(childs(i)) == 0
            queue1(ep) = childs(i);
            ep = ep + 1;
            marked(childs(i)) = 1;
        end
    end
end

if length(find(marked)) < imgNum
    canPano = false;
    fprintf('\n all images can not be stitched to be single pano\n');
else
    canPano = true;
end