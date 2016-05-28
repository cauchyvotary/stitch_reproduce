function [maxTree, treeCenter] = maxSpanningTree(edges, vNum)

disjoints  = [[1:vNum]', [1:vNum]', ones(vNum, 1)];
edge_weights = [];
for i = 1:length(edges)
    edge_weights(i) = edges(i).weight;
end

[~, idxs] = sort(edge_weights, 'descend');
edges = edges(idxs);

tree_edge_cnt = 0;
for i = 1:length(edges)
    src_idx = edges(i).src_idx;
    dst_idx = edges(i).dst_idx;
    
    set1 = findSet(disjoints, src_idx);
    set2 = findSet(disjoints, dst_idx);
    if set1 ~= set2
        tree_edge_cnt = tree_edge_cnt + 1;
        treeEdges(tree_edge_cnt) = edges(i);
        
        %merge set
        if disjoints(set1, 3) > disjoints(set2, 3)
            disjoints(set2, 2) = set1;
            disjoints(set1, 3) = disjoints(set1, 3) + disjoints(set2, 3);
        else
            disjoints(set1, 2) = set2;
            disjoints(set2, 3) = disjoints(set2, 3) + disjoints(set1, 3);
        end
    end
end

assert(length(treeEdges) == vNum -1 )

% constuct a easy-to-use tree, of minimal height
vDegs = zeros(vNum, 1);
adjMat = zeros(vNum, vNum);
for i = 1:length(treeEdges)
    idx1 = treeEdges(i).src_idx;
    idx2 = treeEdges(i).dst_idx;
    adjMat(idx1, idx2) = 1;
    adjMat(idx2, idx1) = 1;
    vDegs(idx1) = vDegs(idx1) + 1;
    vDegs(idx2) = vDegs(idx2) + 1;
end

leafs = find(vDegs == 1);
dist2others = zeros(vNum, 1);
for i = 1:length(leafs)
    idx = leafs(i);
    dist2others = bfsGetDist(adjMat, idx, dist2others)
end
[~, treeCenter] = min(dist2others);
maxTree = genTree(adjMat, treeCenter);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [maxTree] = genTree(adjMat, init_v)

vNum = size(adjMat, 1);
queue1 = zeros(vNum, 1);
queue1(1) = init_v;
adjMat(:, init_v) = 0;
hp = 1;
ep = 2;

for i = 1:vNum
    maxTree(i).childs = [];
end

while (hp < ep)
    v = queue1(hp);
    hp = hp + 1;
   
    childs = find(adjMat(v, :));
    if (length(childs) > 0)
        maxTree(v).childs = childs;
    end
    for i = 1:length(childs)
        child_v = childs(i);
        queue1(ep) = child_v;
        ep = ep + 1;
        adjMat(:, child_v) = 0;
    end
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dist2others] = bfsGetDist(adjMat, init_v, dist2others)
vNum = length(dist2others);
queue1 = zeros(vNum, 2);
hp = 1;
ep = 2;
queue1(1, 1) = init_v;
queue1(1, 2) = 0;

marked = zeros(vNum, 1);
marked(init_v) = 1;
while (hp < ep)
    v = queue1(hp, 1);
    vDist = queue1(hp, 2);
    hp = hp + 1;
    if dist2others(v) < vDist
        dist2others(v) = vDist;
    end
    
    childs = find(adjMat(v, :));
    for i = 1:length(childs)
        child_v = childs(i);
        if marked(child_v) == 0
            marked(child_v) = 1;
            queue1(ep, 1) = child_v;
            queue1(ep, 2) = vDist + 1;
            ep = ep + 1;
            
            if dist2others(child_v) < vDist + 1
                dist2others(child_v) = vDist + 1;
            end
        end
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [set_idx] = findSet(disjoints, idx)
while (idx  ~= disjoints(idx, 2))
    idx = disjoints(idx, 2);
end
set_idx = idx;

