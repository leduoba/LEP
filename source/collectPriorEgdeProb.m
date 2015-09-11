% Collect the prior edge probabilities on all locations
% There are already two pre-computed files storing the results: ep1.mat and ep2.mat

function [bmap1, bmap2] = collectPriorEgdeProb()

iids = dir('D:/mdlseg/mex-train/groundTruth/*.mat');
bmap1 = zeros(481, 321); bmap2 = zeros(321, 481); cnt = 0;
for i = 1 : numel(iids)
    gt = load(['D:/mdlseg/mex-train/groundTruth/', iids(i).name]); gt = gt.groundTruth;
    cnt = cnt + numel(gt);
    for k = 1 : numel(gt)
        seg = gt{1, k}.Segmentation;
        map1 = double(seg(1 : end - 1, 1 : end - 1) ~= seg(2 : end, 2 : end));
        map2 = double(seg(1 : end - 1, 2 : end) ~= seg(2 : end, 1 : end - 1));
        map3 = double(seg(1 : end - 1, :) ~= seg(2 : end, :));
        map4 = double(seg(:, 1 : end - 1) ~= seg(:, 2 : end));

        map11 = imresize(map1, [481, 321]); map21 = imresize(map2, [481, 321]);
        map31 = imresize(map3, [481, 321]); map41 = imresize(map4, [481, 321]);
        map11 = map11 / max(map11(:)); map21 = map21 / max(map21(:));
        map31 = map31 / max(map31(:)); map41 = map41 / max(map41(:));
        bmap1 = bmap1 + map11 + map21 + map31 + map41;

        map12 = imresize(map1, [321, 481]); map22 = imresize(map2, [321, 481]);
        map32 = imresize(map3, [321, 481]); map42 = imresize(map4, [321, 481]);
        map12 = map12 / max(map12(:)); map22 = map22 / max(map22(:));
        map32 = map32 / max(map32(:)); map42 = map42 / max(map42(:));
        bmap2 = bmap2 + map12 + map22 + map32 + map42;
    end
end
bmap1 = bmap1 / cnt / 4; bmap2 = bmap2 / cnt / 4;

end
