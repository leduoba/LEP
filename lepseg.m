function segs = lepseg(imgfn, thrs, ck, tk, aw_c, ew_sf_c, ew_lp_c, ew_p_c)

switch nargin
    case 1
        ck = 128; tk = 256; aw_c = 0.1; ew_sf_c = 8.001; ew_lp_c = 0.01; ew_p_c = 3.501;
        th.thrs = 0.2 + (0 : 199) * 0.1; th.thrn = numel(th.thrs);
    case 2
        ck = 128; tk = 256; aw_c = 0.1; ew_sf_c = 8.001; ew_lp_c = 0.01; ew_p_c = 3.501;
        th.thrs = thrs; th.thrn = numel(thrs);
end
segs = segment(imgfn, ck, tk, aw_c, ew_sf_c, ew_lp_c, ew_p_c, th);

end

function segs = segment(imgfn, ock, otk, aw_c, ew_sf_c, ew_lp_c, ew_p_c, th)

addpath(genpath('./sfedge'));
addpath(genpath('./misc'));
imgpath = './BSDSimages/all/';
thrs = th.thrs; thrn = th.thrn; imgpath = [imgpath, imgfn, '.jpg'];
tic; [h, w, ck, tk, ea, eb, rgb, lbp, ew_lp] = preprocessing(imgpath, ock, otk); toc;
ew_sf = edgeWeight_sf(imgpath); ew = ew_sf_c * ew_sf + ew_lp_c * ew_lp;

amount_v = h * w; amount_e = 4 * w * h - 3 * (w + h) + 2;
aw = aw_c * [0   11.5677   24.0774   24.0774   24.0774   24.0774   27.5677   27.5677   ...
    27.5677   27.5677   27.5677   44.7290   44.7290   44.7290   44.7290   37.1871   29.6387];

load('ep1'); load('ep2'); if h > w, ep = ep1; h0 = 481; w0 = 321; else ep = ep2; h0 = 321; w0 = 481; end;
ind1 = h0 * (w0 - 1); ind2 = (h0 - 1) * w0; ind3 = (h0 - 1) * (w0 - 1);
ep11 = ep(1 : ind1); ep22 = ep(ind1 + 1 : ind1 + ind2);
ep33 = ep(ind1 + ind2 + 1 : ind1 + ind2 + ind3); ep44 = ep(ind1 + ind2 + ind3 + 1 : end);
ep111 = imresize(reshape(ep11, [h0, w0 - 1]), [h, w - 1]);
ep222 = imresize(reshape(ep22, [h0 - 1, w0]), [h - 1, w]);
ep333 = imresize(reshape(ep33, [h0 - 1, w0 - 1]), [h - 1, w - 1]);
ep444 = imresize(reshape(ep44, [h0 - 1, w0 - 1]), [h - 1, w - 1]);
ep111 = ep111'; ep222 = ep222'; ep333 = ep333'; ep444 = ep444';
ep = ew_p_c * [ep111(:); ep222(:); ep333(:); ep444(:)];

tic; segrlt = segment_lep(int32(amount_v), int32(amount_e), int32(ea), int32(eb), double(ew),...
    int32(rgb), int32(ck), int32(lbp), int32(tk), int32(h), int32(w),...
    double(aw), double(ep), double(thrs), int32(thrn));
fprintf('%s in %f seconds\n', imgfn, toc);

segs = genSegCells(segrlt, h, w);

end

function segs = genSegCells(segrlt, h, w)

segn = size(segrlt, 1); segs = cell(1, segn);
for i = 1 : segn
    segs{i} = uint8(mod(reshape(squeeze(segrlt(i, :)), [h, w]) - 1, 255) + 1);
end

end
