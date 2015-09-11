function [emap1,emap2,emap3,emap4] = edgesDetect( I, model )
% Detect edges in image.
%
% For an introductory tutorial please see edgesDemo.m.
%
% The following model params may be altered prior to detecting edges:
%  prm = stride, sharpen, multiscale, nTreesEval, nThreads, nms
% Simply alter model.opts.prm. For example, set model.opts.nms=1 to enable
% non-maximum suppression. See edgesTrain for parameter details.
%
% USAGE
%  [E,O,inds,segs] = edgesDetect( I, model )
%
% INPUTS
%  I          - [h x w x 3] color input image
%  model      - structured edge model trained with edgesTrain
%
% OUTPUTS
%  E          - [h x w] edge probability map
%  O          - [h x w] coarse edge normal orientation (0=left, pi/2=up)
%  inds       - [h/s x w/s x nTreesEval] leaf node indices
%  segs       - [g x g x h/s x w/s x nTreesEval] local segmentations
%
% EXAMPLE
%
% See also edgesDemo, edgesTrain, edgesChns
%
% Structured Edge Detection Toolbox      Version 3.0
% Copyright 2014 Piotr Dollar.  [pdollar-at-microsoft.com]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the MSR-LA Full Rights License [see license.txt]

% get parameters
opts=model.opts; opts.nTreesEval=min(opts.nTreesEval,opts.nTrees);
if(~isfield(opts,'sharpen')), opts.sharpen=0; end
if(~isfield(model,'segs')), model.segs=[]; model.nSegs=[]; end
opts.stride=max(opts.stride,opts.shrink); model.opts=opts;

if( opts.multiscale )
  % if multiscale run edgesDetect multiple times
  ss=2.^(-1:1); k=length(ss); siz=size(I);
  model.opts.multiscale=0; model.opts.nms=0;
  emap1 = 0; emap2 = 0; emap3 = 0; emap4 = 0;
  for i=1:k, s=ss(i); I1=imResample(I,s);
    [emap11,emap21,emap31,emap41]=edgesDetect(I1,model);
    % 分别得到不同尺度上的edge map并且取平均
    emap1=emap1+imResample(emap11,[siz(1),siz(2)-1]);
    emap2=emap2+imResample(emap21,[siz(1)-1,siz(2)]);
    emap3=emap3+imResample(emap31,[siz(1)-1,siz(2)-1]);
    emap4=emap4+imResample(emap41,[siz(1)-1,siz(2)-1]);
    % 由于随机森林是每隔1个像素计算，因此seg map和edge map的尺寸都是当前处理图像I1尺寸的一半
    % segs是当前的16x16的patch的local segmentation labels
    % 对local segs，可以得到四个方向上的edge map，对应八连通的4个map（当然也可以仅取四连通的2个map）
    % 每一对相邻像素，将会在16x16xT/(2x2)个edge maps上出现
    % 对于每个方向上的edge map，计算方式如下：
    % (1) 对于各个scale，由local segs在此方向上的edge maps（约patch大小），求出平均的edge map（scaled图像大小）
    % (2) 由各个scale上的edge map，缩放到原始scale
    % (3) 对各个scale上的edge map求平均做为最终输出
  end;
  emap1=emap1/k;emap2=emap2/k;emap3=emap3/k;emap4=emap4/k;model.opts=opts;
  
else
  % pad image, making divisible by 4
  siz=size(I); r=opts.imWidth/2; p=[r r r r];
  p([2 4])=p([2 4])+mod(4-mod(siz(1:2)+2*r,4),4);
  % 增量: 上p(1)下p(2)左p(3)右p(4)
  % 在这里是上下左右各为：16、19、16、19
  % 也就是说width和height各增加35
  I = imPad(I,p,'symmetric');
  
  % compute features and apply forest to image
  [chnsReg,chnsSim] = edgesChns( I, opts );
  s=opts.sharpen; if(s), I=convTri(rgbConvert(I,'rgb'),1); end
  [emap1,emap2,emap3,emap4,nmap1,nmap2,nmap3,nmap4] = edgesDetectMex(model,I,chnsReg,chnsSim);
  
  % normalize and finalize edge maps
  t=1.66; r=opts.gtWidth/2;
  emap1=emap1(1+r:siz(1)+r,1+r:siz(2)+r-1,:);
  emap2=emap2(1+r:siz(1)+r-1,1+r:siz(2)+r,:);
  emap3=emap3(1+r:siz(1)+r-1,1+r:siz(2)+r-1,:);
  emap4=emap4(1+r:siz(1)+r-1,1+r+1:siz(2)+r,:);
  nmap1=nmap1(1+r:siz(1)+r,1+r:siz(2)+r-1,:);
  nmap2=nmap2(1+r:siz(1)+r-1,1+r:siz(2)+r,:);
  nmap3=nmap3(1+r:siz(1)+r-1,1+r:siz(2)+r-1,:);
  nmap4=nmap4(1+r:siz(1)+r-1,1+r+1:siz(2)+r,:);
  emap1 = emap1 ./ nmap1 * t; emap2 = emap2 ./ nmap2 * t;
  emap3 = emap3 ./ nmap3 * t; emap4 = emap4 ./ nmap4 * t;
  emap1=convTri(emap1,1); emap2=convTri(emap2,1); emap3=convTri(emap3,1); emap4=convTri(emap4,1);
end

% perform nms
if( opts.nms>0 ), emap1=edgesNmsMex(emap1,emap1*0,1,5,1.01,opts.nThreads); end
if( opts.nms>0 ), emap2=edgesNmsMex(emap2,emap2*0+pi*0.5,1,5,1.01,opts.nThreads); end
if( opts.nms>0 ), emap3=edgesNmsMex(emap3,emap3*0+pi*0.75,1,5,1.01,opts.nThreads); end
if( opts.nms>0 ), emap4=edgesNmsMex(emap4,emap4*0+pi*0.25,1,5,1.01,opts.nThreads); end

end
