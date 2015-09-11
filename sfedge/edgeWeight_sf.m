function ew_sf = edgeWeight_sf(imgpath)
%% set opts for training (see edgesTrain.m)
opts=edgesTrain();                % default options (good settings)
opts.modelDir='models/';          % model will be in models/forest
opts.modelFnm='modelBsds';        % model name
opts.nPos=5e5; opts.nNeg=5e5;     % decrease to speedup training
opts.useParfor=0;                 % parallelize if sufficient memory

%% train edge detector (~20m/8Gb per tree, proportional to nPos/nNeg)
model=edgesTrain(opts); % will load model if already trained

%% set detection parameters (can set after training)
model.opts.multiscale=1;          % for top accuracy set multiscale=1
model.opts.sharpen=2;             % for top speed set sharpen=0
model.opts.nTreesEval=4;          % for top speed set nTreesEval=1
model.opts.nThreads=4;            % max number threads for evaluation
model.opts.nms=1;                 % set to true to enable nms

%% detect edge and visualize results
I = imread(imgpath);
tic, [emap1,emap2,emap3,emap4]=edgesDetect(I,model); toc
emap1 = emap1'; emap2 = emap2'; emap3 = emap3'; emap4 = emap4';
ew_sf = single([emap1(:); emap2(:); emap3(:); emap4(:)])';

end