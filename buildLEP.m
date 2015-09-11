function buildLEP()

oldpath = pwd; cd ./source

%% compiling
fprintf('Compiling......\n');
mex segment_lep.cpp -O CXXFLAGS="\$CXXFLAGS -I./ -I/usr/include/opencv -I/usr/include/opencv2 -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -lopencv_calib3d -lopencv_contrib -lopencv_core -lopencv_features2d -lopencv_flann -lopencv_gpu -lopencv_highgui -lopencv_imgproc -lopencv_legacy -lopencv_ml -lopencv_objdetect -lopencv_photo -lopencv_stitching -lopencv_ts -lopencv_video -lopencv_videostab
mex preprocessing.cpp -O CXXFLAGS="\$CXXFLAGS -I./ -I/usr/include/opencv -I/usr/include/opencv2 -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -lopencv_calib3d -lopencv_contrib -lopencv_core -lopencv_features2d -lopencv_flann -lopencv_gpu -lopencv_highgui -lopencv_imgproc -lopencv_legacy -lopencv_ml -lopencv_objdetect -lopencv_photo -lopencv_stitching -lopencv_ts -lopencv_video -lopencv_videostab
fprintf('LEP is built successfully. Please enjoy it.\n');
%% move mex files
movefile('segment_lep.mexa64', '../');
movefile('preprocessing.mexa64', '../');

cd(oldpath);

end
