
%MAIN CODE
clc;
clear all;
close all;
TestPath  = path;
frames = dir(fullfile(TestPath, '*.png'));
c_num = length(frames);

for jj = 1:c_num,
    imgpath = fullfile(TestPath, frames(jj).name);
    I = imread(imgpath);p = 2;q = 1;
            X(jj,:) = feat(I,p,q); %y(jj) = jj;
end

% test part
c=1;
for jj = 1:c_num,
    imgpath = fullfile(TestPath, frames(jj).name);
    I = imnoise(imread(imgpath),'salt & pepper', 0.1);
        for a = .5 : .5 : 1
            for b = .5 :.5 : 1
                y(c,:) = feat(I,a,b); %y(jj) = jj;
                c = c+1;
            end 
        end
        %double for loop runs 4 times and outer loop runs for 16 times.
end
% c becomes 64 here; y is a vector containing 64 elements (as vector of
features)

 for jj = 1 : c_num
    Euc_dist = []; %euc_dist vector made for each of 16 images
    for c = 1 : 64
        temp = 0;
        for i = 1 : 6
           temp = temp + ( X(jj,i) - y(c,i))^2;
        end
        s = sqrt(temp);
        Euc_dist = [Euc_dist s];
    end
    u(jj,:) = Euc_dist;
    [Euc_dist_min , Recognized_index] = min(Euc_dist)
end