clc
clear 
close
img=imread('A.bmp');

%img=imresize(img,[20 20],'nearest'); %binary images-ro make faster for
%slower script

img= double(img)/255; %for binary
imshow(img); 
inLegendre_Cntral_Matrix_c(img);
%inLegendre_Cntral_Matrix(img); %slower