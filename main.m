close all
clear
clc
% load kspace data
load('./data/05_t2_tse_tra_512_s33_3mm_29.mat');
[row, col, coils] = size(ksfull);
% load undersampled mask
mask = imread('./data/mask_random_512_512_SR_20_AC_24.png');
mask = repmat(mask', [1 1 size(ksfull,3)]);
mask = logical(mask);
% unsample
ksdata = mask .* ksfull; 
% save and show reference image
ref = sos(ifft2_3D((ksfull)));
figure;
imshow(ref, []);
% estimate the SPIRiT kernel and Lipschitz constant
ACS = 24;
kernel_size = [5 5];
[Ker, Ker_Tra] = Kernel_Estimation(ksdata, kernel_size, ACS);
Lip_C = Lip_Estimation(ksdata, Ker, kernel_size);
% parameter setting
lambda = 0.055;
max_iter = 50;
gamma = 1.99 / Lip_C;
delta = 1 / gamma;
% calculate and show reconstruction image using ADMM 
disp('ADMM Reconstruction')
tic
[res_SPIRiT_ADMM] = L1_Spirit_ADMM(ksdata,mask,2,lambda,1,max_iter,5,Ker,Ker_Tra);
toc
figure;
imshow(res_SPIRiT_ADMM, []);
% calculate and show reconstruction image using ADMM 
disp('PD3O Reconstruction')
tic
[res_SPIRiT_PD3O] = L1_Spirit_PD3O(ksdata,mask,Lev,lambda,max_iter,Ker,Ker_Tra,gamma,delta);
toc
figure;
imshow(res_SPIRiT_PD3O, []);
