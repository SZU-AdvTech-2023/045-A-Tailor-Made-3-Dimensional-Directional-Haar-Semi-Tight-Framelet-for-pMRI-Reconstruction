% 全采样参考图
figure;
imshow(ref, []);
title('全采样参考图');
% 画方框
pos = [125 310 85 85];
rectangle('Position',pos,'EdgeColor','w');
% 欠采样图像
figure;
imshow(un_ref, []);
title('欠采样图像');
% ADMM重建结果
figure;
imshow(res_SPIRiT_ADMM, []);
title('ADMM重建结果');
% PD3O重建结果
figure;
imshow(res_SPIRiT_PD3O, []);
title('PD3O重建结果');
% 全采样参考图局部放大图
part_ref = imresize(imcrop(ref, pos), [256, 256]);
figure;
imshow(part_ref, []);
title('全采样参考图局部放大图');
% 欠采样图像局部放大图
part_un_ref = imresize(imcrop(un_ref, pos), [256, 256]);
figure;
imshow(part_un_ref, []);
title('欠采样图像局部放大图');
% ADMM重建结果局部放大图
part_res_SPIRiT_ADMM = imresize(imcrop(res_SPIRiT_ADMM, pos), [256, 256]);
figure;
imshow(part_res_SPIRiT_ADMM, []);
title('ADMM重建结果局部放大图');
% PD3O重建结果局部放大图
part_res_SPIRiT_PD3O = imresize(imcrop(res_SPIRiT_PD3O, pos), [256, 256]);
figure;
imshow(part_res_SPIRiT_PD3O, []);
title('PD3O重建结果局部放大图');