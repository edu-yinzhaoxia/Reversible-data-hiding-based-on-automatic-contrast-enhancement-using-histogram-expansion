clear;clc;
N =24;
cover_dir = '..\bossbase620';
% i_img = 1;
for i_img = 1:N
        cover_Path = fullfile([cover_dir,'/',num2str(i_img),'.png']);

        save i_img.mat
        cover = double(imread(cover_Path));
        [m,n] = size(cover);

rng(1) 
payload_length=200000; 
payload=randi([0,1],payload_length,1);
iteration_max =30;

Y = reshape(cover,m*n,1)
[rdh_image,~,~,~,embedding_capacity]=mbp(Y,payload,iteration_max,m,n);
rdh_Y = reshape(rdh_image,m,n);
% rdh_Y = uint8(rdh_Y);

 
% subplot(121),imshow(cover);
% figure,imshow(uint8(cover));

rdh_rce  = RCE(rdh_Y,cover);
rdh_psnr = PSNR(rdh_Y,cover);
rdh_rmbe = RMBE(rdh_Y,cover);
rdh_ssim = SSIM(rdh_Y,cover);
% rdh_brisque =brisque(rdh_Y);
% (100-brisquescore(rdh_Y))/100; brisque
rdh_ree = REE(rdh_Y,cover);
% save rdh_ssim.mat

%avg
Y_rce(i_img,1) = rdh_rce;
Y_psnr(i_img,1) = rdh_psnr;
Y_rmbe(i_img,1) = rdh_rmbe;
Y_ssim(i_img,1) = rdh_ssim;
% Y_brisque(i_img,1) =rdh_brisque;
Y_ree(i_img,1) = rdh_ree;

end
    Y_rceavg= mean(Y_rce);
    Y_psnravg = mean(Y_psnr);
    Y_rmbeavg = mean(Y_rmbe);
    Y_ssimavg = mean(Y_ssim);
    % Y_brisqueavg = mean(Y_brisque);
    Y_reeavg = mean(Y_ree);
	
    save Y_reeavg.mat


