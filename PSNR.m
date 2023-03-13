function psnr =PSNR(I1,I2)
 
% I1与I2之间的PSNR
[h,w]=size(I1);
B=8; % 编码一个像素用多少二进制位
MAX=2^B-1;% 图像有多少灰度级
MES=sum(sum((I1-I2).^2))/(h*w);% 均方差
psnr=20*log10(MAX/sqrt(MES));% 峰值信噪比
end