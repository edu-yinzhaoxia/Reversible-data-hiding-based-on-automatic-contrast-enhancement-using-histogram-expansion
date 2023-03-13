function H_x=info_entropies(I)
% 计算图像的熵值
% I=imread('rice.png');
[C,L]=size(I); %求图像的规格
Img_size=C*L; %图像像素点的总个数
G=256; %图像的灰度级
H_x=0;
nk=zeros(G,1);%产生一个G行1列的全零矩阵
for i=1:C
    for j=1:L
        Img_level=I(i,j)+1; %获取图像的灰度级
        nk(Img_level)=nk(Img_level)+1; %统计每个灰度级像素的点数
    end
end
for k=1:G  %循环
    Ps(k)=nk(k)/Img_size; %计算每一个像素点的概率
    if Ps(k)~=0 %如果像素点的概率不为零
        H_x=-Ps(k)*log2(Ps(k))+H_x; %求熵值的公式
    end
end
