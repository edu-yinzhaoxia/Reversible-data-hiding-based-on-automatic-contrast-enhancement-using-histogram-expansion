function H_x=info_entropies(I)
% ����ͼ�����ֵ
% I=imread('rice.png');
[C,L]=size(I); %��ͼ��Ĺ��
Img_size=C*L; %ͼ�����ص���ܸ���
G=256; %ͼ��ĻҶȼ�
H_x=0;
nk=zeros(G,1);%����һ��G��1�е�ȫ�����
for i=1:C
    for j=1:L
        Img_level=I(i,j)+1; %��ȡͼ��ĻҶȼ�
        nk(Img_level)=nk(Img_level)+1; %ͳ��ÿ���Ҷȼ����صĵ���
    end
end
for k=1:G  %ѭ��
    Ps(k)=nk(k)/Img_size; %����ÿһ�����ص�ĸ���
    if Ps(k)~=0 %������ص�ĸ��ʲ�Ϊ��
        H_x=-Ps(k)*log2(Ps(k))+H_x; %����ֵ�Ĺ�ʽ
    end
end
