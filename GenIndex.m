
ColorFe = zeros(2,27,600);
ShapeFe = zeros(16,600);
TextureFe = zeros(8,600);
for i=1:600
    clear temp_image;
    clear s;
    s = ['(',num2str(i) ,').jpg'];
    temp_image = imread(s);
    temp_CCV=getCCV(temp_image,1,27);
    ColorFe(1,:,i)=temp_CCV(1,:);
    ColorFe(2,:,i)=temp_CCV(2,:);
    ShapeFe(:,i)=getKM(temp_image,0.5,0.5);
    TextureFe(:,i)=getGLCM(temp_image);
    %
end
save colorMom ColorFe;
save shapeMom ShapeFe;
save textureMom TextureFe;
function CCV = getCCV(img,coherentPrec, numberOfColors)
 if ~exist('coherentPrec','var')
     coherentPrec=1;
 end
  if ~exist('coherentPrec','var')
      numberOfColors = 27;
  end
  CCV=zeros(2,numberOfColors);
  %高斯滤波
  Gaus = fspecial('gaussian',[5 5],2);
  img = imfilter(img,Gaus,'same');
  [img, updNumOfPix]= discretizeColors(img,numberOfColors);  %量化
  imgSize = (size(img,1)*size(img,2));
  thresh = int32((coherentPrec/100) *imgSize);
  parfor i=0:updNumOfPix-1
      BW = img==i;
      CC = bwconncomp(BW);
      compsSize = cellfun(@numel,CC.PixelIdxList);
      incoherent = sum(compsSize(compsSize>=thresh));
      CCV(:,i+1) = [incoherent; ...
          sum(compsSize) - incoherent];
  end
end


function [oneChannel, updatedNumC] = discretizeColors(img,numColors)
    width = size(img,2);
    height = size(img,1);
    oneChannel = zeros(height,width);
    numOfBins = floor(pow2(log2(numColors)/3));
    numOfBinsSQ = numOfBins*numOfBins;
    img = floor((img/(256/numOfBins)));
    for i=1:height
        for j=1:width
            oneChannel(i,j) = img(i,j,1)*numOfBinsSQ ...
                + img(i,j,2)*numOfBins + img(i,j,3);
        end
    end
    updatedNumC = power(numOfBins,3);
end

%Krawtchouk


function [Kr]=wkrchkpoly(N,p)
 pc = 1-p;
 pr = pc/p;
 x=zeros(1,N);
 for i=1:N
     x(1,i)=i-1;
 end
 w = zeros(1,N);
    rho = zeros(1,N);
    K = zeros(N,N);
  %temp=zeros(N,N);
    N=N-1;
    w(1,1)=pc^N;
    rho(1,1) = 1;
    K(1,:)=1;
    K(2,:)=1-x/(N*p);
    
     for i =0:N-2
        w(1,i+2) = w(1,i+1)*(N-i)*p/((i+1)*pc);
        rho(1,i+2) = 1/(-1*pr*(i+1)/(-N+i))*rho(i+1);
for j=1:N+1
        K(i+3,j) = ((N*p+(i+1)*(1-2*p)-x(1,j))*K(i+2,j)-(i+1)*(1-p)*K(i+1,j))/(p*(N-i-1));
end
     end
     w(1,N+1)=w(1,N)*p/((N-1+1)*pc);
    rho(1,N+1)=1/(pr*N)*rho(1,N);
   %temp=zeros(N+1,N+1);
    %for i=1:N+1
   % temp(i,:)=sqrt(rho(1,i))*sqrt(w);
    %end
    temp=sqrt(rho).'*sqrt(w);
    Kr=K.*temp;
 
  
end

function [Q,Kr1,Kr2]=wkrchkmoment_single(img,p1,p2)

width = size(img,2);
    height = size(img,1);
   X=double(rgb2gray(img)); 
    Kr1 = wkrchkpoly(width,p1);

    if (height==width && p1==p2)

        Kr2=Kr1;

    else

        Kr2=wkrchkpoly(height,p2);

    end

    Q = Kr2*X*Kr1';
end

function [K]=getKM(img,p1,p2)
K=zeros(1,16);
Q=wkrchkmoment_single(img,p1,p2);
for i=1:4
    K(1,i)=Q(1,i);
end
K(1,5)=Q(2,1);
K(1,6)=Q(2,2);
K(1,7)=Q(2,3);
K(1,8)=Q(2,4);
K(1,9)=Q(3,1);
K(1,10)=Q(3,2);
K(1,11)=Q(3,3);
K(1,12)=Q(3,4);
K(1,13)=Q(4,1);
K(1,14)=Q(4,2);
K(1,15)=Q(4,3);
K(1,16)=Q(4,4);

end

function [GLCM]=getGLCM(img)

Gray = rgb2gray(img);
[M,N] = size(Gray);
% M = 128;
% N = 128;


%将各颜色分量转化为灰度

% Gray = double(0.3*Image(:,:,1)+0.59*Image(:,:,2)+0.11*Image(:,:,3));



%对原始图像灰度级压缩，将Gray量化成16级

for i = 1:M
    for j = 1:N
        for n = 1:256/16
            if (n-1)*16<=Gray(i,j)&&Gray(i,j)<=(n-1)*16+15
                Gray(i,j) = n-1;
            end
        end
    end
end


%计算四个共生矩阵P,取距离为1

P = zeros(16,16,4);
for m = 1:16
    for n = 1:16
        for i = 1:M
            for j = 1:N
                if j<N&&Gray(i,j)==m-1&&Gray(i,j+1)==n-1
                    P(m,n,1) = P(m,n,1)+1;
                    P(n,m,1) = P(m,n,1);
                end
                if i>1&&j<N&&Gray(i,j)==m-1&&Gray(i-1,j+1)==n-1
                    P(m,n,2) = P(m,n,2)+1;
                    P(n,m,2) = P(m,n,2);
                end
                if i<M&&Gray(i,j)==m-1&&Gray(i+1,j)==n-1
                    P(m,n,3) = P(m,n,3)+1;
                    P(n,m,3) = P(m,n,3);
                end
                if i<M&&j<N&&Gray(i,j)==m-1&&Gray(i+1,j+1)==n-1
                    P(m,n,4) = P(m,n,4)+1;
                    P(n,m,4) = P(m,n,4);
                end
            end
        end
        if m==n
            P(m,n,:) = P(m,n,:)*2;
        end
    end
end



% 共生矩阵归一化

for n = 1:4
    P(:,:,n) = P(:,:,n)/sum(sum(P(:,:,n)));
end



%计算能量、熵、对比度、相关性4个纹理参数

H = zeros(1,4);
I = H;
Ux = H;      Uy = H;
deltaX= H;  deltaY = H;
C =H;
E=H;
for n = 1:4
    E(n) = sum(sum(P(:,:,n).^2)); %%能量
    for i = 1:16
        for j = 1:16
            if P(i,j,n)~=0
                H(n) = -P(i,j,n)*log(P(i,j,n))+H(n); %%熵
            end
            I(n) = (i-j)^2*P(i,j,n)+I(n);  %%对比度
          
            Ux(n) = i*P(i,j,n)+Ux(n); %相关性中μx
            Uy(n) = j*P(i,j,n)+Uy(n); %相关性中μy
        end
    end
end
for n = 1:4
    for i = 1:16
        for j = 1:16
            deltaX(n) = (i-Ux(n))^2*P(i,j,n)+deltaX(n); %相关性中σx
            deltaY(n) = (j-Uy(n))^2*P(i,j,n)+deltaY(n); %相关性中σy
            C(n) = i*j*P(i,j,n)+C(n);            
        end
    end
    C(n) = (C(n)-Ux(n)*Uy(n))/deltaX(n)/deltaY(n); %相关性  
end



%能量、熵、对比度、相关的均值和标准差作为最终8维纹理特征

GLCM(1) = mean(E); GLCM(2) = sqrt(cov(E));
GLCM(3) = mean(H); GLCM(4) = sqrt(cov(H));
GLCM(5) = mean(I); GLCM(6) = sqrt(cov(I));
GLCM(7) = mean(C); GLCM(8) = sqrt(cov(C));

end
