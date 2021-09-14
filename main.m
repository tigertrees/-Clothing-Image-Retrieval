
C = load('colorMom.mat');
S = load('shpMom.mat');
T = load('textureMom.mat');
ColorData = C.ColorFe;
ShapeData = S.ShapeFe;
TextureData = T.TextureFe;
order = zeros(1,600);

RetrieveImg = 170;    %´ı¼ìË÷Í¼ÏñµÄ±àºÅ
sum=0;
for i=1:600
    order(1,i)= CalSpDis(ShapeData(:,RetrieveImg),ShapeData(:,i));
    %ÈÚºÏ0.4*0.01*CalSpDis(ShapeData(:,RetrieveImg),ShapeData(:,i))+0.4*0.00001*CalColDis(ColorData(:,:,RetrieveImg),ColorData(:,:,i))+0.2*CalTxtDis(TextureData(:,RetrieveImg),TextureData(:,i));
    %ĞÎ×´CalSpDis(ShapeData(:,RetrieveImg),ShapeData(:,i));
    %ÑÕÉ«CalColDis(ColorData(:,:,RetrieveImg),ColorData(:,:,i));
    %ÎÆÀíCalTxtDis(TextureData(:,RetrieveImg),TextureData(:,i));
end
[sA,index] = sort(order) ;

for i=1:20
    I=imread(['(',num2str(index(i)) ,').jpg']);
    subplot(4,5,i);
    imshow(I);
    if index(i)/100>1&&index(i)/100<=2
        sum=sum+1;
    end
end

sum/20
function  Dis  = CalColDis( arg1,arg2 )

    Dis=0;
    for i=1:27

    C=abs(arg1(1,i)-arg2(1,i));
N=abs(arg1(2,i)-arg2(2,i));
    d=C+N;
    Dis=Dis+d;
    end
end

function  Dis  = CalSpDis( arg1,arg2 )
 Dis=0;

    for i=1:1152

    d1=abs(arg1(i,1)-arg2(i,1));

    Dis=Dis+d1;
    end
 
end

function Dis = CalTxtDis(arg1,arg2)
d=0;

    for i=1:8

    d1=(arg1(i,1)-arg2(i,1))^2;

    d=d+d1;
    end
    Dis=d^0.5;
end