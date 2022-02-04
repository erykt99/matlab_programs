close all;
clear
tic
%% VAR

I=imread('bmpf.bmp');
I=rgb2ycbcr(I);
[width,height,x]=size(I); 
BLOCKSIZE=8;
blockR=zeros(8,8,4096);
blockG=zeros(8,8,4096);
blockB=zeros(8,8,4096);
y=zeros(8,8,4096);
cb=zeros(8,8,4096);
cr=zeros(8,8,4096);
y2=zeros(8,8,4096);
cb2=zeros(8,8,4096);
cr2=zeros(8,8,4096);
tab_kwantyzacji = [16 11 10 16 24 40 51 61;
                   12 12 14 19 26 58 60 55;
                   14 13 16 24 40 57 69 56;
                   14 17 22 29 51 87 80 62;
                   18 22 37 56 68 109 103 77;
                   24 35 55 64 81 104 113 92;
                   49 64 78 87 103 121 120 101;
                   72 92 95 98 112 100 103 99];
%% ARRAYS TO 8X8 CELLS
    n=1;
    for j=1:BLOCKSIZE:height
        if j+BLOCKSIZE>height
             for i = 1:BLOCKSIZE:width
                   if i+BLOCKSIZE>width
                      blockR(:,:,n)=I(i:width,j:height,1);
                      blockG(:,:,n)=I(i:width,j:height,2);
                      blockB(:,:,n)=I(i:width,j:height,3);           
                       break;
                   else 
                      blockR(:,:,n)=I(i:BLOCKSIZE+i-1,j:height,1);
                      blockG(:,:,n)=I(i:BLOCKSIZE+i-1,j:height,2);
                      blockB(:,:,n)=I(i:BLOCKSIZE+i-1,j:height,3);
                   end
                   n=n+1;
             end
             break;
        else
            for i = 1:BLOCKSIZE:width
                   if i+BLOCKSIZE>width
                      blockR(:,:,n)=I(i:width,j:BLOCKSIZE+j-1,1);
                      blockG(:,:,n)=I(i:width,j:BLOCKSIZE+j-1,2);
                      blockB(:,:,n)=I(i:width,j:BLOCKSIZE+j-1,3);
                       n=n+1;
                       break;
                   else 
                       blockR(:,:,n)=I(i:BLOCKSIZE+i-1,j:BLOCKSIZE+j-1,1);
                       blockG(:,:,n)=I(i:BLOCKSIZE+i-1,j:BLOCKSIZE+j-1,2);
                       blockB(:,:,n)=I(i:BLOCKSIZE+i-1,j:BLOCKSIZE+j-1,3);
                       n=n+1;
                   end
           end
        end
    end

numberofBlocks=size(blockR); 
%% DCT
for l=1:1:numberofBlocks(3)
    y(:,:,l)=dct2(blockR(:,:,l));
    cb(:,:,l)=dct2(blockG(:,:,l));
    cr(:,:,l)=dct2(blockB(:,:,l));
end
%% QUANTITY 
for i = 1:1:numberofBlocks(3)
    for j = 1:1:BLOCKSIZE
        for k = 1:1:BLOCKSIZE
            y2(k,j,i) =  round(y(k,j,i)/tab_kwantyzacji(k,j));
            cr2(k,j,i) = round(cr(k,j,i)/tab_kwantyzacji(k,j));
            cb2(k,j,i) = round(cb(k,j,i)/tab_kwantyzacji(k,j));             
        end
    end
end


%% ZIG-ZAG
y3 =zeros(64,4096);
cb3 =zeros(64,4096);
cr3 =zeros(64,4096);

for j = 1:1:numberofBlocks(3)
    y3(:,j) = zigzag(y2(:,:,j));
    cb3(:,j) = zigzag(cb2(:,:,j));
    cr3(:,j) = zigzag(cr2(:,:,j));
end


%% 1 vec

B = [];%zeros(1,262144);
R = [];%zeros(1,262144);
G = [];%zeros(1,262144);

for i =1:1:numberofBlocks(3)
    for k=1:1:BLOCKSIZE^2
      R(end+1)=y3(k,i);
      G(end+1)=cb3(k,i);
      B(end+1)=cr3(k,i);
   end
end

long = [R,G,B];
%% RLE ENCODING
vec = long(1,:);
y_1=[];
c = 1;

for i = 1:length(vec)-1
 if(vec(i) == vec(i+1))
    c=c+1;
 else
    y_1 = [y_1,c,vec(i),];
     c=1;
end
end
y_1 = [y_1,c,vec(length(vec))];

%% HUFFMAN CODING

symbols = unique(y_1);
counts = hist(y_1, symbols);
p = double(counts) ./ sum(counts);
[dict,avglen] = huffmandict(symbols,p); 
comp = huffmanenco(y_1,dict);

%% save

    save('kodek_final_bmpf.jet', 'dict', 'comp');

%% HUFFMAN DECODING
dhsig1 = huffmandeco(comp,dict);


%% RLE DECODING


 y_1= dhsig1;
y_2 = [];
for p = 1:2:length(y_1)
    for q = 1: 1:y_1(p)
        y_2(end+1) = y_1(p+1); 
    end
end

y5 = y_2(1:BLOCKSIZE^2*numberofBlocks(3));
cb5 = y_2(1+BLOCKSIZE^2*numberofBlocks(3):2*BLOCKSIZE^2*numberofBlocks(3));
cr5 = y_2(2*BLOCKSIZE^2*numberofBlocks(3)+1:3*BLOCKSIZE^2*numberofBlocks(3));


y6=reshape(y5,BLOCKSIZE^2,numberofBlocks(3));
cb6=reshape(cb5,BLOCKSIZE^2,numberofBlocks(3));
cr6=reshape(cr5,BLOCKSIZE^2,numberofBlocks(3));
       





%% GIZ-GAZ
unzigx = zeros(8,8,4096);
unzigx1 = zeros(8,8,4096);
unzigx2 = zeros(8,8,4096);


for i = 1:1:numberofBlocks(3)
    x=y6(:,i).';
    unzigx(:,:,i)=izigsc(x,BLOCKSIZE);
    x1=cb6(:,i).';
    unzigx1(:,:,i)=izigsc(x1,BLOCKSIZE);
    x2=cr6(:,i).';
    unzigx2(:,:,i)=izigsc(x2,BLOCKSIZE);
end

%% DECOMPRESION
for i = 1:1:numberofBlocks(3)
    for j = 1:1:BLOCKSIZE
        for k = 1:1:BLOCKSIZE
            y(k,j,i) = round(unzigx(k,j,i)*tab_kwantyzacji(k,j));
            cr(k,j,i) = round(unzigx2(k,j,i)*tab_kwantyzacji(k,j));
            cb(k,j,i) = round(unzigx1(k,j,i)*tab_kwantyzacji(k,j));        
        end
    end
end

for p = 1:1:numberofBlocks(3)
    y2(:,:,p) = idct2(y(:,:,p));
    cb2(:,:,p) = idct2(cb(:,:,p));
    cr2(:,:,p) = idct2(cr(:,:,p));
end


%% CELLS TO ARRAYS
decodedR=zeros(512,512);
decodedG=zeros(512,512);
decodedB=zeros(512,512);
numberofColumns=floor(width/BLOCKSIZE);
numberofRows=floor(height/BLOCKSIZE);

    counter=1;
    for i=0:1:numberofRows-1
        for j=0:1:numberofColumns-1      
            decodedR(j*BLOCKSIZE+1:(j+1)*BLOCKSIZE,1+i*BLOCKSIZE:(i+1)*BLOCKSIZE)=y2(:,:,counter);
            decodedG(j*BLOCKSIZE+1:(j+1)*BLOCKSIZE,1+i*BLOCKSIZE:(i+1)*BLOCKSIZE)=cb2(:,:,counter);
            decodedB(j*BLOCKSIZE+1:(j+1)*BLOCKSIZE,1+i*BLOCKSIZE:(i+1)*BLOCKSIZE)=cr2(:,:,counter);
            counter=counter+1;
        end
    end

img(:,:,1)=(decodedR);
img(:,:,2)=(decodedG);
img(:,:,3)=(decodedB);
%% END
img=uint8(img);
img= ycbcr2rgb(img);
img2 = imbilatfilt(img);
I = ycbcr2rgb(I);
err = immse(I,img);
err2 = immse(I,img2);
%imshow(img2);
montage({I,img,img2})
imwrite(img, 'new_bmp1.jpg')
imwrite(img2, 'new_bmpf2.jpg')
czas = toc;









disp1 = ['time   ',num2str(czas)];
disp2 = ['err without filter   ',num2str(err)];
disp3 = ['err with filter ',num2str(err2)];
disp(disp1)
disp(disp2)
disp(disp3)
%% FUNCTION

function C = zigzag(A)
r = size(A,1);
c = size(A,2); 
kk = 2;
C = [];
while kk <= r+c 
                
                
                
    B = [];
    


for ii = 1:r
    for jj = 1:c
        if ii + jj == kk 
                         
            B = [B,A(ii,jj)];
        end
    end
end
if mod(kk,2) == 0
    C = [C,flip(B)]; 
                    
else
    C = [C,B];
end
kk = kk+1;
end
        
end



function [A] = izigsc(B,dim)
v = ones(1,dim); k = 1;
A = zeros(dim,dim);
for i = 1:2*dim-1
    C1 = diag(v,dim-i);
    C2 = flip(C1(1:dim,1:dim),2);
    C3 = B(k:k+sum(C2(:))-1);
    k = k + sum(C2(:));
    if mod(i,2) == 0
       C3 = flip(C3);
    end
    C4 = zeros(1,dim-size(C3,2));
    if i >= dim
       C5 = cat(2,C4, C3); 
    else       
        C5 = cat(2,C3,C4);
    end
    C6 = C2*diag(C5);
    A = C6 + A;
end
end

