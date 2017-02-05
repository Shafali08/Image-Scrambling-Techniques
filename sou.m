%Encryption Code

function [f_redkey,f_bluekey,f_greenkey] = sou()

global f_redkey; 
global f_bluekey;
global f_greenkey;

plainimage = imread('tree.png'); %source image

figure('name','source');

subplot(1,4,1);
imshow(plainimage);
oneinsource=nnz(plainimage==1); %total number of 1 in source image

keyimage = imread('mandelkey.png'); %key image

subplot(1,4,2); imshow(keyimage);
oneinkey=nnz(keyimage==1); %total number of 1 in key image

p_red=plainimage(:,:,1); %red plan of source Image
k_red=keyimage(:,:,1); %red plan of key Image
sumR=sum(k_red(:));     %sum of red matrix pixels

p_green=plainimage(:,:,2); %green plan of source Image
k_green=keyimage(:,:,2); %green plan of key Image
sumG=sum(k_green(:));     %sum of green matrix pixels

p_blue=plainimage(:,:,3); %blue plan of source Image
k_blue=keyimage(:,:,3); %blue plan of key Image
sumB=sum(k_blue(:));     %sum of blue matrix pixels

%Calculate Arnold cat map parameters
a=13+mod((sumR+oneinsource),29);
b=7+mod((sumG+oneinsource),47);
n1=3+mod((sumB+oneinsource),13);
n2=4+mod(oneinsource,11);
n3=5+mod(oneinkey,9);

[rows,cols,chan]=size(plainimage);
n=n1+n2+n3;
%{
scrambleOrder = (rows*cols);
plainimage = plainimage(scrambleOrder);
plainimage = reshape(plainimage, [rows, cols]);
%}
%scrambleOrder = (rows*cols);
% Extract the individual red, green, and blue color channels.
redChannel = plainimage(:, :, 1);
greenChannel = plainimage(:, :, 2);
blueChannel = plainimage(:, :, 3);
%{
% Scramble according to the scrambling order.
redChannel = redChannel(scrambleOrder);
greenChannel = greenChannel(scrambleOrder);
blueChannel = blueChannel(scrambleOrder);

% Reshape into a 2D image
redChannel = reshape(redChannel, [rows, cols]);
greenChannel = reshape(greenChannel, [rows, cols]);
blueChannel = reshape(blueChannel, [rows, cols]);
%}
%Arnold Scrambling for red channel

iteration = 1;
% Initialize image.
oldredscrambled = redChannel;
% The number of iterations needed to restore the image can be shown never to exceed 3N.
N = rows;
while iteration <= n
	% Scramble the image based on the old image.
	for row = 1 : rows % y
		for col = 1 : cols % x
			c = mod((1296 * col) + (35*row), N) + 1; % x coordinate
			r = mod((37*col) + row, N) + 1; % y coordinate
			% Move the pixel.  Note indexes are (row, column) = (y, x) NOT (x, y)!
			currentredScrambled(row, col, :) = oldredscrambled(r, c, :);
        end
    end
    caption = sprintf('Iteration #%d', iteration);
	fprintf('%s\n', caption);
	%subplot(1,4,3);
     %   imshow(currentredScrambled);
	   
	% Make the current image the prior/old one so we'll operate on that the next iteration.
	oldredscrambled = currentredScrambled;
    
	% Update the iteration counter.
	iteration = iteration+1;
end

%Arnold Scrambling for green channel
iteration = 1;
% Initialize image.
oldgreenScrambled = greenChannel;
% The number of iterations needed to restore the image can be shown never to exceed 3N.
N = rows;
while iteration <= n
	% Scramble the image based on the old image.
	for row = 1 : rows % y
		for col = 1 : cols % x
			c = mod((1296 * col) + (35*row), N) + 1; % x coordinate
			r = mod((37*col) + row, N) + 1; % y coordinate
			% Move the pixel.  Note indexes are (row, column) = (y, x) NOT (x, y)!
			currentgreenScrambled(row, col, :) = oldgreenScrambled(r, c, :);
        end
    end
    caption = sprintf('Iteration #%d', iteration);
	fprintf('%s\n', caption);
	%subplot(1,4,3);
  %      imshow(currentgreenScrambled);
	   
	% Make the current image the prior/old one so we'll operate on that the next iteration.
	oldgreenScrambled = currentgreenScrambled;
    
	% Update the iteration counter.
	iteration = iteration+1;
end

%Arnold Scrambling for blue channel
iteration = 1;
% Initialize image.
oldblueScrambled = blueChannel;
% The number of iterations needed to restore the image can be shown never to exceed 3N.
N = rows;
while iteration <= n
	% Scramble the image based on the old image.
	for row = 1 : rows % y
		for col = 1 : cols % x
			c = mod((1296 * col) + (35*row), N) + 1; % x coordinate
			r = mod((37*col) + row, N) + 1; % y coordinate
			% Move the pixel.  Note indexes are (row, column) = (y, x) NOT (x, y)!
			currentblueScrambled(row, col, :) = oldblueScrambled(r, c, :);
        end
    end
    caption = sprintf('Iteration #%d', iteration);
	fprintf('%s\n', caption);
%	subplot(1,4,3);
 %       imshow(currentblueScrambled);
	   
	% Make the current image the prior/old one so we'll operate on that the next iteration.
	oldblueScrambled = currentblueScrambled;
    
	% Update the iteration counter.
	iteration = iteration+1;
end
currentScrambledImage = cat(3, currentredScrambled, currentgreenScrambled, currentblueScrambled);
subplot(1,4,4);
        imshow(currentScrambledImage);
        imwrite(currentScrambledImage,'C:\Users\nitin\Documents\MATLAB\Enhanced Modulo Image Encryption\currentScrambledImage.png');
        
 [m,n,p]=size(keyimage);
de=40; 
r=zeros(1,256);
for i=1:m
    for j=1:n
        kmax=floor((i-1)/(de+1))+floor((m-i)/(de+1))+1;
        lmax=floor((j-1)/(de+1))+floor((n-j)/(de+1))+1;
        k=i-floor((i-1)/(de+1))*(de+1)+(kmax-1)*(de+1);
        l=j-floor((j-1)/(de+1))*(de+1)+(lmax-1)*(de+1);

        r(i,j)=sqrt((((i-floor((i-1)/(de+1))*(de+1)+(k-1)*(de+1))-i).^2)+((j-floor((j-1)/(de+1))*(de+1)+(l-1)*(de+1))-j).^2);
        
    end
end

%For red layer
for i=1:m
    for j=1:n
        sum4=0;sum1=0;sum2=0;sum3=0;
        for s=i:de:m
            if(s<=m)
            sum4=double(sum4+k_red(s,j));
            end
        end
        for s=i-de:-de:0
            if(s>0)
            sum1=double(sum1+k_red(s,j));
            end
        end
        for t=1:de:n
            if(t<=n)
            sum2=double(sum2+k_red(i,t));
            end
        end
        for t=j-de:-de:0
            if(t>0)
            sum3=double(sum3+k_red(i,t));
            end
        end
%     end  
      d_red(i,j)=sum4+sum1+sum2+sum3;
    end
end
key_red=0; 
f_redkey=zeros(1,256);
 for k=1:m
     for l=1:n
         key_red=key_red+(r(k,l)*d_red(k,l));
         f_redkey(k,l)=key_red;
         key_red=0;
     end
 end      

%For Blue Layer
 for i=1:m
    for j=1:n
        sum4=0;sum1=0;sum2=0;sum3=0;
        for s=i:de:m
            if(s<=m)
            sum4=double(sum4+k_blue(s,j));
            end
        end
        for s=i-de:-de:0
            if(s>0)
            sum1=double(sum1+k_blue(s,j));
            end
        end
        for t=1:de:n
            if(t<=n)
            sum2=double(sum2+k_blue(i,t));
            end
        end
        for t=j-de:-de:0
            if(t>0)
            sum3=double(sum3+k_blue(i,t));
            end
        end
%     end  
      d_blue(i,j)=sum4+sum1+sum2+sum3;
    end
end
key_blue=0; f_bluekey=zeros(1,256);
 for k=1:m
     for l=1:n
         key_blue=key_blue+(r(k,l)*d_blue(k,l));
         f_bluekey(k,l)=key_blue;
         key_blue=0;
     end
 end      

%  for Green Layer
 for i=1:m
    for j=1:n
        sum4=0;sum1=0;sum2=0;sum3=0;
        for s=i:de:m
            if(s<=m)
            sum4=double(sum4+k_green(s,j));
            end
        end
        for s=i-de:-de:0
            if(s>0)
            sum1=double(sum1+k_green(s,j));
            end
        end
        for t=1:de:n
            if(t<=n)
            sum2=double(sum2+k_green(i,t));
            end
        end
        for t=j-de:-de:0
            if(t>0)
            sum3=double(sum3+k_green(i,t));
            end
        end
   %   end  
      d_green(i,j)=sum4+sum1+sum2+sum3;
    end
end
key_green=0; f_greenkey=zeros(1,256);
 for k=1:m
     for l=1:n
         key_green=key_green+(r(k,l)*d_green(k,l));
         f_greenkey(k,l)=key_green;
         key_green=0;
     end
 end     

[enc_red,enc_green,enc_blue] = encryption(f_redkey,f_bluekey,f_greenkey);

end


function [enc_red,enc_green,enc_blue] = encryption(f_redkey,f_bluekey,f_greenkey)

global enc_red;
global enc_green;
global enc_blue;
sourceimage = imread('C:\Users\nitin\Documents\MATLAB\Enhanced Modulo Image Encryption\currentScrambledImage.png');
%figure('name','scrambled');
imshow(sourceimage);
s_red=sourceimage(:,:,1);
s_green=sourceimage(:,:,2);
s_blue=sourceimage(:,:,3);
figure('name','scrambled');
subplot(2,2,1); imagesc(s_red); 
subplot(2,2,2); imagesc(s_green); 
subplot(2,2,3); imagesc(s_blue); 

m=256;n=256; enc_red=zeros(1,256);
for i=1:m
    for j=1:n
        enc_red(i,j) = mod((double(s_red(i,j))+double(f_redkey(i,j))),256);
    end
end
figure('name','encrypted');
subplot(2,2,1);
imagesc(uint8(enc_red));

enc_blue=zeros(1,256);
for i=1:m
    for j=1:n
        enc_blue(i,j) = mod((double(s_blue(i,j))+double(f_bluekey(i,j))),256);
    end
end
%figure;
subplot(2,2,2);
imagesc(enc_blue);

enc_green=zeros(1,256);
for i=1:m
    for j=1:n
        enc_green(i,j) = mod((double(s_green(i,j))+double(f_greenkey(i,j))),256);
    end
end
%figure;
subplot(2,2,3);
imagesc(enc_green);

pic=cat(3, uint8(enc_red), uint8(enc_green), uint8(enc_blue));
  figure; subplot(3,2,1);
  imshow(pic);
  subplot(3,2,2); histogram(pic);
 
end

