
function dec()
clear;
global f_redkey; 
global f_bluekey;
global f_greenkey;
global enc_red;
global enc_green;
global enc_blue;

m=256; n=256; sd_red=zeros(1,256);
    for i=1:m
        for j=1:n
            sd_red(i,j) = mod((double(enc_red(i,j))- double(f_redkey(i,j))),256);
        end
    end
 figure('name','first step decryption');
 subplot(2,2,1);
imagesc(uint8(sd_red));

sd_blue=zeros(1,256);
for i=1:m
    for j=1:n
        sd_blue(i,j) = mod((double(enc_blue(i,j))-double(f_bluekey(i,j))),256);
    end
end
subplot(2,2,2)
 imagesc(uint8(sd_blue));
 
 sd_green=zeros(1,256);
for i=1:m
    for j=1:n
         sd_green(i,j) = mod((double(enc_green(i,j))-double(f_greenkey(i,j))),256);
    end
end
 subplot(2,2,3)
 imagesc(uint8(sd_green));
 
 partialpic=cat(3, uint8(sd_red), uint8(sd_green), uint8(sd_blue));
  figure('name','partial decrypted');
  imshow(partialpic);
 figure; histogram(partialpic);
  imwrite(partialpic,'C:\Users\nitin\Documents\MATLAB\Enhanced Modulo Image Encryption\partialpic.png');
 
[rows,cols,chanel]=size(partialpic);
sd_red=partialpic(:,:,1);
sd_green=partialpic(:,:,2);
sd_blue=partialpic(:,:,3);

 % scrambleOrder = randperm(rows*cols);
  % Recover the image, knowing the sort order

recoverOrder = zeros([rows*cols], 2);
recoverOrder(:, 1) = 1 : (rows*cols);
%recoverOrder(:, 2) = scrambleOrder;
% Sort this to find out where each scrambled location needs to be sent to.
newOrder = sortrows(recoverOrder, 2);
% Extract just column 1, which is the order we need.
newOrder = newOrder(:,1);
% Unscramble according to the recoverOrder order.
%{
redChannel = sd_red(newOrder);
greenChannel = sd_green(newOrder);
blueChannel = sd_blue(newOrder);
% Reshape into a 2D image
redChannel = reshape(redChannel, [rows, cols]);
greenChannel = reshape(greenChannel, [rows, cols]);
blueChannel = reshape(blueChannel, [rows, cols]);
%scrambledImage = cat(3, redChannel, greenChannel, blueChannel);
%}  
 a=37;b=35;n=25;
 
 %Unscramble for Red Color
 iteration = 1;
redpartialdecrypted = sd_red;
% The number of iterations needed to restore the image can be shown never to exceed 3N.
N = rows; 
while iteration <= n
	% Scramble the image based on the old image.
    axis xy;
	for row = 1 : rows % y
		for col = 1 : cols % x
			c = mod(-(35*row)+col, N) +1; % x coordinate
			r = mod((1296*row) - (37*col), N)+1 ; % y coordinate
			% Move the pixel.  Note indexes are (row, column) = (y, x) NOT (x, y)!
			redcurrentdecrypted(row, col, :) = redpartialdecrypted(r, c, :);
        end
    end
    caption = sprintf('Iteration #%d', iteration);
	fprintf('%s\n', caption);
        % Make the current image the prior/old one so we'll operate on that the next iteration.
	redpartialdecrypted = redcurrentdecrypted;
    	% Update the iteration counter.
	iteration = iteration+1;
end
redcurrentdecrypted = redcurrentdecrypted(newOrder);
redcurrentdecrypted = reshape(redcurrentdecrypted, [rows, cols]);

 %Unscramble for green Color
 iteration = 1;
greenpartialdecrypted = sd_green;
% The number of iterations needed to restore the image can be shown never to exceed 3N.
N = rows; 
while iteration <= n
	% Scramble the image based on the old image.
    axis xy;
	for row = 1 : rows % y
		for col = 1 : cols % x
			c = mod(-(35*row)+col, N) +1; % x coordinate
			r = mod((1296*row) - (37*col), N)+1 ; % y coordinate
			% Move the pixel.  Note indexes are (row, column) = (y, x) NOT (x, y)!
			greencurrentdecrypted(row, col, :) = greenpartialdecrypted(r, c, :);
        end
    end
    caption = sprintf('Iteration #%d', iteration);
	fprintf('%s\n', caption);
        % Make the current image the prior/old one so we'll operate on that the next iteration.
	greenpartialdecrypted = greencurrentdecrypted;
    	% Update the iteration counter.
	iteration = iteration+1;
end
greencurrentdecrypted = greencurrentdecrypted(newOrder);
greencurrentdecrypted = reshape(greencurrentdecrypted, [rows, cols]);
 
%Unscramble for blue Color
 iteration = 1;
bluepartialdecrypted = sd_blue;
% The number of iterations needed to restore the image can be shown never to exceed 3N.
N = rows; 
while iteration <= n
	% Scramble the image based on the old image.
    axis xy;
	for row = 1 : rows % y
		for col = 1 : cols % x
			c = mod(-(35*row)+col, N) +1; % x coordinate
			r = mod((1296*row) - (37*col), N)+1 ; % y coordinate
			% Move the pixel.  Note indexes are (row, column) = (y, x) NOT (x, y)!
			bluecurrentdecrypted(row, col, :) = bluepartialdecrypted(r, c, :);
        end
    end
    caption = sprintf('Iteration #%d', iteration);
	fprintf('%s\n', caption);
        % Make the current image the prior/old one so we'll operate on that the next iteration.
	bluepartialdecrypted = bluecurrentdecrypted;
    	% Update the iteration counter.
	iteration = iteration+1;
end
bluecurrentdecrypted = bluecurrentdecrypted(newOrder);
bluecurrentdecrypted = reshape(bluecurrentdecrypted, [rows, cols]);
%{
scrambleOrder = randperm(rows*cols);
  % Recover the image, knowing the sort order
recoverOrder = zeros([rows*cols], 2);
recoverOrder(:, 1) = 1 : (rows*cols);
recoverOrder(:, 2) = scrambleOrder;
% Sort this to find out where each scrambled location needs to be sent to.
newOrder = sortrows(recoverOrder, 2);
% Extract just column 1, which is the order we need.
newOrder = newOrder(:,1);
% Unscramble according to the recoverOrder order.
redChannel = partialpic(newOrder);
%greenChannel = sd_green(newOrder);
%blueChannel = sd_blue(newOrder);
% Reshape into a 2D image
scrambledImage = reshape(redChannel, [rows, cols]);
%blueChannel = reshape(blueChannel, [rows, cols]);
%scrambledImage = cat(3, redChannel, greenChannel, blueChannel);
%}
currentdecryptedImage = cat(3, redcurrentdecrypted, greencurrentdecrypted, bluecurrentdecrypted);
 subplot(1,2,1);
        imshow(currentdecryptedImage);
        imwrite(currentdecryptedImage,'C:\Users\nitin\Documents\MATLAB\Enhanced Modulo Image Encryption\currentdecryptedImage.png');
    im=imtranslate(currentdecryptedImage,[-50,-25]);
    subplot(1,2,2); imshow(im);
end

