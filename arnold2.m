% Arnold transform, v2, 2012-02-18
% Piotr Sklodowski, <piotr.sk@gmail.com>
%
% @in: two dimensional matrix
% @iter: number of iterations
%
function [ out ] = arnold2( in, iter )
in=imread('tree.png');
iter=48;
redim=in(:,:,1);
greenim=in(:,:,2);
blueim=in(:,:,3);
    [m n chan] = size(in);
    if (m ~= n)
        error(['Arnold Transform is defined only for squares. ' ...
        'Please complete empty rows or columns to make the square.']);
    end
    
   %For red color
   
    redout = zeros(m);
    n = n - 1;
    for j=1:iter
        for y=0:n
            for x=0:n
                p = [ 1 1 ; 1 2 ] * [ x ; y ];
                redout(mod(p(2), m)+1, mod(p(1), m)+1) = redim(y+1, x+1);
            end
        end
        redim = redout;
    end
    %redout=uint8(redout);
    subplot(1,4,1);
    imshow(redout);
   % imwrite(out,'C:\Users\nitin\Documents\MATLAB\scrambled.png');
   
   %For green color
   
    greenout = zeros(m);
    n = n - 1;
    for j=1:iter
        for y=0:n
            for x=0:n
                p = [ 1 1 ; 1 2 ] * [ x ; y ];
                greenout(mod(p(2), m)+1, mod(p(1), m)+1) = greenim(y+1, x+1);
            end
        end
        greenim = greenout;
    end
    %greenout=uint8(greenout);
    subplot(1,4,2);
    imshow(greenout);
    
    %For red color
   
    blueout = zeros(m);
    n = n - 1;
    for j=1:iter
        for y=0:n
            for x=0:n
                p = [ 1 1 ; 1 2 ] * [ x ; y ];
                blueout(mod(p(2), m)+1, mod(p(1), m)+1) = blueim(y+1, x+1);
            end
        end
        blueim = blueout;
    end
  %  blueout=uint8(blueout);
    subplot(1,4,3);
    imshow(blueout);
    
    out=cat(3,uint8(redout),uint8(greenout),uint8(blueout));
    subplot(1,4,4);
    imwrite(out,'C:\Users\nitin\Documents\MATLAB\scrambled.png');
    
end