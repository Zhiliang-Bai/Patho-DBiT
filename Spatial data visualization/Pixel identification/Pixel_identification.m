% This matlab script will analyze the tissue image of DBiT-seq samples and
% identify the pixels that are on top tissue. The 'position.txt' generated
% include all the pixels that are on top of a tissue. 

% Enter the Example_Data folder

% image (cropped exactly to the size of working region of DBiT-seq) needs
% to be processed to make it black (background) and white (tissue)
% beforehand, see FFPE-2_BW.jpg for example
I = imread('Mouse_Brain_Reg.9_50um.jpg'); 
I = rgb2gray(I);
BW = imbinarize(I);

% show the figure on the screen
figure
imshow(BW)

%Set pixel number, currently we use 50 rows and 50 columns.
pixel = 50;
pixel_count = 2*pixel-1;
[numRows,numCols] = size(BW);
pixel_w = numCols/pixel_count;
pixel_h = numRows/pixel_count;
str ="";

%Identify the pixel through iteration
for i = 1:50
    y = round(2*(i-1)*pixel_h + 1);
    for j = 1:50
        x = round(2*(j-1)*pixel_w + 1);
        pixel = BW(y:round(y+pixel_h-1),x:round(x+pixel_w-1));
        C = sum(pixel,'all');
        if C > 0
            str = str+','+j+'x'+i;   
        end
    end    
end

%Output the coordinates of pixels that are on top of a tissue.
fid = fopen('position.txt','wt');
fprintf(fid, str);
fclose(fid);