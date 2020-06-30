function [mean_TR, mean_GFP] = rim_int(I_TR,I_GFP)
% The function rim_int takes two same-size grayscale images (I_TR and
% I_GFP) and outputs the mean intensity at the rim of circular objects
% detected in I_TR. In the published work, I_TR is the Texas Red membrane
% dye channel and I_GFP is the channel of the fluorescent LactC2-eGFP
% reporter. 
% Code written by Adriana Serrao and modified by Duco Blanken (2020)

h = [-1 -1 -1; -1 12 -1; -1 -1 -1]./4; % defines a sharpness filter

I_sharp = imfilter(I_TR,h); %applies the sharpness filter

I_filled = imfill(I_sharp); %fill the areas inside liposomes

I_inside = I_filled-I_sharp; % subtract sharpened image from filled image 
                             % to retain the lumen marker 

fgm = zeros(size(I_TR)); %create a matrix the size of I_TR

fgm(I_inside<=100)= 0; %binarize I_inside to create the foreground marker fgm
fgm(I_inside>100) =1;

fgm = imfill(fgm,'holes'); %fill any holes in the binary picture


n_erosion = 5; 
se = strel('disk',n_erosion); 
fgm = imerode(fgm,se);%Erode 5 pixels all around

CC = bwconncomp(fgm); %extract connected components from foreground marker
pixels = CC.PixelIdxList;
props = regionprops(CC,'PixelIdxList','Perimeter','Area','Centroid'); 
%extract properties from each object


% select for circular things, remove all non-spherical things from fgm
for j = 1:length(pixels)
    P(j) = props(j).Perimeter;
    A(j) = props(j).Area;

    C(j) = P(j)^2/(4*pi*A(j)); %circularity criterion
    
    if C(j) >1.05 %not a circle
        fgm(CC.PixelIdxList{j}) = 0;
    elseif C(j) < 0.95 %not a circle
        fgm(CC.PixelIdxList{j}) = 0; 
    end
        %what's left is sufficiently circular for further analysis
end

CC = bwconncomp(fgm); %extract connected components from the marker without 
%non-circular objects
pixels = CC.PixelIdxList;
props = regionprops('table',CC,'Centroid','MajorAxisLength','MinorAxisLength');

centroids = cat(1, props.Centroid);
diameters = mean([props.MajorAxisLength props.MinorAxisLength],2); %determine diameter of objects
radii = diameters/2 + n_erosion; %add erosion disk to get 'real' radius

clear CC
clear pixels
theta_vector = 0:0.1:2*pi; %angles to be looped over
for n = 1:length(centroids) %loop over all detected liposomes
d = 0:0.1:1.5*radii(n); % radial position from 0 (center liposome) to 1.5 * radius 
for i = 1:63 %loop over all angles
theta = theta_vector(i);
for j = 1:length(d) %loop over radial position
    try
    TR_line{n}{i}(j) = double(I_TR(round(centroids(n,2)+sin(theta)*d(j)),round(centroids(n,1)+cos(theta)*d(j))));
    GFP_line{n}{i}(j) = double(I_GFP(round(centroids(n,2)+sin(theta)*d(j)),round(centroids(n,1)+cos(theta)*d(j))));
    catch
    end
    %this loop measures a radial intensity profile along angle theta for TR 
    %and GFP from the liposome center to 1.5x the radius (somewhere beyond
    %the membrane)
end
TR{n}(i) = double(max(TR_line{n}{i}));   %get the peak in the radial profile
GFP{n}(i) = double(max(GFP_line{n}{i})); %that corresponds to the membrane intersection
%TR and GFP are vectors of intensity at the rim at 63 angles
end
mean_TR(n) = mean(TR{n}); %average of the rim intensity for I_TR
mean_GFP(n) = mean(GFP{n});%average of the rim intensity for I_GFP
end

