close all; clear all; clc;

% Specify the input image here.
im = imread('g5.jpg');
% Resize factor should remain 1 for the time being.
resize_factor = 1;
im = imresize(im, 1/resize_factor);

im_hsv = rgb2hsv(im);

im_H = im_hsv(:,:,1);
im_S = im_hsv(:,:,2);
im_V = im_hsv(:,:,3);

im_H_adj = histeq(im_H);
im_S_adj = histeq(im_S);
im_V_adj = histeq(im_V);

im_areamask = imcomplement(im2bw(im,graythresh(rgb2gray(im))));

% Build the global cell edge mask, using the V component
im_edgemask = (im_V_adj > 0.6 & im_V_adj < 0.71);

% Tiny region removal, to suppress extremely small signals
label1 = bwlabel(im_edgemask);
props1 = regionprops(label1, 'Area');
idx = find([props1.Area] > 200);
im_edgemask_clean = ismember(label1, idx);

im_edges = edge(im_edgemask_clean, 'canny');	% Canny edge detection

% Build the circlar masks used for cell detection
% Optimally, this would be done in a separate file and the masks saved and
%  loaded as necessary.

% ************* MASK 1 *************
inner_rad = 15*4/resize_factor;
outer_rad = 20*4/resize_factor;
ec_1 = floor(0.5*(outer_rad+inner_rad)*2*pi);

center = outer_rad+2;
mask1_size = center*2;

circleimg = zeros(mask1_size,mask1_size);
circleimg = MidpointCircle(circleimg, outer_rad, center, center, 1);
% plot of mask 1



cimg2 = zeros(mask1_size,mask1_size);
cimg2 = MidpointCircle(cimg2, inner_rad, center, center, 1);
figure
imshow(cimg2);

cimg3 = imcomplement(or(circleimg, cimg2));
figure
imshow(cimg3);
c_label = bwlabel(cimg3, 4);

cimg4 = zeros(mask1_size,mask1_size);
for i=1:mask1_size
    for j=1:mask1_size
        if (c_label(i,j) == 2) 
            cimg4(i,j) = 1;
        end
    end
end

mask1 = cimg4;
figure
imshow(cimg4);
 
% % ************* MASK 2 *************
% inner_rad = 16*4/resize_factor; 
% outer_rad = 24*4/resize_factor;
% ec_2 = floor(0.5*(outer_rad+inner_rad)*2*pi);
% 
% center = outer_rad+2; 
% mask2_size = center*2;
% 
% circleimg = zeros(mask2_size,mask2_size);
% circleimg = MidpointCircle(circleimg, outer_rad, center, center, 1);
% 
% cimg2 = zeros(mask2_size,mask2_size);
% cimg2 = MidpointCircle(cimg2, inner_rad, center, center, 1);
% 
% cimg3 = imcomplement(or(circleimg, cimg2));
% 
% c_label = bwlabel(cimg3, 4);
% 
% cimg4 = zeros(mask2_size,mask2_size);
% for i=1:mask2_size
%     for j=1:mask2_size
%         if (c_label(i,j) == 2) 
%             cimg4(i,j) = 1;
%         end
%     end
% end
% 
% mask2 = cimg4;
% figure
% imshow(cimg4);
% ************* MASK 3 *************
inner_rad = 18*4/resize_factor;
outer_rad = 26*4/resize_factor;
ec_3 = floor(0.5*(outer_rad+inner_rad)*2*pi);

center = outer_rad+2; 
mask3_size = center*2; 

circleimg = zeros(mask3_size,mask3_size);
circleimg = MidpointCircle(circleimg, outer_rad, center, center, 1);

cimg2 = zeros(mask3_size,mask3_size);
cimg2 = MidpointCircle(cimg2, inner_rad, center, center, 1);

cimg3 = imcomplement(or(circleimg, cimg2));

c_label = bwlabel(cimg3, 4);

cimg4 = zeros(mask3_size,mask3_size);
for i=1:mask3_size
    for j=1:mask3_size
        if (c_label(i,j) == 2) 
            cimg4(i,j) = 1;
        end
    end
end

mask3 = cimg4;

% % ************* MASK 4 *************
% inner_rad = 20*4/resize_factor;
% outer_rad = 28*4/resize_factor;
% ec_4 = floor(0.5*(outer_rad+inner_rad)*2*pi);
% 
% center = outer_rad+2; 
% mask4_size = center*2; 
% 
% circleimg = zeros(mask4_size,mask4_size);
% circleimg = MidpointCircle(circleimg, outer_rad, center, center, 1);
% 
% cimg2 = zeros(mask4_size,mask4_size);
% cimg2 = MidpointCircle(cimg2, inner_rad, center, center, 1);
% 
% cimg3 = imcomplement(or(circleimg, cimg2));
% 
% c_label = bwlabel(cimg3, 4);
% 
% cimg4 = zeros(mask4_size,mask4_size);
% for i=1:mask4_size
%     for j=1:mask4_size
%         if (c_label(i,j) == 2) 
%             cimg4(i,j) = 1;
%         end
%     end
% end
% 
% mask4 = cimg4;

% % ************* MASK 5 *************
% inner_rad = 22*4/resize_factor;
% outer_rad = 30*4/resize_factor;
% ec_5 = floor(0.5*(outer_rad+inner_rad)*2*pi);
% 
% center = outer_rad+2; 
% mask5_size = center*2; 
% 
% circleimg = zeros(mask5_size,mask5_size);
% circleimg = MidpointCircle(circleimg, outer_rad, center, center, 1);
% 
% cimg2 = zeros(mask5_size,mask5_size);
% cimg2 = MidpointCircle(cimg2, inner_rad, center, center, 1);
% 
% cimg3 = imcomplement(or(circleimg, cimg2));
% 
% c_label = bwlabel(cimg3, 4);
% 
% cimg4 = zeros(mask5_size,mask5_size);
% for i=1:mask5_size
%     for j=1:mask5_size
%         if (c_label(i,j) == 2) 
%             cimg4(i,j) = 1;
%         end
%     end
% end
% 
% mask5 = cimg4;

% ************* MASK 6 *************
inner_rad = 24*4/resize_factor;
outer_rad = 32*4/resize_factor;
ec_6 = floor(0.5*(outer_rad+inner_rad)*2*pi);

center = outer_rad+2;
mask6_size = center*2;

circleimg = zeros(mask6_size,mask6_size);
circleimg = MidpointCircle(circleimg, outer_rad, center, center, 1);

cimg2 = zeros(mask6_size,mask6_size);
cimg2 = MidpointCircle(cimg2, inner_rad, center, center, 1);

cimg3 = imcomplement(or(circleimg, cimg2));

c_label = bwlabel(cimg3, 4);

cimg4 = zeros(mask6_size,mask6_size);
for i=1:mask6_size
    for j=1:mask6_size
        if (c_label(i,j) == 2) 
            cimg4(i,j) = 1;
        end
    end
end

mask6 = cimg4;
figure;
imshow(mask6);

% Get masks of candidate nuclei and purple areas.

% Threshold by H and S components (purple areas appear brighter,
% with the dark blotches being a very characteristic light color

% Use S to distinguish nuclei (dark spots). Use H to distinguish purple.
h_median = median(im_H(:));
s_median = median(im_S(:));

a = histeq(im_H);
b = histeq(im_S);
c = histeq(im_V);

im_hsv_mod = im_hsv;
im_hsv_mod(:,:,1) = a;
im_hsv_mod(:,:,2) = b;
im_hsv_mod(:,:,3) = c;


im_nuclei_mask = and(a==1, b>0.98);
im_purple_mask = (a==1);
im_nonnuclei = and(im_purple_mask, imcomplement(im_nuclei_mask));


% Clean up the nuclei and purple area masks (most likely just
% small area removal).
label1 = bwlabel(im_nuclei_mask);
props1 = regionprops(label1, 'Area');
idx = find([props1.Area] > 30);
nuclei_cleaned = ismember(label1, idx);

label1 = bwlabel(im_purple_mask);
props1 = regionprops(label1, 'Area');
idx = find([props1.Area] > 30);
purple_cleaned = ismember(label1, idx);

nonn_cleaned = and(purple_cleaned, imcomplement(nuclei_cleaned));

label1 = bwlabel(nuclei_cleaned);
props1 = regionprops(label1, 'Centroid');
centroids = cat(1, props1.Centroid);
 

% For each successful nuclei try to discern a cell border
% surrounding the cell

num_centroids = size(centroids,1);
final_centers = [];

for i=1:num_centroids
    cur_centroid = centroids(i,:);
    cur_centroid(1) = round(cur_centroid(1));
    cur_centroid(2) = round(cur_centroid(2));
    lbound = (cur_centroid(1)-150/resize_factor);
    rbound = (cur_centroid(1)+150/resize_factor);
    tbound = (cur_centroid(2)-150/resize_factor);
    bbound = (cur_centroid(2)+150/resize_factor);
    
    % For the purposes of the project, skip boundary cases.
    if (lbound < 1 || rbound > size(im,2) || tbound < 1 || bbound > size(im,1))
        fprintf('Skipped %d...\n', i);
        continue;
    end
    
    cur_buf = im(tbound:bbound, lbound:rbound, :);
    
    % Not strictly necessary right now...
    cb_hsv = rgb2hsv(cur_buf);
    cb_gray = histeq(cb_hsv(:,:,3));
    
    cb_edge = im_edges(tbound:bbound, lbound:rbound);
    cb_mask = im_areamask(tbound:bbound, lbound:rbound);
    
    % Run circle detection algorithm
    largestmask = mask6_size;
    cb_big = zeros(size(cur_buf,1)+largestmask, size(cur_buf,2)+largestmask);
    cb_bigmask = 1+cb_big;
    cb_big(largestmask/2:largestmask/2+size(cur_buf,1)-1, ...
        largestmask/2:largestmask/2+size(cur_buf,2)-1) = cb_edge;
    cb_bigmask(largestmask/2:largestmask/2+size(cur_buf,1)-1, ...
        largestmask/2:largestmask/2+size(cur_buf,2)-1) = cb_mask;
    
    centers1 = [];
    centers2 = [];
    centers3 = [];
    centers4 = [];
    centers5 = [];
    centers6 = [];
    step_size = 4*4/resize_factor;
    excl_zone = 6;
    cur_x = 1+largestmask/2;
    cur_y = 1;
    
    x_iters = floor((size(cb_big,1)-largestmask)/step_size);
    y_iters = floor((size(cb_big,2)-largestmask)/step_size);
    
    
    excl_matrix = zeros(x_iters,y_iters);
    
    for ic=1:x_iters
        cur_y = 1+largestmask/2;
        for jc=1:y_iters
            for k=1:6  
                if (excl_matrix(ic,jc) == 0)
                    % Kludge
                    if (k == 1)
                        cur_mask = mask1;
                        mask_size = mask1_size;
                        ec = ec_1;
%                     elseif (k == 2)
%                         cur_mask = mask2;
%                         mask_size = mask2_size;
%                         ec = ec_2;
                      elseif (k == 3)
                          cur_mask = mask3;
                          mask_size = mask3_size;
                          ec = ec_3;
%                     elseif (k == 4)
%                         cur_mask = mask4;
%                         mask_size = mask4_size;
%                         ec = ec_4;
%                     elseif (k == 5)
%                         cur_mask = mask5;
%                         mask_size = mask5_size;
%                         ec = ec_5;
                    elseif (k == 6)
                        cur_mask = mask6;
                        mask_size = mask6_size;
                        ec = ec_6;
                    end
                    segment_thresh = 0.3*2*ec;
                    
                    tempbuf = zeros(mask_size, mask_size);
                    
                    
                    % Copy the image into the temp buffer           
                    tempbuf = cb_big(cur_x-(mask_size/2):cur_x+(mask_size/2)-1, ...
                        cur_y-(mask_size/2):cur_y+(mask_size/2)-1);
                    
                    % 'AND' the contents of the window with the ring-shaped
                    % mask to try and detect a circular outline
                    tempbuf = and(tempbuf, cur_mask);
                    
                    % METHOD: Do image segmentation and then look at the
                    %  size of the longest segment (number of pixels).
                    % (A better idea may be to either take the longest n segments,
                    % or the sum of all segments longer than a certain percentage
                    % of the mask circumference...)
                    [buflabel num] = bwlabel(tempbuf, 8);
                    num_long_segs = 0;      % Alternate metric-count the number of long segments
                    if (num == 1)
                        % If only one segment, just use the simple method
                        val = sum(tempbuf(:));
                    else
                        bufprops = regionprops(buflabel, 'Area'); %, 'Eccentricity');
                        bufareas = cat(1, bufprops.Area);
                        
                        % Find the sum of all segments longer than 0.3x the
                        % expected perimeter.
                        val = 0;
                        for areactr=1:length(bufareas)   
                            % Longer segments are considered more desirable.
                            if (bufareas(areactr) > ec*0.5)
                                val = val+bufareas(areactr)*1.3;
                            elseif (bufareas(areactr) > ec*0.4)
                                val = val+bufareas(areactr)*1.15;
                            elseif (bufareas(areactr) > ec*0.3)
                                val = val+bufareas(areactr)*1;
                            elseif (bufareas(areactr) > ec*0.2)
                                val = val+bufareas(areactr)*0.8;
                            end
                        end
                    end

                    if (val > segment_thresh && cb_bigmask(cur_x,cur_y) == 1)
                        % Store this value in the centers array
                        if (k == 1)
                            centers1 = [centers1; cur_x cur_y];
                        elseif (k == 2)
                            centers2 = [centers2; cur_x cur_y];
                        elseif (k == 3)
                            centers3 = [centers3; cur_x cur_y];
                        elseif (k == 4)
                            centers4 = [centers4; cur_x cur_y];
                        elseif (k == 5)
                            centers5 = [centers5; cur_x cur_y];
                        elseif (k == 6)
                            centers6 = [centers6; cur_x cur_y];
                        end
                        % Mark adjacent values as 'do not check'
                        for ii=-excl_zone:excl_zone
                            for jj=-excl_zone:excl_zone
                                cx = ic+ii;
                                cy = jc+jj;
                                if (cx < 1)
                                    cx = 1;
                                end
                                if (cx > size(excl_matrix,1))
                                    cx = size(excl_matrix,1);
                                end
                                if (cy < 1)
                                    cy = 1;
                                end
                                if (cy > size(excl_matrix,2))
                                    cy = size(excl_matrix,2);
                                end
                                
                                excl_matrix(cx,cy) = -1;
                            end
                        end
                    end
                end
            end
            cur_y = cur_y + step_size;
        end
        cur_x = cur_x + step_size;
    end

    
    % Cleanup
    num_cells_detected = size(centers1,1) + size(centers2,1) + ...
        size(centers3,1) + size(centers4,1) + size(centers5,1) + ...
        size(centers6,1);
    
    all_vec = [centers1 1+zeros(size(centers1,1), 1); ...
        centers2 2+zeros(size(centers2,1), 1); ...
        centers3 3+zeros(size(centers3,1), 1); ...
        centers4 4+zeros(size(centers4,1), 1); ...
        centers5 5+zeros(size(centers5,1), 1); 
        centers6 6+zeros(size(centers6,1), 1)];
    
    % Find the circle center closest to the current centroid
    for z=1:num_cells_detected
        dist = (size(cb_big,1)/2-all_vec(z,1))^2 + (size(cb_big,1)/2-all_vec(z,2))^2;
        if (z==1 || dist < mindist)
            mindist = dist;
            ct1 = [all_vec(z,1) all_vec(z,2)];
            type = all_vec(z,3);
        end
    end
    
    if (num_cells_detected > 0)
        final_centers = [final_centers; ct1 type cur_centroid];
    end
    
    fprintf('Run %d: found %d circles...\n', i, num_cells_detected);
end
