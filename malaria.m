close all; clear all; clc;

% Specify the input image here.
im = imread('g3.jpg');

% Conversion to HSV
im_hsv = rgb2hsv(im);

%Separation of H, S and V components

im_H = im_hsv(:,:,1);
im_S = im_hsv(:,:,2); 
im_V = im_hsv(:,:,3);

% Applying histogram equalization

im_H_adj = histeq(im_H);
im_S_adj = histeq(im_S);
im_V_adj = histeq(im_V);

% Building the area mask

im_areamask = imcomplement(im2bw(im,graythresh(rgb2gray(im))));

% Build the RBC outline mask

im_edgemask = (im_V_adj > 0.6 & im_V_adj < 0.7);

% Tiny region removal

im_edgemask_clean = bwareaopen(im_edgemask, 200);

% Use Canny edge detection to get edge mask

im_edges = edge(im_edgemask_clean, 'canny');

run mask.m 

% Threshold by H and S components.
% Use S to distinguish nuclei. Use H to distinguish purple.

a = histeq(im_H);
b = histeq(im_S);
c = histeq(im_V);

im_nuclei_mask = and(a==1, b>0.98);

% Clean up the nuclei area masks (most likely just
% small area removal).

nuclei_cleaned = bwareaopen(im_nuclei_mask, 30);

label1 = bwlabel(nuclei_cleaned);
props1 = regionprops(label1, 'Centroid');
centroids = cat(1, props1.Centroid);
num_centroids = size(centroids,1);
final_centers = [];

% For each successful nuclei try to discern a cell border
% surrounding the cell

for i=1:num_centroids
    current_centroid = centroids(i,:);
    current_centroid(1) = round(current_centroid(1));
    current_centroid(2) = round(current_centroid(2));
    left_bound = current_centroid(1)-100;
    right_bound = current_centroid(1)+100;
    bottom_bound = current_centroid(2)-100;
    top_bound = current_centroid(2)+100;
    
    % Skipping boundary cases.
    if (left_bound < 1 || right_bound > size(im,2) || bottom_bound < 1 || top_bound > size(im,1))
        continue;
    end

    current_buff = im(bottom_bound:top_bound, left_bound:right_bound, :);
    
    current_edge = im_edges(bottom_bound:top_bound, left_bound:right_bound);
    current_mask = im_areamask(bottom_bound:top_bound, left_bound:right_bound);
    % Run circle detection algorithm
    
    largest = mask3_size;   % size of largest mask
    
    % zeros matrix
    cb_big = zeros(size(current_buff,1)+largest, size(current_buff,2)+largest);
    % ones matrix
    cb_bigmask = 1+cb_big;
    % setting the central regions of the two matrices to current edge and
    % area mask
    cb_big(largest/2:largest/2+size(current_buff,1)-1, ...
        largest/2:largest/2+size(current_buff,2)-1) = current_edge;
    cb_bigmask(largest/2:largest/2+size(current_buff,1)-1, ...
        largest/2:largest/2+size(current_buff,2)-1) = current_mask;
    
    % centers for the three masks
    centers1 = [];
    centers2 = [];
    centers3 = [];
    
    step_size = 3*3;
    excl_zone = 6;
    cur_x = 1+largest/2;
    cur_y = 1;
    
    x_iters = floor((size(cb_big,1)-largest)/step_size);
    y_iters = floor((size(cb_big,2)-largest)/step_size);
    
    excl_matrix = zeros(x_iters,y_iters);
    
    for ic=1:x_iters
        cur_y = 1+largest/2;
        for jc=1:y_iters
            for k=1:3 
                if (excl_matrix(ic,jc) == 0)
                    if (k == 1)
                        cur_mask = mask1;
                        mask_size = mask1_size;
                        ec = ec_1;
                    elseif (k == 2)
                        cur_mask = mask2;
                        mask_size = mask2_size;
                        ec = ec_2;
                    elseif (k == 3)
                        cur_mask = mask3;
                        mask_size = mask3_size;
                        ec = ec_3;
                    end
                    
                    segment_thresh = 0.3*2*ec;
                    
                    temp_buff = zeros(mask_size, mask_size);
                    
                    % Copy the image into the temp buffer           
                    temp_buff = cb_big(cur_x-(mask_size/2):cur_x+(mask_size/2)-1, ...
                        cur_y-(mask_size/2):cur_y+(mask_size/2)-1);
                    
                    % 'AND' the contents of the window with the ring-shaped
                    % mask to try and detect a circular outline
                    temp_buff = and(temp_buff, cur_mask);
                    
                    % METHOD: Do image segmentation and then look at the
                    %  size of the longest segment (number of pixels).
                    % (A better idea may be to either take the longest n segments,
                    % or the sum of all segments longer than a certain percentage
                    % of the mask circumference...)
                    [buflabel num] = bwlabel(temp_buff, 8);
                    num_long_segs = 0;      % Alternate metric-count the number of long segments
                    if (num == 1)
                        % If only one segment, just use the simple method
                        val = sum(temp_buff(:));
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
        size(centers3,1);
    
    all_vec = [centers1 1+zeros(size(centers1,1), 1); ...
        centers2 2+zeros(size(centers2,1), 1); ...
        centers3 3+zeros(size(centers3,1), 1)];
    
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
        final_centers = [final_centers; ct1 type current_centroid];
    end
    
    fprintf('Run %d: found %d circles...\n', i, num_cells_detected);
end
