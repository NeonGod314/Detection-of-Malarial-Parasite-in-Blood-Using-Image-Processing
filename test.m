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