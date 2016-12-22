% Get masks of candidate nuclei and purple areas.
% Build the circlar masks used for cell detection

% Mask 1
inner_rad = 60;
outer_rad = 80;
ec_1 = floor(0.5*(outer_rad+inner_rad)*2*pi);

center = outer_rad+2;
mask1_size = center*2;

circleimg = zeros(mask1_size,mask1_size);
circleimg = MidpointCircle(circleimg, outer_rad, center, center, 1);

cimg2 = zeros(mask1_size,mask1_size);
cimg2 = MidpointCircle(cimg2, inner_rad, center, center, 1);

cimg3 = imcomplement(or(circleimg, cimg2));

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

% ************* MASK 2 *************
inner_rad = 72;
outer_rad = 104;
ec_2 = floor(0.5*(outer_rad+inner_rad)*2*pi);

center = outer_rad+2; 
mask2_size = center*2; 

circleimg = zeros(mask2_size,mask2_size);
circleimg = MidpointCircle(circleimg, outer_rad, center, center, 1);

cimg2 = zeros(mask2_size,mask2_size);
cimg2 = MidpointCircle(cimg2, inner_rad, center, center, 1);

cimg3 = imcomplement(or(circleimg, cimg2));

c_label = bwlabel(cimg3, 4);

cimg4 = zeros(mask2_size,mask2_size);
for i=1:mask2_size
    for j=1:mask2_size
        if (c_label(i,j) == 2) 
            cimg4(i,j) = 1;
        end
    end
end

mask2 = cimg4;

% ************* MASK 3 *************
inner_rad = 96;
outer_rad = 128;
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