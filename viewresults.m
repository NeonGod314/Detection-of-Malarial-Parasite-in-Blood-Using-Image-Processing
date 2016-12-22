% After running parasite1.m, use this file to view results.
cbsize = size(cb_big,1);

imshow(im);
hold on;

for i=1:size(final_centers,1)
    
    actual_x = final_centers(i,2)-cbsize/2;
    actual_x = actual_x + final_centers(i,4);
    actual_y = final_centers(i,1)-cbsize/2;
    actual_y = actual_y + final_centers(i,5);
    
    plot(actual_x, actual_y, 'ro');
    
    if (final_centers(i,3) == 1)
        rad = ec_1;
        fmt = 'b--';
    elseif (final_centers(i,3) == 2)
        rad = ec_2;
        fmt = 'c--';
    elseif (final_centers(i,3) == 3)
        rad = ec_3;
        fmt = 'g--';
    elseif (final_centers(i,3) == 4)
        rad = ec_4;
        fmt = 'm--';
    elseif (final_centers(i,3) == 5)
        rad = ec_5;
        fmt = 'y--';
    elseif (final_centers(i,3) == 6)
        rad = ec_6;
        fmt = 'k--';
    end
    
    rad = rad/(2*pi);  
    circle([actual_x, actual_y], rad, 100, fmt);
end
hold off;
