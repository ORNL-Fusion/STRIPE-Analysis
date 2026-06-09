function [xx_c,xx_s,xx_e,yy_c,yy_s,yy_e] = FindCenter(RawData,ind, pxx_per_cm, pxy_per_cm, rad)
% This funtion locates the center of the plasma by finding the peak
% integrated area over a fixed radius
%  Inputs
%       RawData     3-D array (X,Y,T), either surface temperature or heat flux
%       ind         time index to find center in
%       px_per_cm   pixels per cm number
%       radius      radius to try fit (cm)
%  Output
%       xx_c        x index for center of plasma
%       xx_s        x index for left box of plasma
%       xx_e        x index for right box of plasma
%       yy_c        y index for center of plasma
%       yy_s        y index for left box of plasma
%       yy_e        y index for right box of plasma

    temp = RawData(:,:,ind);
    temp_x = numel(temp(:,1));
    temp_y = numel(temp(1,:));
    xx_size = floor(rad * 2 * pxx_per_cm);
    yy_size = floor(rad * 2 * pxy_per_cm);
    max_mean = 0;
    xx_c = ceil(temp_x/2);
    yy_c = ceil(temp_y/2);
    xx_s = xx_c - ceil(xx_size/2);
    xx_e = xx_c + ceil(xx_size/2);
    yy_s = yy_c - ceil(yy_size/2);
    yy_e = yy_c + ceil(yy_size/2);
    for ii = ceil(xx_size/2)+5:10:temp_x-floor(xx_size/2)-5
        for jj = ceil(yy_size/2)+5:10:temp_y-floor(yy_size/2)-5
            temp_mean = mean(mean(temp(ii-ceil(xx_size/2):ii+floor(xx_size/2), ...
                jj-ceil(yy_size/2):jj+floor(yy_size/2))));
            if temp_mean > max_mean
                max_mean = temp_mean;
                xx_c = ii;
                yy_c = jj;
            	xx_s = xx_c - ceil(xx_size/2);
                xx_e = xx_c + ceil(xx_size/2);
                yy_s = yy_c - ceil(yy_size/2);
                yy_e = yy_c + ceil(yy_size/2);
            end
        end
    end
    temp_xx_c = xx_c;
    temp_yy_c = yy_c;
    for ii = temp_xx_c-4:temp_xx_c+4
        for jj = temp_yy_c-4:temp_yy_c+4
            temp_mean = mean(mean(temp(ii-ceil(xx_size/2):ii+floor(xx_size/2), ...
                jj-ceil(yy_size/2):jj+floor(yy_size/2))));
            if temp_mean > max_mean
                max_mean = temp_mean;
                xx_c = ii;
                yy_c = jj;
            	xx_s = xx_c - ceil(xx_size/2);
                xx_e = xx_c + ceil(xx_size/2);
                yy_s = yy_c - ceil(yy_size/2);
                yy_e = yy_c + ceil(yy_size/2);
            end
        end
    end
end

