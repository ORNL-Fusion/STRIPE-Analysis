function [g] = rotate_IR_data(f,angle)

for ii = 1:size(f,3)
    [g(:,:,ii),~] = rotate_image(angle,f(:,:,ii));
end

end

