function [xx_vec,yy_vec] = Pix2XY(xx_c,yy_c,pix_conv,dim)
%%%returns X and Y vectors from simple pixel to geoemtry ratio, centers
%%%coordinates at xx_c and yy_c
% Inputs
%     xx_c       x location for center of plasma
%     yy_c       y location for center of plasma
%     pix_conv   pixels per cm conversion
%     dim        dimension vector for temperatuer array
% Outputs
%     xx_vec     x vector in real coordinates
%     yy_vec     y vector in real coordinates
xx_vec = linspace(0,dim(1)./pix_conv,dim(1));
yy_vec = linspace(0,dim(2)./pix_conv,dim(2));

xx_vec = (xx_vec-xx_c./pix_conv)/100.;
yy_vec = (yy_vec-yy_c./pix_conv)/100.;
end

