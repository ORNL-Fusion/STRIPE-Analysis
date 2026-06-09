function [xi,yi,zi] = CreateImagePlaneGrid(qx,qy,qz,rotMatrix,L,f,Omega)
% INPUTS:
% qx,qy,qz: [nTheta x nZ] Grid describing 3D object's surface
% focal length
% Reflection factor "Omega"
% L offset vector
% Rotation matrix from constrained view
% 
% OUTPUTS:
% xi, yi, zi: [nTheta x nZ] Grid describing 3D object on image
% plane

% Create 2D grid on image plane:
% =====================================================================
C.f = f;
C.Omega = Omega;
for cc = 1:size(qx,2)
     for rr = 1:size(qx,1)
         qq = [qx(rr,cc),qy(rr,cc),qz(rr,cc)]';
         ri = PinHoleCamera(C,rotMatrix,L,qq);
         xi(rr,cc) = ri(1);
         yi(rr,cc) = ri(2);
         zi(rr,cc) = ri(3);
     end
end

end

