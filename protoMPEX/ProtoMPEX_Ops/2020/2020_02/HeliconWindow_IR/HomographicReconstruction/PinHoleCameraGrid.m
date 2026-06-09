function [r] = PinHoleCameraGrid(C,R,L,q)
% PinHoleCamera: 
% Takes a 3x1 position vector q and maps it to the image plane of a pinhole
% camera given a rotation matrix and relative camera position.
%
% Inputs:
% C.f: Camera's focal length
% C.Omega: Reflection factor, 1: image is inverted just as in a pin-hole
% camera. -1: image is inverted relative to normal output of a pin-hole
% camera.
%
% R: 3x3 Rotation matrix. R is defined as the rotation matrix needed to
% rotate the object's coordinate system to the camera's coordinate system.
% 
% L: 3x1 Vector that defines the location of the camera's pin-hole
% relative to object's datum in terms of the object's coordinate system.
% 
% q: 3x1 Vector that defines an arbitraty point on relative to object's
% datum in terms of the object's coordinate system.
% 
% Outputs:
% r: 3x1 vector that defines the position of q mapped into the image plane
% of the pinhole camera. By construction r(3) = C.f

% Camera geometry information:
f = C.f;
Omega =C.Omega;

% Express Q relative to the camera's pinhole camera and orientation:
Q = R'*(q-L);

% Calculate projection on image plane:
r = -f*Omega*Q/Q(3);
end

