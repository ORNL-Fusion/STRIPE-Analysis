function [a,b,L,rotMatrix] = ConstrainView(m,rA,qA,rB,qB,L)
% INPUTS:
% -------------------------------------------------------------
% m: scalar given as product of Omega*focalLength
% rA: 2x1 vector, image plane conjugate point
% qA: 2x1 vector, Object conjugate point
% rB: 2x1 vector, image plane conjugate point
% qB: 2x1 vector, Object conjugate point
% L : 3x1 vector, offset position vector
%
% OUTPUT:
% -------------------------------------------------------------
% a: scalar, rotation angle
% b: scalar, rotation angle
% L: 3x1 vector, updated offset position vector
% rotMatrix: Rotation matrix of Constrained view
%
% Rotation matrix defined as: [e] = Ryx*[s]
% where [e] is camera's referene frame
% [s] is object's reference frame
%
% a: Rotation about yy to rotate the objects's frame to the camera's frame
% b: Rotation about xx to rotate the objects's frame to the camera's frame

% Conjugate points:
% =============================================================
% Image plane:
r_star.A = rA;
r_star.B = rB;

% Object:
q_star.A = qA;            
q_star.B = qB;

% Point A:
u.A  = q_star.A(1:2);
v.A  = r_star.A(1:2);
Q3.A = q_star.A(3) - L(3);

% Point B:
u.B  = q_star.B(1:2);
v.B  = r_star.B(1:2);
Q3.B = q_star.B(3) - L(3);

% Create rotation matrix:
% =============================================================
% Rotation matrix defined as: [e] = Ryx*[s]
% where [e] is camera's referene frame
% [s] is object's reference frame

% a: Rotation about yy to rotate the objects's frame to the camera's frame
% b: Rotation about xx to rotate the objects's frame to the camera's frame

R{1,1} = @(a,b) +cos(a);
R{1,2} = @(a,b) +sin(a).*sin(b);
R{1,3} = @(a,b) +sin(a).*cos(b);
R{2,1} = @(a,b) 0;
R{2,2} = @(a,b) +cos(b);
R{2,3} = @(a,b) -sin(b);
R{3,1} = @(a,b) -sin(a);
R{3,2} = @(a,b) +cos(a).*sin(b);
R{3,3} = @(a,b) +cos(a).*cos(b);

R1 = @(a,b) [R{1,1}(a,b); R{2,1}(a,b); R{3,1}(a,b)];
R2 = @(a,b) [R{1,2}(a,b); R{2,2}(a,b); R{3,2}(a,b)];
R3 = @(a,b) [R{1,3}(a,b); R{2,3}(a,b); R{3,3}(a,b)];
RR = @(a,b) [R1(a,b), R2(a,b), R3(a,b)];

% M matrix:
% =============================================================
fld = ['A','B'];
for pp  = 1:numel(fld) 
    for rr = 1:numel(u.A)
        for cc = 1:3
            k = fld(pp);
            r_i = r_star.(k)(rr);

            % M_ij       =        r_i*R_j3         + m*R_ji
            M{rr,cc}.(k) = @(a,b) r_i*R{cc,3}(a,b) + m*R{cc,rr}(a,b);

        end
    end
end

% L1 and L2 linear equation:
% =============================================================
for pp = 1:numel(fld)
    % Select conjugate point:
    k = fld(pp);
    % Adjugate matrix of G:
    adjG.(k) = @(a,b) [+M{2,2}.(k)(a,b), -M{1,2}.(k)(a,b);...
                       -M{2,1}.(k)(a,b), +M{1,1}.(k)(a,b)];   
    % Determinant of G:
    detG.(k) = @(a,b)  M{1,1}.(k)(a,b)*M{2,2}.(k)(a,b) -  M{2,1}.(k)(a,b)*M{1,2}.(k)(a,b);                    
    % Factor K
    K.(k) = @(a,b) Q3.(k)./detG.(k)(a,b);
    % Equation for offset vector:
    t.(k) = @(a,b) u.(k) + K.(k)(a,b)*adjG.(k)(a,b)*[M{1,3}.(k)(a,b),M{2,3}.(k)(a,b)]';
end

% Assemble F:
% =============================================================
F = @(a,b) (u.A-u.B) + ...
    K.A(a,b)*adjG.A(a,b)*[M{1,3}.A(a,b);M{2,3}.A(a,b)] - ...
    K.B(a,b)*adjG.B(a,b)*[M{1,3}.B(a,b);M{2,3}.B(a,b)];

% Put function in form needed for non-linear solver:
Fsolve = @(x) F(x(1),x(2));

% Solve non-linear equation:
% =============================================================
tic
[x, ~, ~, exitflag, output, ~] = newtonraphson(Fsolve, [4,-13]*pi/180);
toc

% Assign solution to angles and offset vector:
% =============================================================
% View angles:
a = x(1);
b = x(2);

% rotation matrix:
rotMatrix = RR(a,b);

% Offset vector:
L(1) = [1,0]*t.A(a,b);
L(2) = [0,1]*t.A(a,b);

end

