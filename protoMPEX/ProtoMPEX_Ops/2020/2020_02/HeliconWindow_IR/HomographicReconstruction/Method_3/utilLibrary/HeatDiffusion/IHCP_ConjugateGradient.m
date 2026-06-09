function [q,u,J] = IHCP_ConjugateGradient(P,q,z,dtt,Ni,ExitCond,Val)
% INVERSE HEAT CONDUCTION PROBLEM BASED ON CONJUGATE GRADRIENT MINIMIZATION
% METHOD.

% INPUT DATA:
% P: convolution operator, NxN Toeplitz matrix
% q: initial guess at heat flux, structure, column vector
% z: experimental data, column vector
% dtt: time step size
% Ni: number of iteration
% ExitCond: 1,2. Determines how the loop is terminated
% Val: value associated with ExitCond

% OUTPUT DATA:
% u: calculated temperature, column vector
% J: residual, column vector
% q: minimized heat flux, structure, column vector

for i = 1:Ni
    % Compute temperature fields
    u{i}   = P*q{i}*dtt;
    % Compute residual
    r{i}   = (u{i}-z);
    % Compute value of functional
    J(i)   = 0.5*r{i}'*r{i};
    
    % Compute functional derivative and steepest descent
    if i == 1
        dJ{i}  = P'*r{i};
    else
        % Update the steepest descent and calculate the conjugate gradient
        % on order to speed up the search for the minimum.
        % From reference:
        % J. Wang, A. J. Silva Neto, F. D. Moura Neto, and J. Su,
        % “Function estimation with Alifanov’s iterative regularization
        % method in linear and nonlinear heat conduction problems,”
        % Applied Mathematical Modelling, vol. 26, no. 11, pp. 1093–1111, Nov. 2002.
        % Eq 12b
        
        % Calculate current steepest descent
        dJ{i}  = P'*r{i};
        % Calculate step size
        y(i) = dJ{i}'*dJ{i}/(dJ{i-1}'*dJ{i-1});
        % Improve the search direction
        dJ{i}  = dJ{i} + y(i)*dJ{i-1};
    end
    
    % Having calculated the direction of search dJ{i}, compute new value of heat flux:
    Y{i}   = P*P'*r{i};
    d(i)   = r{i}'*Y{i}/(Y{i}'*Y{i});
    q{i+1} = q{i} - d(i)*dJ{i};
    
%     Exit condition
    switch ExitCond
        case 1
            if abs(J(i))< Val^2 ;% Val = dnoise
                break
            end
        case 2
            if (J(i)/J(1)) < Val % 1e-5
                break
            end
        case 3
    end
end
end

