function [P_RR] = Toeplitz_NM(m_M,n_N,l_L,Lx,Ly,Lz,z_K,t_R)
% Toeplitz_MN computes the Toeplitz matrix associated with the 3D
% transient heat conduction problem. 

M = length(m_M);
N = length(n_N);
L = length(l_L);
K = length(z_K);
R = length(t_R);

%==========================================================================
% Create the 1D Green's function for the 1D finite slab heat conduction
% problem with a surface heat flux.
% =========================================================================

% The normalized 1D finite slab heat conduction problem with an inpulse
% response heat flux (dirac delta function) is the following: 
% PDE: u_t - u_zz = 0 
% Boundary conditions: 
% u_z = -q_0*dirac(t-0) at z = 0 
% u_z = 0               at z = 1
% Initial conditions: 
% u = 0 for all "z" at t = 0;

% We begin by computing the unit step response and then taking the time
% derivative to calculate the impulse response which is the 1D green's
% function

ll = l_L(2:end);
ee = exp( -((pi*ll).^2)*t_R');
dt_R = t_R(2) - t_R(1);
for k = 1:K
    pp{k} = cos(pi*z_K(k)*ll)./((ll*pi).^2);
    StepResponse = t_R + 2*(1-ee')*pp{k};
    % Calculate the Impulse response Kernel:
    ImpulseResponse = diff(StepResponse)/dt_R;
    G_R{k} = [ImpulseResponse;ImpulseResponse(end)];
    
    % The following variable is for diagnositcs purposes and can be plotted
    % to check the quality of the calculation
    G_RK(:,k) = G_R{k};
end

% Testing G
if 0
    figure;
    ribbon(t_R,G_RK,1)
    ylabel('time')
    figure; 
    ribbon(t_R,cumtrapz(G_RK,1))
    ylabel('time')
end

% =========================================================================
% Create Toeplitz matrix
% =========================================================================
% It turns out that in a 3D slab heat conduction problem, the solution can
% be expressed in terms of the 1D Green's function solution "G_R"
% multiplied by some factors that depend on the "m" and "n" modes used on
% the Fourier-cosine representation of the temperature fields.

for k = 1:K
    for m = 1:M
        for n = 1:N
                E_R{n,m} = exp(-((pi*m_M(m)).^2)*(t_R*(Lz/Lx)^2)).*exp(-((pi*n_N(n)).^2)*(t_R*(Lz/Ly)^2));
                K_R{n,m,k} = E_R{n,m}.*G_R{k};
                P_RR{n,m,k} = toeplitz(K_R{n,m,k},[K_R{n,m,k}(1);zeros(size(K_R{n,m,k})-[1,0])])*dt_R;
        end
    end
end
end

