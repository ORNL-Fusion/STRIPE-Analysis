function [params,inputField,outputField] = Slab3dTransientThermalModel(params,inputField)
% Define characteristic scales::
% =========================================================================
params = CharacteristicScales(params);

% Normalize and populate input field:
% =========================================================================
inputField = NormalizeAndPopulateInputField(params,inputField);

% Calculate forward Fourier-cosine transform of input field:
% =========================================================================
inputField = ApplyFourierCosineTransform(params,inputField,'forward');

% Calculate Toeplitz matrix:
% =========================================================================
P_RRMN = CalculateToeplitzMatrix(params,inputField);

% Create output field:
% =========================================================================
outputField = CreateOutputField(inputField);

% Solve heat transfer in Fourier space:
% =========================================================================
switch params.solType
    case 'forward'
         outputField = FourierSpaceForwardSol(params,P_RRMN,inputField,outputField);
    case 'inverse' 
         outputField = FourierSpaceInverseSol(params,P_RRMN,inputField,outputField);
end

% Calculate inverse Fourier-cosine transform of output field:
% =========================================================================
outputField = ApplyFourierCosineTransform(params,outputField,'inverse');

end

% #########################################################################
% Supporting function 1:
% #########################################################################
function params = CharacteristicScales(params)

% Extract basic data from params:
% ===============================
material = params.material;
kt       = params.kt;
rho      = params.rho; 
cp       = params.cp;
f_star   = params.f_star;
f0       = params.f_0;
Lx       = params.Lx;
Ly       = params.Ly;
Lz       = params.Lz;

% Thermal diffusivity:
% ====================
params.d = kt/(rho*cp);

% Characteristic scales:
% ==========================
% time:
params.t_star = (Lz^2)/params.d;

% Heat flux and temperature:
switch params.solType
    case {'forward'}
        params.q_star = f_star;
        params.T_star = f_star*(Lz/kt);
    case {'inverse'}
        params.q_star = f_star/(Lz/kt);
        params.T_star = f_star;
end

% Normalized position where temperarature field is located:
% =========================================================
params.z_u = params.z_T/Lz;

end

% #########################################################################
% Supporting function 2:
% #########################################################################
function inputField = NormalizeAndPopulateInputField(params,inputField)

% Extract values from params:
% ===========================
I      = params.I;
J      = params.J;
R      = params.R;
M      = params.M;
N      = params.N;
Lx     = params.Lx;
Ly     = params.Ly;
Lz     = params.Lz;
t_star = params.t_star;

% Time-dependent field:
% =====================
inputField.f_IJR  = (inputField.f - params.f_0)/params.f_star;
inputField.df_IJR = zeros(size(inputField.f_IJR));

% Time-dependent fourier cosine field:
% ====================================
inputField.F_MNR = zeros(M,N,R);

% Normalized spatial coordinates:
% ===============================
inputField.x_J = inputField.x/Lx + 0.5;
inputField.y_I = inputField.y/Ly + 0.5;

% Temporal coordinate:
% ====================
inputField.t_R  = inputField.t/t_star;
inputField.dt_R = inputField.t_R(2) - inputField.t_R(1); 

% Fourier-cosine eigenfunction:
% =============================
m_M = [0:(M-1)]';
n_N = [0:(N-1)]';
inputField.phi_JM = cos(pi*inputField.x_J*m_M');
inputField.phi_IN = cos(pi*inputField.y_I*n_N');

end

% #########################################################################
% Supporting function 3:
% #########################################################################
function field = ApplyFourierCosineTransform(params,field,operatorType)

% Extract eigenfunction:
% ======================
phi_IN = field.phi_IN;
phi_JM = field.phi_JM;

% Extract number of time steps:
% ============================
R = params.R;

switch operatorType
    
    case {'forward'}
        % Extract field:
        % ==============
        f_IJR = field.f_IJR;
        
        % grid size:
        % ==========
        dx = mean(diff(field.x_J));
        dy = mean(diff(field.y_I));

        % Define Kronecker delta functions:
        % =================================
        d_mm = 0.5*ones(params.M,1)';
        d_mm(1) = 1;

        d_nn = 0.5*ones(params.N,1)';
        d_nn(1) = 1;

        D_MM = diag(1./d_mm);
        D_NN = diag(1./d_nn);

        % Define normalized eigenfunction:
        % ================================
        psi_IN = phi_IN*D_NN;
        psi_JM = phi_JM*D_MM;

        % Apply Fourier-cosine operation: f_IJR -> F_MNR
        % ==============================================
        for r = 1:R
            % Fourier-cosine transform:
            A_MN = psi_JM'*f_IJR(:,:,r)'*psi_IN*dx*dy;

            % Compute the error in approximation with Fourier cosine:
            f_IJR_approx = phi_IN*A_MN'*phi_JM';
            df_IJR = (f_IJR_approx - f_IJR(:,:,r))/max(max(f_IJR(:,:,r)));

            % Assign value to output variable:
            field.F_MNR(:,:,r) = A_MN;
            field.df_IJR(:,:,r) = df_IJR;
        end
        
    case {'inverse'}

        % Apply inverse Fourier-cosine operation:  F_MNR -> f_IJR
        % =======================================================
        for r = 1:R
            F_MN  = field.F_MNR(:,:,r);
            field.f_IJR(:,:,r) = phi_IN*F_MN'*phi_JM';
        end
        
end

end

% #########################################################################
% Supporting function 4:
% #########################################################################
function [P_RRMN] = CalculateToeplitzMatrix(params,inputField)
% Toeplitz_MN computes the Toeplitz matrix associated with the 3D
% transient heat conduction problem. 

% Extract data from inputs:
% =========================
M    = params.M;
N    = params.N;
L    = params.L;
K    = 1;
R    = params.R;
Lx   = params.Lx;
Ly   = params.Ly;
Lz   = params.Lz;
t_R  = inputField.t_R;
dt_R = inputField.dt_R;
z_u  = params.z_u;

% Derived quantities:
% ===================
m_M = [0:(M-1)]';
n_N = [0:(N-1)]';
l_L = [0:(L-1)]';

%==========================================================================
% Create the 1D Green's function for the 1D finite slab heat conduction
% problem with d surface heat flux.
% =========================================================================

% The normalized 1D finite slab heat conduction problem with an inpulse
% response heat flux (dirac delta function) is the following: 
% PDE: u_t - u_zz = 0 
% Boundary conditions: 
% u_z = -q_0*dirac(t-0) at z = 0 
% u_z = 0               at z = 1
% Initial conditions: 
% u = 0 for all "z" at t = 0;


% Calculate fundamental solution:
% ===============================
% We begin by computing the unit step response and then taking the time
% derivative to calculate the impulse response which is the 1D green's
% function

% Eigenvalues:
ll = l_L(2:end);
lambda_l = pi*ll;

% Eigenfunction:
p_l = cos(lambda_l*z_u); 
cosTerm = p_l./(lambda_l.^2);

% Step response:
expTerm = exp( -(lambda_l.^2)*t_R');
StepResponse = t_R + 2*(1-expTerm')*cosTerm;

% Impulse response:
ImpulseResponse = diff(StepResponse)/dt_R;

% Assemble Green's function:
G_R = [ImpulseResponse;ImpulseResponse(end)];

% =========================================================================
% Create Toeplitz matrix
% =========================================================================
% It turns out that in d 3D slab heat conduction problem, the solution can
% be expressed in terms of the 1D Green's function solution "G_R"
% multiplied by some factors that depend on the "m" and "n" modes used on
% the Fourier-cosine representation of the temperature fields.

for m = 1:M
    for n = 1:N
            % Eigenvalues:
            lambda_m = pi*m_M(m)*Lz/Lx;
            lambda_n = pi*n_N(n)*Lz/Ly;

            % Kernel:
            E_R{m,n} = exp(-(lambda_m.^2)*t_R).*exp(-(lambda_n.^2)*t_R);
            K_R{m,n} = E_R{m,n}.*G_R;

            % Toeplitz matrix:
            P_RRMN(:,:,m,n) = toeplitz(K_R{m,n},[K_R{m,n}(1);zeros(size(K_R{m,n})-[1,0])]);
    end
end

end

% #########################################################################
% Supporting function 5:
% #########################################################################
function outputField = CreateOutputField(inputField)

outputField = inputField;

outputField.f      = [];
outputField.f_IJR  = zeros(size(outputField.f_IJR));
outputField.df_IJR = [];
outputField.F_MNR  = zeros(size(outputField.F_MNR));

end

% #########################################################################
% Supporting function 6:
% #########################################################################
function outputField = FourierSpaceForwardSol(params,P_RRMN,inputField,outputField)

M    = params.M;
N    = params.N;
R    = params.R;
dt_R = inputField.dt_R;

for m = 1:M
    for n = 1:N
        % Heat flux Fourier component:
        H_RMN = permute(inputField.F_MNR(m,n,:),[3,1,2]); 
                
        % Temperature fourier component (Convolution operation):
        U_RMN = P_RRMN(:,:,m,n)*H_RMN*dt_R;
        
        % Assign solution to output:
        outputField.F_MNR(m,n,:) = ipermute(U_RMN,[3,1,2]);
    end
end

end

% #########################################################################
% Supporting function 7:
% #########################################################################
function outputField = FourierSpaceInverseSol(params,P_RRMN,inputField,outputField)
M    = params.M;
N    = params.N;
R    = params.R;
dt_R = inputField.dt_R;

% Number of iterations for the inverse method:
num_iterations   = 300; 
% num_iterations = 200; 

% Allocate memory to heat flux variable:
H_RMN    = cell(num_iterations + 1,1);

for m = 1:M
    disp(['m: ',num2str(m)])
    for n = 1:N
        % Temperature fourier component:
        U_RMN = permute(inputField.F_MNR(m,n,:),[3,1,2]); 

        % Initial guess for heat flux in fourier space:
        H_RMN{1} = diff(U_RMN/dt_R);
        H_RMN{1} = [H_RMN{1};H_RMN{1}(end)];
%         H_RMN{1} = zeros(size(U_RMN));

        % Apply inverse method:
        for ii = 1:num_iterations
            
            % Get convolution operator:
            P = P_RRMN(:,:,m,n);
            
            % Inferred temperature:
            U{ii} = P*H_RMN{ii}*dt_R;
            
            % Input temperature data:
            Z = U_RMN;
            
            % Compute residual:
            r{ii}   = U{ii} - Z;
            
            % Compute value of functional:
            J(ii)   = 0.5*r{ii}'*r{ii};
            
            % Calculate search direction based on Conjugate Gradient Method
            % (CGM):
            if ii == 1
                % Gradient based in steepest descent:
                dJ{ii}  = P'*r{ii};
            else
                % Gradient based on CGM from the following reference:
                
                % J. Wang, A. J. Silva Neto, F. D. Moura Neto, and J. Su,
                % “Function estimation with Alifanov’s iterative
                % regularization method in linear and nonlinear heat
                % conduction problems,” Applied Mathematical Modelling,
                % vol. 26, no. 11, pp. 1093–1111, Nov. 2002. Eq 12b

                % Calculate current steepest descent:
                dJ{ii}  = P'*r{ii};
                
                % Calculate step size:
                y(ii) = dJ{ii}'*dJ{ii}/(dJ{ii-1}'*dJ{ii-1});
                
                % Improve the search direction:
                dJ{ii} = dJ{ii} + y(ii)*dJ{ii-1};
            end

            % Compute new value of heat flux based on best search direction:
            Y{ii} = P*P'*r{ii};
            d(ii) = r{ii}'*Y{ii}/(Y{ii}'*Y{ii});
            H_RMN{ii+1} = H_RMN{ii} - d(ii)*dJ{ii};
            
        end 
        
        % Assign solution to outputField:
        outputField.F_MNR(m,n,:) = ipermute(H_RMN{num_iterations+1},[3,1,2]);
    end
end


end

