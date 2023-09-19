function [w_k,a_k,b_k,m_k,P_k,J_k,estQuat] = GP_update(w_kk1,a_kk1,b_kk1,m_kk1,P_kk1,J_kk1,predQuat,Zk,paramGP,inv_P0_extent,basisAngleArray)

addpath('_common');
% Probabilities and detection
p_D_k = 0.95;

% Clutter rate
beta_FA = 20/800;


% Allocate some memory
nZk = size(Zk,2);
w_k = zeros(1,(nZk+1)*J_kk1);
a_k = zeros(1,(nZk+1)*J_kk1);
b_k = zeros(1,(nZk+1)*J_kk1);
m_k = zeros(size(m_kk1,1),(nZk+1)*J_kk1);
P_k = zeros(size(P_kk1,1),size(P_kk1,2),(nZk+1)*J_kk1);

% The expected number of generated measurements
gam(1:J_kk1) = a_kk1(1:J_kk1)./b_kk1(1:J_kk1);   %beta_D;

% Allocate memory for probability of detection
pD = p_D_k*ones(1,J_kk1);
for j = 1:J_kk1
    % Gaussian component corresponding to no detections
    w_k(1,j) = (1-(1-exp(-gam(j)))*pD(j))*w_kk1(1,j);
    a_k(1,j) = a_kk1(1,j);
    b_k(1,j) = b_kk1(1,j);
    m_k(:,j) = m_kk1(:,j);
    P_k(:,:,j) = P_kk1(:,:,j);
end

% l is an index that keeps track of how many measurement cells that have
% been used for updating the filter.
l = 0;
% Number of partitions of Zk
P = numel(Zk); 
% Allocate some memory
% wp = zeros(1,P);
log_wp = zeros(1,P);
W = zeros(1,P);

% Iterate over the partitions of the measurement set Zk
for p = 1:P
    % Number of cells in partition
    W(p) = numel(Zk(p).P);
    % Allocate memory for dw
    dw = zeros(1,W(p));
    % Iterate over the cells w of the partition p
    for w = 1:W(p)
        % Current set of measurements
        pwZ_k = Zk(p).P(w).W;
        % New measurement cell
        l = l+1;
        % Size of current cell
        absW = size(Zk(p).P(w).W,2);
        % Iterate over the prediction components
        for j = 1:J_kk1
            % likelihood function
            
            % update gamma
            a_k(1,J_kk1*l+j) = a_kk1(1,j)+absW;
            b_k(1,J_kk1*l+j) = b_kk1(1,j)+1;
            
            predState = m_kk1(:,j);
            predStateCov = P_kk1(:,:,j);
            
            numStates = size(m_kk1,1);
            predMeas = zeros(absW * 3, 1);                 	% of the form [x1; y1; z1;...; xn; yn; zn]
            measCov = zeros(absW*3);                         	% measurement noise covariance matrix
            linMeasMatrix = zeros(absW * 3, numStates);       % linearized measurement model
            
            for nz =1:absW          %Iterate over the measurements one by one
                [iMeasPredicted,  iMeasCovariance]  = compute_meas_prediction_and_covariance...
                    (pwZ_k(:,nz), predState, predQuat, paramGP, basisAngleArray, inv_P0_extent);
                % Log these variables to later utilize in the update mechanism
                predMeas(1+(nz-1)*3 : nz*3) = iMeasPredicted;
                measCov(((nz-1)*3+1):nz*3, ((nz-1)*3+1):nz*3) = iMeasCovariance;

                % Obtain linearized measurement model
                iLinMeasMatrix = compute_linearized_measurement_matrix(pwZ_k(:,nz), predState...
                    , predQuat, paramGP, basisAngleArray, inv_P0_extent, 1e-6);

                % Linearized measurement model for the current measurement
                % linMeasMatrix = [dh/dc, dh/dv = zero, dh/dq, dh/dw = zero, dh/dextent]
                linMeasMatrix(1+(nz-1)*3:nz*3, :) = iLinMeasMatrix;
            end

            % Put the measurements in a column the form: [x1; y1; z1;...; xn; yn; zn]
            tempArray = pwZ_k;
            curMeasArrayColumn = tempArray(:);

            % Realize measurement update
            S = linMeasMatrix*predStateCov*linMeasMatrix' + measCov;
            kalmanGain = predStateCov * linMeasMatrix'...
                / S;
            m_k(:,J_kk1*l+j) = predState + kalmanGain * (curMeasArrayColumn - predMeas);
            P_k(:,:,J_kk1*l+j) = (eye(numStates) - kalmanGain*linMeasMatrix) * predStateCov;
            P_k(:,:,J_kk1*l+j) = 0.5*(P_k(:,:,J_kk1*l+j)+P_k(:,:,J_kk1*l+j)');

            % Likelihood
            PL = mvnpdf(curMeasArrayColumn,predMeas,S);
            w_k(1,J_kk1*l+j) = exp(-gam(j)+log(pD(j))+absW*log(gam(j))+log(PL)-absW*log(beta_FA))*w_kk1(1,j); 
            w_k(1,J_kk1*l+j) = w_k(1,J_kk1*l+j)*updateGammaweightLog(a_k(1,J_kk1*l+j),a_kk1(1,j),b_k(1,J_kk1*l+j),b_kk1(1,j));
%           w_k(1,J_kk1*l+j) = w_k(1,J_kk1*l+j)*nbinpdf(absW,a_kk1(1,j),b_kk1(1,j)/(b_kk1(1,j)+1));
            
            % Quaternion Treatment
            estQuatDev = m_k(7:9,J_kk1*l+j);
            estQuat = apply_quat_deviation(estQuatDev, predQuat);
            m_k(7:9,J_kk1*l+j) = zeros(3,1);                         % Reset the orientation deviation

        end
        % Compute dw
        dw(w) = deltaFunction(absW,1)+sum(w_k(1,J_kk1*l+(1:J_kk1)));
        % Divide the weights by dw
        w_k(1,J_kk1*l+(1:J_kk1)) = w_k(1,J_kk1*l+(1:J_kk1))/dw(w);
    end
    % Replace inf with realmax in dw
    idx = dw==Inf;
    dw(idx) = realmax;
    % Take sum of log(dw) instead of prod(dw) to avoid numerical problems.
    log_wp(p) = sum(log(dw));
end
J_k = J_kk1*(l+1);

% Different normalisation since log(wp) are kept instead of wp.
% See http://en.wikipedia.org/wiki/List_of_logarithmic_identities
if P == 0
    wp = [];
elseif P == 1
    wp = 1;
else
    log_wp = log_wp-(log_wp(1)+log(1+sum(exp(log_wp(2:end)-log_wp(1)))));
    wp = exp(log_wp);
end


prev = J_kk1;
for p = 1:length(wp)
    w_k(1,prev+(1:J_kk1*W(p))) = w_k(1,prev+(1:J_kk1*W(p)))*wp(p);  
    prev = prev + J_kk1*W(p);
end

w_k = w_k(1,1:J_k);
a_k = a_k(1,1:J_k);
b_k = b_k(1,1:J_k);
m_k = m_k(:,1:J_k);
P_k = P_k(:,:,1:J_k);

end


function [measPredicted,   measCovariance] = compute_meas_prediction_and_covariance...
    (meas, state, quat, paramGP, basisAngleArray, inv_P0_extent)

meanGP = paramGP{7};    % mean of the GP

% Extract relevant information from the predicted state
center = state(1:3);
extent = state(13:end);

% Compute the following variables to exploit in measurement model
[p, H_f, covMeasBasis, angle_L] = compute_measurement_model_matrices...
    (meas, state, quat, paramGP, basisAngleArray, inv_P0_extent);

% Obtain the measurement prediction by the original nonlinear model
measPredicted = center + p * H_f * extent + p*meanGP*(1-H_f*ones(size(extent,1),1));

% Definition: iCovMeas = k(uk,uk), uk: argument of the current measurement
covGP = compute_GP_covariance(angle_L, angle_L, paramGP);
% Definition: iR = k(uk,uk) - k(uk, uf) * inv(K(uf,uf)) * k(uf, uk)
iR_f = covGP - covMeasBasis * inv_P0_extent * covMeasBasis';
% Obtain the covariance of the measurement
stdMeasGP = paramGP{6};
measCovariance = stdMeasGP^2 * eye(3) + p * iR_f * p';
end

function [ linearMeasMat ] = compute_linearized_measurement_matrix( meas, state...
    , quat, paramGP, basisAngleArray, inv_P0_extent, eps)

meanGP = paramGP{7};    % mean of the GP

%% Obtain the measurement prediction by the original nonlinear meas model
[p, H_f, ~, ~] = compute_measurement_model_matrices(meas, state, quat, paramGP...
    , basisAngleArray, inv_P0_extent);
center = state(1:3);
extent = state(13:end);
measPredicted = center + p * H_f * extent + p*meanGP*(1-H_f*ones(size(extent,1),1));

%% Calculate partial derivative wrt target extent
H_extent = p * H_f;     % This partial derivative is taken analytically.

%% Calculate partial derivative wrt linear positions
H_center = zeros(3,3);
for i = 1:3
    stateIncremented = state;   % Initialize the vector
    stateIncremented(i) = stateIncremented(i) + eps;    % increment the state
    centerIncremented = stateIncremented(1:3);          % incremented center
    
    % Compute the output of original measurement model for the incremented state vector
    [p, H_f, ~, ~] = compute_measurement_model_matrices(meas, stateIncremented...
        , quat, paramGP, basisAngleArray, inv_P0_extent);
    measPredictedForIncState = centerIncremented + p * H_f * extent + p*meanGP*(1-H_f*ones(size(extent,1),1));
    
    % Compute numeric derivative
    H_center(:, i) = (measPredictedForIncState - measPredicted) / eps;
end

%% Calculate partial derivative wrt linear velocities
H_linVel = zeros(3,3); % Measurement model does not depend on linear velocities

%% Calculate partial derivative wrt quaternion deviations
H_quatDev = zeros(3,3);
for i = (7:9)
    stateIncremented = state;   % Initialize the vector
    stateIncremented(i) = stateIncremented(i) + eps;
    
    % Compute the output of original measurement model for the incremented state vector
    [p, H_f, ~, ~] = compute_measurement_model_matrices(meas, stateIncremented...
        , quat, paramGP, basisAngleArray, inv_P0_extent);
    measPredictedForIncState = center + p * H_f * extent + p*meanGP*(1-H_f*ones(size(extent,1),1));
    
    % Compute numeric derivative
    H_quatDev(:, i-6) = (measPredictedForIncState - measPredicted) / eps;
end

%% Calculate partial derivative wrt angular velocities
H_angVel = zeros(3,3); % Measurement model does not depend on linear velocities

%% Form the linearized measurement matrix
linearMeasMat = [H_center H_linVel H_quatDev H_angVel H_extent];

end

function [diffUnitVector_G, H_f, covMeasBasis, measAngle_L] = compute_measurement_model_matrices...
    (meas, state, quatPrev, paramGP, basisAngleArray, inv_P0_extent)

% Extract relevant information from the state
center = state(1:3);
quatDev = state(7:9);
quat = apply_quat_deviation(quatDev, quatPrev); % Apply predicted orientation deviation

% Compute the following variables to exploit in measurement model
diffVector_G = meas - center;   % the vector from the center to the measurement
diffVectorMag = norm(diffVector_G);
if diffVectorMag == 0
    diffVectorMag = eps;
end
diffUnitVector_G = diffVector_G / diffVectorMag;

% Express the diffVectorGlobal in Local frame
R_from_G_to_L = rotation_matrix_from_global_to_local(quat);
diffVector_L = R_from_G_to_L * diffUnitVector_G;

% Find the Spherical angle in Local frame corresponding to this vector
[azi_L, elv_L, ~] = cart2sph(diffVector_L(1), diffVector_L(2), diffVector_L(3));
azi_L = mod(azi_L, 2*pi);       % to keep it between 0 and 2*pi
angle_L = [azi_L elv_L];        % angle in the local coordinate frame

% Compute the covariance matrix component from the GP model
covMeasBasis = compute_GP_covariance(angle_L, basisAngleArray, paramGP);

H_f = covMeasBasis * inv_P0_extent; % GP model relating the extent and the radial
% function value evaluated at the current measurement angle

measAngle_L = angle_L;  % the angle of the measurement in local coordinate frame

end

function [F, Q] = compute_rotational_dynamics(w, std, T)
% Inputs:
%              w:        Angular rate
%              std:     Standard deviation of the angular velocity
%              T:         Sampling time

% Dummy variables
wNorm = norm(w);
if wNorm == 0
    wNorm = 1e-3;
end

S = skew_symmetric_matrix(-w);
c = cos(1/2*T*wNorm);
s = sin(1/2*T*wNorm);
expMat = eye(3,3) + s/wNorm*S + (1-c)/wNorm^2*S^2;

% Construct the state transition matrix
F = [expMat  T*expMat ; zeros(3,3)  eye(3,3)];

% Construct the process noise covariance matrix
G11 = T*eye(3,3) + 2/wNorm^2*(1-c)*S + 1/wNorm^2*(T-2/wNorm*s)*S^2;
G12 = 1/2*T^2*eye(3,3) + 1/wNorm^2*(4/wNorm*s-2*T*c)*S + 1/wNorm^2*(1/2*T^2+2/wNorm*T*s+4/wNorm^2*(c-1))*S^2;
G21 = zeros(3,3);
G22 = T * eye(3,3);
G = [G11 G12; G21 G22];

B = G * [zeros(3,3); eye(3,3)];

% cov = std^2 * diag([0 0 1]);
cov = std^2 * eye(3,3);
Q = B * cov * B';
end


function [qOut] = apply_quat_deviation(a, qIn)
qa = 1/ sqrt(4+norm(a)) * [a; 2];   % It depends on the Rodrigues parametrization
qOut = quat_product(qa, qIn);

qOut = qOut/ norm(qOut);    % Normalize due to numeric inprecisions
end

function [out] = quat_product(q1, q2)
q1V = q1(1:3);
q1S = q1(4);
q2V = q2(1:3);
q2S = q2(4);

out = [q1S*q2V + q2S*q1V - cross(q1V,q2V);...
    q1S*q2S - q1V'*q2V];
end

function [df] = deltaFunction(x,y)
% this is just an implementation of Diracs delta...
if x == y
    df = 1;
else
    df = 0;
end
end