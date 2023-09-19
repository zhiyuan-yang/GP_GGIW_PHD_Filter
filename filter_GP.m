function [est_x_GP,est_X_GP,est_n_GP,OSPA_GP,RMSEx,RMSEv,IoU] = filter_GP(meas,groundTruth)
    
addpath("_common\");
addpath("func_GP\");
addpath("spheretri-master\");

% Simulation Parameters
expDuration = 30;   % (in seconds)
T = 0.1;            % sampling time (in seconds)
eps = 1e-6;         % a small scalar
numInstants = ceil(expDuration/ T);


% Gaussian Process (GP) Parameters
minNumBasisAngles = 200;        % minimum number of basis angles on which the extent is maintained
meanGP = 0;                     % mean of the GP
stdPriorGP = 1;                 % prior standard deviation of the GP
stdRadiusGP = 0.2;              % standard deviation of the radius
scaleLengthGP = pi/8;           % lengthscale
kernelTypeGP = 1;               % 1:Squared exponential , 2: Matern3_2, 3: Matern5_2
isObjSymmetric = 1;             % it is used to adjust the kernel function accordingly
stdMeasGP = 0.1;                % standard deviation of the measurements (used in the GP model)
paramGP = {kernelTypeGP, stdPriorGP, stdRadiusGP, scaleLengthGP, isObjSymmetric, stdMeasGP, meanGP};


% Determine the Basis Angles
[basisVertices, ~] = spheretri(minNumBasisAngles);      % produces evenly spaced points on the sphere
[azimuthBasisArray, elevationBasisArray, ~] = cart2sph(basisVertices(:,1), basisVertices(:,2)...
    , basisVertices(:,3));
azimuthBasisArray = mod(azimuthBasisArray, 2*pi);       % to make it consistent with the convention
numBasisAngles = size(basisVertices, 1);
% Arrange the basis array so that azimuth and elevation are in ascending order
basisAngleArray = [azimuthBasisArray elevationBasisArray];
basisAngleArray = sortrows(basisAngleArray, 2);
basisAngleArray = sortrows(basisAngleArray, 1);

% Gaussian P0
P0_extent = compute_GP_covariance(basisAngleArray, basisAngleArray, paramGP);
P0_extent = P0_extent + eps*eye(size(basisAngleArray,1));       % to prevent numerical errors thrown in matrix inversion
chol_P0_extent = chol(P0_extent);
inv_chol_P0_extent = inv(chol_P0_extent);
inv_P0_extent = inv_chol_P0_extent * inv_chol_P0_extent';

% initial value
J =2;
w = [0.1 0.1];
b = [5 5];                                                                  
a = [100 100]; 

dev0 = [0 0 0];
q0 = [0,0,0,1]';
omega0 = [0 0 0];
f0 = meanGP * ones(1, numBasisAngles);          % initial extent (is initialized regarding the GP model)
m = [0 0 1 0 1 0 dev0 omega0 f0;
     15 0 1 0.5 0 0 dev0 omega0 f0]';
% m = [0 0 1 1 0 0 dev0 omega0 f0;
%      0 -20 1 2 0 0 dev0 omega0 f0]';
Quat = q0;

Pc = 0.1 * eye(3);
Pv = 0.1 * eye(3);
P_dev = 1e-5 * eye(3,3);
P_omega = 1e-0 * eye(3,3);
Pe  = compute_GP_covariance(basisAngleArray, basisAngleArray, paramGP);
P = blkdiag(Pc,Pv,P_dev,P_omega,Pe);
P = cat(3,P,P);


% partition
max_dist_prob = 0.80;
min_dist_prob = 0.40;
obs.r = 2;
max_distance = obs.r*chi2inv(max_dist_prob,3); % 置信度为0.8时的距离
min_distance = obs.r*chi2inv(min_dist_prob,3); % 置信度为0.4时的距离
N_low = 2;

% merge&reduction
threshold = 1e-3;
Jmax = 100;
U=30;

est_x_GP = cell(numInstants,1);
est_X_GP = cell(numInstants,1);
est_n_GP = zeros(numInstants,1);
OSPA_GP = zeros(numInstants,1);
RMSEx = zeros(1,length(meas));
RMSEv = zeros(1,length(meas));
estQuat = zeros(4,length(meas));
IoU = zeros(1,length(meas));

for k=1:numInstants
%     fprintf("Current Time Step: %d\n",k)
    %Predict
    [w,a,b,m,P,J] = GP_predict(w,a,b,m,P,J,numBasisAngles);

    %当前帧量测
    [Zkp,~] = distance_partition_with_sub_partition(meas{k},max_distance,min_distance,N_low);

    %updte
    [w,a,b,m,P,J,Quat] = GP_update(w,a,b,m,P,J,Quat,Zkp,paramGP,inv_P0_extent,basisAngleArray);

    %假设删减
    [w,a,b,m,P,J] = GP_reduction(w,a,b,m,P,J,threshold,U,Jmax);

    %状态提取
    [~,est_x_GP{k},est_X_GP{k}] = GP_stateExtraction(w,a,b,m,P,J);
    est_n_GP(k) = size(est_x_GP{k},2);

    OSPA_GP(k) = ospa_dist(groundTruth{k}(1:3,:), est_x_GP{k}(1:3,:),15,1);

    RMSEx(1,k) = cal_RMSE(groundTruth{k}(1:3,:),est_x_GP{k}(1:3,:));
    RMSEv(1,k) = cal_RMSE(groundTruth{k}(4:6,:),est_x_GP{k}(4:6,:));
    
    %IoU(k) = cal_IoU_GP(groundTruth{k},est_x_GP{k},basisAngleArray);
end

end