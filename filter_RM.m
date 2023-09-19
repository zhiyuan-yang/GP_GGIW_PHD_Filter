function [est_x_RM,est_X_RM,est_n_RM,OSPA_RM,RMSEx,RMSEv,IoU] = filter_RM(meas,groundTruth)

addpath("_common\");
addpath("func_RM\");

est_x_RM = cell(length(meas),1);
est_X_RM = cell(length(meas),1);
est_n_RM = zeros(length(meas),1);
OSPA_RM = zeros(length(meas),1);
RMSEx = zeros(1,length(meas));
RMSEv = zeros(1,length(meas));
IoU = zeros(1,length(meas));

%initial
w = [0.1 0.1];
a = [100 100];
b = [5 5];
m = [0,0,1,0,1,0;
    15,0,1,0.5,0,0]';
% m = [0 0 1 1 0 0;
%      0 -20 1 2 0 0]';
P = 0.1*eye(2);
P = cat(3,P,P);
nu = [10,10];
V = eye(3);
V = cat(3,V,V);
J = 2;

% partition
max_dist_prob = 0.80;
min_dist_prob = 0.40;
obs.r = 2;
max_distance = obs.r*chi2inv(max_dist_prob,3); % 置信度为0.8时的距离
min_distance = obs.r*chi2inv(min_dist_prob,3); % 置信度为0.4时的距离
N_low = 2;

% merge&reduction
threshold = 1e-5;
Jmax = 100;
U_G = 30;
U_N = 30;
U_IW = 30;


for k=1:length(meas)
    [w,a,b,m,P,nu,V,J] = GGIW_PHD_predict(w,a,b,m,P,nu,V,J);

    [Zkp,~] = distance_partition_with_sub_partition(meas{k},max_distance,min_distance,N_low);

    %updte
    [w,a,b,m,P,nu,V,J] = GGIW_PHD_update(w,a,b,m,P,nu,V,J,Zkp);

    %假设删减
    [w,a,b,m,P,nu,V,J] = GGIW_reduction(w,a,b,m,P,nu,V,J,threshold,Jmax,U_G,U_N,U_IW);

    %状态提取
    [~,est_x_RM{k},~,est_X_RM{k}] = GGIW_stateExtraction(w,a,b,m,P,nu,V,J);
    est_n_RM(k) = size(est_x_RM{k},2);

    OSPA_RM(k) = ospa_dist(groundTruth{k}(1:3,:),est_x_RM{k}(1:3,:),15,1); 
    
    RMSEx(1,k) = cal_RMSE(groundTruth{k}(1:3,:),est_x_RM{k}(1:3,:));
    RMSEv(1,k) = cal_RMSE(groundTruth{k}(4:6,:),est_x_RM{k}(4:6,:));

    %IoU(k) = cal_IoU_RM(groundTruth{k},est_x_RM{k}(1:3,:),est_X_RM{k});
end

end