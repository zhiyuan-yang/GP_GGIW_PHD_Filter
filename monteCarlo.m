clear
clc   
tic

numMonteCarlo = 50;
expDuration = 30;   % (in seconds)
T = 0.1;            % sampling time (in seconds)
numInstants = ceil(expDuration/T);
lambda = 30;        % poisson rate for meas
N_FA = 10;          % poisson rate for clutter
objType = [1,2];        % 1 for cylinder  2 for cube 3 for sphere
OSPA_RM = zeros(numMonteCarlo,numInstants);
RMSE_x_RM = zeros(numMonteCarlo,numInstants);
RMSE_v_RM = zeros(numMonteCarlo,numInstants);
OSPA_GP = zeros(numMonteCarlo,numInstants);
RMSE_x_GP = zeros(numMonteCarlo,numInstants);
RMSE_v_GP = zeros(numMonteCarlo,numInstants);
IoURM = zeros(numMonteCarlo,numInstants);
IoUGP = zeros(numMonteCarlo,numInstants);

%%
parpool(3);
parfor i=1:numMonteCarlo
    %fprintf("Current Monte Carlo Turn: %d\n", i)
    [meas,groundTruth] = get_measurements(objType,lambda);
    [~,~,~,OSPA_RM(i,:),~,~,IoURM(i,:)] = filter_RM(meas,groundTruth);
    %[~,~,~,OSPA_GP(i,:),~,~,IoUGP(i,:)] = filter_GP(meas,groundTruth);
end
toc