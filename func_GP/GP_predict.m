function [w_k,a_k,b_k,m_k,P_k,J_k] = GP_predict(w_k1,alpha_k1,beta_k1,m_k1,P_k1,J_k1,numBasisAngles)
addpath('_common\');

T = 0.1;
Fk = kron([1 T; 0 1], eye(3));
Fe = eye(numBasisAngles);
angVelEst = m_k1(10:12);
stdAngVel = 1e-1;
[FRot, QRot] = compute_rotational_dynamics(angVelEst, stdAngVel, T);
F = blkdiag(Fk,FRot,Fe);
stdCenter = 1e-1;
Qk = kron([T^3/3 T^2/2; T^2/2 T], stdCenter^2*eye(3));
Qe = zeros(numBasisAngles);
Q = blkdiag(Qk, QRot, Qe);
eta_k = 25/24;
lambda = 0.99;
p_S_k = 0.99;


% The spontaneous birth distributions
% w_gam_k = model.birth.w;
% a_gam_k = model.birth.a;
% b_gam_k = model.birth.b;
% m_gam_k = model.birth.m_RHM;
% P_gam_k = model.birth.P_RHM;
J_gam_k = 0;

% The spawn distributions
J_beta_k = 0;% no spawn consider

%% 给需要预测的参数预分配内存
w_k = zeros(1,J_gam_k+J_beta_k*J_k1+J_k1);
a_k = zeros(1,J_gam_k+J_beta_k*J_k1+J_k1);
b_k = zeros(1,J_gam_k+J_beta_k*J_k1+J_k1);
m_k = zeros(size(m_k1,1),J_gam_k+J_beta_k*J_k1+J_k1);
P_k = zeros(size(P_k1,1),size(P_k1,2),J_gam_k+J_beta_k*J_k1+J_k1);

%% 新生目标
i = 0;                                                                     % 所有的预测项数 = 新生 + 存活
% for j = 1:J_gam_k
%     i = i+1;
%     w_k(1,i) = w_gam_k(1,j);
%     a_k(1,i) = a_gam_k;
%     b_k(1,i) = b_gam_k;
%     m_k(:,i) = m_gam_k(:,j);
%     P_k(:,:,i) = P_gam_k;
% end

%% 存活目标预测
for j = 1:J_k1
    i = i+1;
    w_k(1,i) = p_S_k*w_k1(1,j);
    a_k(1,i) = alpha_k1(1,j)/eta_k;                       
    b_k(1,i) = beta_k1(1,j)/eta_k;
    m_k(:,i) = F*m_k1(:,j); 
    P_k(:,:,i) = Q + F*P_k1(:,:,j)*F';       
    P_k(:,:,i) = 0.5*(P_k(:,:,i)+P_k(:,:,i)');                             % 使矩阵对称
    P_k(13:end,13:end,i) = 1/lambda * P_k(13:end,13:end,i);  
end
J_k = i;

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