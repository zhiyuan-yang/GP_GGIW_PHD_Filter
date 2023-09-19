function [w_k,a_k,b_k,m_k,P_k,v_k,V_k,J_k] = GGIW_PHD_predict(w_k1,alpha_k1,beta_k1,m_k1,P_k1,v_k1,V_k1,J_k1)

% Sampling time and temporal decay
Ts = 0.1;
tau = 5;
sigma_v = 1e-1;
F_k1 = [1,Ts;0,1];
Q_k1 = [Ts^3/3,Ts^2/2;Ts^2/2,Ts]*sigma_v^2;
v_min = 8;
eta_k = 25/24;
d = 3;

%% 参数声明
% Probabilities of survival and detection
p_S_k = 0.99;
% The spontaneous birth distributions
% w_gam_k = [0.01 0.01];
% a_gam_k = [100 100];
% b_gam_k = [50 50];
% m_gam_k = ;
% P_gam_k = model.birth.P_GGIW;
% v_gam_k = model.birth.v;
% V_gam_k = model.birth.V;
J_gam_k =  0;

% The spawn distributions
J_beta_k = 0;% no spawn consider





%% 给需要预测的参数预分配内存
w_k = zeros(1,J_gam_k+J_beta_k*J_k1+J_k1);
a_k = zeros(1,J_gam_k+J_beta_k*J_k1+J_k1);
b_k = zeros(1,J_gam_k+J_beta_k*J_k1+J_k1);
m_k = zeros(size(m_k1,1),J_gam_k+J_beta_k*J_k1+J_k1);
P_k = zeros(size(P_k1,1),size(P_k1,2),J_gam_k+J_beta_k*J_k1+J_k1);
v_k = zeros(1,J_gam_k+J_beta_k*J_k1+J_k1);
V_k = zeros(size(V_k1,1),size(V_k1,2),J_gam_k+J_beta_k*J_k1+J_k1);

%% 新生目标
i = 0;                                                                     % 所有的预测项数 = 新生 + 存活
% for j = 1:J_gam_k
%     i = i+1;
%     w_k(1,i) = w_gam_k(1,j);
%     a_k(1,i) = a_gam_k;
%     b_k(1,i) = b_gam_k;
%     m_k(:,i) = m_gam_k(:,j);
%     P_k(:,:,i) = P_gam_k;
%     v_k(1,i) = v_gam_k;
%     V_k(:,:,i) = V_gam_k;
% end

%% 存活目标预测
for j = 1:J_k1
    i = i+1;
    w_k(1,i) = p_S_k*w_k1(1,j);
    a_k(1,i) = alpha_k1(1,j)/eta_k;                       
    b_k(1,i) = beta_k1(1,j)/eta_k;
    m_k(:,i) = kron(F_k1,eye(d))*m_k1(:,j);               
    P_k(:,:,i) = Q_k1 + F_k1*P_k1(:,:,j)*F_k1';       
    P_k(:,:,i) = 0.5*(P_k(:,:,i)+P_k(:,:,i)');                             % 使矩阵对称
    v_k(1,i) = max(exp(-Ts/tau)*v_k1(j),v_min);                            % Note the numerical hack here. We must have nu>nu_min for the first two moments to be well defined.
    V_k(:,:,i) = ((v_k(1,i)-d-1)/(v_k1(j)-d-1))*V_k1(:,:,j);
    V_k(:,:,i) = 0.5*(V_k(:,:,i)+V_k(:,:,i)');                             % 使矩阵对称
end
J_k = i;

end