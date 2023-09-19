function [ghat_k,Xhat_k,Phat_k] = GP_stateExtraction(w_k,a_k,b_k,m_k,P_k,J_k)
% Function that extracts state estimates from the PHD filter.
%ghat_k 量测率
%Xhat_k 运动状态
%Phat_k 状态协方差

    
    ghat_k = zeros(1,0);
    Xhat_k = zeros(size(m_k,1),0);
    Phat_k = zeros(size(m_k,1),size(m_k,1),0);
    counter = 0;
    for i = 1:J_k
        if w_k(1,i) > 0.4
            counter = counter+1;
            ghat_k = [ghat_k a_k(i)/b_k(i)];
            Xhat_k = [Xhat_k , m_k(:,i)];
            Phat_k(:,:,counter) = P_k(:,:,i);
        end
    end
end