function RMSE = cal_RMSE(X,Y)
    n = size(X,2);
    m = size(Y,2);
    if m == n
        XX= repmat(X,[1 m]);
        YY= reshape(repmat(Y,[n 1]),[size(Y,1) n*m]);
        D = sqrt(sum((XX-YY).^2,1));
        RMSE = min(D(1:n)) + min(D(n+1:n*m));
    else 
        RMSE = 1;
    end
end