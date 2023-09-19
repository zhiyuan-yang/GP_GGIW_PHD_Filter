function IoU = cal_IoU_RM(groundTruth,x,X)
IoU = zeros(1,size(x,2));
for i=1:size(x,2)
    erro = groundTruth(1:3,:) - x(:,i);  
    erro = vecnorm(erro);
    [~,idx] = min(erro);
    truth = groundTruth(1:3,idx);
    objType = groundTruth(7,idx);
    %egiValue = sort(eig(X(:,:,idx)),'descend');
    theta = atan2(groundTruth(5,i),groundTruth(4,i));
    transformMatrix = [cos(theta),-sin(theta),0;
                          sin(theta),cos(theta),0;
                          0,0,1];
    if objType == 1
        %maxTruth = truth + [1;1;1];
        %minTruth = truth - [1;1;1];
        maxTruth = truth + [2;2;2];
        minTruth = truth - [2;2;2];
    elseif objType == 2
        %maxTruth = truth + transformMatrix*[1;2;1];
        %minTruth = truth - transformMatrix*[1;2;1];
        maxTruth = truth + [2;2;2];
        minTruth = truth - [2;2;2];
    end
    %maxestimate = x(:,i) + transformMatrix * egiValue;
    %minestimate = x(:,i) - transformMatrix * egiValue;
    %maxlength = max([maxTruth,maxestimate],[],2);
    %minlength = min([minTruth,minestimate],[],2);
    maxlength = maxTruth;
    minlength = minTruth;
    xlength = linspace(minlength(1),maxlength(1),50);
    ylength = linspace(minlength(2),maxlength(2),50);
    zlength = linspace(minlength(3),maxlength(3),50);
    inBoth = 0;
    inEither = 0;
    inellis = 0;
    incylinder = 0;
    for j1 = 1:50
        for j2=1:50
            for j3=1:50
                inTruth = 0;
                inEst = 0;
                unitTruth = [xlength(j1);ylength(j2);zlength(j3)] - truth;
                unitTruth = transformMatrix * unitTruth;
                unitEst = [xlength(j1);ylength(j2);zlength(j3)] - x(:,i);
                if objType == 1
                    if unitTruth(1)^2 + unitTruth(2)^2 < 1 && abs(unitTruth(3))<1
                        inTruth = 1;
                        incylinder = incylinder + 1;
                    end
                    if unitEst'*inv(X(:,:,i))*unitEst<2
                        inEst = 1;
                        inellis = inellis + 1;
                    end
                    if inTruth ==1 && inEst == 1
                        inBoth = inBoth +1;
                    elseif inTruth == 1 || inEst == 1
                        inEither = inEither + 1;
                    end
                elseif objType == 2
                    if abs(unitTruth(1)) < 1 && abs(unitTruth(2)) < 2 && abs(unitTruth(3))<1
                        inTruth = 1;
                    end
                    if unitEst'*inv(X(:,:,i))*unitEst<1
                        inEst = 1;
                    end
                    if inTruth ==1 && inEst == 1
                        inBoth = inBoth +1;
                    elseif inTruth == 1 || inEst == 1
                        inEither = inEither + 1;
                    end
                end
            end
        end
    end
    IoU(i) = inBoth/inEither ;
end
IoU = mean(IoU);
end