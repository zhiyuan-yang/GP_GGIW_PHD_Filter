function IoU = cal_IoU_GP(groundTruth,x,basisAngleArray)
    IoU = zeros(1,size(x,2));
    for i=1:size(x,2)
        erro = groundTruth(1:3,:) - x(1:3,i);
        erro = vecnorm(erro);
        [~,idx] = min(erro);
        objType = groundTruth(7,idx);
        inboth = 0;
        ineither = 0;
        if objType == 1
            volume = linspace(0,sqrt(2),200);
            for r = volume
                for k = 1:642
                    intruth = 0;
                    inest = 0;
                    unitz = r*sin(basisAngleArray(k,2));
                    unitx = r*cos(basisAngleArray(k,2))*cos(basisAngleArray(k,1));
                    unity = r*cos(basisAngleArray(k,2))*sin(basisAngleArray(k,1));
                    if x(k+12,i) > r
                        inest = 1;
                    end
                    if unitx^2+unity^2 <1 && abs(unitz)<1
                        intruth = 1;
                    end
                    if intruth == 1 && inest == 1
                        inboth = inboth + 1;
                    elseif intruth == 1 ||  inest == 1
                        ineither = ineither + 1;
                    end
                end
            end
        elseif objType == 2
            volume = linspace(0,sqrt(6),200);
            for r = volume
                for k = 1:642
                    intruth = 0;
                    inest = 0;
                    unitz = r*sin(basisAngleArray(k,2));
                    unitx = r*cos(basisAngleArray(k,2))*cos(basisAngleArray(k,1));
                    unity = r*cos(basisAngleArray(k,2))*sin(basisAngleArray(k,1));
                    if x(k+12,i) > r
                        inest = 1;
                    end
                    if abs(unitx)<1 && abs(unity)<2 && abs(unitz)<1
                        intruth = 1;
                    end
                    if intruth == 1 && inest == 1
                        inboth = inboth + 1;
                    elseif intruth == 1 ||  inest == 1
                        ineither = ineither + 1;
                    end
                end
            end
        end
        IoU(i) = inboth/(inboth + ineither);
    end
IoU = mean(IoU);
end