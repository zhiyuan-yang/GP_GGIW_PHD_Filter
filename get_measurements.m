function [meas,groundTruth] = get_measurements(objType,lambda)
numInstants = 300;
T = 0.1;
groundTruth = cell(numInstants,1);
groundTruth1 = zeros(7,numInstants);
groundTruth2 = zeros(7,numInstants);
groundTruth1(7,:) = objType(1);
groundTruth2(7,:) = objType(2);
meas = cell(numInstants,1);
stdMeas = 0.1;
lambda_c = 2;
obsArea = [-30 30;
           -30 30;
           -5  5];

%% A. Trajectoray 1
initialPosition = [0,0,1;
                   15,0,1]';
groundTruth1(1:3,:) = initialPosition(:,1) + [zeros(1,300); T*(1:300); zeros(1,300)]; %沿y轴做1m/s匀速直线运动
groundTruth1(4:6,:) = repmat([1 0 0]',1,numInstants);

%z轴无速度无运动
groundTruth2(3,:) = initialPosition(3,2);
groundTruth2(6,:) = 0;

groundTruth2(1:2,1:100) = initialPosition(1:2,2) + [0.5*T*(1:100); zeros(1,100)];%先沿x轴做0.5m/s匀速直线运
groundTruth2(4:5,1:100) = repmat([0.5,0]',1,100);

omega = pi/20;
deltaOmega = omega*T;
F = [1,0,sin(deltaOmega)/omega,(cos(deltaOmega)-1)/omega;
     0,1,(1-cos(deltaOmega))/omega,sin(deltaOmega)/omega;
     0,0,cos(deltaOmega),-sin(deltaOmega);
     0,0,sin(deltaOmega),cos(deltaOmega)];
for i=101:200
     temp = F*[groundTruth2(1:2,i-1);groundTruth2(4:5,i-1)];
     groundTruth2(1:2,i)=temp(1:2);
     groundTruth2(4:5,i) = temp(3:4);
end
F = kron([1,T;0,1],eye(2));
for i=201:300
    temp = F*[groundTruth2(1:2,i-1);groundTruth2(4:5,i-1)];
    groundTruth2(1:2,i)=temp(1:2);
    groundTruth2(4:5,i) = temp(3:4);
end

%% B. Trajectory 2
% groundTruth1(3,:) = 1;
% groundTruth2(3,:) = 1;  %只在xy平面运动
% omega1 = pi/90;
% omega2 = -pi/100;
% F1 = get_motion_matrix(omega1,T);
% F2 = get_motion_matrix(omega2,T);
% position_k1 = [0;0;1;0];
% position_k2 = [0;-20;2;0];
% for i=1:300
%     position_k1 = F1 * position_k1;
%     position_k2 = F2 * position_k2;
%     groundTruth1(1:2,i) = position_k1(1:2);
%     groundTruth1(4:5,i) = position_k1(3:4);
%     groundTruth2(1:2,i) = position_k2(1:2);
%     groundTruth2(4:5,i) = position_k2(3:4);
% end


%% Generate Measurements
for i=1:numInstants
    groundTruth{i} = [groundTruth1(:,i),groundTruth2(:,i)];
    for k = 1:2
        numMeasurements = random('Poisson',lambda);
        currGroundTruth = groundTruth{i}(:,k); 
        temp_meas = produce_surface_measurements(currGroundTruth,numMeasurements,stdMeas,objType(k));
        meas{i} = [meas{i} temp_meas];
    end
    numClutter = random('Poisson',lambda_c); %add clutter
    clutter = repmat(obsArea(:,1),[1,numClutter]) + diag(obsArea*[-1;1])*rand(3,numClutter);
    meas{i} = [meas{i} clutter];
end

%%
figure(1);
%plot3(groundTruth1(1,:),groundTruth1(2,:),groundTruth1(3,:),'red');
hold on;
plot3(groundTruth2(1,:),groundTruth2(2,:),groundTruth2(3,:),'blue');
 %scatter3(groundTruth1(1,:),groundTruth1(2,:),groundTruth1(3,:));
 %scatter3(groundTruth1(1,:),groundTruth1(2,:),groundTruth1(3,:));
grid on;
xlabel('X(m)');
ylabel('Y(m)');
%zlabel('Z/m');
for k=1:30:300
    s1 = scatter3(meas{k}(1,:),meas{k}(2,:),meas{k}(3,:),'.b');
    s2 = scatter3(groundTruth1(1,k),groundTruth1(2,k),groundTruth1(3,k),'xr');
    s3 = scatter3(groundTruth2(1,k),groundTruth2(2,k),groundTruth2(3,k),'xr');
end
legend([s1,s2],{'Measurements','Groundtruth'});
hold off;

end

function temp_meas = produce_surface_measurements(currGroundTruth,numMeasurements,stdMeas,objType)
    xTrue_L = zeros(1,numMeasurements);
    yTrue_L = zeros(1,numMeasurements);
    zTrue_L = zeros(1,numMeasurements);
    theta = atan2(currGroundTruth(5),currGroundTruth(4));
    transformMatix = [cos(theta),-sin(theta),0;
                      sin(theta),cos(theta),0;
                      0,0,1];
    switch objType 
        case 1 %cylinder r=1 d=2
            r = 1;
            d = 2;
            theta = 2*pi*rand(1,numMeasurements);
            xTrue_L = r*cos(theta);
            yTrue_L = r*sin(theta);
            zTrue_L = -1 + d*rand(1,numMeasurements);
        case 2 % Box 
            h = 2;
            l = 4;
            w = 2;
            faceArray = randi([0 5], numMeasurements, 1);
            pointer = 1;
            for iFace = 0:5
                num = sum(faceArray == iFace);
                switch iFace
                    case 0 % on X-Y plane, Z is positive
                        xTrue_L(pointer:(pointer+num-1)) = -l/2 + l*rand(1,num);
                        yTrue_L(pointer:(pointer+num-1)) = -w/2 + w*rand(1,num);
                        zTrue_L(pointer:(pointer+num-1)) = h/2*ones(1,num);
                    case 1 % on X-Y plane, Z is negative
                        xTrue_L(pointer:(pointer+num-1)) = -l/2 + l*rand(1,num);
                        yTrue_L(pointer:(pointer+num-1)) = -w/2 + w*rand(1,num);
                        zTrue_L(pointer:(pointer+num-1)) = -h/2*ones(1,num);
                    case 2 % on X-Z plane, Y is positive
                        xTrue_L(pointer:(pointer+num-1)) = -l/2 + l*rand(1,num);
                        yTrue_L(pointer:(pointer+num-1)) = w/2*ones(1,num);
                        zTrue_L(pointer:(pointer+num-1)) = -h/2 + h*rand(1,num);
                    case 3 % on X-Z plane, Y is negative
                        xTrue_L(pointer:(pointer+num-1)) = -l/2 + l*rand(1,num);
                        yTrue_L(pointer:(pointer+num-1)) = -w/2*ones(1,num);
                        zTrue_L(pointer:(pointer+num-1)) = -h/2 + h*rand(1,num);
                    case 4 % on Y-Z plane, X is positive
                        xTrue_L(pointer:(pointer+num-1)) = l/2*ones(1,num);
                        yTrue_L(pointer:(pointer+num-1)) = -w/2 + w*rand(1,num);
                        zTrue_L(pointer:(pointer+num-1)) = -h/2 + h*rand(1,num);
                    case 5 % on Y-Z plane, X is negative
                        xTrue_L(pointer:(pointer+num-1)) = -l/2*ones(1,num);
                        yTrue_L(pointer:(pointer+num-1)) = -w/2 + w*rand(1,num);
                        zTrue_L(pointer:(pointer+num-1)) = -h/2 + h*rand(1,num);
                end
                pointer = pointer + num;
            end
        case 3  %sphere
            radius = 2;
            sampleArray = randn(3, numMeasurements);
            sampleArray = normc(sampleArray);
            sampleArray = radius .* sampleArray;
            xTrue_L = sampleArray(1,:);
            yTrue_L = sampleArray(2,:);
            zTrue_L = sampleArray(3,:);
    end
    rotationTrue_L = transformMatix*[xTrue_L;yTrue_L;zTrue_L];
    xTrue_L = rotationTrue_L(1,:);
    yTrue_L = rotationTrue_L(2,:);
    zTrue_L = rotationTrue_L(3,:);
    xL = xTrue_L + stdMeas * randn(1, numMeasurements);
    yL = yTrue_L + stdMeas * randn(1, numMeasurements);
    zL = zTrue_L + stdMeas * randn(1, numMeasurements);
    temp_meas = [xL + currGroundTruth(1); yL + currGroundTruth(2); zL + currGroundTruth(3)];
end


function F = get_motion_matrix(omega,T)
deltaOmega = omega*T;
F = [1,0,sin(deltaOmega)/omega,(cos(deltaOmega)-1)/omega;
     0,1,(1-cos(deltaOmega))/omega,sin(deltaOmega)/omega;
     0,0,cos(deltaOmega),-sin(deltaOmega);
     0,0,sin(deltaOmega),cos(deltaOmega)];
end