function [meas] = test(groundTruth,numInstants,objType,lambda)

initialQuaternion = [0 ; 0; 0; 1];
initialPosition = [0; 0; 0];
linearVelocityMagnitude = 10; % in m/sec
angularVelocityArray = zeros(numInstants, 3);
angularVelocityPhi = 0;         % in rad/sec, around X axis of the local frame
angularVelocityTheta = 0;   % in rad/sec, around Y axis of the local frame
angularVelocityPsi = 0;    % in rad/sec, around Z axis of the local frame
angularVelocityArray(:, 1) = angularVelocityPhi; 
angularVelocityArray(:, 2) = angularVelocityTheta; 
angularVelocityArray(31:230, 3) = angularVelocityPsi; 

posArray = zeros(numInstants, 3); % : [xArray yArray zArray]
quatArray = zeros(numInstants, 4); % : [q0Array q1Array q2Array q3Array]
velArray = zeros(numInstants, 3); 
groundTruth = zeros(numInstants, 14); % : [timeStamp centerPosition' quaternions' angularRate' linearVel']


% produce groundTruth
for i= 1:numInstants

end



% % produce measurements
% for i = 1:numInstants
%     for k = 1:2
%         numMeasurements = random('Poisson',lambda);
%         currGroundTruth = groundTruth{i}(:,k); 
%         temp_meas = produce_surface_measurements(currGroundTruth,numMeasurements,stdMeas,objType);
%         meas{i} = [meas{i} temp_meas];
%     end
%     numClutter = random('Poisson',lambda_c); %add clutter
%     clutter = repmat(obsArea(:,1),[1,numClutter]) + diag(obsArea*[-1;1])*rand(3,numClutter);
%     meas{i} = [meas{i} clutter];
% end



end