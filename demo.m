clear
clc
tic
lambda = 20;        % poisson rate for meas
objType = [3,3];        % 1 for cylinder 2 for box 3 for sphere
rng(2022);

[meas,groundTruth] = get_measurements(objType,lambda);

[est_x_RM,est_X_RM,est_n_RM,OSPA_RM,~,~,IoURM] = filter_RM(meas,groundTruth);

[est_x_GP,est_X_GP,est_n_GP,OSPA_GP,~,~,IoUGP] = filter_GP(meas,groundTruth);
toc

%%
figure(1);
groundTruth1 = zeros(3,300);
groundTruth2 = zeros(3,300);
for k=1:300
    groundTruth1(:,k) = groundTruth{k}(1:3,1);
    groundTruth2(:,k) = groundTruth{k}(1:3,2);
end
plot3(groundTruth1(1,:),groundTruth1(2,:),groundTruth1(3,:),'red');
hold on;
plot3(groundTruth2(1,:),groundTruth2(2,:),groundTruth2(3,:),'blue');
grid on;
xlabel('X/m');
ylabel('Y/m');
zlabel('Z/m');
for k=1:20:300
    scatter3(meas{k}(1,:),meas{k}(2,:),meas{k}(3,:),'+k');
end
hold off;

%%
figure(2);
RM = plot(0.1 * (1:300),OSPA_RM,'red');
hold on;
GP = plot(0.1 * (1:300),OSPA_GP,'blue');
%legend([RM,GP],{'RM','GP'});
hold off;

%%
figure(3);
for k=1:300
    GP = scatter3(est_x_GP{k}(1,:),est_x_GP{k}(2,:),est_x_GP{k}(3,:),'r.');
    hold on;
end
xlabel('X/m');
ylabel('Y/m');
zlabel('Z/m');
for i=1:300
    RM = scatter3(est_x_RM{i}(1,:),est_x_RM{i}(2,:),est_x_RM{i}(3,:),'b.');
end
plot3(groundTruth1(1,:),groundTruth1(2,:),groundTruth1(3,:),'black');
truth = plot3(groundTruth2(1,:),groundTruth2(2,:),groundTruth2(3,:),'black');
legend([GP,RM,truth],{'GP','RM','Truth'},'Location','best');
zlim([-3,3]);
hold off;

%% 
figure(4)
addpath("_common\");
axis equal;
xlabel('X(m)');
ylabel('Y(m)');
zlabel('Z(m)');
grid on;
hold on;
[xSphere, ySphere, zSphere] = sphere(20);
minNumBasisAngles = 200;        % minimum number of basis angles on which the extent is maintained
% Determine the Basis Angles
[basisVertices, ~] = spheretri(minNumBasisAngles);      % produces evenly spaced points on the sphere
[azimuthBasisArray, elevationBasisArray, ~] = cart2sph(basisVertices(:,1), basisVertices(:,2)...
    , basisVertices(:,3));
azimuthBasisArray = mod(azimuthBasisArray, 2*pi);       % to make it consistent with the convention
numBasisAngles = size(basisVertices, 1);
% Arrange the basis array so that azimuth and elevation are in ascending order
basisAngleArray = [azimuthBasisArray elevationBasisArray];
basisAngleArray = sortrows(basisAngleArray, 2);
basisAngleArray = sortrows(basisAngleArray, 1);

for k=120:60:240
    s1 = scatter3(meas{k}(1,:),meas{k}(2,:),meas{k}(3,:),'+b',LineWidth=1.5);
    for i=1:2
        % plot groundtruth
        theta = atan2(groundTruth{k}(5,1),groundTruth{k}(4,1));
        transformMatrix = [cos(theta),-sin(theta),0;
                          sin(theta),cos(theta),0;
                          0,0,1];
        switch objType(i)
            case 1  %Cylinder
                r = 1;
                d = 2;
                [xCylinder, yCylinder, zCylinder] = cylinder(r);
                x = xCylinder + groundTruth{k}(1,i);
                y = yCylinder + groundTruth{k}(2,i);
                z = d * zCylinder + groundTruth{k}(3,i) - 1;
                s2 = surf(x,y,z,'facecolor','blue','facealpha',0,'edgecolor','black','linewidth',0.1); % surface plot of the sphere
            case 2  %box
                height = 2;
                length = 4;
                width = 2;
                vertices_L = [length width -height;
                    -length width -height;
                    -length width height;
                    length width height;
                    -length -width height;
                    length -width height;
                    length -width -height;
                    -length -width -height]* 0.5;
                vertices_L = vertices_L * transformMatrix';
                fac = [1 2 3 4;
                    4 3 5 6;
                    6 7 8 5;
                    1 2 8 7;
                    6 7 1 4;
                    2 3 5 8];
                box = [vertices_L(:,1)+groundTruth{k}(1,i),vertices_L(:,2)+groundTruth{k}(2,i),vertices_L(:,3)+groundTruth{k}(3,i)];
                s2 = patch('Faces',fac,'Vertices',box, 'EdgeColor'...
                         , [0.3 0.3 0.3], 'LineWidth', 1,'FaceAlpha', 0);
            case 3  %sphere
                r = 2;
                x = r * xSphere + groundTruth{k}(1,i);
                y = r * ySphere + groundTruth{k}(2,i);
                z = r * zSphere + groundTruth{k}(3,i);
                s2 = surf(x,y,z,'facecolor','blue','facealpha',0,'edgecolor','black','linewidth',0.1);
        end

        % plot random matrix
        extend = est_X_RM{k}(:,:,i);
        aix_length = eig(extend);
        aix_length = sort(aix_length,'descend');
        thetaRM = atan2(est_x_RM{k}(5,i),est_x_RM{k}(4,i));
        transformMatrixRM = [cos(thetaRM),-sin(thetaRM),0;
                          sin(thetaRM),cos(thetaRM),0;
                          0,0,1];
        for j=1:size(xSphere,1)
            rotation = transformMatrixRM * [xSphere(j,:);ySphere(j,:);zSphere(j,:)];
            xSphere(j,:) = rotation(1,:);
            ySphere(j,:) = rotation(2,:);
            zSphere(j,:) = rotation(3,:);
        end
        aix_length = transformMatrixRM * aix_length;
        x = 2*aix_length(1)*xSphere + est_x_RM{k}(1,i);
        y = 2*aix_length(2)*ySphere + est_x_RM{k}(2,i);
        z = 2*aix_length(3)*zSphere + est_x_RM{k}(3,i);
        s3 = surf(x,y,z,'facecolor',[0.3 0.75 0.93],'facealpha',0.6,'edgecolor','none','linewidth',0.1); % surface plot of the sphere

        % plot GP
        estCenter = est_x_GP{k}(1:3,i);
        estExtend = est_x_GP{k}(13:end,i);
        [x_L, y_L, z_L] = sph2cart(basisAngleArray(:,1), basisAngleArray(:,2), estExtend);
        cart_L = [x_L'; y_L'; z_L'];
        
        % Transform the extent to the Global frame
        extentG = estCenter + cart_L;
        xExtent = transpose(extentG(1, :));
        yExtent = transpose(extentG(2, :));
        zExtent = transpose(extentG(3, :));
        triangles = boundary(extentG', 0.5);
        s4 = trisurf(triangles, xExtent, yExtent, zExtent, 'FaceAlpha', 0.1, 'FaceColor',  [1 0.9 0]); % Plot the triangles
        s4.EdgeColor = 'none';
    end
end

legend([s1,s2,s3,s4],{'Measurements','Truth','GGIW-PHD','GP-PHD'},'Location','northeast');
view([-37.5,30]);
xlim([-5 30]);
ylim([0 30]);
zlim([-2 8]);