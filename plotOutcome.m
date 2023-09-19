clear
clc

addpath("data\");
addpath("_common\");
addpath("spheretri-master\");
load("groundTruth.mat");
load("meas.mat");
load('traject1.mat');
load("traject2.mat");
load("est_x_GP.mat");
load("est_x_RM.mat");
load("OSPA_monte_GP.mat");
load("OSPA_monte_RM.mat");
load("RMSE_v_monte_GP.mat");
load("RMSE_v_monte_RM.mat");
load("RMSE_x_monte_GP.mat");
load("RMSE_x_monte_RM.mat");
load("est_extend_RM.mat");
load("meas.mat");

%%
figure(1);
truth = plot3(trajec1.dataLog(:,2),trajec1.dataLog(:,3),trajec1.dataLog(:,4),'k');
hold on;
plot3(trajec2.dataLog(:,2),trajec2.dataLog(:,3),trajec2.dataLog(:,4),'k');
for i=1:300
    GP = scatter3(est_x_GP{i}(1,:),est_x_GP{i}(2,:),est_x_GP{i}(3,:),'r.');
    RM = scatter3(est_x_RM{i}(1,:),est_x_RM{i}(2,:),est_x_RM{i}(3,:),'b.');
    Meas = scatter3(meas{i}(1,:),meas{i}(2,:),meas{i}(3,:),'r+');
end

grid on;
xlabel('X/m');
ylabel('Y/m');
zlabel('Z/m');
legend([truth,GP,RM,Meas],{'真实轨迹','GP-PHD','GGIW-PHD','量测'},'Location','northeast');
hold off;

%%
figure(2);
t = (1:300)*0.1;
plot(t,OSPA_monte_GP,'r',t,OSPA_monte_RM,'b');
legend('GP','GGIW');
xlabel('Time/s');
ylabel('OSPA/m');

figure(3);
plot(t,RMSE_x_monte_GP,'r',t,RMSE_x_monte_RM,'b');
legend('GP','GGIW');


xlabel('Time/s');
ylabel('RMSE/m');

figure(4);
plot(t,RMSE_v_monte_GP,'r',t,RMSE_v_monte_RM,'b');
legend('GP','GGIW');
xlabel('Time/s');
ylabel('RMSE/(m/s)');

%%
figure(5)

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

axis equal;
view(45, 25);
xlabel('X/m');
ylabel('Y/m');
zlabel('Z/m');
grid on;
hold on;

for k=120:60:120
    s1 = scatter3(meas{k}(1,:),meas{k}(2,:),meas{k}(3,:),'+b',LineWidth=1.5);
    for i=1:2
        % plot groundtruth
        x = 2*xSphere + groundTruth{k}(2,i);
        y = 2*ySphere + groundTruth{k}(3,i);
        z = 2*zSphere + groundTruth{k}(4,i);
        s2 = surf(x,y,z,'facecolor','blue','facealpha',0,'edgecolor','black','linewidth',0.1); % surface plot of the sphere       

        % plot random matrix
        extend = est_X_RM{k}(:,:,i);
        aix_length = eig(extend);
        x = aix_length(1)*xSphere + est_x_RM{k}(1,i);
        y = aix_length(2)*ySphere + est_x_RM{k}(2,i);
        z = aix_length(3)*zSphere + est_x_RM{k}(3,i);
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
legend([s1,s2,s3,s4],{'量测','真实扩展目标','GGIW','GP'},'Location','northeast');


