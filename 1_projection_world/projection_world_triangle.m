%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%相机参数
%相机的平移向量
t = [6;6;2.5];
%相机的旋转角度
th_x = -pi*3/5;
th_y = 0;
th_z = pi*2/4; 
%相机的旋转矩阵
R_x = [1 0 0; 0 cos(th_x) -sin(th_x); 0 sin(th_x) cos(th_x)];
R_y = [cos(th_y) 0 sin(th_y); 0 1 0; -sin(th_y) 0 cos(th_y)];
R_z = [cos(th_z) -sin(th_z) 0; sin(th_z) cos(th_z) 0; 0 0 1];
R = R_z * R_y * R_x;
%相机的变换矩阵；已知相机坐标，通过相机的旋转平移Tcw，可得世界坐标：Pw = Tcw.Pc
Tcw = [R t; 0 0 0 1]; 
%已知世界坐标，求相机坐标：Pc = Tcw^(-1)Pw
T = inv(Tcw);

%内参矩阵
fx = 518.0;
fy = 519.0;
cx = 325.5; %横向651
cy = 253.5; %纵向507
K = [fx 0 cx; 0 fy cy; 0 0 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%三维空间

%每1个单位取ngap个点，两点最小距离为gap
ngap = 50;
gap = 1/ngap;

%房间长宽高
a = 10;
b = 10;
h = 3;

%房间
axe_x = [linspace(0,a,ngap*a) linspace(0,a,ngap*a) linspace(0,a,ngap*a)...
    linspace(0,a,ngap*a) zeros(1,ngap*b) zeros(1,ngap*b)...
    ones(1,ngap*b)*a ones(1,ngap*b)*a zeros(1,ngap*h)...
    ones(1,ngap*h)*a zeros(1,ngap*h) ones(1,ngap*h)*a];

axe_y = [zeros(1,ngap*a) zeros(1,ngap*a) ones(1,ngap*a)*b...
    ones(1,ngap*a)*b linspace(0,b,ngap*b) linspace(0,b,ngap*b)...
    linspace(0,b,ngap*b) linspace(0,b,ngap*b) zeros(1,ngap*h)...
    zeros(1,ngap*h) ones(1,ngap*h)*b ones(1,ngap*h)*b];

axe_z = [zeros(1,ngap*a) ones(1,ngap*a)*h zeros(1,ngap*a)...
    ones(1,ngap*a)*h zeros(1,ngap*b) ones(1,ngap*b)*h...
    zeros(1,ngap*b) ones(1,ngap*b)*h linspace(0,h,ngap*h)...
    linspace(0,h,ngap*h) linspace(0,h,ngap*h) linspace(0,h,ngap*h)];

%物体六个面
%第一个面
S1x = ones(1,ngap*ngap) + 2;
S1y = repmat(linspace(0,1,ngap), 1, ngap) + 4; %向右复制ngap次
S1z = repmat(linspace(0,1,ngap), ngap, 1); %先向下复制ngap次
S1z = S1z(:)' ;%转成列向量再转置->行向量
%第二个面
S2x = zeros(1,ngap*ngap) + 2;
S2y = repmat(linspace(0,1,ngap), 1, ngap) + 4;
S2z = repmat(linspace(0,1,ngap), ngap, 1);
S2z = S2z(:)';
%第三个面
S3x = repmat(linspace(0,1,ngap), 1, ngap) + 2;
S3y = zeros(1,ngap*ngap) + 4;
S3z = repmat(linspace(0,1,ngap), ngap, 1);
S3z = S3z(:)';
%第四个面
S4x = repmat(linspace(0,1,ngap), 1, ngap) + 2;
S4y = ones(1,ngap*ngap) + 4;
S4z = repmat(linspace(0,1,ngap), ngap, 1);
S4z = S4z(:)';
%第五个面
S5x = repmat(linspace(0,1,ngap), ngap, 1);
S5x = S5x(:)' + 2;
S5y = repmat(linspace(0,1,ngap), 1, ngap) + 4;
S5z = ones(1,ngap*ngap);
%第六个面
S6x = repmat(linspace(0,1,ngap), ngap, 1);
S6x = S6x(:)' + 2;
S6y = repmat(linspace(0,1,ngap), 1, ngap) + 4;
S6z = zeros(1,ngap*ngap);

%物体全部
obj_x = [S1x S2x S3x S4x S5x S6x];
obj_y = [S1y S2y S3y S4y S5y S6y];
obj_z = [S1z S2z S3z S4z S5z S6z];

%三维空间模拟图
% figure(1)
% scatter3(axe_x,axe_y,axe_z,'.')
% hold on;
% scatter3(S1x,S1y,S1z,'.');
% hold on;
% scatter3(S2x,S2y,S2z,'.');
% hold on;
% scatter3(S3x,S3y,S3z,'.');
% hold on;
% scatter3(S4x,S4y,S4z,'.');
% hold on;
% scatter3(S5x,S5y,S5z,'.');
% hold on;
% scatter3(S6x,S6y,S6z,'.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%转换到相机坐标系

%世界坐标系下的全部（物体+房间）
x = [obj_x axe_x];
y = [obj_y axe_y];
z = [obj_z axe_z];
Pw = [x; y; z]; %普通三维坐标(xw,yw,zw)
Pw = [Pw; ones(1, length(x))]; %转换成四维归一化坐标(xw,yw,zw,1)

%相机坐标系下的全部
Pc = T * Pw; %四维归一化坐标(x,y,z,1)（<-转换到相机坐标系中）
Pc = Pc(1:3, :); %三维普通坐标(x,y,z)（<-去掉一维1）

%八个顶点(列向量)
A = Pc(1:3,3*ngap^2);
B = Pc(1:3,ngap^2);
C = Pc(1:3,2*ngap^2);
D = Pc(1:3,4*ngap^2+1);
E = Pc(1:3,1);
F = Pc(1:3,6*ngap^2);
G = Pc(1:3,3*ngap^2+1);
H = Pc(1:3,ngap^2+1);
%顶点集合 (xyz, 8顶点)
Pks = [A B C D E F G H];

%三角形集合 (xyz, 3顶点, 12三角形)
Trs = ones(3, 3, 12);
Trs(:,:,1) = [B C A];
Trs(:,:,2) = [D A C];
Trs(:,:,3) = [F E G];
Trs(:,:,4) = [H G E];
Trs(:,:,5) = [B A F];
Trs(:,:,6) = [E F A];
Trs(:,:,7) = [C G D];
Trs(:,:,8) = [H D G];
Trs(:,:,9) = [B F C];
Trs(:,:,10) = [G C F];
Trs(:,:,11) = [A D E];
Trs(:,:,12) = [H E D];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%处理遮盖问题
Pcc = [Pc; Pc(3,:)]; %四维坐标(x,y,z,z)
for i = 1:length(Pcc(1,:))
    if(Pcc(4,i)>0) %无需判断不入镜的点
        for j = 1:12
            %三角形法向量:N (列向量)
            N = cross(Trs(:,2,j)-Trs(:,1,j),Trs(:,3,j)-Trs(:,1,j));
            N = N/norm(N);
            %三角形平面方程常数项
            Dp = -N'*Trs(:,1,j);
            %物体点方向向量:R (列向量)
            R = Pc(:,i)/norm(Pc(:,i));
            %平视(N'*R=0)：物体点不被遮挡
            %不能被同一平面三角形遮挡
            if(N'*R < 0 && N'*Pc(:,i)+Dp < 0-gap*0.01) 
                t = -Dp/(N'*R);
                I = R*t; %物体点在三角形平面上的投影点
                %三角形三边 (列向量)
                Edg1 = Trs(:,2,j)-Trs(:,1,j);
                Edg2 = Trs(:,3,j)-Trs(:,2,j);
                Edg3 = Trs(:,1,j)-Trs(:,3,j);
                %三角形三顶点到交点I (列向量)
                SI1 = I - Trs(:,1,j);
                SI2 = I - Trs(:,2,j);
                SI3 = I - Trs(:,3,j);
                %被遮挡的点第四维置为-1
                if(N'*cross(Edg1, SI1) > 0-gap*0.01 ...
                    && N'*cross(Edg2, SI2) > 0-gap*0.01 ...
                    && N'*cross(Edg3, SI3) > 0-gap*0.01)
                    Pcc(4,i) = -1;
                end
            end
        end
    end
end


%分割获取6个面
Ptest = [Pc; Pcc(4,:)];
S1 = Ptest(:,1:ngap^2);
S2 = Ptest(:,ngap^2+1:2*ngap^2);
S3 = Ptest(:,2*ngap^2+1:3*ngap^2);
S4 = Ptest(:,3*ngap^2+1:4*ngap^2);
S5 = Ptest(:,4*ngap^2+1:5*ngap^2);
S6 = Ptest(:,5*ngap^2+1:6*ngap^2);
axe = Ptest(:,6*ngap^2+1:length(Ptest(1,:)));

%去除被遮挡的点
S1 = S1(:,S1(4,:)>0);
S2 = S2(:,S2(4,:)>0);
S3 = S3(:,S3(4,:)>0);
S4 = S4(:,S4(4,:)>0);
S5 = S5(:,S5(4,:)>0);
S6 = S6(:,S6(4,:)>0);
axe = axe(:,axe(4,:)>0);

%三维空间模拟图（像机坐标系中）
figure(2)
% scatter3(axe(1,:),axe(2,:),axe(3,:),'.')
% hold on;
scatter3(S1(1,:),S1(2,:),S1(3,:),'.');
hold on;
scatter3(S2(1,:),S2(2,:),S2(3,:),'.');
hold on;
scatter3(S3(1,:),S3(2,:),S3(3,:),'.');
hold on;
scatter3(S4(1,:),S4(2,:),S4(3,:),'.');
hold on;
scatter3(S5(1,:),S5(2,:),S5(3,:),'.');
hold on;
scatter3(S6(1,:),S6(2,:),S6(3,:),'.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%投影成像
%三维归一化坐标(x,y,1)
% Pc = Pc./[Pc(3,:); Pc(3,:); Pc(3,:)]; 
S1 = S1(1:3,:)./[S1(3,:); S1(3,:); S1(3,:)];
S2 = S2(1:3,:)./[S2(3,:); S2(3,:); S2(3,:)];
S3 = S3(1:3,:)./[S3(3,:); S3(3,:); S3(3,:)];
S4 = S4(1:3,:)./[S4(3,:); S4(3,:); S4(3,:)];
S5 = S5(1:3,:)./[S5(3,:); S5(3,:); S5(3,:)];
S6 = S6(1:3,:)./[S6(3,:); S6(3,:); S6(3,:)];
axe = axe(1:3,:)./[axe(3,:); axe(3,:); axe(3,:)];
%像素坐标系的照片(u,v,1)
S1 = K * S1; 
S2 = K * S2;
S3 = K * S3;
S4 = K * S4;
S5 = K * S5;
S6 = K * S6;
axe = K * axe;
%在视野范围内的像素平面普通坐标(u,v)（不便展示）
S1 = S1(1:2, S1(1,:)>0 & S1(1,:)<2*cx & S1(2,:)>0 & S1(2,:)<2*cy);
S2 = S2(1:2, S2(1,:)>0 & S2(1,:)<2*cx & S2(2,:)>0 & S2(2,:)<2*cy);
S3 = S3(1:2, S3(1,:)>0 & S3(1,:)<2*cx & S3(2,:)>0 & S3(2,:)<2*cy);
S4 = S4(1:2, S4(1,:)>0 & S4(1,:)<2*cx & S4(2,:)>0 & S4(2,:)<2*cy);
S5 = S5(1:2, S5(1,:)>0 & S5(1,:)<2*cx & S5(2,:)>0 & S5(2,:)<2*cy);
S6 = S6(1:2, S6(1,:)>0 & S6(1,:)<2*cx & S6(2,:)>0 & S6(2,:)<2*cy);
axe = axe(1:2, axe(1,:)>0 & axe(1,:)<2*cx & axe(2,:)>0 & axe(2,:)<2*cy);
%方便可视化展示：翻转y轴
S1 = [S1(1,:); S1(2,:)*(-1)+2*cy];
S2 = [S2(1,:); S2(2,:)*(-1)+2*cy];
S3 = [S3(1,:); S3(2,:)*(-1)+2*cy];
S4 = [S4(1,:); S4(2,:)*(-1)+2*cy];
S5 = [S5(1,:); S5(2,:)*(-1)+2*cy];
S6 = [S6(1,:); S6(2,:)*(-1)+2*cy];
axe = [axe(1,:); axe(2,:)*(-1)+2*cy];

figure(3)
scatter(S1(1,:), S1(2,:), '.');
hold on;
scatter(S2(1,:), S2(2,:), '.');
hold on;
scatter(S3(1,:), S3(2,:), '.');
hold on;
scatter(S4(1,:), S4(2,:), '.');
hold on;
scatter(S5(1,:), S5(2,:), '.');
hold on;
scatter(S6(1,:), S6(2,:), '.');
hold on;
scatter(axe(1,:), axe(2,:), '.');
axis([0 2*cx 0 2*cy]);
















