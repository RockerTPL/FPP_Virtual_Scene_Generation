%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%相机参数
%相机的平移向量
t = [6;2;2.5];
%相机的旋转角度
th_x = -pi*3/5;
th_y = 0;
th_z = pi*1/4; 
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
global fx;
global fy;
global cx;
global cy;
global K;
fx = 518.0;
fy = 519.0;
cx = 325.5; %横向651 % c代表光心和像素坐标系之间的偏移量，2c即为像素平面尺寸
cy = 253.5; %纵向507
% cx = 200;
% cy = 200;
K = [fx 0 cx; 0 fy cy; 0 0 1];

global crit; %极小值判断
global Crit; %极大值判断
crit = -0.0001;
Crit = 10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 三维空间

%房间长宽高
a = 10;
b = 10;
h = 3;
%正方体棱长
c = 1;
%正方体平移向量
tcb = [2; 4; 0];


%房间顶点 (3个xyz坐标 * 8个顶点)
Rmp = [[0; 0; 0] [0; 0; h] ...
        [0; b; 0] [0; b; h]...
        [a; 0; 0] [a; 0; h]...
        [a; b; 0] [a; b; h]];

%正方体顶点 (3个xyz坐标 * 8个顶点)
Cbp = [[0; 0; 0] [0; 0; c] ...
        [0; c; 0] [0; c; c]...
        [c; 0; 0] [c; 0; c]...
        [c; c; 0] [c; c; c]]...
        + tcb; 

    
% 所有房间三角形向后排列 (3个xyz坐标 * 3个顶点 * 10个三角形)
Rm = Rmp2Tr(Rmp);
% 所有物体三角形向后排列 (3个xyz坐标 * 3个顶点 * 12个三角形)
Cb = Cbp2Tr(Cbp);


% 房间颜色 (3个RGB * 10个三角形)
Rmcl = zeros(3,10);
for i = 1:5
    Rmcl(:, 2*i-1) = [1-i*0.1; 1-i*0.1; 1-i*0.1];
    Rmcl(:, 2*i) = Rmcl(:, 2*i-1);
end
% 物体颜色 (3个RGB * 12个三角形)
Cbcl = zeros(3,12);
for i = 1:6
    Cbcl(:, 2*i-1) = [1-i*0.15; 0.7-i*0.1; i*0.15];
    Cbcl(:, 2*i) = Cbcl(:, 2*i-1);
end
% 全部22个三角形的颜色集合 (3个RGB * 22个三角形)
Allcl = [Rmcl, Cbcl];

% 所有房间三角形各取一点向右排列 (3个点的x/y/z坐标 * 10个三角形)
Rmx = zeros(3,10);
Rmy = zeros(3,10);
Rmz = zeros(3,10);
for i = 1:10
    Rmx(:, i) = Rm(1, :, i)';
    Rmy(:, i) = Rm(2, :, i)';
    Rmz(:, i) = Rm(3, :, i)';
end
% 所有房间三角形各取一点向右排列 (3个点的x/y/z坐标 * 12个三角形)
Cbx = zeros(3,12);
Cby = zeros(3,12);
Cbz = zeros(3,12);
for i = 1:12
    Cbx(:, i) = Cb(1, :, i)';
    Cby(:, i) = Cb(2, :, i)';
    Cbz(:, i) = Cb(3, :, i)';
end


%世界坐标三维空间模拟图
figure(1)
for i = 1:10
    fill3(Rmx(:, i), Rmy(:, i), Rmz(:, i), Rmcl(:, i)');
    hold on;
end
for i = 1:12
    fill3(Cbx(:, i), Cby(:, i), Cbz(:, i), Cbcl(:, i)');
    hold on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 转换到相机坐标

%世界坐标系下的全部
Allp = [Rmp Cbp]; % (3个xyz坐标 * 16个顶点)
Allp = [Allp; ones(1, length(Allp(1,:)))]; %扩张一维 (4个xyz1坐标 * 16个顶点)
%相机坐标系下的全部
Allp = T * Allp; % 转换到相机坐标系中 (4个xyz1坐标 * 16个顶点)
Allp = Allp(1:3, :); % 去掉一维1 (3个xyz坐标 * 16个顶点)
%相机坐标系下三角形
Rm = Rmp2Tr(Allp(:,1:8)); % (3个xyz坐标 * 3个顶点 * 10个三角形)
Cb = Cbp2Tr(Allp(:,9:16)); % (3个xyz坐标 * 3个顶点 * 12个三角形)
All = cat(3, Rm, Cb); % (3个xyz坐标 * 3个顶点 * 22个三角形)


%最终图片 (2cx * 2cy * 3个RGB)
Im = zeros(2*cx, 2*cy, 3);

%遍历 2cx*2cy 个像素点
for i = 1:2*cy %450%
    for j = 1:2*cx %550%
        
        %像素对应方向向量
        Dir = [j; i; 1];
        Dir = K\Dir;
        Dir = Dir/norm(Dir);
        %用于存放与三角形的交点信息
        temp = [0; Crit];
        
        % 遍历22个三角形
        for k = 1:size(All, 3)
            
            %三角形三顶点
            A = [All(1,1,k); All(2,1,k); All(3,1,k)];
            B = [All(1,2,k); All(2,2,k); All(3,2,k)];
            C = [All(1,3,k); All(2,3,k); All(3,3,k)];
            %三角形三边
            Edg1 = B-A;
            Edg2 = C-B;
            Edg3 = A-C;
            %三角形法向量:N (列向量)
            N = cross(Edg1, Edg2);
            N = N/norm(N);
            %三角形平面方程常数项
            Dp = dot(N,A)*(-1);
            %光路与三角形平面交点
            NDir = dot(N, Dir);
            if(NDir ~= 0)
                I = Dir*(-Dp/NDir);
            else
                I = [0;0;-1];
            end
            
            % 交点需在视野范围内 (z>0)
            if(I(3) > crit)
                %三角形三顶点到交点 (列向量)
                AI = I - A;
                BI = I - B;
                CI = I - C;
                %交点需在三角形内 && 交点需更近 (z较小)
                if(dot(N, cross(Edg1, AI)) > crit ...
                    && dot(N, cross(Edg2, BI)) > crit ...
                    && dot(N, cross(Edg3, CI)) > crit ...
                    && I(3) < temp(2))
                    % 记录列向量 (三角形编号, 交点深度)
                    temp = [k; I(3)];
                end
            end
        end
        
        % 三角形编号
        num = temp(1);
        % 给图像上色
        if(num ~= 0)
            Im(i, j, 1) = Allcl(1, num);
            Im(i, j, 2) = Allcl(2, num);
            Im(i, j, 3) = Allcl(3, num);
        end
        
    end
end

% 最终照片
figure(2)
image(Im);
axis([1 2*cx 1 2*cy]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%自定义函数

%将房间的顶点转化为10个三角形 
%(3个xyz坐标 * 8个顶点) -> (3个xyz坐标 * 3个顶点 * 10个三角形)
function Rm = Rmp2Tr(Rmp)
    Rm = zeros(3,3,10);
    Rm(:,:,1) = [Rmp(:,8) Rmp(:,7) Rmp(:,6)];
    Rm(:,:,2) = [Rmp(:,5) Rmp(:,6) Rmp(:,7)];
    Rm(:,:,3) = [Rmp(:,4) Rmp(:,2) Rmp(:,3)];
    Rm(:,:,4) = [Rmp(:,1) Rmp(:,3) Rmp(:,2)];
    Rm(:,:,5) = [Rmp(:,8) Rmp(:,4) Rmp(:,7)];
    Rm(:,:,6) = [Rmp(:,3) Rmp(:,7) Rmp(:,4)];
    Rm(:,:,7) = [Rmp(:,6) Rmp(:,5) Rmp(:,2)];
    Rm(:,:,8) = [Rmp(:,1) Rmp(:,2) Rmp(:,5)];
    Rm(:,:,9) = [Rmp(:,7) Rmp(:,3) Rmp(:,5)];
    Rm(:,:,10) = [Rmp(:,1) Rmp(:,5) Rmp(:,3)];
end

%将正方体的顶点转化为12个三角形 
%(3个xyz坐标 * 8个顶点) -> (3个xyz坐标 * 3个顶点 * 12个三角形)
function Cb = Cbp2Tr(Cbp)
    Cb = zeros(3,3,10);
    Cb(:,:,1) = [Cbp(:,8) Cbp(:,6) Cbp(:,7)];
    Cb(:,:,2) = [Cbp(:,5) Cbp(:,7) Cbp(:,6)];
    Cb(:,:,3) = [Cbp(:,4) Cbp(:,3) Cbp(:,2)];
    Cb(:,:,4) = [Cbp(:,1) Cbp(:,2) Cbp(:,3)];
    Cb(:,:,5) = [Cbp(:,8) Cbp(:,7) Cbp(:,4)];
    Cb(:,:,6) = [Cbp(:,3) Cbp(:,4) Cbp(:,7)];
    Cb(:,:,7) = [Cbp(:,6) Cbp(:,2) Cbp(:,5)];
    Cb(:,:,8) = [Cbp(:,1) Cbp(:,5) Cbp(:,2)];
    Cb(:,:,9) = [Cbp(:,7) Cbp(:,5) Cbp(:,3)];
    Cb(:,:,10) = [Cbp(:,1) Cbp(:,3) Cbp(:,5)];
    Cb(:,:,11) = [Cbp(:,8) Cbp(:,4) Cbp(:,6)];
    Cb(:,:,12) = [Cbp(:,2) Cbp(:,6) Cbp(:,4)];
end
