%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%�������
%�����ƽ������
t = [6;6;2.5];
%�������ת�Ƕ�
th_x = -pi*3/5;
th_y = 0;
th_z = pi*2/4; 
%�������ת����
R_x = [1 0 0; 0 cos(th_x) -sin(th_x); 0 sin(th_x) cos(th_x)];
R_y = [cos(th_y) 0 sin(th_y); 0 1 0; -sin(th_y) 0 cos(th_y)];
R_z = [cos(th_z) -sin(th_z) 0; sin(th_z) cos(th_z) 0; 0 0 1];
R = R_z * R_y * R_x;
%����ı任������֪������꣬ͨ���������תƽ��Tcw���ɵ��������꣺Pw = Tcw.Pc
Tcw = [R t; 0 0 0 1]; 
%��֪�������꣬��������꣺Pc = Tcw^(-1)Pw
T = inv(Tcw);

%�ڲξ���
fx = 518.0;
fy = 519.0;
cx = 325.5; %����651
cy = 253.5; %����507
K = [fx 0 cx; 0 fy cy; 0 0 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%��ά�ռ�

%ÿ1����λȡngap���㣬������С����Ϊgap
ngap = 50;
gap = 1/ngap;

%���䳤���
a = 10;
b = 10;
h = 3;

%����
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

%����������
%��һ����
S1x = ones(1,ngap*ngap) + 2;
S1y = repmat(linspace(0,1,ngap), 1, ngap) + 4; %���Ҹ���ngap��
S1z = repmat(linspace(0,1,ngap), ngap, 1); %�����¸���ngap��
S1z = S1z(:)' ;%ת����������ת��->������
%�ڶ�����
S2x = zeros(1,ngap*ngap) + 2;
S2y = repmat(linspace(0,1,ngap), 1, ngap) + 4;
S2z = repmat(linspace(0,1,ngap), ngap, 1);
S2z = S2z(:)';
%��������
S3x = repmat(linspace(0,1,ngap), 1, ngap) + 2;
S3y = zeros(1,ngap*ngap) + 4;
S3z = repmat(linspace(0,1,ngap), ngap, 1);
S3z = S3z(:)';
%���ĸ���
S4x = repmat(linspace(0,1,ngap), 1, ngap) + 2;
S4y = ones(1,ngap*ngap) + 4;
S4z = repmat(linspace(0,1,ngap), ngap, 1);
S4z = S4z(:)';
%�������
S5x = repmat(linspace(0,1,ngap), ngap, 1);
S5x = S5x(:)' + 2;
S5y = repmat(linspace(0,1,ngap), 1, ngap) + 4;
S5z = ones(1,ngap*ngap);
%��������
S6x = repmat(linspace(0,1,ngap), ngap, 1);
S6x = S6x(:)' + 2;
S6y = repmat(linspace(0,1,ngap), 1, ngap) + 4;
S6z = zeros(1,ngap*ngap);

%����ȫ��
obj_x = [S1x S2x S3x S4x S5x S6x];
obj_y = [S1y S2y S3y S4y S5y S6y];
obj_z = [S1z S2z S3z S4z S5z S6z];

%��ά�ռ�ģ��ͼ
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

%ת�����������ϵ

%��������ϵ�µ�ȫ��������+���䣩
x = [obj_x axe_x];
y = [obj_y axe_y];
z = [obj_z axe_z];
Pw = [x; y; z]; %��ͨ��ά����(xw,yw,zw)
Pw = [Pw; ones(1, length(x))]; %ת������ά��һ������(xw,yw,zw,1)

%�������ϵ�µ�ȫ��
Pc = T * Pw; %��ά��һ������(x,y,z,1)��<-ת�����������ϵ�У�
Pc = Pc(1:3, :); %��ά��ͨ����(x,y,z)��<-ȥ��һά1��

%�˸�����(������)
A = Pc(1:3,3*ngap^2);
B = Pc(1:3,ngap^2);
C = Pc(1:3,2*ngap^2);
D = Pc(1:3,4*ngap^2+1);
E = Pc(1:3,1);
F = Pc(1:3,6*ngap^2);
G = Pc(1:3,3*ngap^2+1);
H = Pc(1:3,ngap^2+1);
%���㼯�� (xyz, 8����)
Pks = [A B C D E F G H];

%�����μ��� (xyz, 3����, 12������)
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

%�����ڸ�����
Pcc = [Pc; Pc(3,:)]; %��ά����(x,y,z,z)
for i = 1:length(Pcc(1,:))
    if(Pcc(4,i)>0) %�����жϲ��뾵�ĵ�
        for j = 1:12
            %�����η�����:N (������)
            N = cross(Trs(:,2,j)-Trs(:,1,j),Trs(:,3,j)-Trs(:,1,j));
            N = N/norm(N);
            %������ƽ�淽�̳�����
            Dp = -N'*Trs(:,1,j);
            %����㷽������:R (������)
            R = Pc(:,i)/norm(Pc(:,i));
            %ƽ��(N'*R=0)������㲻���ڵ�
            %���ܱ�ͬһƽ���������ڵ�
            if(N'*R < 0 && N'*Pc(:,i)+Dp < 0-gap*0.01) 
                t = -Dp/(N'*R);
                I = R*t; %�������������ƽ���ϵ�ͶӰ��
                %���������� (������)
                Edg1 = Trs(:,2,j)-Trs(:,1,j);
                Edg2 = Trs(:,3,j)-Trs(:,2,j);
                Edg3 = Trs(:,1,j)-Trs(:,3,j);
                %�����������㵽����I (������)
                SI1 = I - Trs(:,1,j);
                SI2 = I - Trs(:,2,j);
                SI3 = I - Trs(:,3,j);
                %���ڵ��ĵ����ά��Ϊ-1
                if(N'*cross(Edg1, SI1) > 0-gap*0.01 ...
                    && N'*cross(Edg2, SI2) > 0-gap*0.01 ...
                    && N'*cross(Edg3, SI3) > 0-gap*0.01)
                    Pcc(4,i) = -1;
                end
            end
        end
    end
end


%�ָ��ȡ6����
Ptest = [Pc; Pcc(4,:)];
S1 = Ptest(:,1:ngap^2);
S2 = Ptest(:,ngap^2+1:2*ngap^2);
S3 = Ptest(:,2*ngap^2+1:3*ngap^2);
S4 = Ptest(:,3*ngap^2+1:4*ngap^2);
S5 = Ptest(:,4*ngap^2+1:5*ngap^2);
S6 = Ptest(:,5*ngap^2+1:6*ngap^2);
axe = Ptest(:,6*ngap^2+1:length(Ptest(1,:)));

%ȥ�����ڵ��ĵ�
S1 = S1(:,S1(4,:)>0);
S2 = S2(:,S2(4,:)>0);
S3 = S3(:,S3(4,:)>0);
S4 = S4(:,S4(4,:)>0);
S5 = S5(:,S5(4,:)>0);
S6 = S6(:,S6(4,:)>0);
axe = axe(:,axe(4,:)>0);

%��ά�ռ�ģ��ͼ���������ϵ�У�
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

%ͶӰ����
%��ά��һ������(x,y,1)
% Pc = Pc./[Pc(3,:); Pc(3,:); Pc(3,:)]; 
S1 = S1(1:3,:)./[S1(3,:); S1(3,:); S1(3,:)];
S2 = S2(1:3,:)./[S2(3,:); S2(3,:); S2(3,:)];
S3 = S3(1:3,:)./[S3(3,:); S3(3,:); S3(3,:)];
S4 = S4(1:3,:)./[S4(3,:); S4(3,:); S4(3,:)];
S5 = S5(1:3,:)./[S5(3,:); S5(3,:); S5(3,:)];
S6 = S6(1:3,:)./[S6(3,:); S6(3,:); S6(3,:)];
axe = axe(1:3,:)./[axe(3,:); axe(3,:); axe(3,:)];
%��������ϵ����Ƭ(u,v,1)
S1 = K * S1; 
S2 = K * S2;
S3 = K * S3;
S4 = K * S4;
S5 = K * S5;
S6 = K * S6;
axe = K * axe;
%����Ұ��Χ�ڵ�����ƽ����ͨ����(u,v)������չʾ��
S1 = S1(1:2, S1(1,:)>0 & S1(1,:)<2*cx & S1(2,:)>0 & S1(2,:)<2*cy);
S2 = S2(1:2, S2(1,:)>0 & S2(1,:)<2*cx & S2(2,:)>0 & S2(2,:)<2*cy);
S3 = S3(1:2, S3(1,:)>0 & S3(1,:)<2*cx & S3(2,:)>0 & S3(2,:)<2*cy);
S4 = S4(1:2, S4(1,:)>0 & S4(1,:)<2*cx & S4(2,:)>0 & S4(2,:)<2*cy);
S5 = S5(1:2, S5(1,:)>0 & S5(1,:)<2*cx & S5(2,:)>0 & S5(2,:)<2*cy);
S6 = S6(1:2, S6(1,:)>0 & S6(1,:)<2*cx & S6(2,:)>0 & S6(2,:)<2*cy);
axe = axe(1:2, axe(1,:)>0 & axe(1,:)<2*cx & axe(2,:)>0 & axe(2,:)<2*cy);
%������ӻ�չʾ����תy��
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
















