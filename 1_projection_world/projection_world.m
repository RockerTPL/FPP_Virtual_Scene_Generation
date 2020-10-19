%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%�������
%�����ƽ������
t = [6;2;2.5];
%�������ת�Ƕ�
th_x = -pi*3/5;
th_y = 0;
th_z = pi*1/4; 
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
cy = 253.5; %����307
K = [fx 0 cx; 0 fy cy; 0 0 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%��ά�ռ�

%ÿ1����λȡngap���㣬������С����Ϊgap
ngap = 40;
gap = 1/ngap;

%���䳤���
a = 10;
b = 10;
h = 3;

%����
axe_x = [linspace(0,a,ngap*a) linspace(0,a,ngap*a) linspace(0,a,ngap*a) linspace(0,a,ngap*a)...
    zeros(1,ngap*b) zeros(1,ngap*b) ones(1,ngap*b)*a ones(1,ngap*b)*a ...
    zeros(1,ngap*h) ones(1,ngap*h)*a zeros(1,ngap*h) ones(1,ngap*h)*a];

axe_y = [zeros(1,ngap*a) zeros(1,ngap*a) ones(1,ngap*a)*b ones(1,ngap*a)*b ...
    linspace(0,b,ngap*b) linspace(0,b,ngap*b) linspace(0,b,ngap*b) linspace(0,b,ngap*b)...
    zeros(1,ngap*h) zeros(1,ngap*h) ones(1,ngap*h)*b ones(1,ngap*h)*b];

axe_z = [zeros(1,ngap*a) ones(1,ngap*a)*h zeros(1,ngap*a) ones(1,ngap*a)*h ...
    zeros(1,ngap*b) ones(1,ngap*b)*h zeros(1,ngap*b) ones(1,ngap*b)*h ...
    linspace(0,h,ngap*h) linspace(0,h,ngap*h) linspace(0,h,ngap*h) linspace(0,h,ngap*h)];

%���������棨ƽ��ǰ��
%��һ����
S1x = ones(1,ngap*ngap);
S1y = repmat(linspace(0,1,ngap), 1, ngap); %���Ҹ���ngap��
S1z = repmat(linspace(0,1,ngap), ngap, 1); %�����¸���ngap��
S1z = S1z(:)' ;%ת����������ת��->������
%�ڶ�����
S2x = zeros(1,ngap*ngap);
S2y = repmat(linspace(0,1,ngap), 1, ngap);
S2z = repmat(linspace(0,1,ngap), ngap, 1);
S2z = S2z(:)';
%��������
S3x = repmat(linspace(0,1,ngap), 1, ngap);
S3y = zeros(1,ngap*ngap);
S3z = repmat(linspace(0,1,ngap), ngap, 1);
S3z = S3z(:)';
%���ĸ���
S4x = repmat(linspace(0,1,ngap), 1, ngap);
S4y = ones(1,ngap*ngap);
S4z = repmat(linspace(0,1,ngap), ngap, 1);
S4z = S4z(:)';
%�������
S5x = repmat(linspace(0,1,ngap), ngap, 1);
S5x = S5x(:)';
S5y = repmat(linspace(0,1,ngap), 1, ngap);
S5z = ones(1,ngap*ngap);
%��������
S6x = repmat(linspace(0,1,ngap), ngap, 1);
S6x = S6x(:)';
S6y = repmat(linspace(0,1,ngap), 1, ngap);
S6z = zeros(1,ngap*ngap);

%����ȫ����ƽ�ƺ�
obj_x = [S1x S2x S3x S4x S5x S6x] + 2;
obj_y = [S1y S2y S3y S4y S5y S6y] + 4;
obj_z = [S1z S2z S3z S4z S5z S6z];

%���������棨ƽ�ƺ�
S1x = obj_x(1:ngap^2);
S1y = obj_y(1:ngap^2);
S1z = obj_z(1:ngap^2);
S2x = obj_x(ngap^2+1:2*ngap^2);
S2y = obj_y(ngap^2+1:2*ngap^2);
S2z = obj_z(ngap^2+1:2*ngap^2);
S3x = obj_x(2*ngap^2+1:3*ngap^2);
S3y = obj_y(2*ngap^2+1:3*ngap^2);
S3z = obj_z(2*ngap^2+1:3*ngap^2);
S4x = obj_x(3*ngap^2+1:4*ngap^2);
S4y = obj_y(3*ngap^2+1:4*ngap^2);
S4z = obj_z(3*ngap^2+1:4*ngap^2);
S5x = obj_x(4*ngap^2+1:5*ngap^2);
S5y = obj_y(4*ngap^2+1:5*ngap^2);
S5z = obj_z(4*ngap^2+1:5*ngap^2);
S6x = obj_x(5*ngap^2+1:6*ngap^2);
S6y = obj_y(5*ngap^2+1:6*ngap^2);
S6z = obj_z(5*ngap^2+1:6*ngap^2);

%��ά�ռ�ģ��ͼ
figure(1)
scatter3(axe_x,axe_y,axe_z,'.')
hold on;
scatter3(S1x,S1y,S1z,'.');
hold on;
scatter3(S2x,S2y,S2z,'.');
hold on;
scatter3(S3x,S3y,S3z,'.');
hold on;
scatter3(S4x,S4y,S4z,'.');
hold on;
scatter3(S5x,S5y,S5z,'.');
hold on;
scatter3(S6x,S6y,S6z,'.');

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
% % Pc = Pc(:, Pc(3, :)>0); %���뾵����ά��ͨ����(x,y,z)��Z>0��

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%�����ڸ�����
Pcc = [Pc; Pc(3,:)]; %��ά����(x,y,z,z)
Pcc(1:3,:) = Pcc(1:3,:)./[Pcc(3,:); Pcc(3,:); Pcc(3,:)]; %��ά����(x/z, y/z, 1, z)

%�����Լ��ڵ��ĵ�ĵ���ά���-1

% %ԭʼ������ÿ���㶼��ȫ����������
% for i = 1:length(Pcc(1,:))
%     if(Pcc(4,i)>0)
%         for j = 1:length(Pcc(1,:))
%             if(max(abs(Pcc(1,i)-Pcc(1,j)),abs(Pcc(2,i)-Pcc(2,j))) ...
%                < gap*0.71/Pcc(4,i)  &&  Pcc(4,j)-Pcc(4,i) > gap*5)
%                 Pcc(4,j) = -1;
%             end
%         end
%     end
% end

%�Ż���ͬһ�����ϵĵ㲻�ụ���ڵ������ٱ���ͬ��ĵ㣩
for i = 1:length(obj_x(1,:)) %����
    if(Pcc(4,i)>0)
        for j = [1:i-rem(i,length(S1x)) ...
                length(S1x)*(floor(i/length(S1x))+1)+1:length(Pcc(1,:))]
            if(max(abs(Pcc(1,i)-Pcc(1,j)),abs(Pcc(2,i)-Pcc(2,j))) ...
               < gap*0.6/Pcc(4,i)  &&  Pcc(4,j)-Pcc(4,i) > gap*2)
                Pcc(4,j) = -1;
            end
        end
    end
end
for i = length(obj_x(1,:))+1:length(Pcc(1,:)) %����
    if(Pcc(4,i)>0)
        for j = 1:length(obj_x(1,:))
            if(max(abs(Pcc(1,i)-Pcc(1,j)),abs(Pcc(2,i)-Pcc(2,j))) ...
               < gap*0.6/Pcc(4,i)  &&  Pcc(4,j)-Pcc(4,i) > gap*2)
                Pcc(4,j) = -1;
            end
        end
    end
end


%�ָ��ȡ6����
Ptest = [Pc; Pcc(4,:)];
S1x = Ptest(1,1:ngap^2);
S1y = Ptest(2,1:ngap^2);
S1z = Ptest(3,1:ngap^2);
S1j = Ptest(4,1:ngap^2);
S2x = Ptest(1,ngap^2+1:2*ngap^2);
S2y = Ptest(2,ngap^2+1:2*ngap^2);
S2z = Ptest(3,ngap^2+1:2*ngap^2);
S2j = Ptest(4,ngap^2+1:2*ngap^2);
S3x = Ptest(1,2*ngap^2+1:3*ngap^2);
S3y = Ptest(2,2*ngap^2+1:3*ngap^2);
S3z = Ptest(3,2*ngap^2+1:3*ngap^2);
S3j = Ptest(4,2*ngap^2+1:3*ngap^2);
S4x = Ptest(1,3*ngap^2+1:4*ngap^2);
S4y = Ptest(2,3*ngap^2+1:4*ngap^2);
S4z = Ptest(3,3*ngap^2+1:4*ngap^2);
S4j = Ptest(4,3*ngap^2+1:4*ngap^2);
S5x = Ptest(1,4*ngap^2+1:5*ngap^2);
S5y = Ptest(2,4*ngap^2+1:5*ngap^2);
S5z = Ptest(3,4*ngap^2+1:5*ngap^2);
S5j = Ptest(4,4*ngap^2+1:5*ngap^2);
S6x = Ptest(1,5*ngap^2+1:6*ngap^2);
S6y = Ptest(2,5*ngap^2+1:6*ngap^2);
S6z = Ptest(3,5*ngap^2+1:6*ngap^2);
S6j = Ptest(4,5*ngap^2+1:6*ngap^2);
axe_x = Ptest(1,6*ngap^2+1:length(Ptest(1,:)));
axe_y = Ptest(2,6*ngap^2+1:length(Ptest(1,:)));
axe_z = Ptest(3,6*ngap^2+1:length(Ptest(1,:)));
axe_j = Ptest(4,6*ngap^2+1:length(Ptest(1,:)));

%ȥ�����ڵ��ĵ�
S1x = S1x(S1j>0);
S1y = S1y(S1j>0);
S1z = S1z(S1j>0);
S2x = S2x(S2j>0);
S2y = S2y(S2j>0);
S2z = S2z(S2j>0);
S3x = S3x(S3j>0);
S3y = S3y(S3j>0);
S3z = S3z(S3j>0);
S4x = S4x(S4j>0);
S4y = S4y(S4j>0);
S4z = S4z(S4j>0);
S5x = S5x(S5j>0);
S5y = S5y(S5j>0);
S5z = S5z(S5j>0);
S6x = S6x(S6j>0);
S6y = S6y(S6j>0);
S6z = S6z(S6j>0);
axe_x = axe_x(axe_j>0);
axe_y = axe_y(axe_j>0);
axe_z = axe_z(axe_j>0);

%��ά�ռ�ģ��ͼ���������ϵ�У�
figure(2)
% scatter3(axe_x,axe_y,axe_z,'.')
% hold on;
scatter3(S1x,S1y,S1z,'.');
hold on;
scatter3(S2x,S2y,S2z,'.');
hold on;
scatter3(S3x,S3y,S3z,'.');
hold on;
scatter3(S4x,S4y,S4z,'.');
hold on;
scatter3(S5x,S5y,S5z,'.');
hold on;
scatter3(S6x,S6y,S6z,'.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ͶӰ����
% Pc = Pc./[Pc(3,:); Pc(3,:); Pc(3,:)]; %��ά��һ������
%��ά��ͨ����
S1 = [S1x; S1y; S1z];
S2 = [S2x; S2y; S2z];
S3 = [S3x; S3y; S3z];
S4 = [S4x; S4y; S4z];
S5 = [S5x; S5y; S5z];
S6 = [S6x; S6y; S6z];
axe = [axe_x; axe_y; axe_z];
%��ά��һ������(x,y,1)
S1 = S1./[S1(3,:); S1(3,:); S1(3,:)];
S2 = S2./[S2(3,:); S2(3,:); S2(3,:)];
S3 = S3./[S3(3,:); S3(3,:); S3(3,:)];
S4 = S4./[S4(3,:); S4(3,:); S4(3,:)];
S5 = S5./[S5(3,:); S5(3,:); S5(3,:)];
S6 = S6./[S6(3,:); S6(3,:); S6(3,:)];
axe = axe./[axe(3,:); axe(3,:); axe(3,:)];
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




