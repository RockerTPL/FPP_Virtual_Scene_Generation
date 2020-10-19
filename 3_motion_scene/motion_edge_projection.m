%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%�ڲξ���
global fx;
global fy;
global cx;
global cy;
global ci;
global cj;
global K;
global invK;
fx = 518.0;
fy = 519.0;
cx = 200;  % 325.5; %����651 % c������ĺ���������ϵ֮���ƫ������2c��Ϊ����ƽ��ߴ�
cy = 155;  % 253.5; %����507
ci = 2*cy;
cj = 2*cx;
K = [fx 0 cx; 0 fy cy; 0 0 1];
invK = inv(K);

global crit; %��Сֵ�ж�
global Crit; %����ֵ�ж�
crit = -0.0001;
Crit = 10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%�������
%�������ת�Ƕ�
th_x = -pi*3/5;
th_y = 0;
th_z = pi*0.6/4;%(-0.3)/4;%1.2/4; 
%�������ת����
R_x = [1 0 0; 0 cos(th_x) -sin(th_x); 0 sin(th_x) cos(th_x)];
R_y = [cos(th_y) 0 sin(th_y); 0 1 0; -sin(th_y) 0 cos(th_y)];
R_z = [cos(th_z) -sin(th_z) 0; sin(th_z) cos(th_z) 0; 0 0 1];
R = R_z * R_y * R_x;
%%% �����ƽ������
t = [4.8; 1.3; 2.5];%[7; 1.5; 2.5] + [-0.4; 0.2; 0]*5.5;
%����ı任������֪������꣬ͨ���������תƽ��Tcw���ɵ��������꣺Pw = Tcw.Pc
Tcw = [R t; 0 0 0 1]; 
%��֪�������꣬��������꣺Pc = Tcw^(-1)Pw
T = inv(Tcw);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ��ά�ռ�

% ���䳤���
a = 10;
b = 10;
h = 3;
% ���䶥�� (3��xyz���� * 8������)
Rmp = [[0; 0; 0] [0; 0; h] ...
        [0; b; 0] [0; b; h]...
        [a; 0; 0] [a; 0; h]...
        [a; b; 0] [a; b; h]];
% ���з���������������� (3��xyz���� * 3������ * 10��������)
Rm = Rmp2Tr(Rmp);
    

% ������ƽ������
tcb = [2.5; 5; 0];
% �������ⳤ
c = 1;
% �����嶥�� (3��xyz���� * 8������)
Cbp = [[0; 0; 0] [0; 0; c] ...
        [0; c; 0] [0; c; c]...
        [c; 0; 0] [c; 0; c]...
        [c; c; 0] [c; c; c]]...
        + tcb;% + [0.3; 0.15; 0]*5.5;
% ����������������������� (3��xyz���� * 3������ * 12��������)
Cb = Cbp2Tr(Cbp);

    
% ������ƽ������
tte = [5; 7.5; 0];
% �������ⳤ
c2 = 1.5;
% �����嶥�㣨3��xyz���� * 4�����㣩
Tep = [[0; 0; 0] [0; c2; 0] ...
       [c2/2*sqrt(3); c2/2; 0] ...
       [c/2/sqrt(3); c2/2; sqrt(c2^2 - c2^2/3)]] ...
       + tte;
% ����������������������У�3��xyz���� * 3������ * 4�������Σ�
Te = Tep2Tr(Tep);


% �ѵ���ƽ������
tct = [2.5; 7; 0];
% �ѵ����ⳤ
c3 = 1.2;
% �ѵ��嶥�㣨3��xyz���� * 7�����㣩
Ctp = [[c3/2*sqrt(3); 0; 0] [0; c3/2; 0] [c3/2*sqrt(3); c3; 0] ...
       [c3/2*sqrt(3); 0; c3] [0; c3/2; c3] [c3/2*sqrt(3); c3; c3] ...
       [c3/sqrt(3); c3/sqrt(3); 1.5*c3]] ...
       + tct;
% ���жѵ���������������У�3��xyz���� * 3������ * 10�������Σ�
Ct = Ctp2Tr(Ctp);


% ������ɫ (3��RGB * 10��������)
Rmcl = zeros(3,10);
for i = 1:5
    Rmcl(:, 2*i-1) = [1-i*0.1; 1-i*0.1; 1-i*0.1];
    Rmcl(:, 2*i) = Rmcl(:, 2*i-1);
end
% ��������ɫ (3��RGB * 12��������)
Cbcl = zeros(3,12);
for i = 1:6
    Cbcl(:, 2*i-1) = [1-i*0.15; 0.7-i*0.1; i*0.15];
    Cbcl(:, 2*i) = Cbcl(:, 2*i-1);
end
%��������ɫ��3��RGB * 4�������Σ�
Tecl = zeros(3,4);
for i = 1:4
    Tecl(:, i) = [1-i*0.15; 0.7-i*0.1; i*0.15];
end
%�ѵ�����ɫ��3��RGB * 10�������Σ�
Ctcl = zeros(3,10);
for i = 1:3
    Ctcl(:, 2*i-1) = [1-i*0.15; 0.7-i*0.1; i*0.15];
    Ctcl(:, 2*i) = Ctcl(:, 2*i-1);
end
for i = 7:10
    Ctcl(:, i) = [i*0.07; 0.7-i*0.06; 1-i*0.07];
end



% ȫ�������ε���ɫ���� (3��RGB * nTr��������)
Allcl = [Rmcl, Cbcl, Tecl, Ctcl];

% ���з��������θ�ȡһ���������� (3�����x/y/z���� * 10��������)
% Rmx = zeros(3,10);
Rmx = getxyz(Rm, 10, 'x');
Rmy = getxyz(Rm, 10, 'y');
Rmz = getxyz(Rm, 10, 'z');
% ���������������θ�ȡһ���������� (3�����x/y/z���� * 12��������)
Cbx = getxyz(Cb, 12, 'x');
Cby = getxyz(Cb, 12, 'y');
Cbz = getxyz(Cb, 12, 'z');
% ���������������θ�ȡһ���������� (3�����x/y/z���� * 4��������)
Tex = getxyz(Te, 4, 'x');
Tey = getxyz(Te, 4, 'y');
Tez = getxyz(Te, 4, 'z');
% ���жѵ��������θ�ȡһ���������� (3�����x/y/z���� * 10��������)
Ctx = getxyz(Ct, 10, 'x');
Cty = getxyz(Ct, 10, 'y');
Ctz = getxyz(Ct, 10, 'z');


% ����������ά�ռ�ģ��ͼ
% figure(1)
% for i = 1:10
%     fill3(Rmx(:, i), Rmy(:, i), Rmz(:, i), Rmcl(:, i)');
%     hold on;
% end
% for i = 1:12
%     fill3(Cbx(:, i), Cby(:, i), Cbz(:, i), Cbcl(:, i)');
%     hold on;
% end
% for i = 1:4
%     fill3(Tex(:, i), Tey(:, i), Tez(:, i), Tecl(:, i)');
%     hold on;
% end
% for i = 1:10
%     fill3(Ctx(:, i), Cty(:, i), Ctz(:, i), Ctcl(:, i)');
%     hold on;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ת�����������

%��������ϵ�µ�ȫ��
Allp0 = [Rmp Cbp Tep Ctp]; % (3��xyz���� * nSm������)
Allp0 = [Allp0; ones(1, length(Allp0(1,:)))]; %����һά (4��xyz1���� * nSm������)
%%% �������ϵ�µ�ȫ��
Allp = T * Allp0; % ת�����������ϵ�� (4��xyz1���� * nSm������)
Allp = Allp(1:3, :); % ȥ��һά1 (3��xyz���� * nSm������)
%�������ϵ��������
Rm = Rmp2Tr(Allp(:,1:8));     % (3��xyz���� * 3������ * 10��������)
Cb = Cbp2Tr(Allp(:,9:16));    % (3��xyz���� * 3������ * 12��������)
Te = Tep2Tr(Allp(:,17:20));   % (3��xyz���� * 3������ * 4��������)
Ct = Ctp2Tr(Allp(:,21:27));   % (3��xyz���� * 3������ * 10��������)
All = cat(3, Rm, Cb, Te, Ct); % (3��xyz���� * 3������ * nTr��������)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ���ɳ�ʼͼƬ

%���������θ���
global nTr;
nTr = size(All, 3);
%���������ζ������
global nSm;
nSm = length(Allp0(1,:));
%����ֱ�Ǳߡ�������ı� ����
global nEg;
nEg = 42;


%����ͼƬ (2cx * 2cy * 3��RGB)
Im = zeros(ci, cj, 3);

%ÿ�����ص��Ӧ������
Map = zeros(ci, cj);

%�������ص��Ӧ������
Cover = zeros(ci, cj, nTr);

%������������
List = zeros(1, 1, nTr);


%���� ci * cj �����ص�
for i = 1:ci
    for j = 1:cj
        
        %���ض�Ӧ��������
        Dir = [j; i; 1];
        Dir = invK * Dir;
        Dir = Dir/norm(Dir);
        %���ڴ���������εĽ�����Ϣ
        temp = [0; Crit];
        List = linspace(1, nTr, nTr);
        
        % ����nTr��������
        travTr(i, j, List, All, Cover, temp, Dir, crit, nTr);
        
        % �����α��
        num = temp(1);
        Map(i,j) = num;
        % ��ͼ����ɫ
        if(num ~= 0)
            Im(i, j, 1) = Allcl(1, num);
            Im(i, j, 2) = Allcl(2, num);
            Im(i, j, 3) = Allcl(3, num);
        end
    end
end


% ������Ƭ
figure(2)
image(Im);
axis([1 cj 1 ci]);

















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










% �˶�

gap = 6;%6;
timegap = 0.05;
% ���ƽ���ٶȡ�С��
vc = 0.4;  %0.4
dt = [-vc; vc/2; 0]*timegap;
wc = 0.13;
dth = wc * timegap;
% ������ƽ���ٶȡ�С��
vo = 0.65;  %0.2;
dCbp = [vo; 0; 0]*timegap;


count = 0;

filename = 'move.gif';

for time = 0:timegap:7%3.9
    
%     count = count+1;
    
    %�������ת�Ƕ�
    th_x = -pi*3/5;
    th_y = 0;
    th_z = th_z - dth;%1.2/4; 
    %�������ת����
    R_x = [1 0 0; 0 cos(th_x) -sin(th_x); 0 sin(th_x) cos(th_x)];
    R_y = [cos(th_y) 0 sin(th_y); 0 1 0; -sin(th_y) 0 cos(th_y)];
    R_z = [cos(th_z) -sin(th_z) 0; sin(th_z) cos(th_z) 0; 0 0 1];
    R = R_z * R_y * R_x;
%     %%% �����ƽ������
%     t = t - dt;
    %����ı任������֪������꣬ͨ���������תƽ��Tcw���ɵ��������꣺Pw = Tcw.Pc
    Tcw = [R t; 0 0 0 1]; 
    %��֪�������꣬��������꣺Pc = Tcw^(-1)Pw
    T = inv(Tcw);
    
  
    %��������ϵ�µ�ȫ��
    Cbp = Cbp + dCbp;
    Allp0 = [Rmp Cbp Tep Ctp; ones(1, nSm)]; %����һά (4��xyz1���� * nSm������)
    %%% �������ϵ�µ�ȫ��
    Allp = T * Allp0; % ת�����������ϵ�� (4��xyz1���� * 16������)
    Allp = Allp(1:3, :); % ȥ��һά1 (3��xyz���� * 16������)
    %�������ϵ��������
    Rm = Rmp2Tr(Allp(:,1:8));   % (3��xyz���� * 3������ * 10��������)
    Cb = Cbp2Tr(Allp(:,9:16));  % (3��xyz���� * 3������ * 12��������)
    Te = Tep2Tr(Allp(:,17:20)); % (3��xyz���� * 3������ * 4��������)
    Ct = Ctp2Tr(Allp(:,21:27)); % (3��xyz���� * 3������ * 10��������)
    All = cat(3, Rm, Cb, Te, Ct); % (3��xyz���� * 3������ * nTr��������)

    
    %�߿� 
    frm = zeros(2*gap*(ci+cj),2);
    kk = 1;
    for i = [1:gap, ci-gap+1:ci]
        for j = 1:cj
            frm(kk, :) = [i j];
            kk = kk + 1;
        end
    end
    for j = [1:gap, cj-gap+1:cj]
        for i = 1:ci
            frm(kk,:) = [i j];
            kk = kk + 1;
        end
    end
    
    % �����߿����ص�
    for kkk = 1:kk-1
        i = frm(kkk, 1);
        j = frm(kkk, 2);

        % ���ض�Ӧ��������
        Dir = [j; i; 1];
        Dir = invK * Dir;
        Dir = Dir/(Dir(1)^2 + Dir(2)^2 + Dir(3)^2)^0.5;
        % ���ڴ���������εĽ�����Ϣ
        temp = [0; Crit];
        List = linspace(1, nTr, nTr);
        
        % ��������������
        travTr(i, j, List, All, Cover, temp, Dir, crit, nTr);

        % �����α��
        num = temp(1);
        Map(i,j) = num;
        % ��ͼ����ɫ
        if(num ~= 0)
            Im(i, j, 1) = Allcl(1, num);
            Im(i, j, 2) = Allcl(2, num);
            Im(i, j, 3) = Allcl(3, num);
        end
    end

    
    % Ѱ�ұ�Ե
    ij = [];
    % ������ֱ�Ǳ� (6��xyzxyz���� * nEg)
    Alleg = [Allp(:,1), Allp(:,3), Allp(:,7), Allp(:,5), Allp(:,2), Allp(:,4), ...
             Allp(:,8), Allp(:,6), Allp(:,1), Allp(:,3), Allp(:,5), Allp(:,7), ...
             ...
             Allp(:,9), Allp(:,11), Allp(:,15), Allp(:,13), Allp(:,10), Allp(:,12), ...
             Allp(:,16), Allp(:,14), Allp(:,9), Allp(:,11), Allp(:,13), Allp(:,15), ...
             ...
             Allp(:,17), Allp(:,18), Allp(:,19), Allp(:,17), Allp(:,18), Allp(:,19), ...
             ...
             Allp(:,22), Allp(:,22), Allp(:,23), Allp(:,25), Allp(:,25), Allp(:,26),...
             Allp(:,25), Allp(:,26), Allp(:,24), Allp(:,27), Allp(:,27), Allp(:,27);
             ...
             ...
             Allp(:,3), Allp(:,7), Allp(:,5), Allp(:,1), Allp(:,4), Allp(:,8), ...
             Allp(:,6), Allp(:,2), Allp(:,2), Allp(:,4), Allp(:,6), Allp(:,8), ...
             ...
             Allp(:,11), Allp(:,15), Allp(:,13), Allp(:,9), Allp(:,12), Allp(:,16), ...
             Allp(:,14), Allp(:,10), Allp(:,10), Allp(:,12), Allp(:,14), Allp(:,16), ...
             ...
             Allp(:,18), Allp(:,19), Allp(:,17), Allp(:,20), Allp(:,20), Allp(:,20),...
             ...
             Allp(:,21), Allp(:,23), Allp(:,21), Allp(:,24), Allp(:,26), Allp(:,24),...
             Allp(:,22), Allp(:,23), Allp(:,21), Allp(:,25), Allp(:,26), Allp(:,24)];

         
    for e = 1:nEg
        pt1x = Alleg(1,e);
        pt1y = Alleg(2,e);
        pt1z = Alleg(3,e);
        pt2x = Alleg(4,e);
        pt2y = Alleg(5,e);
        pt2z = Alleg(6,e);
        pt1 = [pt1x; pt1y; pt1z];
        pt2 = [pt2x; pt2y; pt2z];


        %��������ǰ��
        if(pt1z > 0 && pt2z > 0)
            pt1 = pt1 / pt1(3);
            pt1uv = K * pt1;
            pt1uv = pt1uv(1:2);
            pt1uv = round(pt1uv);
            
            pt2 = pt2 / pt2(3);
            pt2uv = K * pt2;
            pt2uv = pt2uv(1:2);
            pt2uv = round(pt2uv);
            
            gpi = abs(pt1uv(2) - pt2uv(2));
            gpj = abs(pt1uv(1) - pt2uv(1));
            gp = max(gpi, gpj);
            idi = round(linspace(pt1uv(2), pt2uv(2), gp))';
            idj = round(linspace(pt1uv(1), pt2uv(1), gp))';
            id0 = [idi, idj];
            id0 = id0(id0(:,1)>=1 & id0(:,1)<=ci & id0(:,2)>=1 & id0(:,2)<=cj, :);
            
            ij = [ij; id0];


        %һ���ں�
        elseif(pt1z > 0 && pt2z <= 0)
            ptdir = pt1 - pt2;
            ptw = - pt2z / (pt1z - pt2z);
            pt2 = pt2 + ptdir * (ptw - crit);

            pt1 = pt1 / pt1(3);
            pt1uv = K * pt1;
            pt1uv = pt1uv(1:2);
            pt1uv = round(pt1uv);
            pt2 = pt2 / pt2(3);
            pt2uv = K * pt2;
            pt2uv = pt2uv(1:2);
            pt2uv = round(pt2uv);
            gpi = abs(pt1uv(2) - pt2uv(2));
            gpj = abs(pt1uv(1) - pt2uv(1));
            gp = max(gpi, gpj);
            idi = round(linspace(pt1uv(2), pt2uv(2), gp))';
            idj = round(linspace(pt1uv(1), pt2uv(1), gp))';
            id0 = [idi, idj];
            id0 = id0(id0(:,1)>=1 & id0(:,1)<=ci & id0(:,2)>=1 & id0(:,2)<=cj, :);
            
            ij = [ij; id0];

        elseif(pt1z <= 0 && pt2z > 0)
            ptdir = pt2 - pt1;
            ptw = - pt1z / (pt2z - pt1z);
            pt1 = pt1 + ptdir * (ptw - crit);

            pt1 = pt1 / pt1(3);
            pt1uv = K * pt1;
            pt1uv = pt1uv(1:2);
            pt1uv = round(pt1uv);
            pt2 = pt2 / pt2(3);
            pt2uv = K * pt2;
            pt2uv = pt2uv(1:2);
            pt2uv = round(pt2uv);
            gpi = abs(pt1uv(2) - pt2uv(2));
            gpj = abs(pt1uv(1) - pt2uv(1));
            gp = max(gpj, gpi);
            idi = round(linspace(pt1uv(2), pt2uv(2), gp))';
            idj = round(linspace(pt1uv(1), pt2uv(1), gp))';
            id0 = [idi, idj];
            id0 = id0(id0(:,1)>=1 & id0(:,1)<=ci & id0(:,2)>=1 & id0(:,2)<=cj, :);
            
            ij = [ij; id0];
            
        end
    end

%     ij = getEg(Alleg, K, nEg, ci, cj, crit, gap);

    % ���� distance-map
%     Dmp = ones(ci, cj) * max(ci, cj) * 2;
%     idmp = sub2ind(size(Dmp), ij(:,1), ij(:,2));
%     Dmp(idmp) = 0;
%     for i = 2:ci
%         for j = 2:cj
%             Dmp(i,j) = min(Dmp(i,j), min(Dmp(i-1,j), Dmp(i, j-1))+1);
%         end
%     end
%     for j = linspace(cj-1, 1, cj-1)
%         Dmp(ci,j) = min(Dmp(ci,j), Dmp(ci, j+1)+1);
%     end
%     for i = linspace(ci-1, 1, ci-1)
%         Dmp(i,cj) = min(Dmp(i,cj), Dmp(i+1, cj)+1);
%     end
%     for i = linspace(ci-1, 1, ci-1)
%         for j = linspace(cj-1, 1, cj-1)
%             Dmp(i,j) = min(Dmp(i,j), min(Dmp(i+1,j), Dmp(i, j+1))+1);
%         end
%     end
%     
%     for i = 1:ci
%         for j = 1:cj
%             if Dmp(i,j) <= gap
%                 ij = [ij; [i j]];
%             end
%         end
%     end

    sizeij = size(ij, 1);
    ij = bigEg(ij, sizeij, ci, cj, gap);

    ij = unique(ij, 'rows');
    
    
% %     % ��Ե��Ϊ��ɫ��ʾ
% %     Imtst = Im;
% %     idr = sub2ind(size(Imtst), ij(:,1), ij(:,2), ones(size(ij(:,1))));
% %     idg = sub2ind(size(Imtst), ij(:,1), ij(:,2), ones(size(ij(:,1)))*2);
% %     idb = sub2ind(size(Imtst), ij(:,1), ij(:,2), ones(size(ij(:,1)))*3);
% %     idrgb = [idr idg idb];
% %     Imtst(idrgb) = 0;
% %     figure(4)
% %     image(Imtst);
% %     axis([1 cj 1 ci]);
    
    
    % ������Ե���ص�
    kk = size(ij,1); 
    for kkk = 1:kk
        i = ij(kkk, 1);
        j = ij(kkk, 2);
            
        % ���ڴ�����������α��
%         List = zeros(1, nTr);
        List = Cover(i, j, :);

        % �����Լ������������α�� 
        getNeibTr(i, j, Map, Cover, List, ci, cj, gap, nTr);
        
        
        Cover(i, j, :) = 0;
%         List = linspace(1,nTr,nTr);

        % ���ض�Ӧ��������
        Dir = [j; i; 1];
        Dir = K\Dir;
        Dir = Dir/(Dir(1)^2 + Dir(2)^2 + Dir(3)^2)^0.5;
        % ���ڴ���������εĽ�����Ϣ
        temp = [0; Crit];
        
        % ��������������
        travTr(i, j, List, All, Cover, temp, Dir, crit, nTr);

        % �����α��
        num = temp(1);
        Map(i,j) = num;
        % ��ͼ����ɫ
        if(num ~= 0)
            Im(i, j, 1) = Allcl(1, num);
            Im(i, j, 2) = Allcl(2, num);
            Im(i, j, 3) = Allcl(3, num);
        end
    end
    
    pause(0.0001)
    image(Im);
    
    
%     drawnow;
%     F = getframe(gcf);
%     I = frame2im(F);
%     [I,map] = rgb2ind(I,256);
%     if time == 0
%         imwrite(I,map,'test.gif','gif', 'Loopcount',inf,'DelayTime',0.08);
%     else
%         imwrite(I,map,'test.gif','gif','WriteMode','append','DelayTime',0.08);
%     end

    
end
    






















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%�Զ��庯��

%������Ķ���ת��Ϊ10�������� 
%(3��xyz���� * 8������) -> (3��xyz���� * 3������ * 10��������)
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

%��������Ķ���ת��Ϊ12�������� 
%(3��xyz���� * 8������) -> (3��xyz���� * 3������ * 12��������)
function Cb = Cbp2Tr(Cbp)
    Cb = zeros(3,3,12);
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

%��������Ķ���ת��Ϊ4�������� 
%(3��xyz���� * 8������) -> (3��xyz���� * 3������ * 4��������)
function Te = Tep2Tr(Tep)
    Te = zeros(3,3,4);
    Te(:,:,1) = [Tep(:,1) Tep(:,2) Tep(:,3)];
    Te(:,:,2) = [Tep(:,1) Tep(:,2) Tep(:,4)];
    Te(:,:,3) = [Tep(:,2) Tep(:,3) Tep(:,4)];
    Te(:,:,4) = [Tep(:,3) Tep(:,1) Tep(:,4)];
end

%���ѵ���Ķ���ת��Ϊ10�������� 
%(3��xyz���� * 8������) -> (3��xyz���� * 3������ * 10��������)
function Ct = Ctp2Tr(Ctp)
    Ct = zeros(3,3,10);
    Ct(:,:,1) = [Ctp(:,3) Ctp(:,6) Ctp(:,1)];
    Ct(:,:,2) = [Ctp(:,4) Ctp(:,1) Ctp(:,6)];
    Ct(:,:,3) = [Ctp(:,3) Ctp(:,2) Ctp(:,6)];
    Ct(:,:,4) = [Ctp(:,5) Ctp(:,6) Ctp(:,2)];
    Ct(:,:,5) = [Ctp(:,4) Ctp(:,5) Ctp(:,1)];
    Ct(:,:,6) = [Ctp(:,2) Ctp(:,1) Ctp(:,5)];
    Ct(:,:,7) = [Ctp(:,1) Ctp(:,2) Ctp(:,3)];
    Ct(:,:,8) = [Ctp(:,7) Ctp(:,5) Ctp(:,4)];
    Ct(:,:,9) = [Ctp(:,7) Ctp(:,4) Ctp(:,6)];
    Ct(:,:,10) = [Ctp(:,7) Ctp(:,6) Ctp(:,5)];
end

% �����������μ��ϻ�ȡ������x/y/z���꼯
% (3�����x/y/z���� * numTr��������)
function Objxyz = getxyz(Obj,numTr, xyz)
    Objxyz = zeros(3,numTr);
    switch xyz
        case 'x'
            for i = 1:numTr
                Objxyz(:, i) = Obj(1, :, i)';
            end
        case 'y'
            for i = 1:numTr
                Objxyz(:, i) = Obj(2, :, i)';
            end
        case 'z'
            for i = 1:numTr
                Objxyz(:, i) = Obj(3, :, i)';
            end
    end
end



