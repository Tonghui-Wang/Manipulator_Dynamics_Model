clc
clear

syms a1 a2 a3  %关节DH参数a
syms d1 d2 d3  %关节DH参数d
syms y1 y2 yc1 yc2 yc3  %关节1y方向偏移及重心y方向偏移
syms ac1 ac2 ac3  %质心DH参数a  相对于J1、J2、J3
syms dc1 dc2 dc3  %质心DH参数d  相对于J1、J2、J3
syms q1 q2 q3 real  %关节DH参数q
syms dq1 dq2 dq3 real   %关节速度
syms ddq1 ddq2 ddq3 real  %关节加速度
syms m1 m2 m3 real	  %杆质量
syms Ix1 Iy1 Iz1 Ix2 Iy2 Iz2 Ix3 Iy3 Iz3 g0 real  %质心的转动惯量及重力加速度

%三连杆变换矩阵
c1=cos(q1);
s1=sin(q1);
c2=cos(q2);
s2=sin(q2);
c3=cos(q3);
s3=sin(q3);

%三连杆模型
A1=[c1 0 s1 a1*c1-y1*s1
    s1 0 -c1 a1*s1+y1*c1;
    0  1   0    d1;
    0  0   0     1];
A2=[c2 -s2 0 a2*c2-y2*s2;
    s2 c2 0 a2*s2+y2*c2;
    0  0   1   d2;
    0  0   0   1];
A3=[c3 -s3 0 a3*c3;
    s3 c3 0 a3*s3;
    0 0 1 d3;
    0 0 0 1];

%连杆2及连杆3末端变换矩阵
A12=A1*A2;
A123=A1*A2*A3;

%三连杆位姿变换矩阵
R1=A1(1:3,1:3);
R2=A2(1:3,1:3);
R3=A3(1:3,1:3);
R12=R1*R2;
R13=R1*R2*R3;

%基座标矩阵及向量g
z0=[0 0 1]';
z1=R1*z0;
z2=R1*R2*z0;

%基座标系、关节1末端、关节2末端、关节3末端位置
P0=[0 0 0]';
P1=A1(1:3,4);
P2=A12(1:3,4);
P3=A123(1:3,4);

%关节臂1质心、关节臂2质心、关节臂3质心位置
Ac1=[c1 0 s1 ac1*c1-yc1*s1
    s1 0 -c1 ac1*s1+yc1*c1;
    0  1   0    dc1;
    0  0   0     1];
Ac2=[c2 -s2 0 ac2*c2-yc2*s2;
    s2 c2 0 ac2*s2+yc2*c2;
    0  0   1   dc2;
    0  0   0   1];
Ac3=[c3 -s3 0 ac3*c3-yc3*s3;
    s3 c3 0 ac3*s3+yc3*c3;
    0  0   1   dc3;
    0  0   0   1];
Ac12=A1*Ac2;
Ac123=A1*A2*Ac3;
Pc1=Ac1(1:3,4);
Pc2=Ac12(1:3,4);
Pc3=Ac123(1:3,4);

%定义矩阵a为零向量，填充雅克比矩阵中的零项
a=[0;0;0];

%计算雅克比矩阵
J3=[cross(z0,(Pc3-P0)) cross(z1,(Pc3-P1)) cross(z2,(Pc3-P2)); z0 z1 z2];
J3=simplify(J3);
J2=[cross(z0,(Pc2-P0)) cross(z1,(Pc2-P1)) cross(z1,(P2-P2)); z0 z1 a];
J2=simplify(J2);
J1=[cross(z0,(Pc1-P0)) cross(z1,(P2-P2)) cross(z1,(P1-P1)); z0 a a];
J1=simplify(J1);

%雅克比速度矩阵及加速度矩阵
Jp3=J3(1:3,:);
Jp2=J2(1:3,:);
Jp1=J1(1:3,:);
Jo3=J3(4:6,:);
Jo2=J2(4:6,:);
Jo1=J1(4:6,:);

%关节臂1、2及3的转动惯量
I1=[Ix1 0 0; 0 Iy1 0; 0 0 Iz1];
I2=[Ix2 0 0; 0 Iy2 0; 0 0 Iz2];
I3=[Ix3 0 0; 0 Iy3 0; 0 0 Iz3];

%关节1、关节2、关节3的位置、速度、加速度
q=[q1;q2;q3];
dq=[dq1;dq2;dq3];
ddq=[ddq1;ddq2;ddq3];

%关节臂1、关节臂2、关节臂3的质量及转动惯量矩阵，即惯性矩阵
Dq1=simplify(m1*Jp1'*Jp1+Jo1'*R1*I1*R1'*Jo1);
Dq2=simplify(m2*Jp2'*Jp2+Jo2'*R12*I2*R12'*Jo2);
Dq3=simplify(m3*Jp3'*Jp3+Jo3'*R13*I3*R13'*Jo3);
Dq=simplify(Dq1+Dq2+Dq3);

d11=simplify(Dq(1,1));
d12=simplify(Dq(1,2));
d13=simplify(Dq(1,3));
d21=simplify(Dq(2,1));
d22=simplify(Dq(2,2));
d23=simplify(Dq(2,3));
d31=simplify(Dq(3,1));
d32=simplify(Dq(3,2));
d33=simplify(Dq(3,3));

%动能
Ke=simplify(1/2*dq'*Dq*dq);

%势能
g=[0 0 -g0]';
Po1=-m1*g'*Pc1;
Po2=-m2*g'*Pc2;
Po3=-m3*g'*Pc3;
Po=simplify(Po1+Po2+Po3);

%仅涉及q而与导数无关的项-重力矩
gq1=diff(Po,q1);
gq2=diff(Po,q2);
gq3=diff(Po,q3);
gq=[gq1;gq2;gq3];

%计算(第一类)Christoffel符号cijk
c111=simplify((diff(d11,q1)+diff(d11,q1)-diff(d11,q1))/2);
c121=simplify((diff(d12,q1)+diff(d11,q2)-diff(d12,q1))/2);
c131=simplify((diff(d13,q1)+diff(d11,q3)-diff(d13,q1))/2);
c211=simplify((diff(d11,q2)+diff(d12,q1)-diff(d21,q1))/2);
c221=simplify((diff(d12,q2)+diff(d12,q2)-diff(d22,q1))/2);
c231=simplify((diff(d13,q2)+diff(d12,q3)-diff(d23,q1))/2);
c311=simplify((diff(d13,q1)+diff(d11,q3)-diff(d31,q1))/2);
c321=simplify((diff(d12,q3)+diff(d13,q2)-diff(d32,q1))/2);
c331=simplify((diff(d13,q3)+diff(d13,q3)-diff(d33,q1))/2);
c112=simplify((diff(d21,q1)+diff(d21,q1)-diff(d11,q2))/2);
c122=simplify((diff(d22,q1)+diff(d21,q2)-diff(d12,q2))/2);
c132=simplify((diff(d23,q1)+diff(d21,q3)-diff(d13,q2))/2);
c212=simplify((diff(d21,q2)+diff(d22,q1)-diff(d21,q2))/2);
c222=simplify((diff(d22,q2)+diff(d22,q2)-diff(d22,q2))/2);
c232=simplify((diff(d23,q2)+diff(d22,q3)-diff(d23,q2))/2);
c312=simplify((diff(d21,q3)+diff(d23,q1)-diff(d31,q2))/2);
c322=simplify((diff(d22,q3)+diff(d23,q2)-diff(d32,q2))/2);
c332=simplify((diff(d23,q3)+diff(d23,q3)-diff(d33,q2))/2);
c113=simplify((diff(d31,q1)+diff(d31,q1)-diff(d11,q3))/2);
c123=simplify((diff(d32,q1)+diff(d31,q2)-diff(d12,q3))/2);
c133=simplify((diff(d33,q1)+diff(d31,q3)-diff(d13,q3))/2);
c213=simplify((diff(d31,q2)+diff(d32,q1)-diff(d21,q3))/2);
c223=simplify((diff(d32,q2)+diff(d32,q2)-diff(d22,q3))/2);
c233=simplify((diff(d33,q2)+diff(d32,q3)-diff(d23,q3))/2);
c313=simplify((diff(d31,q3)+diff(d33,q1)-diff(d31,q3))/2);
c323=simplify((diff(d32,q3)+diff(d33,q2)-diff(d32,q3))/2);
c333=simplify((diff(d33,q3)+diff(d33,q3)-diff(d33,q3))/2);

C11=c111*dq1+c121*dq2+c131*dq3;
C12=c211*dq1+c221*dq2+c231*dq3;
C13=c311*dq1+c321*dq2+c331*dq3;
C21=c112*dq1+c122*dq2+c132*dq3;
C22=c212*dq1+c222*dq2+c232*dq3;
C23=c312*dq1+c322*dq2+c332*dq3;
C31=c113*dq1+c123*dq2+c133*dq3;
C32=c213*dq1+c223*dq2+c233*dq3;
C33=c313*dq1+c323*dq2+c333*dq3;

C=[ C11 C12 C13;
    C21 C22 C23; 
    C31 C32 C33];

%动力学方程
tau11=simplify(Dq*ddq);
tau12=simplify(C*dq);
tau13=simplify(gq);

tau=simplify(Dq*ddq+C*dq+gq);
tau1=simplify(tau(1))
tau2=simplify(tau(2))
tau3=simplify(tau(3))
