%%
clc;
clear;
close all;

% Settling time for x and theta of less than 5 seconds
% Rise time for x of less than 0.5 seconds
% Pendulum angle theta never more than 20 degrees (0.35 radians) from the vertical
% Steady-state error of less than 2% for x and theta

%%  system parameter
M = 0.5 * 1.1;
m = 0.2* 1.1;
b = 0.1* 1.1;
I = 0.006* 1.1;
g = 9.8* 1.1;
l = 0.3* 1.1;
q = (M+m)*(I+m*l^2)-(m*l)^2;
s = tf('s');

p = I*(M+m)+M*m*l^2; %denominator for the A and B matrices

F = [0      1              0           0;
     0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0;
     0      0              0           1;
     0 -(m*l*b)/p       m*g*l*(M+m)/p  0];
G = [     0;
     (I+m*l^2)/p;
          0;
        m*l/p];
    
H = [1 0 0 0;
     0 0 1 0];
J = [0;
     0];

 H_pos = H(1,:);
 H_angle = H(2,:);
 
states = {'x' 'x_dot' 'phi' 'phi_dot'};
inputs = {'u'};
outputs = {'x'; 'phi'};

sys = ss(F,G,H,J,'statename',states,'inputname',inputs,'outputname',outputs);


poles = eig(F);        % pole
transFc = tf(sys);           % system 

zero_Group = transFc.Numerator;
pole_Group = transFc.Denominator;

x_trans_num = conv(zero_Group{1,1},zero_Group{1,1});
x_trans_den = conv(pole_Group{1,1},pole_Group{1,1});
sys_x = tf(x_trans_num,x_trans_den);

angle_trans_num = conv(zero_Group{2,1},zero_Group{2,1});
angle_trans_den = conv(pole_Group{2,1},pole_Group{2,1});
sys_angle = tf(angle_trans_num,angle_trans_den);

Sc = ctrb(sys);        %controlalbe Matrix   ---> 랭크 4  컨트롤러블 
So = obsv(sys);        %Observable Matrix    ---> 랭크 3  옵버버블 하지 않다. Cuz 랭크가 4보다 작기 때문이다.(angle)
                                            %---> 랭크 4  옵버버블 하다(x)
rank(Sc)
rank(So)


%% SRL(LQR) design to determine K

% R = 1 and Q = C'C


Q = H'*H;

%      1     0     0     0      ---> Weight's for Car's Position 
%      0     0     0     0      
%      0     0     1     0      ---> Weight's Car's Angle
%      0     0     0     0

% Q의 가중치값
Q(1,1) = 5000;
Q(3,3) = 100;

R = 1;
K_lqr = lqr(F,G,Q,R)



% Q와 R의 비율이 중요하다!   --- >p 값임  putting more weight on the errors at the cost of increased control effort
% 그래서 개쩌는 K값을 실험하면서 찾아보자는 것!

%% close Loop representation

Fc = [(F-G*K_lqr)];
Gc = [G];
Hc = [H];
Jc = [J];


states = {'x' 'x_dot' 'phi' 'phi_dot'};
inputs = {'r'};
outputs = {'x'; 'phi'};

sys_cl = ss(Fc,Gc,Hc,Jc,'statename',states,'inputname',inputs,'outputname',outputs);


% visualization
t = 0:0.01:5;
r =0.2*ones(size(t));
[y,t,x]=lsim(sys_cl,r,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','cart position (m)')
set(get(AX(2),'Ylabel'),'String','pendulum angle (radians)')
title('Step Response with LQR Control')


Hn = [1 0 0 0];                 % 레퍼런스가 오직 Car postion만 있을경우 (angle을 빼버리고 생각합시다)
sys_ss = ss(F,G,Hn,0);
Nbar = rscale(sys_ss,K_lqr);


%%
sys_cl = ss(Fc,Gc*Nbar,Hc,Jc,'statename',states,'inputname',inputs,'outputname',outputs);

t = 0:0.01:5;
r =15*ones(size(t));
[y,t,x]=lsim(sys_cl,r,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','cart position (m)')
set(get(AX(2),'Ylabel'),'String','pendulum angle (radians)')
title('Step Response with Precompensation and LQR Control')


poles = eig(Fc);



%%
s1 = -1.44 - 1.92i;
s2 = -1.44 + 1.92i;
s3 = -1.44*10;
s4 = -1.44*15;

S = [s1 ; s2 ; s3 ; s4];
P = S*10;
L = place(F',H',P)';

P_lqr = [-40 -41 -42 -43]
L_lqr = place(F',H',P)';
% 
% Faa = 0;
% Fab = [ 1 0 0 ];
% Fba = [ 0; 0; 0];
% Fbb = [ -0.1818 2.6727 0;
%          0     0   1 ;
%          -0.4545 31.1818 0];
% Ga = 0;
% Gb = [ 1.8182 ; 0 ; 4.5455];
% pe = [s1 s2 s3]

Faa = [0 0;0 2.6727];
Fab = [ 1 0 ;-0.1818 0];
Fba = [ 0 0; 0 31.1818];
Fbb = [ 0 1;
        -0.4545 0];    
Ga = [0; 0];
Gb = [ 1.8182 ; 4.5455];
pe = 10*[s1 s2]
Lt_re = place(Fbb', Fab', pe)
L_re = Lt_re'

save('result.m');



%% 흐지 코드

NF = [0 [1 0 0 0] ;[0 ;0 ;0 ;0] F];
NG = [0 ;G];
   NH=[0 1 0 0 0 ];
%        0 0 0 0 1 0] ;
% NH = H
NJ = [0]


%use LQR Poles

Npoles = eig(Fc)'
% inputpole;
%  NPole2 = [-1.4 ; -14 ; -28;-56 ]
 NPole1 = [ -8.4910 - 7.9283i  -8.4910 + 7.9283i  -4.7592 - 0.8309i  -4.7592 + 0.8309i -47 ]
% K_itg = acker(NF,NG, [Npoles 47.59 47.59*1.2])
 K_itg = acker(NF,NG, NPole1);

K_itg0 = [K_itg(1)]
K_itg1 = [K_itg(2) K_itg(3) K_itg(4) K_itg(5) ]

