clc
Q_s2o = [1,0,0;0,-1,0;0,0,-1]
Q_s2e = [0,-1,0; -1,0,0; 0,0,1]
Q_e2s = [0,-1,0; -1,0,0; 0,0,-1]
Q_e2o = [0,-1,0; 1,0,0; 0,0,1]

disp("Q_e2o = Q_s2o * Q_e2s = ")
disp(num2str(Q_s2o*Q_e2s))

g_old = [-0.3118    0.8183    0.4828;...
   -0.8783   -0.4421    0.1822;...
    0.3625   -0.3672    0.8566];

g_noflag = [0.8183    0.3118    0.4828;...
   -0.4421    0.8783    0.1822;...
   -0.3672   -0.3625    0.8566];

g_new = [0.8183   -0.3118   -0.4828;...
   -0.4421   -0.8783   -0.1822;...
   -0.3672    0.3625   -0.8566];

disp(g_noflag'*g_old)

disp(g_noflag'*g_old * Q_e2s')

K = [0,0,0; 0,0,-1; 0,1,0];

% angle = pi;
% R = eye(3) + sin(angle)*K + (1-cos(angle))*K*K;
% disp("R(pi) = ")
% disp(R)
% 
% angle = -pi;
% R = eye(3) + sin(angle)*K + (1-cos(angle))*K*K;
% disp(" ")
% disp("R(-pi) = ")
% disp(R)
% 
% 
% angle = -pi;
% R = eye(3) + sin(angle)*K + (1-cos(angle))*K*K;
% disp(" ")
% disp("R(-pi) = ")
% disp(R)

%%
sampkleTilt = 10;