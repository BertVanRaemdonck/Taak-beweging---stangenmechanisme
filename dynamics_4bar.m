%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kinematica en werkuigendynamica.
%
% Voorbeeldanalyse van een vierstangenmechanisme.
%
% Bram Demeulenaere <bram.demeulenaere@mech.kuleuven.be>
% Maarten De Munck <maarten.demunck@mech.kuleuven.be>
% Johan Rutgeerts <johan.rutgeerts@mech.kuleuven.be>
% Wim Meeussen <wim.meeussen@mech.kuleuven.be>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [F12x, F12y, F23x, F23y, F212x, F212y, F34x, F34y, F14x, F14y, F45, F56x, F56y, F67x, F67y, ...
    F68x, F68y, F17x, F17y, F89x, F89y, F810x, F810y, F19, F1011x, F1011y, F1112x, F1112y, F111, ...
    M12, M19, M111, M45] ...
    = dynamics_4bar(phi2,  phi3,  phi4,  x5,  phi6,  phi7,  phi8,  x9,  phi10,  x11,  phi12, ...
                    dphi2, dphi3, dphi4, dx5, dphi6, dphi7, dphi8, dx9, dphi10, dx11, dphi12, ...
                    ddphi2,ddphi3,ddphi4,ddx5,ddphi6,ddphi7,ddphi8,ddx9,ddphi10,ddx11,ddphi12, ...
                    r2l, r2k, r3, a, b, r6l, r6k, r7, r8l, r8k, r10, r11, r12, x4, y4, x7, y7, y9, ...
                    m2,m3,ma,mb,m4,m5,m6k,m6l,m6,m7,m8k,m8l,m8,m9,m10,m11,m12, mpiston1, mpiston2,...
                    X2,X3,X4,X5,X6k,X6l,X6,X7,X8k,X8l,X8,X9,X10,X11,X12, ...
                    Y2,Y3,Y4,Y5,Y6k,Y6l,Y6,Y7,Y8k,Y8l,Y8,Y9,Y10,Y11,Y12, ...
                    J2,J3,J4,J5,J6k,J6l,J6,J7,J8k,J8l,J8,J9,J10,J11,J12, t,fig_dyn_4bar);


% a lot of definitions to make the matrix A and B a bit clear.
% skip the definitions for now (move down to "force analysis")
% and check them when you need them.

% Declaratie nieuwe variabelen:
g = 9.81;           % valversnelling


% cogi_P_x, cogn_P_y = vector from the centre of gravity of bar i to point P
cog2_23_x = -r2k*cos(phi2-pi/2);        % voor stang 2 gerekend vanaf vaste punt
cog2_23_y = r2k*sin(phi2-pi/2);
cog2_212_x = r2l*cos(phi2);
cog2_212_y = -r2l*sin(phi2);

cog3_23_x = X3*cos(phi3);               % voor stang 3 gerekend vanaf massacentrum
cog3_23_y = X3*sin(phi3);
cog3_34_x = -(r3-X3)*cos(phi3);
cog3_34_y = -(r3-X3)*sin(phi3);

proj_45_x = cos(phi4-pi/2);             % projectie van kracht F45
proj_45_y = sin(phi4-pi/2);

cog4_45 = x5;                           % voor stang 4 gerekend vanaf vaste punt
cog4_34_x = (b*cos(phi4)) + (a*cos(phi4+pi/2));
cog4_34_y = (b*sin(phi4)) - (a*cos(phi4+pi/2));

cog6_56_x = -(X6-r6k)*cos(phi6);        % voor stang 6 gerekend vanaf massacentrum
cog6_56_y = -(X6-r6k)*sin(phi6);
cog6_67_x = -X6*cos(phi6);
cog6_67_y = -X6*sin(phi6);
cog6_68_x = (r6k + r6l - X6)*cos(phi6);
cog6_68_y = (r6k + r6l - X6)*sin(phi6);

cog7_17_x = -X7*cos(phi7);              % voor stang 7 gerekend vanaf massacentrum
cog7_17_y = X7*sin(phi7);
cog7_67_x = (r7-X7)*cos(phi7);
cog7_67_y = -(r7-X7)*sin(phi7);

cog8_68_x = (r8-X8)*cos(phi8);          % voor stang 8 gerekend vanaf massacentrum
cog8_68_y = (r8-X8)*sin(phi8);
cog8_89_x = X9*cos(phi8);
cog8_89_y = X9*sin(phi8);
cog8_810_x = -X8*cos(phi8);
cog8_810_y = -X8*sin(phi8);

cog10_810_x = (r10-X10)*cos(phi10);     % voor stang 10 gerekend vanaf massacenrum
cog10_810_y = -(r10-X10)*sin(phi10);
cog10_1011_x = -X10*cos(phi10);
cog10_1011_y = X10*sin(phi10);

cog11_111_y = X11;                      % voor stang 11 gerekend vanaf massacentrum
cog11_1011_y = -(r11-X11);

cog12_212_x = -X12*cos(phi12);          % voor stang 12 gerekend vanaf massacentrum
cog12_212_y = -X12*sin(phi12);
cog12_1112_x = (r12-X12)*cos(phi12);
cog12_1112_y = (r12-X12)*sin(phi12);


% Voorbeeldcode:
% cog2_P_x = -X2*cos(phi2)-Y2*cos(phi2+pi/2);
% cog2_P_y = -X2*sin(phi2)-Y2*sin(phi2+pi/2);
% cog2_Q_x = (r2-X2)*cos(phi2)-Y2*cos(phi2+pi/2);
% cog2_Q_y = (r2-X2)*sin(phi2)-Y2*sin(phi2+pi/2);
% cog3_Q_x = -Y3*cos(phi3+pi/2)-X3*cos(phi3);
% cog3_Q_y = -Y3*sin(phi3+pi/2)-X3*sin(phi3);
% cog3_R_x = (r3-X3)*cos(phi3)-Y3*cos(phi3+pi/2);
% cog3_R_y = (r3-X3)*sin(phi3)-Y3*sin(phi3+pi/2);
% cog4_S_x = -X4*cos(phi4)-Y4*cos(phi4+pi/2);
% cog4_S_y = -X4*sin(phi4)-Y4*sin(phi4+pi/2);
% cog4_R_x = (r4-X4)*cos(phi4)-Y4*cos(phi4+pi/2);
% cog4_R_y = (r4-X4)*sin(phi4)-Y4*sin(phi4+pi/2);

% 3D omega (dphi) and alpha (ddphi) vectors)    NOG AAN TE PASSEN!!!!
omega2 = [zeros(size(phi2)) zeros(size(phi2)) dphi2];
omega3 = [zeros(size(phi2)) zeros(size(phi2)) dphi3];
omega4 = [zeros(size(phi2)) zeros(size(phi2)) dphi4];
alpha2 = [zeros(size(phi2)) zeros(size(phi2)) ddphi2];
alpha3 = [zeros(size(phi2)) zeros(size(phi2)) ddphi3];
alpha4 = [zeros(size(phi2)) zeros(size(phi2)) ddphi4];

% 3D model vectors    NOG AAN TE PASSEN!!!!
P_cog2_vec = [-cog2_P_x    -cog2_P_y    zeros(size(phi2))];
Q_cog3_vec = [-cog3_Q_x    -cog3_Q_y    zeros(size(phi2))];
S_cog4_vec = [-cog4_S_x    -cog4_S_y    zeros(size(phi2))];
PQ_vec = [r2*cos(phi2) r2*sin(phi2) zeros(size(phi2))];

% acceleration vectors    NOG AAN TE PASSEN!!!
acc_2 =       cross(omega2,cross(omega2,P_cog2_vec))+cross(alpha2,P_cog2_vec);
acc_Q =       cross(omega2,cross(omega2,PQ_vec    ))+cross(alpha2,PQ_vec    );
acc_3 = acc_Q+cross(omega3,cross(omega3,Q_cog3_vec))+cross(alpha3,Q_cog3_vec);
acc_4 =       cross(omega4,cross(omega4,S_cog4_vec))+cross(alpha4,S_cog4_vec);
acc_2x = acc_2(:,1);
acc_2y = acc_2(:,2);
acc_3x = acc_3(:,1);
acc_3y = acc_3(:,2);
acc_4x = acc_4(:,1);
acc_4y = acc_4(:,2);


% **********************
% *** force analysis ***
% **********************

% allocate matrices for force (F) and moment (M)
F12x = zeros(size(phi2));
F12y = zeros(size(phi2));
F23x = zeros(size(phi2));
F23y = zeros(size(phi2));
F212x = zeros(size(phi2));
F67y = zeros(size(phi2));
F68x = zeros(size(phi2));
F212y = zeros(size(phi2));
F34x = zeros(size(phi2));
F34y = zeros(size(phi2));
F14x = zeros(size(phi2));
F14y = zeros(size(phi2));
F45 = zeros(size(phi2));
F56x = zeros(size(phi2));
F56y = zeros(size(phi2));
F67x = zeros(size(phi2));
F68y = zeros(size(phi2));
F17x = zeros(size(phi2));
F17y = zeros(size(phi2));
F89x = zeros(size(phi2));
F89y = zeros(size(phi2));
F810x = zeros(size(phi2));
F810y = zeros(size(phi2));
F19 = zeros(size(phi2));
F1011x = zeros(size(phi2));
F1011y = zeros(size(phi2));
F1112x = zeros(size(phi2));
F1112y = zeros(size(phi2));
F111 = zeros(size(phi2));
M12 = zeros(size(phi2));
M19 = zeros(size(phi2));
M111 = zeros(size(phi2));
M45 = zeros(size(phi2));


% calculate dynamics for each time step
t_size = size(t,1);    % number of simulation steps
for k=1:t_size
  A = [ 1           0            1             0             1              0            0             0                  0          0           0            0             0            0             0             0             0            0           0           0           0           0           0           0           0           0           0           0           0           0           0            0           0;
        0           1            0             1             0              1            0             0                  0          0           0            0             0            0             0             0             0            0           0           0           0           0           0           0           0           0           0           0           0           0           0            0           0;
        0           0            cog2_23_y(k) -cog2_23_x(k)  cog2_212_y(k) -cog2_212_x   0             0                  0          0           0            0             0            0             0             0             0            0           0           0           0           0           0           0           0           0           0           0           0           1           0            0           0;
        0           0            0             0             0              0            0             0                  0          0           0            0             0            0             0             0             0            0           0           0           0           0           0           0           0           0           0           0           0           0           0            0           0;
        0           0           -1             0             0              0           -1             0                  0          0           0            0             0            0             0             0             0            0           0           0           0           0           0           0           0           0           0           0           0           0           0            0           0;
        0           0            0            -1             0              0            0            -1                  0          0           0            0             0            0             0             0             0            0           0           0           0           0           0           0           0           0           0           0           0           0           0            0           0;
        0           0            cog3_23_y(k) -cog3_23_x(k)  0              0            cog3_34_y(k) -cog3_34_x(k)       0          0           0            0             0            0             0             0             0            0           0           0           0           0           0           0           0           0           0           0           0           0           0            0           0;
        0           0            0             0             0              0            1             0                  1          0           proj_45_x(k) 0             0            0             0             0             0            0           0           0           0           0           0           0           0           0           0           0           0           0           0            0           0;
        0           0            0             0             0              0            0             1                  0          1           proj_45_y(k) 0             0            0             0             0             0            0           0           0           0           0           0           0           0           0           0           0           0           0           0            0           0;
        0           0            0             0             0              0            cog4_34_y(k) -cog4_34_x(k)       0          0           cog4_45(k)   0             0            0             0             0             0            0           0           0           0           0           0           0           0           0           0           0           0           0           0            0           1;
        0           0            0             0             0              0            0             0                  0          0          -proj_45_x(k) 1             0            0             0             0             0            0           0           0           0           0           0           0           0           0           0           0           0           0           0            0           0;
        0           0            0             0             0              0            0             0                  0          0          -proj_45_y(k) 0             1            0             0             0             0            0           0           0           0           0           0           0           0           0           0           0           0           0           0            0           0;
        0           0            0             0             0              0            0             0                  0          0           0            0             0            0             0             0             0            0           0           0           0           0           0           0           0           0           0           0           0           0           0            0           1;
        0           0            0             0             0              0            0             0                  0          0           0           -1             0           -1             0            -1             0            0           0           0           0           0           0           0           0           0           0           0           0           0           0            0           0;
        0           0            0             0             0              0            0             0                  0          0           0            0            -1            0            -1             0            -1            0           0           0           0           0           0           0           0           0           0           0           0           0           0            0           0;
        0           0            0             0             0              0            0             0                  0          0           0            cog6_56_y(k) -cog6_56_x(k) cog6_67_y(k) -cog6_67_x(k)  cog6_68_y(k) -cog6_68_x(k) 0           0           0           0           0           0           0           0           0           0           0           0           0           0           0            0           0;
        0           0            0             0             0              0            0             0                  0          0           0            0             0            1             0             0             0            1           0           0           0           0           0           0           0           0           0           0           0           0           0            0           0;
        0           0            0             0             0              0            0             0                  0          0           0            0             0            0             1             0             0            0           1           0           0           0           0           0           0           0           0           0           0           0           0            0           0;
        0           0            0             0             0              0            0             0                  0          0           0            0             0            0             0             0             0            0           0           0           0           0           0           0           0           0           0           0           0           0           0            0           0;
        
        
        0           0            0             0             0              0            0             0                  0          0           0            0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0           0            0           0];
    
  B = [ m2*acc_2x(k);
        m2*(acc_2y(k)+g);
        J2*ddphi2(k);
        m3*acc_3x(k);
        m3*(acc_3y(k)+g);
        J3*ddphi3(k);
        m4*acc_4x(k);
        m4*(acc_4y(k)+g);
        J4*ddphi4(k);
        m5*acc_5x(k);
        m5*(acc_5y(k)+g);
        J5*ddphi4(k);
        m6*acc_6x(k);
        m6*(acc_6y(k)+g);
        J6*ddphi6(k);
        m7*acc_7x(k);
        m7*(acc_7y(k)+g);
        J7*ddphi7(k);
        m8*acc_8x(k);
        m8*(acc_8y(k)+g);
        J8*ddphi8(k);
        m9*acc_9x(k);
        m9*(acc_9y(k)+g);
        0;
        m10*acc_10x(k);
        m10*(acc_10y(k)+g);
        J10*ddphi10(k);
        m11*acc_11x(k);
        m11*(acc_11y(k)+g);
        0;
        m12*acc_12x(k);
        m12*(acc_12y(k)+g);
        J12*ddphi12(k)];
        
    
    x = A\B;
    
    % save results
    F12x(k) = x(1);
    F12y(k) = x(2);
    F23x(k) = x(3);
    F23y(k) = x(4);
    F212x(k) = x(5);
    F67y(k) = x(6);
    F68x(k) = x(7);
    F212y(k) = x(8);
    F34x(k) = x(9);
    F34y(k) = x(10);
    F14x(k) = x(11);
    F14y(k) = x(12);
    F45(k) = x(13);
    F56x(k) = x(14);
    F56y(k) = x(15);
    F67x(k) = x(16);
    F68y(k) = x(17);
    F17x(k) = x(18);
    F17y(k) = x(19);
    F89x(k) = x(20);
    F89y(k) = x(21);
    F810x(k) = x(22);
    F810y(k) = x(23);
    F19(k) = x(24);    
    F1011x(k) = x(25);
    F1011y(k) = x(26);
    F1112x(k) = x(27);
    F1112y(k) = x(28);
    F111(k) = x(29);
    M12(k) = x(30);
    M19(k) = x(31);
    M111(k) = x(32);
    M45(k) = x(33);
end



% **********************
% *** plot figures ***
% **********************

% NOG AAN TE PASSEN !!!!


% if fig_dyn_4bar
%     
%     figure
%     subplot(221)
%     plot(F_P_x,F_P_y),grid
%     xlabel('F_P_x [N]')
%     ylabel('F_P_y [N]')
%     axis tight
%     subplot(222)
%     plot(F_Q_x,F_Q_y),grid
%     xlabel('F_Q_x [N]')
%     ylabel('F_Q_y [N]')
%     axis tight
%     subplot(223)
%     plot(F_R_x,F_R_y),grid
%     xlabel('F_R_x [N]')
%     ylabel('F_R_y [N]')
%     axis tight
%     subplot(224)
%     plot(F_S_x,F_S_y),grid
%     xlabel('F_S_x [N]')
%     ylabel('F_S_y [N]')
%     axis tight
%     
%     figure
%     plot(t,M_P)
%     ylabel('M_P [N-m]')
%     xlabel('t [s]')
    
end


