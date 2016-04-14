function [phi2,  phi3,  phi4,  x5,  phi6,  phi7,  phi8,  x9,  phi10,  x11,  phi12, ...
          dphi2, dphi3, dphi4, dx5, dphi6, dphi7, dphi8, dx9, dphi10, dx11, dphi12, ...
          ddphi2,ddphi3,ddphi4,ddx5,ddphi6,ddphi7,ddphi8,ddx9,ddphi10,ddx11,ddphi12, ...
          F12x, F12y, F23x, F23y, F212x, F212y, F34x, F34y, F14x, F14y, F45, F56x, F56y, F67x, F67y, ...
          F68x, F68y, F17x, F17y, F89x, F89y, F810x, F810y, F19, F1011x, F1011y, F1112x, F1112y, F111, ...
          M19, M111, M45]...
          =forward_dynamics(M12, phi2_init, dphi2_init, ... 
                            phi3_init, phi4_init, x5_init, phi6_init, phi7_init, phi8_init, x9_init, phi10_init, x11_init, phi12_init, ...
                            r2l, r2k, r3, r4l, r4k, r6l, r6k, r7, r8l, r8k, r10, r11, r12, x4, y4, x7, y7, y9, L9, ...
                            m2,m3,ma,mb,m4,m5,m6k,m6l,m6,m7,m8k,m8l,m8,m9,m10,m11,m12, mpiston1, mpiston2,...
                            X2,X3,X4,X5,X6k,X6l,X6,X7,X8k,X8l,X8,X9,X10,X11,X12, ...
                            Y2,Y3,Y4,Y5,Y6k,Y6l,Y6,Y7,Y8k,Y8l,Y8,Y9,Y10,Y11,Y12, ...
                            J2,J3,J4,J5,J6k,J6l,J6,J7,J8k,J8l,J8,J9,J10,J11,J12, t, fig_forward_dyn)

g=9.81;

%  initialization**

% allocation of the result vectors (this results in better performance because we don't have to reallocate and
% copy the vector each time we add an element.
phi2    = zeros(size(M12)); 
phi3    = zeros(size(M12));
phi4    = zeros(size(M12));
x5      = zeros(size(M12));
phi6    = zeros(size(M12));
phi7    = zeros(size(M12));
phi8    = zeros(size(M12));
x9      = zeros(size(M12));
phi10   = zeros(size(M12));
x11     = zeros(size(M12));
phi12   = zeros(size(M12));

dphi2    = zeros(size(M12));
dphi3    = zeros(size(M12));
dphi4    = zeros(size(M12));
dx5      = zeros(size(M12));
dphi6    = zeros(size(M12));
dphi7    = zeros(size(M12));
dphi8    = zeros(size(M12));
dx9      = zeros(size(M12));
dphi10   = zeros(size(M12));
dx11     = zeros(size(M12));
dphi12   = zeros(size(M12));

ddphi2    = zeros(size(M12));
ddphi3    = zeros(size(M12));
ddphi4    = zeros(size(M12));
ddx5      = zeros(size(M12));
ddphi6    = zeros(size(M12));
ddphi7    = zeros(size(M12));
ddphi8    = zeros(size(M12));
ddx9      = zeros(size(M12));
ddphi10   = zeros(size(M12));
ddx11     = zeros(size(M12));
ddphi12   = zeros(size(M12));

F12x   = zeros(size(M12));
F12y   = zeros(size(M12));
F23x   = zeros(size(M12));
F23y   = zeros(size(M12));
F212x  = zeros(size(M12));
F212y  = zeros(size(M12));
F34x   = zeros(size(M12));
F34y   = zeros(size(M12));
F14x   = zeros(size(M12));
F14y   = zeros(size(M12));
F45    = zeros(size(M12));
F56x   = zeros(size(M12));
F56y   = zeros(size(M12));
F67x   = zeros(size(M12));
F67y   = zeros(size(M12));
F68x   = zeros(size(M12));
F68y   = zeros(size(M12));
F17x   = zeros(size(M12));
F17y   = zeros(size(M12));
F89x   = zeros(size(M12));
F89y   = zeros(size(M12));
F810x  = zeros(size(M12));
F810y  = zeros(size(M12));
F19    = zeros(size(M12));
F1011x = zeros(size(M12));
F1011y = zeros(size(M12));
F1112x = zeros(size(M12));
F1112y = zeros(size(M12));
F111   = zeros(size(M12));
M19    = zeros(size(M12));
M111   = zeros(size(M12));
M45    = zeros(size(M12));


cog2_23_x = zeros(size(M12));
cog2_23_y = zeros(size(M12));
cog2_212_x = zeros(size(M12));
cog2_212_y = zeros(size(M12));

cog3_23_x = zeros(size(M12));
cog3_23_y = zeros(size(M12));
cog3_34_x = zeros(size(M12));
cog3_34_y = zeros(size(M12));

proj_45_x = zeros(size(M12));
proj_45_y = zeros(size(M12));

vp4_45 = zeros(size(M12));
vp4_34_x = zeros(size(M12));
vp4_34_y = zeros(size(M12));

cog6_56_x = zeros(size(M12));
cog6_56_y = zeros(size(M12));
cog6_67_x = zeros(size(M12));
cog6_67_y = zeros(size(M12));
cog6_68_x = zeros(size(M12));
cog6_68_y = zeros(size(M12));

cog7_17_x = zeros(size(M12));
cog7_17_y = zeros(size(M12));             
cog7_67_x = zeros(size(M12));
cog7_67_y = zeros(size(M12));

cog8_68_x = zeros(size(M12));
cog8_68_y = zeros(size(M12));
cog8_89_x = zeros(size(M12));             
cog8_89_y = zeros(size(M12));
cog8_810_x = zeros(size(M12));
cog8_810_y = zeros(size(M12));

cog10_810_x = zeros(size(M12));
cog10_810_y = zeros(size(M12));
cog10_1011_x = zeros(size(M12));
cog10_1011_y = zeros(size(M12));

cog12_212_x = zeros(size(M12));
cog12_212_y = zeros(size(M12));
cog12_1112_x = zeros(size(M12));
cog12_1112_y = zeros(size(M12));

  
omega2 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
omega3 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
omega4 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
omega6 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
omega7 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
omega8 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
omega10 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
omega12 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 

vec_23_cog3 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
vec_212_cog12 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
vec_vp7_cog7 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
vec_67_cog6 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
vec_vp7_67 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
vec_68_cog8 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
vec_67_68 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
vec_67_56 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
vec_68_810 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 

vec_vp2_cog2 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
vec_212_cog2 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
vec_23_cog2 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
vec_23_34 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
vec_vp4_cog4 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
vec_vp4_cog5 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
vec_68_89 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
vec_1011_cog10 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
vec_212_1112 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
vec_1011_810 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
vec_34_cog4 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))];

vec_23_cog3x = zeros(size(M12));
vec_23_cog3y = zeros(size(M12));
vec_212_cog12x = zeros(size(M12));
vec_212_cog12y = zeros(size(M12));
vec_vp7_cog7x = zeros(size(M12));
vec_vp7_cog7y = zeros(size(M12));
vec_67_cog6x = zeros(size(M12));
vec_67_cog6y = zeros(size(M12));
vec_vp7_67x = zeros(size(M12));
vec_vp7_67y = zeros(size(M12));
vec_68_cog8x = zeros(size(M12));
vec_68_cog8y = zeros(size(M12));
vec_67_68x = zeros(size(M12));
vec_67_68y = zeros(size(M12));
vec_67_56x = zeros(size(M12));
vec_67_56y = zeros(size(M12));
vec_vp2_cog2x = zeros(size(M12));
vec_vp2_cog2y = zeros(size(M12));
vec_cog2_212 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))];
vec_vp2_212 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))];
vec_vp2_212x = zeros(size(M12));
vec_vp2_212y = zeros(size(M12));
vec_cog2_23 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))];
vec_vp2_23 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))];
vec_vp2_23x = zeros(size(M12));
vec_vp2_23y = zeros(size(M12));
vec_vp4_cog4x = zeros(size(M12));
vec_vp4_cog4y = zeros(size(M12));
vec_1011_cog10x = zeros(size(M12));
vec_1011_cog10y = zeros(size(M12));

acc_2 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
acc_2_12 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
acc_2_3 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 

acc_3 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
acc_3_4 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 

acc_4 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
acc_4_check = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 

acc_7 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
acc_7_6 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 

acc_6 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
acc_6_8 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 

acc_5 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 

acc_8 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
acc_8_9 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
acc_8_10 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 


acc_12 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
acc_11_12 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 


acc_10 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 
acc_10_8 = [zeros(size(M12)) zeros(size(M12)) zeros(size(M12))]; 


acc_2x  = zeros(size(M12));
acc_2y = zeros(size(M12));
acc_3x = zeros(size(M12));
acc_3y = zeros(size(M12));
acc_4x = zeros(size(M12));
acc_4y = zeros(size(M12));
acc_5x = zeros(size(M12));
acc_5y = zeros(size(M12));
acc_6x = zeros(size(M12));
acc_6y = zeros(size(M12));
acc_7x = zeros(size(M12));
acc_7y = zeros(size(M12));
acc_8x = zeros(size(M12));
acc_8y = zeros(size(M12));
acc_9x = zeros(size(M12));
acc_9y = zeros(size(M12));
acc_10x = zeros(size(M12));
acc_10y = zeros(size(M12));
acc_11x = zeros(size(M12));
acc_11y = zeros(size(M12));
acc_12x = zeros(size(M12));
acc_12y = zeros(size(M12));

cog11_111_y = Y11*ones(size(M12));                      % voor stang 11 gerekend vanaf massacentrum, mag geen scalar zijn, maar vector met dezelfde waarde overal
cog11_1011_y = -(r11-Y11)*ones(size(M12));              % X11 vanaf scharnier 1,11, mag geen scalar zijn, maar vector met dezelfde waarde overal

optim_options = optimset('Display','off');
Ts = t(2) - t(1);      % timestep
t_size = size(t,1);    % number of simulation steps
for k=1:t_size
    
    [x, fval, exitflag] = fsolve('loop_closure_eqs', ...
                                 [phi3_init, phi4_init, x5_init, phi6_init, phi7_init, phi8_init, x9_init, phi10_init, x11_init, phi12_init]', ...
                                 optim_options, ...
                                 phi2(k), ...
                                 r2l, r2k, r3, r4l, r4k, r6l, r6k, r7, r8l, r8k, r10, r11, r12, x4, y4, x7, y7, y9);
    if (exitflag ~= 1)
        display 'The fsolve exit flag was not 1, probably no convergence!'
        exitflag
    end
    % Declaratie nieuwe variabele
    r6 = r6k + r6l;
    
    
    % save results of fsolve
    phi3(k) =   x(1);
    phi4(k) =   x(2);
    x5(k)  =    x(3);
    phi6(k) =   x(4);
    phi7(k) =   x(5);
    phi8(k) =   x(6);
    x9(k) =     x(7);
    phi10(k) =  x(8);
    x11(k) =    x(9);
    phi12(k) =  x(10);
    
    
    
    A = [0,                 0,                  0,                  0,                  0,                  0,                  0,                  0,                  1,                  -r12*sin(phi12(k));
         0,                 0,                  0,                  0,                  0,                  0,                  0,                  0,                  0,                  r12*cos(phi12(k));
         0,                 0,                  0,                  0,                  0,                  -r8l*sin(phi8(k)),  -1,                 -r10*sin(phi10(k)), -1,                 0;
         0,                 0,                  0,                  0,                  0,                  r8l*cos(phi8(k)),   0,                  r10*cos(phi10(k)),  0,                  0;
         0,                 -x5(k)*sin(phi4(k)),cos(phi4(k)),       -r6k*sin(phi6(k)),  r7*sin(phi7(k)),    0,                  0,                  0,                  0,                  0;
         0,                 x5(k)*cos(phi4(k)), sin(phi4(k)),       r6k*cos(phi6(k)),   -r7*cos(phi7(k)),   0,                  0,                  0,                  0,                  0;
         0,                 0,                  0,                  -r6*sin(phi6(k)),   r7*sin(phi7(k)),    r8k*sin(phi8(k)),  -1,                 0,                  0,                  0;
         0,                 0,                  0,                  r6*cos(phi6(k)),    -r7*cos(phi7(k)),   -r8k*cos(phi8(k)),  0,                  0,                  0,                  0;
         -r3*sin(phi3(k)),  -r4l*sin(phi4(k))+r4k*cos(phi4(k)),  0,     0,                  0,                  0,                  0,                  0,                  0,                  0;
         r3*cos(phi3(k)),   r4l*cos(phi4(k))+r4k*sin(phi4(k)), 0,       0,                  0,                  0,                  0,                  0,                  0,                  0];
     
    B = [-r2l*sin(phi2(k))*dphi2(k);
         r2l*cos(phi2(k))*dphi2(k);
         0;
         0;
         0;
         0;
         0;
         0;
         -r2k*cos(phi2(k))*dphi2(k);
         -r2k*sin(phi2(k))*dphi2(k)];
        
    x = A\B;
    
    % save results
    dphi3(k)  = x(1);
    dphi4(k)  = x(2);
    dx5(k)    = x(3);
    dphi6(k)  = x(4);
    dphi7(k)  = x(5);
    dphi8(k)  = x(6);
    dx9(k)    = x(7);
    dphi10(k) = x(8);
    dx11(k)   = x(9);
    dphi12(k) = x(10);
    
cog2_23_x(k) = r2k*cos(phi2(k)-pi/2);        % voor stang 2 gerekend vanaf vaste punt
cog2_23_y(k) = r2k*sin(phi2(k)-pi/2);
cog2_212_x(k) = -r2l*cos(phi2(k));
cog2_212_y(k) = -r2l*sin(phi2(k));

cog3_23_x(k) = -X3*cos(phi3(k));               % voor stang 3 gerekend vanaf massacentrum
cog3_23_y(k) = -X3*sin(phi3(k));               % X3 vanaf scharnier 1,3
cog3_34_x(k) = (r3-X3)*cos(phi3(k));
cog3_34_y(k) = (r3-X3)*sin(phi3(k));

proj_45_x(k) = cos(phi4(k)-pi/2);             % projectie van kracht F45
proj_45_y(k) = sin(phi4(k)-pi/2);

vp4_45(k) = x5(k);                           % voor stang 4 gerekend vanaf vaste punt
vp4_34_x(k) = (-r4l*cos(phi4(k))) + (r4k*cos(phi4(k)+pi/2));      % Fout ontdekt: a is de lange zijde, b is de korte zijde
vp4_34_y(k) = (-r4l*sin(phi4(k))) + (r4k*sin(phi4(k)+pi/2));     % Inderdaad fouten in de mintekens, zou zou het volgens mij moeten kloppen

cog6_56_x(k) = -(X6-r6k)*cos(phi6(k));        % voor stang 6 gerekend vanaf massacentrum
cog6_56_y(k) = -(X6-r6k)*sin(phi6(k));        % X6 vanaf scharnier 6,7
cog6_67_x(k) = -X6*cos(phi6(k));
cog6_67_y(k) = -X6*sin(phi6(k));
cog6_68_x(k) = (r6k + r6l - X6)*cos(phi6(k));
cog6_68_y(k) = (r6k + r6l - X6)*sin(phi6(k));

cog7_17_x(k) = X7*cos(phi7(k));               % voor stang 7 gerekend vanaf massacentrum
cog7_17_y(k) = X7*sin(phi7(k));               % X7 vanaf scharnier 1,7
cog7_67_x(k) = -(r7-X7)*cos(phi7(k));
cog7_67_y(k) = -(r7-X7)*sin(phi7(k));

cog8_68_x(k) = (r8k + r8l -X8)*cos(phi8(k));          % voor stang 8 gerekend vanaf massacentrum
cog8_68_y(k) = (r8k + r8l -X8)*sin(phi8(k));          % X8 vanaf scharnier 8,10
cog8_89_x(k) = L9*cos(phi8(k));               % L9 de afstand van cog tot aangrijping F89
cog8_89_y(k) = L9*sin(phi8(k));
cog8_810_x(k) = -X8*cos(phi8(k));
cog8_810_y(k) = -X8*sin(phi8(k));


cog10_810_x(k) = (r10-X10)*cos(phi10(k));     % voor stang 10 gerekend vanaf massacenrum
cog10_810_y(k) = (r10-X10)*sin(phi10(k));    % X10 vanaf scharnier 10,11
cog10_1011_x(k) = -X10*cos(phi10(k));
cog10_1011_y(k) = -X10*sin(phi10(k));



cog12_212_x(k) = -X12*cos(phi12(k));          % voor stang 12 gerekend vanaf massacentrum
cog12_212_y(k) = -X12*sin(phi12(k));          % X12 vanaf scharnier 2,12
cog12_1112_x(k) = (r12-X12)*cos(phi12(k));
cog12_1112_y(k) = (r12-X12)*sin(phi12(k));



% Declaratie 3D omega en alpha vectoren ***

% 3D omega (dphi) and alpha (ddphi) vectors)    
omega2(k,:) = [0 0 dphi2(k)]; 
omega3(k,:) = [0 0 dphi3(k)]; 
omega4(k,:) = [0 0 dphi4(k)]; 
omega6(k,:) = [0 0 dphi6(k)]; 
omega7(k,:) = [0 0 dphi7(k)]; 
omega8(k,:) = [0 0 dphi8(k)]; 
omega10(k,:) = [0 0 dphi10(k)]; 
omega12(k,:) = [0 0 dphi12(k)]; 




%  3D vectoren voor versnellingen ***

% 3D model vectors:

vec_23_cog3(k,:) = [-cog3_23_x(k)    -cog3_23_y(k)   0];
vec_23_cog3x(k) = vec_23_cog3(k,1);
vec_23_cog3y(k) = vec_23_cog3(k,2);
vec_212_cog12(k,:) = [-cog12_212_x(k)    -cog12_212_y(k)   0];
vec_212_cog12x(k) = vec_212_cog12(k,1);
vec_212_cog12y(k) = vec_212_cog12(k,2);
vec_vp7_cog7(k,:) = [-cog7_17_x(k)    -cog7_17_y(k)   0];
vec_vp7_cog7x(k) = vec_vp7_cog7(k,1);
vec_vp7_cog7y(k) = vec_vp7_cog7(k,2);
vec_67_cog6(k,:) = [-cog6_67_x(k)    -cog6_67_y(k)   0];
vec_67_cog6x(k) = vec_67_cog6(k,1);
vec_67_cog6y(k) = vec_67_cog6(k,2);
vec_vp7_67(k,:) = [-r7*cos(phi7(k))   -r7*sin(phi7(k)) 0];
vec_vp7_67x(k) = vec_vp7_67(k,1);
vec_vp7_67y(k) = vec_vp7_67(k,2);
vec_68_cog8(k,:) = [-cog8_68_x(k)    -cog8_68_y(k)   0];
vec_68_cog8x(k) = vec_68_cog8(k,1);
vec_68_cog8y(k) = vec_68_cog8(k,2);
vec_67_68(k,:) = [(r6k+r6l)*cos(phi6(k))   (r6k+r6l)*sin(phi6(k))  0];
vec_67_68x(k) = vec_67_68(k,1);
vec_67_68y(k) = vec_67_68(k,2);
vec_67_56(k,:) = [r6k*cos(phi6(k))   r6k*sin(phi6(k))   0];
vec_67_56x(k) = vec_67_56(k,1);
vec_67_56y(k) = vec_67_56(k,2);
vec_68_810(k,:) = [-(r8k+r8l)*cos(phi8(k))   -(r8k+r8l)*sin(phi8(k))  0];

%Extra vectoren nodig:
vec_vp2_cog2(k,:) = [(X2*cos(phi2(k)-pi/2))-(Y2*cos(phi2(k))) (X2*sin(phi2(k)-pi/2))-(Y2*sin(phi2(k))) 0];
vec_vp2_cog2x(k) = vec_vp2_cog2(k,1);
vec_vp2_cog2y(k) = vec_vp2_cog2(k,2);
vec_212_cog2(k,:) = [-cog2_212_x(k)    -cog2_212_y(k)   0];
vec_cog2_212 = -vec_212_cog2;
vec_vp2_212 = vec_vp2_cog2 + vec_cog2_212;
vec_vp2_212x(k) = vec_vp2_212(k,1);
vec_vp2_212y(k) = vec_vp2_212(k,2);
vec_23_cog2(k,:) = [-cog2_23_x(k)    -cog2_23_y(k)   0];
vec_cog2_23 = -vec_23_cog2;
vec_vp2_23 = vec_vp2_cog2 + vec_cog2_23;
vec_vp2_23x(k) = vec_vp2_23(k,1);
vec_vp2_23y(k) = vec_vp2_23(k,2);
vec_23_34(k,:) = [r3*cos(phi3(k)) r3*sin(phi3(k)) 0];
vec_vp4_cog4(k,:) = [(-Y4*cos(phi4(k)))+(X4*cos(phi4(k)+pi/2))  (-Y4*sin(phi4(k)))+(X4*sin(phi4(k)+pi/2)) 0];       % moet hetzelfde zijn als vp4_34 maar Y4 ipv r2k en X4 ipv r2l                         
vec_vp4_cog4x(k) = vec_vp4_cog4(k,1);
vec_vp4_cog4y(k) = vec_vp4_cog4(k,2);
vec_vp4_cog5(k,:) = [-1*x5(k)*cos(phi4(k))  -1*x5(k)*sin(phi4(k))  0];       % Moet kleine x5 zijn, want deze positie varieert continu + elk element apart vermenigvuldigen

vec_68_89(k,:) = [-r8k*cos(phi8(k))  -r8k*sin(phi8(k))  0];
vec_1011_cog10(k,:) = [X10*cos(phi10(k))  X10*sin(phi10(k))  0];
vec_1011_cog10x(k) = vec_1011_cog10(k,1);
vec_1011_cog10y(k) = vec_1011_cog10(k,2);
vec_212_1112(k,:) = [r12*cos(phi12(k))  r12*sin(phi12(k))   0];

vec_1011_810(k,:) = [r10*cos(phi10(k))  r10*sin(phi10(k))  0];

vec_34_cog4(k,:) = [(-(r4k-X4)*cos(phi4(k)+pi/2))+((r4l-Y4)*cos(phi4(k)))  (-(r4k-X4)*sin(phi4(k)+pi/2))+((r4l-Y4)*sin(phi4(k))) 0];       % moet hetzelfde zijn als vp4_34 maar (a-Y4) ipv a en (b-X4) ipv b + tegengesteld gericht




%%%%% 3D versnellingsvectoren voor matrix B ***

% acceleration vectors
% bv. acc_2 = versnelling massacentrum van stang 2
%     acc_2_12 = versnelling van scharnierpunt 2,12

acc_2(k,:) =                   cross(omega2(k,:),cross(omega2(k,:),vec_vp2_cog2(k,:)));  % normaal gezien nul
acc_2_12(k,:) = acc_2(k,:) +     cross(omega2(k,:),cross(omega2(k,:),-vec_212_cog2(k,:)));
acc_2_3(k,:) = acc_2(k,:) +      cross(omega2(k,:),cross(omega2(k,:),-vec_23_cog2(k,:)));

acc_3(k,:) = acc_2_3(k,:) +      cross(omega3(k,:),cross(omega3(k,:),vec_23_cog3(k,:)));
acc_3_4(k,:) = acc_2_3(k,:) +    cross(omega3(k,:),cross(omega3(k,:),vec_23_34(k,:)));   % tot besef gekomen dat we deze niet nodig hebben, tenzij voor controle

acc_4(k,:) =                   cross(omega4(k,:),cross(omega4(k,:),vec_vp4_cog4(k,:)));
acc_4_check(k,:) = acc_3_4(k,:) +cross(omega4(k,:),cross(omega4(k,:),vec_34_cog4(k,:)));
%Deze controle klopt!

acc_7(k,:) =                cross(omega7(k,:),cross(omega7(k,:),vec_vp7_cog7(k,:))) ;
acc_7_6(k,:) =              cross(omega7(k,:),cross(omega7(k,:),vec_vp7_67(k,:)));

acc_6(k,:) = acc_7_6(k,:) +      cross(omega6(k,:),cross(omega6(k,:),vec_67_cog6(k,:)));
acc_6_8(k,:) = acc_7_6(k,:) +    cross(omega6(k,:),cross(omega6(k,:),vec_67_68(k,:)));

%acc_5_lin = [-1*times(ddx5,cos(phi4))  -1*times(ddx5,sin(phi4))  zeros(size(phi2))];            % times(A,B) geeft de vector terug waar elk element van A vermenigvuldigd wordt met het overeenkomstige element van B
%acc_5 = acc_5_lin +    cross(omega4,cross(omega4,vec_vp4_cog5)) +    cross(alpha4,vec_vp4_cog5);
% => Dit klopt niet: er moet nog een coriolisvernelling in. Gemakkelijker
% is de versnelling via stang 6 te definiëren
acc_5(k,:) = acc_7_6(k,:) +      cross(omega6(k,:),cross(omega6(k,:),vec_67_56(k,:))) ;

acc_8(k,:) = acc_6_8(k,:) +      cross(omega8(k,:),cross(omega8(k,:),vec_68_cog8(k,:))) ;
acc_8_9(k,:) = acc_6_8(k,:) +    cross(omega8(k,:),cross(omega8(k,:),vec_68_89(k,:))) ;
acc_8_10(k,:) = acc_6_8(k,:)+    cross(omega8(k,:),cross(omega8(k,:),vec_68_810(k,:))) ;  % tot besef gekomen dat we deze niet nodig hebben, tenzij voor controle


acc_12(k,:) = acc_2_12(k,:) +    cross(omega12(k,:),cross(omega12(k,:),vec_212_cog12(k,:)));
acc_11_12(k,:) = acc_2_12(k,:) + cross(omega12(k,:),cross(omega12(k,:),vec_212_1112(k,:)));



acc_10(k,:) =       cross(omega10(k,:),cross(omega10(k,:),vec_1011_cog10(k,:)));
acc_10_8(k,:) =     cross(omega10(k,:),cross(omega10(k,:),vec_1011_810(k,:)));

% Declaratie van de deelversnellingen voor matrix B:
acc_2x(k) = acc_2(k,1);
acc_2y(k) = acc_2(k,2);
acc_3x(k) = acc_3(k,1);
acc_3y(k) = acc_3(k,2);
acc_4x(k) = acc_4(k,1);
acc_4y(k) = acc_4(k,2);
acc_5x(k) = acc_5(k,1);
acc_5y(k) = acc_5(k,2);
acc_6x(k) = acc_6(k,1);
acc_6y(k) = acc_6(k,2);
acc_7x(k) = acc_7(k,1);
acc_7y(k) = acc_7(k,2);
acc_8x(k) = acc_8(k,1);
acc_8y(k) = acc_8(k,2);
acc_9x(k) = acc_8_9(k,1);
acc_9y(k) = acc_8_9(k,2);
acc_10x(k) = acc_10(k,1);
acc_10y(k) = acc_10(k,2);
acc_11x(k) = acc_11_12(k,1);
acc_11y(k) = acc_11_12(k,2);
acc_12x(k) = acc_12(k,1);
acc_12y(k) = acc_12(k,2);

%         F12x        F12y         F23x          F23y          F212x           F212y           F34x          F34y               F14x       F14y         F45          F56x         F56y         F67x          F67y          F68x          F68y          F17x          F17y          F89x          F89y          F810x           F810y          F19         F1011x           F1011y          F1112x           F1112y          F111                   M19          M111        M45   ddphi2               ddphi3              ddphi4                             ddx5              ddphi6              ddphi7               ddphi8              ddx9  ddphi10                 ddphi11 ddphi12
    A = [ 1           0            1             0             1               0               0             0                  0          0           0            0             0            0             0             0             0             0             0             0             0             0               0              0           0                0               0                0               0                      0            0           0     m2*vec_vp2_cog2y(k)  0                   0                                  0                 0                   0                    0                   0     0                       0       0;
          0           1            0             1             0               1               0             0                  0          0           0            0             0            0             0             0             0             0             0             0             0             0               0              0           0                0               0                0               0                      0            0           0     -m2*vec_vp2_cog2x(k) 0                   0                                  0                 0                   0                    0                   0     0                       0       0;
          0           0            -cog2_23_y(k) cog2_23_x(k) -cog2_212_y(k)   cog2_212_x(k)   0             0                  0          0           0            0             0            0             0             0             0             0             0             0             0             0               0              0           0                0               0                0               0                      0            0           0     -J2                  0                   0                                  0                 0                   0                    0                   0     0                       0       0;
          0           0           -1             0             0               0              -1             0                  0          0           0            0             0            0             0             0             0             0             0             0             0             0               0              0           0                0               0                0               0                      0            0           0     m3*vec_vp2_23y(k)    m3*vec_23_cog3y(k)  0                                  0                 0                   0                    0                   0     0                       0       0;
          0           0            0            -1             0               0               0            -1                  0          0           0            0             0            0             0             0             0             0             0             0             0             0               0              0           0                0               0                0               0                      0            0           0     -m3*vec_vp2_23x(k)   -m3*vec_23_cog3x(k) 0                                  0                 0                   0                    0                   0     0                       0       0;
          0           0            cog3_23_y(k) -cog3_23_x(k)  0               0               cog3_34_y(k) -cog3_34_x(k)       0          0           0            0             0            0             0             0             0             0             0             0             0             0               0              0           0                0               0                0               0                      0            0           0     0                    -J3                 0                                  0                 0                   0                    0                   0     0                       0       0;
          0           0            0             0             0               0               1             0                  1          0           proj_45_x(k) 0             0            0             0             0             0             0             0             0             0             0               0              0           0                0               0                0               0                      0            0           0     0                    0                   m4*vec_vp4_cog4y(k)                0                 0                   0                    0                   0     0                       0       0;
          0           0            0             0             0               0               0             1                  0          1           proj_45_y(k) 0             0            0             0             0             0             0             0             0             0             0               0              0           0                0               0                0               0                      0            0           0     0                    0                   -m4*vec_vp4_cog4x(k)               0                 0                   0                    0                   0     0                       0       0;
          0           0            0             0             0               0              -vp4_34_y(k)   vp4_34_x(k)        0          0           vp4_45(k)    0             0            0             0             0             0             0             0             0             0             0               0              0           0                0               0                0               0                      0            0           1     0                    0                   -(J4+(X4^2+Y4^2)*m4)               0                 0                   0                    0                   0     0                       0       0;
          0           0            0             0             0               0               0             0                  0          0          -proj_45_x(k) 1             0            0             0             0             0             0             0             0             0             0               0              0           0                0               0                0               0                      0            0           0     0                    0                   0                                  0                 m5*vec_67_56y(k)    m5*vec_vp7_67y(k)    0                   0     0                       0       0;
          0           0            0             0             0               0               0             0                  0          0          -proj_45_y(k) 0             1            0             0             0             0             0             0             0             0             0               0              0           0                0               0                0               0                      0            0           0     0                    0                   0                                  0                 -m5*vec_67_56x(k)   -m5*vec_vp7_67x(k)   0                   0     0                       0       0;
          0           0            0             0             0               0               0             0                  0          0           0            0             0            0             0             0             0             0             0             0             0             0               0              0           0                0               0                0               0                      0            0          -1     0                    0                   -J5                                0                 0                   0                    0                   0     0                       0       0;
          0           0            0             0             0               0               0             0                  0          0           0           -1             0           -1             0            -1             0             0             0             0             0             0               0              0           0                0               0                0               0                      0            0           0     0                    0                   0                                  0                 m6*vec_67_cog6y(k)  m6*vec_vp7_67y(k)    0                   0     0                       0       0;
          0           0            0             0             0               0               0             0                  0          0           0            0            -1            0            -1             0            -1             0             0             0             0             0               0              0           0                0               0                0               0                      0            0           0     0                    0                   0                                  0                 -m6*vec_67_cog6x(k) -m6*vec_vp7_67x(k)   0                   0     0                       0       0;
          0           0            0             0             0               0               0             0                  0          0           0            cog6_56_y(k) -cog6_56_x(k) cog6_67_y(k) -cog6_67_x(k)  cog6_68_y(k) -cog6_68_x(k)  0             0             0             0             0               0              0           0                0               0                0               0                      0            0           0     0                    0                   0                                  0                 -J6                 0                    0                   0     0                       0       0;
          0           0            0             0             0               0               0             0                  0          0           0            0             0            1             0             0             0             1             0             0             0             0               0              0           0                0               0                0               0                      0            0           0     0                    0                   0                                  0                 0                   m7*vec_vp7_cog7y(k)  0                   0     0                       0       0;
          0           0            0             0             0               0               0             0                  0          0           0            0             0            0             1             0             0             0             1             0             0             0               0              0           0                0               0                0               0                      0            0           0     0                    0                   0                                  0                 0                   -m7*vec_vp7_cog7x(k) 0                   0     0                       0       0;
          0           0            0             0             0               0               0             0                  0          0           0            0             0           -cog7_67_y(k)  cog7_67_x(k)  0             0            -cog7_17_y(k)  cog7_17_x(k)  0             0             0               0              0           0                0               0                0               0                      0            0           0     0                    0                   0                                  0                 0                   -J7                  0                   0     0                       0       0;
          0           0            0             0             0               0               0             0                  0          0           0            0             0            0             0             1             0             0             0             1             0             1               0              0           0                0               0                0               0                      0            0           0     0                    0                   0                                  0                 m8*vec_67_68y(k)    m8*vec_vp7_67y(k)    m8*vec_68_cog8y(k)  0     0                       0       0;
          0           0            0             0             0               0               0             0                  0          0           0            0             0            0             0             0             1             0             0             0             1             0               1              0           0                0               0                0               0                      0            0           0     0                    0                   0                                  0                 -m8*vec_67_68x(k)   -m8*vec_vp7_67x(k)   -m8*vec_68_cog8x(k) 0     0                       0       0;
          0           0            0             0             0               0               0             0                  0          0           0            0             0            0             0            -cog8_68_y(k)  cog8_68_x(k)  0             0            -cog8_89_y(k)  cog8_89_x(k) -cog8_810_y(k)   cog8_810_x(k)  0           0                0               0                0               0                      0            0           0     0                    0                   0                                  0                 0                   0                    -J8                 0     0                       0       0;
          0           0            0             0             0               0               0             0                  0          0           0            0             0            0             0             0             0             0             0            -1             0             0               0              0           0                0               0                0               0                      0            0           0     0                    0                   0                                  0                 0                   0                    0                   -m9   0                       0       0;
          0           0            0             0             0               0               0             0                  0          0           0            0             0            0             0             0             0             0             0             0            -1             0               0             -1           0                0               0                0               0                      0            0           0     0                    0                   0                                  0                 0                   0                    0                   0     0                       0       0;
          0           0            0             0             0               0               0             0                  0          0           0            0             0            0             0             0             0             0             0             0             0             0               0              0           0                0               0                0               0                     -1            0           0     0                    0                   0                                  0                 0                   0                    0                   0     0                       0       0;
          0           0            0             0             0               0               0             0                  0          0           0            0             0            0             0             0             0             0             0             0             0            -1               0              0          -1                0               0                0               0                      0            0           0     0                    0                   0                                  0                 0                   0                    0                   0     m10*vec_1011_cog10y(k)  m10     0;    
          0           0            0             0             0               0               0             0                  0          0           0            0             0            0             0             0             0             0             0             0             0             0              -1              0           0               -1               0                0               0                      0            0           0     0                    0                   0                                  0                 0                   0                    0                   0     -m10*vec_1011_cog10x(k) 0       0;                                                                                                            
          0           0            0             0             0               0               0             0                  0          0           0            0             0            0             0             0             0             0             0             0             0             cog10_810_y(k) -cog10_810_x(k) 0           cog10_1011_y(k) -cog10_1011_x(k) 0                0               0                      0            0           0     0                    0                   0                                  0                 0                   0                    0                   0     -J10                    0       0;
          0           0            0             0             0               0               0             0                  0          0           0            0             0            0             0             0             0             0             0             0             0             0               0              0           1                0               1                0               0                      0            0           0     0                    0                   0                                  0                 0                   0                    0                   0     0                       m11     0;     
          0           0            0             0             0               0               0             0                  0          0           0            0             0            0             0             0             0             0             0             0             0             0               0              0           0                1               0                1               1                      0            0           0     0                    0                   0                                  0                 0                   0                    0                   0     0                       0       0;
          0           0            0             0             0               0               0             0                  0          0           0            0             0            0             0             0             0             0             0             0             0             0               0              0          -cog11_1011_y(k)  0              -cog11_111_y(k)   0               0                      0            1           0     0                    0                   0                                  0                 0                   0                    0                   0     0                       0       0;
          0           0            0             0             -1              0               0             0                  0          0           0            0             0            0             0             0             0             0             0             0             0             0               0              0           0                0              -1                0               0                      0            0           0     m12*vec_vp2_212y(k)  0                   0                                  0                 0                   0                    0                   0     0                       0       m12*vec_212_cog12y(k);               
          0           0            0             0             0               -1              0             0                  0          0           0            0             0            0             0             0             0             0             0             0             0             0               0              0           0                0               0               -1               0                      0            0           0     -m12*vec_vp2_212x(k) 0                   0                                  0                 0                   0                    0                   0     0                       0       -m12*vec_212_cog12x(k);               
          0           0            0             0             cog12_212_y(k) -cog12_212_x(k)  0             0                  0          0           0            0             0            0             0             0             0             0             0             0             0             0               0              0           0                0               cog12_1112_y(k) -cog12_1112_x(k) 0                      0            0           0     0                    0                   0                                  0                 0                   0                    0                   0     0                       0       -J12;
          0           0            0             0             0               0               0             0                  0          0           0            0             0            0             0             0             0             0             0             0             0             0               0              0           0                0               0                0               0                      0            0           0     r2l*sin(phi2(k))     0                   0                                  0                 0                   0                    0                   0     0                       1       -r12*sin(phi12(k));
          0           0            0             0             0               0               0             0                  0          0           0            0             0            0             0             0             0             0             0             0             0             0               0              0           0                0               0                0               0                      0            0           0     -r2l*cos(phi2(k))    0                   0                                  0                 0                   0                    0                   0     0                       0       r12*cos(phi12(k));
          0           0            0             0             0               0               0             0                  0          0           0            0             0            0             0             0             0             0             0             0             0             0               0              0           0                0               0                0               0                      0            0           0      0                   0                   0                                  0                 0                   0                    -r8l*sin(phi8(k))   -1    -r10*sin(phi10(k))      -1      0;
          0           0            0             0             0               0               0             0                  0          0           0            0             0            0             0             0             0             0             0             0             0             0               0              0           0                0               0                0               0                      0            0           0      0                   0                   0                                  0                 0                   0                    r8l*cos(phi8(k))    0     r10*cos(phi10(k))       0       0;
          0           0            0             0             0               0               0             0                  0          0           0            0             0            0             0             0             0             0             0             0             0             0               0              0           0                0               0                0               0                      0            0           0      0                   0                   -x5(k)*sin(phi4(k))                cos(phi4(k))      -r6k*sin(phi6(k))   r7*sin(phi7(k))      0                   0     0                       0       0;
          0           0            0             0             0               0               0             0                  0          0           0            0             0            0             0             0             0             0             0             0             0             0               0              0           0                0               0                0               0                      0            0           0      0                   0                   x5(k)*cos(phi4(k))                 sin(phi4(k))      r6k*cos(phi6(k))    -r7*cos(phi7(k))     0                   0     0                       0       0;
          0           0            0             0             0               0               0             0                  0          0           0            0             0            0             0             0             0             0             0             0             0             0               0              0           0                0               0                0               0                      0            0           0      0                   0                   0                                  0                 -r6*sin(phi6(k))    r7*sin(phi7(k))      r8k*sin(phi8(k))    -1    0                       0       0;
          0           0            0             0             0               0               0             0                  0          0           0            0             0            0             0             0             0             0             0             0             0             0               0              0           0                0               0                0               0                      0            0           0      0                   0                   0                                  0                 r6*cos(phi6(k))     -r7*cos(phi7(k))     -r8k*cos(phi8(k))   0     0                       0       0;
          0           0            0             0             0               0               0             0                  0          0           0            0             0            0             0             0             0             0             0             0             0             0               0              0           0                0               0                0               0                      0            0           0      r2k*cos(phi2(k))    -r3*sin(phi3(k))    -r4l*sin(phi4(k))+r4k*cos(phi4(k)) 0                 0                   0                    0                   0     0                       0       0;
          0           0            0             0             0               0               0             0                  0          0           0            0             0            0             0             0             0             0             0             0             0             0               0              0           0                0               0                0               0                      0            0           0      r2k*sin(phi2(k))    r3*cos(phi3(k))     r4l*cos(phi4(k))+r4k*sin(phi4(k))  0                 0                   0                    0                   0     0                       0       0];
%         F12x        F12y         F23x          F23y          F212x           F212y           F34x          F34y               F14x       F14y        F45          F56x          F56y         F67x          F67y          F68x          F68y          F17x          F17y          F89x          F89y          F810x           F810y          F19         F1011x           F1011y          F1112x           F1112y          F111                   M19          M111        M45

   B = [m2*acc_2x(k);
        m2*(acc_2y(k)+g);
        -M12(k);                    % Was positief in A dynamica, moet nu van teken veranderen
        m3*acc_3x(k);
        m3*(acc_3y(k)+g);
        0;
        m4*acc_4x(k);
        m4*(acc_4y(k)+g);
        m4*g*vec_vp4_cog4(k,1);     % Aangepast
        m5*acc_5x(k);
        m5*(acc_5y(k)+g);
        0;
        m6*acc_6x(k);
        m6*(acc_6y(k)+g);
        0;
        m7*acc_7x(k);
        m7*(acc_7y(k)+g);
        0;
        m8*acc_8x(k);
        m8*(acc_8y(k)+g);
        0;
        m9*acc_9x(k);
        m9*(acc_9y(k)+g);
        0;
        m10*acc_10x(k);
        m10*(acc_10y(k)+g);
        0;
        m11*acc_11x(k);
        m11*(acc_11y(k)+g);
        0;
        m12*acc_12x(k);
        m12*(acc_12y(k)+g);
        0;
        -r2l*cos(phi2(k))*dphi2(k)^2+r12*cos(phi12(k))*dphi12(k)^2;
        -r2l*sin(phi2(k))*dphi2(k)^2+r12*sin(phi12(k))*dphi12(k)^2;
        r10*cos(phi10(k))*dphi10(k)^2+r8l*cos(phi8(k))*dphi8(k)^2;
        r10*sin(phi10(k))*dphi10(k)^2+r8l*sin(phi8(k))*dphi8(k)^2;
        -r7*cos(phi7(k))*dphi7(k)^2+r6k*cos(phi6(k))*dphi6(k)^2+2*sin(phi4(k))*dphi4(k)*dx5(k)+x5(k)*cos(phi4(k))*dphi4(k)^2;
        -r7*sin(phi7(k))*dphi7(k)^2+r6k*sin(phi6(k))*dphi6(k)^2-2*cos(phi4(k))*dphi4(k)*dx5(k)+x5(k)*sin(phi4(k))*dphi4(k)^2;
        -r7*cos(phi7(k))*dphi7(k)^2+r6*cos(phi6(k))*dphi6(k)^2-r8k*cos(phi8(k))*dphi8(k)^2;
        -r7*sin(phi7(k))*dphi7(k)^2+r6*sin(phi6(k))*dphi6(k)^2-r8k*sin(phi8(k))*dphi8(k)^2;
        r2k*sin(phi2(k))*dphi2(k)^2+r3*cos(phi3(k))*dphi3(k)^2+(r4l*cos(phi4(k))+r4k*sin(phi4(k)))*dphi4(k)^2;
        -r2k*cos(phi2(k))*dphi2(k)^2+r3*sin(phi3(k))*dphi3(k)^2+(r4l*sin(phi4(k))-r4k*cos(phi4(k)))*dphi4(k)^2];
     
   
    
     x=A\B;
     
 % save results
    F12x(k) = x(1);
    F12y(k) = x(2);
    F23x(k) = x(3);
    F23y(k) = x(4);
    F212x(k) = x(5);
    F212y(k) = x(6);
    F34x(k) = x(7);
    F34y(k) = x(8);
    F14x(k) = x(9);
    F14y(k) = x(10);
    F45(k) = x(11);
    F56x(k) = x(12);
    F56y(k) = x(13);
    F67x(k) = x(14);
    F67y(k) = x(15);
    F68x(k) = x(16);
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
    M19(k) = x(30);
    M111(k) = x(31);
    M45(k) = x(32);
    ddphi2(k) = x(33);
    ddphi3(k) = x(34);
    ddphi4(k) = x(35);
    ddx5(k) = x(36);
    ddphi6(k) = x(37);
    ddphi7(k) = x(38);
    ddphi8(k) = x(39);
    ddx9(k) = x(40);
    ddphi10(k) = x(41);
    ddx11(k) = x(42);
    ddphi12(k) = x(43);
    
    phi2(k+1)    = phi2(k) + Ts * dphi2(k);
    dphi2(k+1)    = dphi2(k) + Ts * ddphi2(k);

end

phi2 = phi2(1:k,:);
dphi2 = dphi2(1:k,:);

if fig_forward_dyn
       
    screen_size = get(groot, 'ScreenSize');
    figure('Name', 'Controle voorwaartse dynamica', 'NumberTitle', 'off', ...
           'Position', [screen_size(3)/3 screen_size(4)/6 screen_size(3)/3 screen_size(4)/1.5])
    subplot(3,1,1)
    plot(t, phi2)
    ylabel('\theta_2 [rad]')   
    
    subplot(3,1,2)
    plot(t, dphi2)
    ylabel('d\theta_2 [rad/s]')
    
    subplot(3,1,3)
    plot(t, ddphi2)
    ylabel('dd\theta_2 [rad/s²]')
    
    set(gcf,'NextPlot','add');
    axes;
    h = title({'Controle voorwaartse dynamica'; ''});
    set(gca,'Visible','off');
    set(h,'Visible','on')
end
