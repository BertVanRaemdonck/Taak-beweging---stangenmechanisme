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

function [  phi3,   phi4,   x5,     phi6,   phi7,   phi8,   x9,     phi10,      x11,    phi12, ... 
            dphi3,  dphi4,  dx5,    dphi6,  dphi7,  dphi8,  dx9,    dphi10,     dx11,   dphi12, ...
            ddphi3, ddphi4, ddx5,   ddphi6, ddphi7, ddphi8, ddx9,   ddphi10,    ddx11,  ddphi12   ] = ...
            kinematics_4bar(r2l, r2k, r3, a, b, r6l, r6k, r7, r8l, r8k, r10, r11, r12, x4, y4, x7, y7, y9, ...
                    phi1, phi2, dphi2, ddphi2, omega, alpha, ...
                    phi3_init, phi4_init, x5_init, phi6_init, phi7_init, phi8_init, x9_init, phi10_init, x11_init, phi12_init, ...
                    t, fig_kin_4bar)

% allocation of the result vectors (this results in better performance because we don't have to reallocate and
% copy the vector each time we add an element.
phi3    = zeros(size(t));
phi4    = zeros(size(t));
x5      = zeros(size(t));
phi6    = zeros(size(t));
phi7    = zeros(size(t));
phi8    = zeros(size(t));
x9      = zeros(size(t));
phi10   = zeros(size(t));
x11     = zeros(size(t));
phi12   = zeros(size(t));

dphi3    = zeros(size(t));
dphi4    = zeros(size(t));
dx5      = zeros(size(t));
dphi6    = zeros(size(t));
dphi7    = zeros(size(t));
dphi8    = zeros(size(t));
dx9      = zeros(size(t));
dphi10   = zeros(size(t));
dx11     = zeros(size(t));
dphi12   = zeros(size(t));

ddphi3    = zeros(size(t));
ddphi4    = zeros(size(t));
ddx5      = zeros(size(t));
ddphi6    = zeros(size(t));
ddphi7    = zeros(size(t));
ddphi8    = zeros(size(t));
ddx9      = zeros(size(t));
ddphi10   = zeros(size(t));
ddx11     = zeros(size(t));
ddphi12   = zeros(size(t));

% fsolve options (help fsolve, help optimset)
optim_options = optimset('Display','off');

% *** loop over positions ***
Ts = t(2) - t(1);      % timestep
t_size = size(t,1);    % number of simulation steps
for k=1:t_size
    
    % *** position analysis ***
    
    % fsolve solves the non-linear set of equations
    % loop closure equations: see loop_closure_eqs.m
    % argument loop_closure_eqs: file containing closure equations
    % argument [..]': initial values of unknown angles phi3 and phi4
    % argument optim options: parameters for fsolve
    % argument phi2(k): input angle phi2 for which we want to calculate the unknown angles phi3 and phi4
    % argument a1 ... phi1: constants
    % return value x: solution for the unknown angles phi3 and phi4
    % return exitflag: indicates convergence of algorithm
    [x, fval, exitflag] = fsolve('loop_closure_eqs', ...
                                 [phi3_init, phi4_init, x5_init, phi6_init, phi7_init, phi8_init, x9_init, phi10_init, x11_init, phi12_init]', ...
                                 optim_options, ...
                                 phi2(k), ...
                                 r2l, r2k, r3, a, b, r6l, r6k, r7, r8l, r8k, r10, r11, r12, x4, y4, x7, y7, y9);
    if (exitflag ~= 1)
        display 'The fsolve exit flag was not 1, probably no convergence!'
        exitflag
    end
    
    % Declaratie nieuwe variabele
    r6 = r6k + r6l;
    dphi2(k)= omega;
    ddphi2(k) = alpha;
    
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
    
    % *** velocity analysis ***
    A = [0,                 0,                  0,                  0,                  0,                  0,                  0,                  0,                  1,                  -r12*sin(phi12(k));
         0,                 0,                  0,                  0,                  0,                  0,                  0,                  0,                  0,                  r12*cos(phi12(k));
         0,                 0,                  0,                  0,                  0,                  -r8l*sin(phi8(k)),  -1,                 r10*sin(phi10(k)),  1,                  0;
         0,                 0,                  0,                  0,                  0,                  r8l*cos(phi8(k)),  0,                  -r10*cos(phi10(k)), 0,                  0;
         0,                 -x5(k)*sin(phi4(k)),cos(phi4(k)),       -r6k*sin(phi6(k)),  r7*sin(phi7(k)),    0,                  0,                  0,                  0,                  0;
         0,                 x5(k)*cos(phi4(k)), sin(phi4(k)),       r6k*cos(phi6(k)),   -r7*cos(phi7(k)),   0,                  0,                  0,                  0,                  0;
         0,                 0,                  0,                  -r6*sin(phi6(k)),   r7*sin(phi7(k)),    -r8k*sin(phi8(k)),  -1,                 0,                  0,                  0;
         0,                 0,                  0,                  r6*cos(phi6(k)),    -r7*cos(phi7(k)),   -r8k*cos(phi8(k)),  0,                  0,                  0,                  0;
         r3*sin(phi3(k)),   -a*sin(phi4(k))+b*cos(phi4(k)),  0,      0,                  0,                  0,                  0,                  0,                  0,                  0;
         -r3*cos(phi3(k)),  a*cos(phi4(k))+b*sin(phi4(k)), 0,      0,                  0,                  0,                  0,                  0,                  0,                  0];
     
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
    
    % *** acceleration analysis ***
    
    % A = zelfde als daarnet
    
    B = [-r2l*sin(phi2(k))*ddphi2(k) - r2l*cos(phi2(k))*dphi2(k)^2 + r12*cos(phi12(k))*dphi12(k)^2;
         r2l*cos(phi2(k))*ddphi2(k) - r2l*sin(phi2(k))*dphi2(k)^2 + r12*sin(phi12(k))*dphi12(k)^2;
         -r10*cos(phi10(k))*dphi10(k)^2 + r8l*cos(phi8(k))*dphi8(k)^2;
         -r10*sin(phi10(k))*dphi10(k)^2 + r8l*sin(phi8(k))*dphi8(k)^2;
         -r7*cos(phi7(k))*dphi7(k)^2 + r6k*cos(phi6(k))*dphi6(k)^2 + 2*sin(phi4(k))*dphi4(k)*dx5(k) + x5(k)*cos(phi4(k))*dphi4(k)^2;
         -r7*sin(phi7(k))*dphi7(k)^2 + r6k*sin(phi6(k))*dphi6(k)^2 - 2*cos(phi4(k))*dphi4(k)*dx5(k) + x5(k)*sin(phi4(k))*dphi4(k)^2;
         -r7*cos(phi7(k))*dphi7(k)^2 + r6*cos(phi6(k))*dphi6(k)^2 + r8k*cos(phi8(k))*dphi8(k)^2;
         -r7*sin(phi7(k))*dphi7(k)^2 + r6*sin(phi6(k))*dphi6(k)^2 - r8k*sin(phi8(k))*dphi8(k)^2;
         -r2k*cos(phi2(k))*ddphi2(k) + r2k*sin(phi2(k))*dphi2(k)^2 - r3*cos(phi3(k))*dphi3(k)^2 + (a*cos(phi4(k)) + b*sin(phi4(k)))*dphi4(k)^2;
         -r2k*sin(phi2(k))*ddphi2(k) - r2k*cos(phi2(k))*dphi2(k)^2 - r3*sin(phi3(k))*dphi3(k)^2 + (a*sin(phi4(k)) - b*cos(phi4(k)))*dphi4(k)^2];
    
    x = A\B;
    
    % save results
    ddphi3(k)  = x(1);
    ddphi4(k)  = x(2);
    ddx5(k)    = x(3);
    ddphi6(k)  = x(4);
    ddphi7(k)  = x(5);
    ddphi8(k)  = x(6);
    ddx9(k)    = x(7);
    ddphi10(k) = x(8);
    ddx11(k)   = x(9);
    ddphi12(k) = x(10);
    
    
    % *** calculate initial values for next iteration step ***
    phi3_init  = phi3(k) + Ts*dphi3(k);
    phi4_init  = phi4(k) + Ts*dphi4(k);
    x5_init    = x5(k) + Ts*dx5(k);
    phi6_init  = phi6(k) + Ts*dphi6(k);
    phi7_init  = phi7(k) + Ts*dphi7(k);
    phi8_init  = phi8(k) + Ts*dphi8(k);
    x9_init    = x9(k) + Ts*dx9(k);
    phi10_init = phi10(k) + Ts*dphi10(k);
    x11_init   = x11(k) + Ts*dx11(k);
    phi12_init = phi12(k) + Ts*dphi12(k);
    
    
end % loop over positions



% *** create movie ***

% startpunt
% nomenclatuur: scharnier tussen staaf x en staaf y = sch_x_y
sch_1_2 = 0;
% andere vaste punten (coördinaten = complexe getallen)
sch_1_7 = (x7 + j*y7)*exp(j*phi1);
sch_1_4 = (x4 + j*y4)*exp(j*phi1);

% define which positions we want as frames in our movie
frames = 40;    % number of frames in movie
delta = floor(t_size/frames); % time between frames
index_vec = [1:delta:t_size]';

% Create a window large enough for the whole mechanisme in all positions, to prevent scrolling.
% This is done by plotting a diagonal from (x_left, y_bottom) to (x_right, y_top), setting the
% axes equal and saving the axes into "movie_axes", so that "movie_axes" can be used for further
% plots.
x_left = -1.5*r2l;
y_bottom = -1.5*max(r2l, r11);
x_right = x7 + 1.5*(r7 + r6);
y_top = 1.2*max(y4,max(y7,y9)); 

figure(10)
hold on
plot([x_left, x_right], [y_bottom, y_top]);
axis equal;
movie_axes = axis;   %save current axes into movie_axes

% draw and save movie frame
for m=1:length(index_vec)
    index = index_vec(m);
    
    % Volgens mij moet je ze in de volgorde van de loop definiëren ten
    % opzichte van elkaar.
    
    sch_2_12  = sch_1_2 + r2l*exp(j*(phi2(index) + pi));
    sch_2_3   = sch_1_2 + r2k*exp(j*(phi2(index) - pi/2));
    sch_3_4a  = sch_2_3 + r3*exp(j*(phi3(index) + pi));
    
    sch_11_12 = sch_2_12 + r12*exp(j*phi12(index));
    sch_10_11 = sch_11_12 - j*r11;
    sch_8_101 = sch_10_11 + r10*exp(j*(phi10(index) + pi));
        
    sch_8_9   = r2l + r12 + x9(index) + j*y9; % De bovenste van de twee zuigers
        
    sch_6_7   = sch_1_7 + r7*exp(j*(phi7(index) + pi));
    sch_5_6   = sch_6_7 + r6k*exp(j*(phi6(index))); % Het blokje dat over de roterende staaf glijdt
    sch_6_8   = sch_5_6 + r6l*exp(j*phi6(index));
    sch_8_102  = sch_6_8 + (r8l + r8k)*exp(j*(phi8(index) + pi));
    
    hoekpunt_4 = sch_1_4 + a*exp(j*(phi4(index) + pi)); % Het hoekpunt van staaf 4, staat ook naar 'beneden' gericht + verkeerde afstand genomen
    sch_3_4b   = hoekpunt_4 + b*exp(j*(phi4(index) + pi/2));
    
    
    staaf2 = [sch_1_2 sch_2_3 sch_2_12 sch_1_2]; % De driehoekige staaf 2
    loop1  = [sch_2_12 sch_11_12 sch_10_11 sch_8_101];
    loop2  = [sch_1_7 sch_6_7 sch_5_6 sch_6_8 sch_8_9 sch_8_102]; 
    staaf4 = [sch_1_4 hoekpunt_4 sch_3_4b];
    staaf3  = [sch_2_3 sch_3_4a sch_3_4b];
    knoop810 = [sch_8_101 sch_8_102];
      
    
    figure(10)
    clf
    hold on
    plot(real(staaf2), imag(staaf2), '-o')
    plot(real(loop1), imag(loop1), '-o')
    plot(real(loop2), imag(loop2), '-o')
    plot(real(staaf4), imag(staaf4), '-o')
    plot(real(staaf3), imag(staaf3), '-o')
    plot(real(knoop810), imag(knoop810), '-o')
       
    
    axis(movie_axes);     % set axes as in movie_axes
    Movie(m) = getframe;  % save frame to a variable Film
end

% save movie
save fourbar_movie Movie
close(10)


% *** plot figures ***

if fig_kin_4bar
    
    %plot assembly at a certain timestep 
    index = 1; %select 1st timestep
    sch_1_2 = 0;
    sch_1_7 = (x7 + j*y7)*exp(j*phi1);
    sch_1_4 = (x4 + j*y4)*exp(j*phi1);
    
    % Volgens mij moet je ze in de volgorde van de loop definiëren ten
    % opzichte van elkaar.
    
    sch_2_12  = sch_1_2 + r2l*exp(j*(phi2(index) + pi));
    sch_2_3   = sch_1_2 + r2k*exp(j*(phi2(index) - pi/2)); 
    sch_3_4a  = sch_2_3 + r3*exp(j*(phi3(index) + pi));
    
    sch_11_12 = sch_2_12 + r12*exp(j*phi12(index));
    sch_10_11 = sch_11_12 - j*r11;
    sch_8_101 = sch_10_11 + r10*exp(j*(phi10(index) + pi));
        
    sch_8_9   = r2l + r12 + x9(index) + j*y9; % De bovenste van de twee zuigers
        
    sch_6_7   = sch_1_7 + r7*exp(j*(phi7(index) + pi));
    sch_5_6   = sch_6_7 + r6k*exp(j*(phi6(index))); % Het blokje dat over de roterende staaf glijdt
    sch_6_8   = sch_5_6 + r6l*exp(j*phi6(index));
    sch_8_102  = sch_6_8 + (r8l + r8k)*exp(j*(phi8(index) + pi));
    
    hoekpunt_4 = sch_1_4 + a*exp(j*(phi4(index) + pi)); % Het hoekpunt van staaf 4, staat ook naar 'beneden' gericht + verkeerde afstand genomen
    sch_3_4b   = hoekpunt_4 + b*exp(j*(phi4(index) + pi/2));
    
        
    figure
    
    staaf2 = [sch_1_2 sch_2_3 sch_2_12 sch_1_2]; % De driehoekige staaf 2
    loop1  = [sch_2_12 sch_11_12 sch_10_11 sch_8_101];
    loop2  = [sch_1_7 sch_6_7 sch_5_6 sch_6_8 sch_8_9 sch_8_102]; 
    staaf4 = [sch_1_4 hoekpunt_4 sch_3_4b];
    staaf3  = [sch_2_3 sch_3_4a sch_3_4b];
    knoop810 = [sch_8_101 sch_8_102];
    
    
    plot(real(knoop810), imag(knoop810),'ro-')
    
    hold on;
    plot(real(loop1), imag(loop1), 'ro-')
    plot(real(loop2), imag(loop2), 'ro-')
    plot(real(staaf2), imag(staaf2), 'ro-')
    plot(real(staaf4), imag(staaf4), 'ro-')
    plot(real(staaf3), imag(staaf3), 'ro-')
    xlabel('[m]')
    ylabel('[m]')
    title('assembly')
    axis equal
    
    figure
    subplot(311)
    plot(t,phi2)
    ylabel('\phi_2 [rad]')
    subplot(312)
    plot(t,phi3)
    ylabel('\phi_3 [rad]')
    subplot(313)
    plot(t,phi4)
    ylabel('\phi_4 [rad]')
    xlabel('t [s]')
    
    figure
    subplot(311)
    plot(t,dphi2)
    ylabel('d\phi_2 [rad/s]')
    subplot(312)
    plot(t,dphi3)
    ylabel('d\phi_3 [rad/s]')
    subplot(313)
    plot(t,dphi4)
    ylabel('d\phi_4 [rad/s]')
    xlabel('t [s]')
    
    figure
    subplot(311)
    plot(t,ddphi2)
    ylabel('dd\phi_2 [rad/s^2]')
    subplot(312)
    plot(t,ddphi3)
    ylabel('dd\phi_3 [rad/s^2]')
    subplot(313)
    plot(t,ddphi4)
    ylabel('dd\phi_4 [rad/s^2]')
    xlabel('t [s]')
end