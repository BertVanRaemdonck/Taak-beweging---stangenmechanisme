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

clear
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data initialization (all data is converted to SI units)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% program data
fig_kin_4bar = 1;        % draw figures of kinematic analysis if 1
fig_dyn_4bar = 1;        % draw figures of dynamic analysis if 1

% kinematic parameters (link lengths)
r2k = 2.406;                % eccentric crank circle diameter
r2l = 2.0;                  % 'langere' zijde van de eccentric crank
r3 = 9.750;                 % eccentric rod
a = 1.844;                  % vervang lengte van link crank (dg) (verticaal)
b = 0.454;                  % vervang lengte van link crank (dg) (horizontaal)
r6k = 1.625;                % radius rod extension
r6l = 9.125;                % radius rod
r7 = 1.625;                 % lifting link
r8k = 0.336;                % combination lever upper
r8l = 3.664;                % combination lever lower
r10 = 2.922;                % union link
r11 = 1;                    % drop link to cross head vert
r12 = 16.5;                 % main rod

x4 = 10.330;                % zelf berekend
y4 = 3.538;                 % zelf berekend
x7 = 8.442;                 % zelf berekend
y7 = 3.474;                 % zelf berekend
y9 = 3.789;                 % zelf berekend


% dynamic parameters, defined in a local frame on each of the bars.
% NOG EENS GOED NAKIJKEN, IN PRINCIPE DAN ENKEL X COORDINAAT voor stangen?
X3 = r3/2;                  % zwaartepunt
X4 = a/2;
X6 = (r6k + r6l)/2;
X7 = r7/2;
X8 = (r8k + r8l)/2;
X10 = r10/2;
X11 = r11/2;
X12 = r12/2;

% zwaartepunt van driehoek constructie 2 nog berekenen voor algemene geval
X2 = r2k/3;
Y2 = r2l/3;
% zwaartepunten van pistons ook bepalen? (dan wel nieuw symbool voor 11
%       nodig want r11 is al in gebruik)

Y3 = 0;                     % Y coordinates of cog
Y4 = b/2;
Y6 = 0;
Y7 = 0;
Y8 = 0;
Y10 = 0;
Y11 = 0;
Y12 = 0;

%Vanaf hieronder nog aanpassen:
m2 = r2*1.76;
m3 = r3*1.76;
m4 = r4*0.54;

J2 = m2*r2^2/12;
J3 = m3*r3^2/12;
J4 = m4*r4^2/12;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1. Determination of Kinematics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% position analysis

% initial condition for first step of position analysis with fsolve (phi3 and phi4)
% VERY IMPORTANT because it determines which branch of the mechanism you're in
phi_init = [Pi; 2*Pi/3 ; 1.25 ; Pi/12 ; 2*Pi/3 ; 7*Pi/12 ; 2 ; 13*Pi/12 ; 1 ; Pi/12]    
        %phi3=phi_init(1); %phi4=phi_init(2); %x5=phi_init(3); %phi6=phi_init(4);
        %phi7=phi_init(5); %phi8=phi_init(6); %x9=phi_init(7); %phi10=phi_init(8);
        %x11=phi_init(9);  %phi12=phi_init(10);


t_begin = 0;                   % start time of simulation
t_end = 10;                    % end time of simulation
Ts = 0.05;                     % time step of simulation
t = [t_begin:Ts:t_end]';       % time vector

% initialization of driver
omega = 0.5;
A = 1;
phi2=1+A*sin(omega*t);
dphi2=omega*A*cos(omega*t);
ddphi2=-omega^2*A*sin(omega*t);

% calculation of the kinematics (see kin_4bar.m)
[phi3,phi4,dphi3,dphi4,ddphi3,ddphi4] = kinematics_4bar(r1,r2,r3,r4,phi1,phi2,dphi2,ddphi2,phi3_init,phi4_init,t,fig_kin_4bar);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2. Dynamics Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculation of the dynamics (see dyn_4bar.m)
[F_P_x,F_Q_x,F_R_x,F_S_x,F_P_y,F_Q_y,F_R_y,F_S_y,M_P] = dynamics_4bar(phi2,phi3,phi4,dphi2,dphi3,dphi4,ddphi2,ddphi3,ddphi4,r2,r3,r4, ...
  m2,m3,m4,X2,X3,X4,Y2,Y3,Y4,J2,J3,J4,t,fig_dyn_4bar);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3. Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
load fourbar_movie Movie
movie(Movie)

