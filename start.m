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
conv = 0.19;                % conversie factor reële waardes en lengtes simulatie programma

r2k = 1.203 * conv;         % eccentric crank circle diameter
r2l = 2.0 * conv;           % 'langere' zijde van de eccentric crank
r3 = 9.750 * conv;          % eccentric rod
a = 3.021 * conv;           % vervang lengte van link crank (dg) (verticaal) = link crank vert
b = 0.801 * conv;           % vervang lengte van link crank (dg) (horizontaal) = link crank back set
r6k = 1.625 * conv;         % radius rod extension
r6l = 9.125 * conv;         % radius rod
r7 = 1.625 * conv;          % lifting link
r8k = 0.336 * conv;         % combination lever upper
r8l = 3.664 * conv;         % combination lever lower
r10 = 2.922 * conv;         % union link
r11 = 1 * conv;             % drop link to cross head vert
r12 = 16.5 * conv;          % main rod

x4 = 10.297 * conv;         % link center pivot horizontal
y4 = 3.500 * conv;          % link center pivot vertical
x7 = 8.439 * conv;          % gemeten in simulatie programma
y7 = 3.472 * conv;          % gemeten in simulatie programma
y9 = 3.747 * conv;          % gemeten in simulatie programma

phi1 = 0;                   % omdat dat ook in de voorbeelden zo staat


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


% massa's (nog na te kijken)
rho_l1 = 14,72;             % massa per lengte stang, alles behalve drijfstang
rho_l2 = 71.76;             % massa per lengte drijfstang
rho_A = 10;                 % massa per oppervlakte van het element (bij stang 2 driehoek), voorlopig genomen als waarde

m2 = ((r2k*r2l)/2) * rho_A;     % totale massa van stang 2 aangezien driehoekige stang
m3 = r3 * rho_l1;
m4 = (a+b) * rho_l1;         % totale massa van stang 4, ma en mb zijn de massa's van de aparte delen
ma = a * rho_l1;
mb = b * rho_l1;
m6 = (r6k+r6l) * rho_l1;     % totale massa van stang 6, m6k en m6l zijn de massa's van de aparte delen
m6k = r6k * rho_l1;
m6l = r6l * rho_l1;
m7 = r7 * rho_l1;
m8 = (r8k +r8l) * rho_l1;    % totale massa van stang 8, m8k en m8l zijn de massa's van de aparte dele
m8k = r8k * rho_l1;
m8l = r8l * rho_l1;
m10 = r10 * rho_l1;
m11 = r11 * rho_l1;
m12 = r12 * rho_l2;

m_piston_1 = 20;            % voorlopig gekozen
m_piston_2 = 15;            % voorlopig gekozen


%Vanaf hieronder nog aanpassen:

%J2 = m2*r2^2/12;
%J3 = m3*r3^2/12;
%J4 = m4*r4^2/12;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1. Determination of Kinematics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% position analysis

% initial condition for first step of position analysis with fsolve (phi3 and phi4)
% VERY IMPORTANT because it determines which branch of the mechanism you're in
phi_init = [pi; 2*pi/3 ; 1.25 ; pi/12 ; 2*pi/3 ; 7*pi/12 ; 2 ; 13*pi/12 ; 1 ; pi/12];    
        
phi3_init =  phi_init(1);
phi4_init =  phi_init(2);
x5_init =    phi_init(3);
phi6_init =  phi_init(4);
phi7_init =  phi_init(5);
phi8_init =  phi_init(6);
x9_init =    phi_init(7);
phi10_init = phi_init(8);
x11_init =   phi_init(9);
phi12_init = phi_init(10);


t_begin = 0;                   % start time of simulation
t_end = 10;                    % end time of simulation
Ts = 0.05;                     % time step of simulation
t = [t_begin:Ts:t_end]';       % time vector

% initialization of driver
omega = 0.5;
alpha = 0;
phi2 = omega*t + pi/2;
dphi2 = omega;
ddphi2 = alpha;

% calculation of the kinematics (see kin_4bar.m)
[   phi3,   phi4,   x5,     phi6,   phi7,   phi8,   x9,     phi10,      x11,    phi12, ... 
    dphi3,  dphi4,  dx5,    dphi6,  dphi7,  dphi8,  dx9,    dphi10,     dx11,   dphi12, ...
    ddphi3, ddphi4, ddx5,   ddphi6, ddphi7, ddphi8, ddx9,   ddphi10,    ddx11,  ddphi12   ] = ...
    kinematics_4bar(r2l, r2k, r3, a, b, r6l, r6k, r7, r8l, r8k, r10, r11, r12, x4, y4, x7, y7, y9, ...
                    phi1, phi2, dphi2, ddphi2, omega, alpha, ...
                    phi3_init, phi4_init, x5_init, phi6_init, phi7_init, phi8_init, x9_init, phi10_init, x11_init, phi12_init, ...
                    t, fig_kin_4bar);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2. Dynamics Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculation of the dynamics (see dyn_4bar.m)

% Uitgecommentarieerd zodat dit geen error geeft, aangezien we hier nog
% niet mee bezig zijn.

% [F_P_x,F_Q_x,F_R_x,F_S_x,F_P_y,F_Q_y,F_R_y,F_S_y,M_P] = dynamics_4bar(phi2,phi3,phi4,dphi2,dphi3,dphi4,ddphi2,ddphi3,ddphi4,r2,r3,r4, ...
 % m2,m3,m4,X2,X3,X4,Y2,Y3,Y4,J2,J3,J4,t,fig_dyn_4bar);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3. Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
load fourbar_movie Movie
movie(Movie)

