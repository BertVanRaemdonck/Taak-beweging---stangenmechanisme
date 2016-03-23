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
rev_arm_angle = asin(1.30/2.344);                 % hoek tussen horizontale en lifting arm in de uiterste stand
x7 = (6.164 + 2.75 * cos(rev_arm_angle)) * conv;  % horizontale positie scharnier 1-7 = reverse arm pivot hor + lifting arm length (met hoekcorrectie)
y7 = (5.031 - 2.75 * sin(rev_arm_angle)) * conv;  % verticale positie scharnier 1-7 = reverse arm pivot vert + ligting arm length (met hoekcorrectie)
y9 = 3.5 * conv;            % valve CL to cyl CL

phi1 = 0;                   % omdat dat ook in de voorbeelden zo staat



% dynamic parameters, defined in a local frame on each of the bars.
% NOG EENS GOED NAKIJKEN, IN PRINCIPE DAN ENKEL X COORDINAAT voor stangen?
X3 = r3/2;                  % zwaartepunt
X4 = a/2;
X5 = 0;                     % in lokaal assenstelsel is cog schuifscharnier in scharnierpunt
X6k = r6k/2;
X6l = r6l/2;
X6 = (r6k + r6l)/2;
X7 = r7/2;
X8k = r8k/2;
X8l = r8l/2;
X8 = (r8k + r8l)/2;
X9 = 0;                  
X10 = r10/2;
X11 = r11/2;                % NOG AAN TE PASSEN, IN VERGELIJKINGEN IS DAT HET ZWAARTEPUNT VAN STANG 11 + DE PISTON ZIJN!!!!
X12 = r12/2;

% driehoek is geen driehoek maar cirkel dus cog in scharnier (kan eventueel
% verwijderd worden)
X2 = 0;
Y2 = 0;

% zwaartepunten van pistons ook bepalen? (dan wel nieuw symbool voor 11
%       nodig want r11 is al in gebruik)

Y3 = 0;                     % Y coordinates of cog
Y4 = b/2;
Y5 = 0;                     % in lokaal assenstelsel is cog schuifscharnier in scharnierpunt
Y6k = 0;
Y6l = 0;
Y6 = 0;
Y7 = 0;
Y8k = 0;
Y8l = 0;
Y8 = 0;
Y9 = 0;
Y10 = 0;
Y11 = 0;
Y12 = 0;


% Definitie extra parameters:   (voorlopig gekozen)
breedte5 = 0.05;            % breedte schuifscharnier
hoogte5 = 0.05;             % hoogte schuifscharnier
breedte9 = 0.2;             % breedte piston1, stang 9
hoogte9 = 0.1;              % hoogte piston1, stang 9
breedte11 = 0.3;            % breedte piston2, stang 11
hoogte11 = 0.2;             % hoogte piston2, stang11

L9 = r8l - X8;              % Afstand tussen scharnierpunt 8,9 en het zwaartepunt van stang 8

% massa's (nog na te kijken)
rho_l1 = 14.72;             % massa per lengte stang, alles behalve drijfstang
rho_l2 = 71.76;             % massa per lengte drijfstang
rho_A = 10;                 % massa per oppervlakte van het element (bij stang 2 driehoek), voorlopig genomen als waarde

m2 = ((r2k*r2l)/2) * rho_A;     % totale massa van stang 2 aangezien driehoekige stang
m3 = r3 * rho_l1;
m4 = (a+b) * rho_l1;         % totale massa van stang 4, ma en mb zijn de massa's van de aparte delen
ma = a * rho_l1;
mb = b * rho_l1;
m5 = 1;                      % massa van schuifscharnier
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

m_piston_1 = 15;            % voorlopig gekozen, stang 9
m_piston_2 = 20;            % voorlopig gekozen, stang 11

m9 = m_piston_1;            % extra nodig 

%Vanaf hieronder nog aanpassen en controleren:

J2 = (r2l*(r2k^3))/12;   % door het vaste punt, zie "https://en.wikipedia.org/wiki/List_of_area_moments_of_inertia"
J3 = m3*r3^2/12;
J4 = m4*a^2/12;         % NU FOUTIEF; Moeten we eens over nadenken hoe we deze opstellen, moet door vaste punt, zie eventueel site hierboven
J5 = m5*((hoogte5^2)+(breedte5^2)) / 12 ;     % te benaderen als gevulde balk?  zie "https://en.wikipedia.org/wiki/List_of_moments_of_inertia"
J6k = m6k*r6k^2/12;
J6l = m6l*r6l^2/12;
J6 = m6*(r6k+r6l)^2/12;
J7 = m7*r7^2/12;
J8k = m8k*r8k^2/12;
J8l = m8l*r8l^2/12;
J8 = m8*(r8k+r8l)^2/12;
J9 = m_piston_1*((hoogte9^2)+(breedte9^2)) / 12 ;    % benaderd als volle balk
J10 = m10*r10^2/12;
J11 = m_piston_2*((hoogte11^2)+(breedte11^2)) / 12;  % benaderd als volle balk
J12 = m12*r12^2/12;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1. Determination of Kinematics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% position analysis

% initial condition for first step of position analysis with fsolve (phi3 and phi4)
% VERY IMPORTANT because it determines which branch of the mechanism you're in
phi_init = [0; 2*pi/3 ; 1.25 ; pi/12 ; 2*pi/3 ; 7*pi/12 ; 2 ; pi/12 ; 1 ; pi/12];    
        
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
t_end = 100;                   % end time of simulation
Ts = 0.05;                     % time step of simulation
t = [t_begin:Ts:t_end]';       % time vector

% initialization of driver
omega = -0.1;
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


[F12x, F12y, F23x, F23y, F212x, F212y, F34x, F34y, F14x, F14y, F45, F56x, F56y, F67x, F67y, ...
   F68x, F68y, F17x, F17y, F89x, F89y, F810x, F810y, F19, F1011x, F1011y, F1112x, F1112y, F111, ...
   M12, M19, M111, M45] ...
   = dynamics_4bar(phi2,  phi3,  phi4,  x5,  phi6,  phi7,  phi8,  x9,  phi10,  x11,  phi12, ...
                   dphi2, dphi3, dphi4, dx5, dphi6, dphi7, dphi8, dx9, dphi10, dx11, dphi12, ...
                   ddphi2,ddphi3,ddphi4,ddx5,ddphi6,ddphi7,ddphi8,ddx9,ddphi10,ddx11,ddphi12, ...
                   r2l, r2k, r3, a, b, r6l, r6k, r7, r8l, r8k, r10, r11, r12, x4, y4, x7, y7, y9, L9, ...
                   m2,m3,ma,mb,m4,m5,m6k,m6l,m6,m7,m8k,m8l,m8,m9,m10,m11,m12, m_piston_1, m_piston_2,...
                   X2,X3,X4,X5,X6k,X6l,X6,X7,X8k,X8l,X8,X9,X10,X11,X12, ...
                   Y2,Y3,Y4,Y5,Y6k,Y6l,Y6,Y7,Y8k,Y8l,Y8,Y9,Y10,Y11,Y12, ...
                   J2,J3,J4,J5,J6k,J6l,J6,J7,J8k,J8l,J8,J9,J10,J11,J12, t,fig_dyn_4bar);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3. Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
load fourbar_movie Movie
movie(Movie)
