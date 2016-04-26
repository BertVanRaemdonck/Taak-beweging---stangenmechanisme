clear;
close all;

%% Declaration of variables and loading in matcam

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT
%
% To get this program to run smoothly, please fill in the location of the
% motion law and external load files on your computer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mot_law_location = 'C:\Users\Bert\School\Beweging en trillingen\Code\nok_hefwet.mot'; 
%mot_law_location = 'E:\Data\KULeuven\3de bachelor\2de semester\Beweging\Taak nokken\Taak-beweging-stangenmechanisme\nok_hefwet.mot';
ext_load_location = 'C:\Users\Bert\School\Beweging en trillingen\Code\nok_externe_krachten.exl';
%ext_load_location = 'E:\Data\KULeuven\3de bachelor\2de semester\Beweging\Taak nokken\Taak-beweging-stangenmechanisme\nok_externe_krachten.exl';

% assigned values
m_follower = 20;
zeta = 0.077;
T_cycle = 0.5;
omega = 2*pi/T_cycle;
cam_rpm = 60/T_cycle;

R_0 = 80;
R_v = 20;

matcam()

% setting the additional parameters for the analysis
global genmotlaw_startangle_edit genmotlaw_endangle_edit
global genmotlaw_startlift_edit genmotlaw_endlift_edit
global bcr_edit rof_edit exc_edit contourgrad_edit          
global mass_edit spring_edit sprload_edit rpm_edit
set(genmotlaw_startangle_edit, 'string', '330');            % not relevant for calculations, just keep matcam from throwing errors
set(genmotlaw_endangle_edit, 'string', '360');              % not relevant for calculations, just keep matcam from throwing errors
set(genmotlaw_startlift_edit, 'string', '0');               % not relevant for calculations, just keep matcam from throwing errors
set(genmotlaw_endlift_edit, 'string', '0');                 % not relevant for calculations, just keep matcam from throwing errors
set(bcr_edit, 'string', num2str(R_0));
set(rof_edit, 'string', num2str(R_v));
set(exc_edit, 'string', '0');                               % Analysis first without an eccentric follower
set(contourgrad_edit, 'string', '0');
set(mass_edit, 'string', num2str(m_follower));
set(spring_edit, 'string', '0');                            % Analysis without a spring first, see section 3) for calculation of the spring
set(sprload_edit, 'string', '0');                           % Analysis without a spring first, see section 3) for calculation of the spring
set(rpm_edit, 'string', num2str(cam_rpm));

% loading the files with the motion law and load profile
if exist(mot_law_location, 'file') == 2
    matcam('genmotlawload', mot_law_location)
else
    matcam('genmotlawload')
end

if exist(ext_load_location, 'file') == 2
    matcam('genextloadload', ext_load_location)
else
    matcam('genextloadload')
end

matcam('calc')
close all

%% Making matcam variables locally accessible

global Stot Vtot Atot;    
S = Stot;                                       % Lift in mm
V = Vtot;                                       % Velocity in mm/degree
A = Atot;                                       % Acceleration in mm/degree^2
global tetatot tetatotrad rpm;
theta = tetatot;                                % Cam angle in degrees
theta_rad = tetatotrad;                         % Cam angle in radians
global alfa;
alpha = alfa;                                   % Pressure angle
global force_spring force_load;
F_spring = force_spring;                        % Spring force
F_load = force_load;                            % External load
global force_acc force_tot;
F_acc = force_acc;                              % Inertial force
F_tot = force_tot;                              % Total force
global force_x force_y;
F_x = force_x;                                  % Total force in x direction
F_y = force_y;                                  % Total force in y direction
global roc;                                     % Radius curvature
radius_of_curvature = roc;
global Rtotrad

% test to see if the loading of the variables workes
%figure()
%plot(alpha)

%% 2) Calculating motion laws

figure()
plot(S)                                     % plotting the lift
title('Graph of lift')

figure()
plot(V)                                     % plotting the velocity
title('Graph of velocity')

figure()
plot(A)                                     % plotting the acceleration
title('Graph of acceleration')



%% 2) Calculating cam dimensions:

% Without eccentricity

figure()
plot(radius_of_curvature)                       % plotting radius of curvature
axis([0 360 0 500])
title('Graph of radius of curvature without eccentricity')

figure()
plot(alpha)                                     % plotting pressure angle
title('Graph of pressure angle without eccentricity')

% With eccentricity

eccentricity = 10;                              % Chosen by hand


if 1 == 1                                       % Recalibration matcam
    
    matcam()
    % setting the additional parameters for the analysis
    set(genmotlaw_startangle_edit, 'string', '330');            % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_endangle_edit, 'string', '360');              % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_startlift_edit, 'string', '0');               % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_endlift_edit, 'string', '0');                 % not relevant for calculations, just keep matcam from throwing errors
    set(bcr_edit, 'string', num2str(R_0));
    set(rof_edit, 'string', num2str(R_v));
    set(exc_edit, 'string', num2str(eccentricity));             % Analysis first without an eccentric follower
    set(contourgrad_edit, 'string', '0');
    set(mass_edit, 'string', num2str(m_follower));
    set(spring_edit, 'string', '0');                            % Analysis without a spring first, see section 3) for calculation of the spring
    set(sprload_edit, 'string', '0');                           % Analysis without a spring first, see section 3) for calculation of the spring
    set(rpm_edit, 'string', num2str(cam_rpm));

    % loading the files with the motion law and load profile
    if exist(mot_law_location, 'file') == 2
        matcam('genmotlawload', mot_law_location)
    else
        matcam('genmotlawload')
    end

    if exist(ext_load_location, 'file') == 2
        matcam('genextloadload', ext_load_location)
    else
        matcam('genextloadload')
    end

    matcam('calc')
    % Making matcam variables locally accessible

    S = Stot;                                       % Lift in mm
    V = Vtot;                                       % Velocity in mm/degree
    A = Atot;                                       % Acceleration in mm/degree^2

    theta = tetatot;                                % Cam angle in degrees
    theta_rad = tetatotrad;                         % Cam angle in radians

    alpha = alfa;                                   % Pressure angle

    F_spring = force_spring;                        % Spring force
    F_load = force_load;                            % External load

    F_acc = force_acc;                              % Inertial force
    F_tot = force_tot;                              % Total force

    F_x = force_x;                                  % Total force in x direction
    F_y = force_y;                                  % Total force in y direction

    radius_of_curvature = roc;                      % Radius curvature

    % Close matcam window
    matcam_figure = gcf;
    close(matcam_figure.Number)
end 



figure()
plot(radius_of_curvature)                       % plotting radius of curvature
axis([0 360 -500 500])
title('Graph of radius of curvature with eccentricity')

figure()
plot(alpha)                                     % plotting pressure angle
title('Graph of pressure angle with eccentricity')


%% 3) Calculation spring constant:

figure()
plot(F_tot)                                     % plotting total contact force without a spring
title('Graph of total contact force without a spring')


F_v0 = 20;                                      % Spring preload [N]   (chosen by hand)
g = 9.81;                                       % Gravitational acceleration [m/s^2]
gamma = 0;                                      % Angle between follower and vertical axis [°]

conversie_factor = (180^2)/(pi^2*1000);         % Conversion factor of converting [(kg*rad^2)/(s^2*°^2)] to [(N*rad^2)/mm]

% Min teken moet erbij, maar doet dan ambetant
k = max((-F_load - F_v0*ones(size(S)) - m_follower.*omega.*omega.*A.*conversie_factor )./S)               % - follower_mass*g*cos(gamma)*ones(size(S))  weggedaan omdat zwaartekracht te verwaarlozen is

k = 5*ceil(k/5);                                % Rounds k up to the next integer wich is a multiple of 5, in order to have a strong enough spring


if 1 == 1                                       % Recalibration matcam
    
    matcam()
    % setting the additional parameters for the analysis
    set(genmotlaw_startangle_edit, 'string', '330');            % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_endangle_edit, 'string', '360');              % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_startlift_edit, 'string', '0');               % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_endlift_edit, 'string', '0');                 % not relevant for calculations, just keep matcam from throwing errors
    set(bcr_edit, 'string', num2str(R_0));
    set(rof_edit, 'string', num2str(R_v));
    set(exc_edit, 'string', num2str(eccentricity));             % Analysis with an eccentric follower
    set(contourgrad_edit, 'string', '0');
    set(mass_edit, 'string', num2str(m_follower));
    set(spring_edit, 'string', num2str(k));                     % Analysis with a spring 
    set(sprload_edit, 'string', num2str(F_v0));                 % Analysis with a spring
    set(rpm_edit, 'string', num2str(cam_rpm));

    % loading the files with the motion law and load profile
    if exist(mot_law_location, 'file') == 2
        matcam('genmotlawload', mot_law_location)
    else
        matcam('genmotlawload')
    end

    if exist(ext_load_location, 'file') == 2
        matcam('genextloadload', ext_load_location)
    else
        matcam('genextloadload')
    end

    matcam('calc')
    % Making matcam variables locally accessible

    S = Stot;                                       % Lift in mm
    V = Vtot;                                       % Velocity in mm/degree
    A = Atot;                                       % Acceleration in mm/degree^2

    theta = tetatot;                                % Cam angle in degrees
    theta_rad = tetatotrad;                         % Cam angle in radians

    alpha = alfa;                                   % Pressure angle

    F_spring = force_spring;                        % Spring force
    F_load = force_load;                            % External load

    F_acc = force_acc;                              % Inertial force
    F_tot = force_tot;                              % Total force

    F_x = force_x;                                  % Total force in x direction
    F_y = force_y;                                  % Total force in y direction

    radius_of_curvature = roc;                      % Radius curvature

    % Close matcam window
    matcam_figure = gcf;
    close(matcam_figure.Number)
end 


figure()
plot(F_tot)                                     % plotting total contact force with a spring
title('Graph of total contact force with a spring')


% Doubling the rotation speed, same spring:

double_rpm = 2*cam_rpm;

if 1 == 1                                       % Recalibration matcam
    
    matcam()
    % setting the additional parameters for the analysis
    set(genmotlaw_startangle_edit, 'string', '330');            % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_endangle_edit, 'string', '360');              % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_startlift_edit, 'string', '0');               % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_endlift_edit, 'string', '0');                 % not relevant for calculations, just keep matcam from throwing errors
    set(bcr_edit, 'string', num2str(R_0));
    set(rof_edit, 'string', num2str(R_v));
    set(exc_edit, 'string', num2str(eccentricity));             % Analysis with an eccentric follower
    set(contourgrad_edit, 'string', '0');
    set(mass_edit, 'string', num2str(m_follower));
    set(spring_edit, 'string', num2str(k));                     % Analysis with same spring 
    set(sprload_edit, 'string', num2str(F_v0));                 % Analysis with same spring
    set(rpm_edit, 'string', num2str(double_rpm));               % Analysis with double rotation speed

    % loading the files with the motion law and load profile
    if exist(mot_law_location, 'file') == 2
        matcam('genmotlawload', mot_law_location)
    else
        matcam('genmotlawload')
    end

    if exist(ext_load_location, 'file') == 2
        matcam('genextloadload', ext_load_location)
    else
        matcam('genextloadload')
    end

    matcam('calc')
    % Making matcam variables locally accessible

    S = Stot;                                       % Lift in mm
    V = Vtot;                                       % Velocity in mm/degree
    A = Atot;                                       % Acceleration in mm/degree^2

    theta = tetatot;                                % Cam angle in degrees
    theta_rad = tetatotrad;                         % Cam angle in radians

    alpha = alfa;                                   % Pressure angle

    F_spring = force_spring;                        % Spring force
    F_load = force_load;                            % External load

    F_acc = force_acc;                              % Inertial force
    F_tot = force_tot;                              % Total force

    F_x = force_x;                                  % Total force in x direction
    F_y = force_y;                                  % Total force in y direction

    radius_of_curvature = roc;                      % Radius curvature

    % Close matcam window
    matcam_figure = gcf;
    close(matcam_figure.Number)
end 

figure()
plot(F_tot)                                     % plotting total contact force with the same spring and double the rotation speed
title('Graph of total contact force, double rotation speed')

% Doubling the rotation speed, different spring

F_v0_double = 150;                               % Spring preload [N]   (chosen by hand)

double_omega = (2*pi*double_rpm)/60;

% Min teken moet erbij, maar doet dan ambetant
k_double = max((-F_load - F_v0_double*ones(size(S)) - m_follower.*double_omega.*double_omega.*A.*conversie_factor )./S)   % - follower_mass*g*cos(gamma)*ones(size(S)) weggedaan omdat zwaartekracht te verwaarlozen is, geeft ook te zwakke veer

k_double = 5*ceil(k_double/5);                  % Rounds k_double up to the next integer wich is a multiple of 5, in order to have a strong enough spring

if 1 == 1                                       % Recalibration matcam
    
    matcam()
    % setting the additional parameters for the analysis
    set(genmotlaw_startangle_edit, 'string', '330');            % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_endangle_edit, 'string', '360');              % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_startlift_edit, 'string', '0');               % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_endlift_edit, 'string', '0');                 % not relevant for calculations, just keep matcam from throwing errors
    set(bcr_edit, 'string', num2str(R_0));
    set(rof_edit, 'string', num2str(R_v));
    set(exc_edit, 'string', num2str(eccentricity));             % Analysis with an eccentric follower
    set(contourgrad_edit, 'string', '0');
    set(mass_edit, 'string', num2str(m_follower));
    set(spring_edit, 'string', num2str(k_double));              % Analysis with different spring 
    set(sprload_edit, 'string', num2str(F_v0_double));          % Analysis with different spring
    set(rpm_edit, 'string', num2str(double_rpm));               % Analysis with double rotation speed

    % loading the files with the motion law and load profile
    if exist(mot_law_location, 'file') == 2
        matcam('genmotlawload', mot_law_location)
    else
        matcam('genmotlawload')
    end

    if exist(ext_load_location, 'file') == 2
        matcam('genextloadload', ext_load_location)
    else
        matcam('genextloadload')
    end

    matcam('calc')
    % Making matcam variables locally accessible

    S = Stot;                                       % Lift in mm
    V = Vtot;                                       % Velocity in mm/degree
    A = Atot;                                       % Acceleration in mm/degree^2

    theta = tetatot;                                % Cam angle in degrees
    theta_rad = tetatotrad;                         % Cam angle in radians

    alpha = alfa;                                   % Pressure angle

    F_spring = force_spring;                        % Spring force
    F_load = force_load;                            % External load

    F_acc = force_acc;                              % Inertial force
    F_tot = force_tot;                              % Total force

    F_x = force_x;                                  % Total force in x direction
    F_y = force_y;                                  % Total force in y direction

    radius_of_curvature = roc;                      % Radius curvature

    % Close matcam window
    matcam_figure = gcf;
    close(matcam_figure.Number)
end 

figure()
plot(F_tot)                                     % plotting total contact force with the same spring and double the rotation speed
title('Graph of total contact force, double rotation speed, different spring')

%% 3) Calculation instantaneous power

if 1 == 1                                       % Recalibration matcam
    
    matcam()
    % setting the additional parameters for the analysis
    set(genmotlaw_startangle_edit, 'string', '330');            % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_endangle_edit, 'string', '360');              % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_startlift_edit, 'string', '0');               % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_endlift_edit, 'string', '0');                 % not relevant for calculations, just keep matcam from throwing errors
    set(bcr_edit, 'string', num2str(R_0));
    set(rof_edit, 'string', num2str(R_v));
    set(exc_edit, 'string', '0');                               % Analysis without eccentricity
    set(contourgrad_edit, 'string', '0');
    set(mass_edit, 'string', num2str(m_follower));
    set(spring_edit, 'string', num2str(k));                     % Analysis with a spring 
    set(sprload_edit, 'string', num2str(F_v0));                 % Analysis with a spring
    set(rpm_edit, 'string', num2str(cam_rpm));

    % loading the files with the motion law and load profile
    if exist(mot_law_location, 'file') == 2
        matcam('genmotlawload', mot_law_location)
    else
        matcam('genmotlawload')
    end

    if exist(ext_load_location, 'file') == 2
        matcam('genextloadload', ext_load_location)
    else
        matcam('genextloadload')
    end

    matcam('calc')
    % Making matcam variables locally accessible

    S = Stot;                                       % Lift in mm
    V = Vtot;                                       % Velocity in mm/degree
    A = Atot;                                       % Acceleration in mm/degree^2

    theta = tetatot;                                % Cam angle in degrees
    theta_rad = tetatotrad;                         % Cam angle in radians

    alpha = alfa;                                   % Pressure angle

    F_spring = force_spring;                        % Spring force
    F_load = force_load;                            % External load

    F_acc = force_acc;                              % Inertial force
    F_tot = force_tot;                              % Total force

    F_x = force_x;                                  % Total force in x direction
    F_y = force_y;                                  % Total force in y direction

    radius_of_curvature = roc;                      % Radius curvature
    
    R_tot = R_0 + R_v ;                             % Total distance between center of rotation and center of follower
    R_tot_rad = Rtotrad ;
    
    
    % Close matcam window
    matcam_figure = gcf;
    close(matcam_figure.Number)
end 

P1 = F_tot .*(omega*ones(size(S))) .*(sind(alpha)) .* ((R_tot + S).*0.001);      % instantaneous power

P1_2 = F_tot .*(omega*ones(size(S))) .*(sind(alpha)) .* ((R_tot_rad).*0.001);    % instantaneous power calculated by values of matcam

figure()
plot(P1)                                     % plotting total contact force with the same spring and double the rotation speed
title('Graph of instantaneous power 1')

if 1 == 1                                       % Recalibration matcam
    
    matcam()
    % setting the additional parameters for the analysis
    set(genmotlaw_startangle_edit, 'string', '330');            % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_endangle_edit, 'string', '360');              % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_startlift_edit, 'string', '0');               % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_endlift_edit, 'string', '0');                 % not relevant for calculations, just keep matcam from throwing errors
    set(bcr_edit, 'string', num2str(R_0));
    set(rof_edit, 'string', num2str(R_v));
    set(exc_edit, 'string', num2str(eccentricity));             % Analysis with eccentricity
    set(contourgrad_edit, 'string', '0');
    set(mass_edit, 'string', num2str(m_follower));
    set(spring_edit, 'string', num2str(k));                     % Analysis with a spring 
    set(sprload_edit, 'string', num2str(F_v0));                 % Analysis with a spring
    set(rpm_edit, 'string', num2str(cam_rpm));

    % loading the files with the motion law and load profile
    if exist(mot_law_location, 'file') == 2
        matcam('genmotlawload', mot_law_location)
    else
        matcam('genmotlawload')
    end

    if exist(ext_load_location, 'file') == 2
        matcam('genextloadload', ext_load_location)
    else
        matcam('genextloadload')
    end

    matcam('calc')
    % Making matcam variables locally accessible

    S = Stot;                                       % Lift in mm
    V = Vtot;                                       % Velocity in mm/degree
    A = Atot;                                       % Acceleration in mm/degree^2

    theta = tetatot;                                % Cam angle in degrees
    theta_rad = tetatotrad;                         % Cam angle in radians

    alpha = alfa;                                   % Pressure angle

    F_spring = force_spring;                        % Spring force
    F_load = force_load;                            % External load

    F_acc = force_acc;                              % Inertial force
    F_tot = force_tot;                              % Total force

    F_x = force_x;                                  % Total force in x direction
    F_y = force_y;                                  % Total force in y direction

    radius_of_curvature = roc;                      % Radius curvature
    
    R_tot = R_0 + R_v;                              % Total distance between center of rotation and center of follower
    R_tot_rad = Rtotrad;
    
     % Close matcam window
%     matcam_figure = gcf;
%     close(matcam_figure.Number)
end 

P2 = F_tot .*(omega*ones(size(S))) .*( ( ( (sqrt(R_tot^2 - eccentricity^2) + S).* 0.001).*sind(alpha)) + (eccentricity*0.001.*cosd(alpha)) );         % instantaneous power
% Wel hetzelfde, formule zou nu moeten kloppen 

P2_2 = F_tot .*(omega*ones(size(S))) .*( ( ( (sqrt(R_tot_rad.^2 - eccentricity^2)).* 0.001).*sind(alpha)) + (eccentricity*0.001.*cosd(alpha)) );      % instantaneous power calculated by values of matcam


figure()
plot(P2)                                     % plotting total contact force with the same spring and double the rotation speed
title('Graph of instantaneous power 2')      

figure()
plot(P1-P2)
title('Graph of difference P1 and P2')

figure()
plot(P1-P1_2)
title('Graph of difference P1 and P1_2')

figure()
plot(P2-P2_2)
title('Graph of difference P2 and P2_2')


%% 3) Calculation average power usage

% Defining integration boundaries:
theta_min = 1;
theta_max = length(S);                      % Normally 361 values

P_tot1 = zeros(size(S));
P_tot2 = zeros(size(S));

i = theta_min;
while i < theta_max
    P_tot1(i+1) = P_tot1(i) + (P1(i) + P1(i+1))/2;
    P_tot2(i+1) = P_tot2(i) + (P2(i) + P2(i+1))/2;
    
    i = i+1;

end 

P_average1 = (sum(P_tot1))/(theta_max-theta_min)
P_average2 = (sum(P_tot2))/(theta_max-theta_min)



%% 4) Dynamics of a flexible follower

% Values critical rise/fall
beta_min = 330-265;                 % minimal angle of a rise in degrees
t_min = beta_min*(pi/180) / omega;  % minimal time of a rise in seconds

% Calculate minimal follower stiffness
k_follower = m_follower*(0.75*2*pi/(zeta*t_min))^2-k*1000;  % Spring constant of the follower in N/m
k_follower = k_follower/1000                                % Spring constant of the follower in N/mm

% Numerical simulation cam/follower system
omega_n = sqrt(1000*(k_follower+k)/m_follower);
t_n = 2*pi/omega_n;
lambda = t_min / t_n;

numerator = (2*pi*lambda)^2;
denominator = [1, 2*zeta*(2*pi*lambda), (2*pi*lambda)^2];
sys = tf(numerator, denominator);

T_s = 1/66;
tau = 0:T_s:10;
crit_rise_input = zeros(size(tau));
crit_rise_input(1:66) = S(266:331) / (max(S(266:331))-min(S(266:331)));

init_rise = 1;
init_vel = 0;
[A,B,C,D] = tf2ss(numerator, denominator);
X0 = [1/C(2)*init_vel; 1/C(2)*init_rise];
lsim(A,B,C,D, crit_rise_input, tau, X0);            % picture
gamma = lsim(A,B,C,D, crit_rise_input, tau, X0);    % saving data

% Approximate analysis
