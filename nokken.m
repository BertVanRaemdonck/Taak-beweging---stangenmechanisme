clear;
close all;

%% Declaration of variables and loading in matcam

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT
%
% To get this program to run smoothly, please fill in the location of the
% motion law and external load files on your computer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mot_law_location = 'E:\Data\KULeuven\3de bachelor\2de semester\Beweging\Taak nokken\Taak-beweging-stangenmechanisme\nok_hefwet.mot';
ext_load_location = 'E:\Data\KULeuven\3de bachelor\2de semester\Beweging\Taak nokken\Taak-beweging-stangenmechanisme\nok_externe_krachten.exl';

% assigned values
follower_mass = 20;
zeta = 0.077;
T_cycle = 0.5;
omega = 2*pi/T_cycle;
cam_rpm = 60/T_cycle;

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
set(bcr_edit, 'string', '80');
set(rof_edit, 'string', '20');
set(exc_edit, 'string', '0');                               % Analysis first without an eccentric follower
set(contourgrad_edit, 'string', '0');
set(mass_edit, 'string', num2str(follower_mass));
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
S = Stot;                                       % Lift
V = Vtot;                                       % Velocity
A = Atot;                                       % Acceleration
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
% test to see if the loading of the variables workes
%figure()
%plot(alpha)

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
    set(bcr_edit, 'string', '80');
    set(rof_edit, 'string', '20');
    set(exc_edit, 'string', num2str(eccentricity));             % Analysis first without an eccentric follower
    set(contourgrad_edit, 'string', '0');
    set(mass_edit, 'string', num2str(follower_mass));
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

    S = Stot;                                       % Lift
    V = Vtot;                                       % Velocity
    A = Atot;                                       % Acceleration

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

% Min teken moet erbij, maar doet dan ambetant
k = max((-F_load - F_v0*ones(size(S)) - follower_mass.*omega.*omega.*A)./S)               % - follower_mass*g*cos(gamma)*ones(size(S))  weggedaan omdat zwaartekracht te verwaarlozen is

k = ceil(k);                                    % Rounds k up to the next integer


if 1 == 1                                       % Recalibration matcam
    
    matcam()
    % setting the additional parameters for the analysis
    set(genmotlaw_startangle_edit, 'string', '330');            % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_endangle_edit, 'string', '360');              % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_startlift_edit, 'string', '0');               % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_endlift_edit, 'string', '0');                 % not relevant for calculations, just keep matcam from throwing errors
    set(bcr_edit, 'string', '80');
    set(rof_edit, 'string', '20');
    set(exc_edit, 'string', num2str(eccentricity));             % Analysis with an eccentric follower
    set(contourgrad_edit, 'string', '0');
    set(mass_edit, 'string', num2str(follower_mass));
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

    S = Stot;                                       % Lift
    V = Vtot;                                       % Velocity
    A = Atot;                                       % Acceleration

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
    set(bcr_edit, 'string', '80');
    set(rof_edit, 'string', '20');
    set(exc_edit, 'string', num2str(eccentricity));             % Analysis with an eccentric follower
    set(contourgrad_edit, 'string', '0');
    set(mass_edit, 'string', num2str(follower_mass));
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

    S = Stot;                                       % Lift
    V = Vtot;                                       % Velocity
    A = Atot;                                       % Acceleration

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
k_double = max((-F_load - F_v0_double*ones(size(S)) - follower_mass.*double_omega.*double_omega.*A )./S)   % - follower_mass*g*cos(gamma)*ones(size(S))  weggedaan omdat zwaartekracht te verwaarlozen is

k_double = ceil(k_double);

if 1 == 1                                       % Recalibration matcam
    
    matcam()
    % setting the additional parameters for the analysis
    set(genmotlaw_startangle_edit, 'string', '330');            % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_endangle_edit, 'string', '360');              % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_startlift_edit, 'string', '0');               % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_endlift_edit, 'string', '0');                 % not relevant for calculations, just keep matcam from throwing errors
    set(bcr_edit, 'string', '80');
    set(rof_edit, 'string', '20');
    set(exc_edit, 'string', num2str(eccentricity));             % Analysis with an eccentric follower
    set(contourgrad_edit, 'string', '0');
    set(mass_edit, 'string', num2str(follower_mass));
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

    S = Stot;                                       % Lift
    V = Vtot;                                       % Velocity
    A = Atot;                                       % Acceleration

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
    set(bcr_edit, 'string', '80');
    set(rof_edit, 'string', '20');
    set(exc_edit, 'string', '0');                               % Analysis without eccentricity
    set(contourgrad_edit, 'string', '0');
    set(mass_edit, 'string', num2str(follower_mass));
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

    S = Stot;                                       % Lift
    V = Vtot;                                       % Velocity
    A = Atot;                                       % Acceleration

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
    
    global Rtotrad
    R_tot = Rtotrad;                                % Total distance between center of rotation and center of follower
    
    % Close matcam window
    matcam_figure = gcf;
    close(matcam_figure.Number)
end 

P1 = F_tot .*(omega*ones(size(S))) .*(sin(alpha)) .* R_tot;      % instantaneous power

if 1 == 1                                       % Recalibration matcam
    
    matcam()
    % setting the additional parameters for the analysis
    set(genmotlaw_startangle_edit, 'string', '330');            % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_endangle_edit, 'string', '360');              % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_startlift_edit, 'string', '0');               % not relevant for calculations, just keep matcam from throwing errors
    set(genmotlaw_endlift_edit, 'string', '0');                 % not relevant for calculations, just keep matcam from throwing errors
    set(bcr_edit, 'string', '80');
    set(rof_edit, 'string', '20');
    set(exc_edit, 'string', num2str(eccentricity));             % Analysis with eccentricity
    set(contourgrad_edit, 'string', '0');
    set(mass_edit, 'string', num2str(follower_mass));
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

    S = Stot;                                       % Lift
    V = Vtot;                                       % Velocity
    A = Atot;                                       % Acceleration

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
    
    global Rtotrad
    R_tot = Rtotrad;                                % Total distance between center of rotation and center of follower
    
    % Close matcam window
    matcam_figure = gcf;
    close(matcam_figure.Number)
end 

P2 = F_tot .*(omega*ones(size(S))) .*( ( ( ((R_tot.^2)-((eccentricity*ones(size(S))).^2)).^(-1/2)).*sin(alpha)) + (eccentricity.*cos(alpha)) ); 
% Niet hetzelfde, dus afgeleide formule werkt niet...

figure()
plot(P1)                                     % plotting total contact force with the same spring and double the rotation speed
title('Graph of instantaneous power 1')

figure()
plot(P2)                                     % plotting total contact force with the same spring and double the rotation speed
title('Graph of instantaneous power 2')




%% 4) Dynamics of a flexible follower

beta_min = 330-265;                 % minimal angle of a rise in degrees
t_min = beta_min*(pi/180) / omega;  % minimal time of a rise in seconds