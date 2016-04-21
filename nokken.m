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
ext_load_location = 'C:\Users\Bert\School\Beweging en trillingen\Code\nok_externe_krachten.exl';

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
set(exc_edit, 'string', '0');
set(contourgrad_edit, 'string', '0');
set(mass_edit, 'string', num2str(follower_mass));
set(spring_edit, 'string', '2');
set(sprload_edit, 'string', '10');
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

% test to see if the loading of the variables workes
figure()
plot(alpha)

%% 4) Dynamics of a flexible follower

beta_min = 330-265;                 % minimal angle of a rise in degrees
t_min = beta_min*(pi/180) / omega;  % minimal time of a rise in seconds