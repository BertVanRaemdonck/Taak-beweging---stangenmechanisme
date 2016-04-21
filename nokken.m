clear;
close all;

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
follower_c = 0;
T_cycle = 0.5;
omega = 2*pi/T_cycle;
cam_rpm = 60/T_cycle;

matcam();
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

matcam('genmotlawcalc')

%% inladen matcam variabelen

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

close all;
figure()
plot(alpha)