clear;
close all;

% assigned values
follower_mass = 20;
follower_c = 0;
T_cycle = 0.5;
omega = 2*pi/T_cycle;
cam_rpm = 60/T_cycle;

matcam();
% loading the files with the motion law and load profile
matcam('genmotlawload')
matcam('genextloadload')

% setting the additional parameters for the analysis
global bcr_edit rof_edit exc_edit contourgrad_edit          
global mass_edit spring_edit sprload_edit rpm_edit          
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

figure()
plot(F_load)


% global beta0 beta1
% global lambda_plot lambda
% global roc_plot roc gamma

% global RRtotgamma VVtotgamma AAtotgamma
% global Rtotrad 
% global rof
% global extload extloadtot
% global betaload betaloadtot
% global betaload0 betaload1
% global x_pitch y_pitch
% global x_pitch_plot y_pitch_plot
% global x_contour y_contour
% global x_contour_plot y_contour_plot
% global x_rolfol y_rolfol
% global x_rod y_rod x_rod_ref y_rod_ref
% global x_center y_center
% global x_pressangle y_pressangle
% global x_refangle y_refangle
% global x_boreleft x_boreright y_bore
% global currentloadfile currentliftfile
% global bcr
% global rof