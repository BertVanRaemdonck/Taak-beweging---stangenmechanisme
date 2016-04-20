matcam('genmotlawload');

%% inladen matcam variabelen

global Stot Vtot Atot;    
S = Stot;                                       % Lift
V = Vtot;                                       % Velocity
A = Atot;                                       % Acceleration
global tetatot tetatotrad rpm;
theta = tetatot;                                % Cam angle in degrees
theta_rad = tetatotrad;                         % Cam angle in radians
rpm_cam = rpm;                                  % Angular velocity of the cam
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