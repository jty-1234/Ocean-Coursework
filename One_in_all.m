%% --------------------------------------------------------------
%   Wind (Frøya) + Wave Particle Velocity Profiles (Chart Datum)
%   - Input: U_extreme (10-min extreme wind at z_CD_meas above CD)
%            Wave height H, period T, water depth h
%   - Output: Combined plot of wind profile + wave particle velocity
%             with vertical axis referenced to Chart Datum (CD)
%             and "normal" orientation (no YDir reverse).
% --------------------------------------------------------------

clear; clc; close all;

%% --------------------------------------------------------------
%  USER INPUTS
% --------------------------------------------------------------

% --- Wind inputs ---
U_extreme  = 22.82;     % extreme 10-min wind speed at measurement height (m/s)
z_CD_meas  = 20;     % measurement height above Chart Datum (m)

MSL_above_CD = 0;  % mean sea level above Chart Datum (m)
Tide         = 7.37;  % tide elevation above MSL at the time (m)


z_min_sea_wind = 0;      % minimum height above sea surface (m)
z_max_sea_wind = 50;    % maximum height above sea surface (m)
Nz_wind        = 200;


H_wave = 4.7;        % wave height H (m)
T_wave = 8.5;        % wave period T (s)
h_water = 40;       % water depth below still water level (m) from digimap
Nz_wave = 200;      % number of vertical points for wave profile

g = 9.81;           % gravity (m/s^2)

%% --------------------------------------------------------------
%  STEP 0 – Sea surface elevation relative to CD
% --------------------------------------------------------------
% eta = sea surface elevation above Chart Datum (can be positive or negative)
eta = MSL_above_CD + Tide;      

% Convert measurement height from CD -> height above sea surface
z_meas_sea = z_CD_meas - eta;   

fprintf('Sea surface elevation above CD (eta)         = %.3f m\n', eta);
fprintf('Measurement height above sea surface (z_meas) = %.3f m\n', z_meas_sea);

if z_meas_sea <= 0
    warning('Measurement point is at or below sea surface (z_meas = %.2f m).', z_meas_sea);
end

H  = 10;          % reference height for Frøya (m)
T0 = 1.0;         % hours
T_meas = 10/60;   % 10 minutes in hours

U_T_meas = U_extreme;    % your 10-min mean at z_meas_sea

% Initial guess for U0 = 1-hour mean at 10 m
U0 = U_T_meas;

for i = 1:100
    % Coefficients based on current U0
    C  = 5.73e-2 * sqrt(1 + 0.148 * U0);
    Iv = 0.06 * (1 + 0.043 * U0) * (z_meas_sea / H)^(-0.22);

    % Frøya prediction at (T_meas, z_meas)
    factor = (1 + C * log(z_meas_sea / H)) * ...
             (1 - 0.41 * Iv * log(T_meas / T0));

    % Fixed-point update for U0
    U0_new = U_T_meas / factor;

    if abs(U0_new - U0) < 1e-6
        break;
    end

    U0 = U0_new;
end
fprintf('1-hour mean wind speed at 10 m above sea level (U0) = %.3f m/s\n', U0);

% Height of that 10 m point relative to Chart Datum
z10_CD = H + eta;
%% --------------------------------------------------------------
%  STEP 3 – Generate Frøya wind profile (above sea surface)
% --------------------------------------------------------------
z_sea_wind = linspace(z_min_sea_wind+0.1, z_max_sea_wind, Nz_wind);  % heights above sea surface

C = 0.0573 * sqrt(1 + 0.148 * U0);
U_wind =U0 .* (1 + C .* log(z_sea_wind ./ (z_CD_meas-Tide)))*(1 - 0.41 * Iv * log(T_meas / T0));

% Convert wind profile heights to Chart Datum reference
z_CD_wind = z_sea_wind + eta;  % height above CD


% The value used for wave prediction 
u_wind_10m = U0 .* (1 + C .* log(10 ./ (z_CD_meas-Tide)))*(1 - 0.41 * Iv * log(T_meas / T0))
%% --------------------------------------------------------------
%  STEP 4 – Wave particle horizontal velocity profile (linear theory)
% --------------------------------------------------------------
% z_sea_wave: from seabed (-h) to still water level (0), relative to SWL
z_sea_wave = linspace(-h_water-eta, 0, Nz_wave);
d= h_water+eta;
% Solve dispersion relation: omega^2 = g k tanh(k h)
L = g*T_wave^2/(2*pi);   % initial guess
for i = 1:100
    L = g*T_wave^2 / (2*pi)*tanh(2*pi*d/L);
end

a = H_wave / 2;    % wave amplitude

% Maximum horizontal particle velocity:

u_wave = a * g*T_wave/L.* cosh(2*pi * (z_sea_wave + d)/L) ./ cosh(2*pi*d/L);

% Particle acceleration

a_wave = g*pi*H_wave/L.* cosh(2*pi * (z_sea_wave + d)/L) ./ cosh(2*pi*d/L)

% Convert wave depths to Chart Datum reference
z_CD_wave = z_sea_wave + eta;   % (z relative SWL) + eta -> relative to CD
%  STEP 5 – PLOT wind induced profile

vc_wind_0= U0*0.03
z_wind_current = linspace(-d, 0, Nz_wave);
vc_wind = vc_wind_0*((50+z_wind_current)/50)
z_CD_wind_current = z_wind_current + eta

%  STEP 6 – PLOT tide current profile

vc_tide_0 = 1.5/((d-18)/d)^(1/7)
z_tide_current = linspace(-d, 0, Nz_wave);
vc_tide = vc_tide_0*(((d+z_tide_current)/d).^(1/7))
z_CD_tide_current = z_tide_current + eta
%% --------------------------------------------------------------
%  STEP 7 – PLOT combined profiles 
% --------------------------------------------------------------
figure;
hold on;
     % wind profile
plot(u_wave, z_CD_wave, 'LineWidth', 2);       % wave particle velocity profile
plot(vc_wind, z_CD_wind_current, 'LineWidth', 2);
plot(vc_tide, z_CD_tide_current, 'LineWidth', 2);
plot(a_wave, z_CD_wave, 'LineWidth', 2);
hold off;

xlabel('Speed (m/s)');
ylabel('Vertical coordinate relative to Chart Datum (m)');
title('Velocity Profile (Chart Datum reference)');

% IMPORTANT: no YDir reverse here – heights increase upward,
% so this is "upside down" compared with the earlier reversed plots.
% (CD = 0; sea surface ≈ eta; above that wind, below that wave motion.)

grid on;
legend({'Maximum Wave particle velocity', 'Wind induced current veloctiy','tide current velocity', 'Maximum wave particle acceleration'}, ...
       'Location', 'best');
set(gca, 'FontSize', 9);
grid minor;
Total_current = vc_tide + vc_wind;
figure;
plot(Total_current, z_CD_tide_current, 'LineWidth', 2);
xlabel('Speed (m/s)');
ylabel('Vertical coordinate relative to Chart Datum (m)');
title('Total Current Profile (Chart Datum reference)');
grid on;
grid minor;
figure;
hold on
plot(U_wind, z_CD_wind, 'LineWidth', 2);
xlabel('Speed (m/s)');
ylabel('Vertical coordinate relative to Chart Datum (m)');
title('Wind velocity Profile (Chart Datum reference)');
grid on;
grid minor;
%% --------------------------------------------------------------
%  STEP 8 – Hydrodynamic forces with Morison equation
%           (wave + current, with varying diameter)
% --------------------------------------------------------------

% --- Hydrodynamic parameters (edit as needed) ---
rho = 1025;     % water density (kg/m^3)
CM  = 2.0;      % inertia coefficient  (typical order)
CD  = 0.7;      % drag coefficient     (typical order)

% --- Structural geometry ---
D_base   = 40;       % [m] base diameter at seabed (your unknown D – set value)
D_shaft  = 15;       % [m] shaft diameter
H_struct = h_water;  % [m] seabed to CD (H in sketch)
h_plinth = 5;        % [m] height of constant-D base block
theta_cone = 60;     % [deg] cone angle to horizontal

z_seabed     = -H_struct;                     % seabed level in CD coordinates
z_plinth_top = z_seabed + h_plinth;           % top of plinth (start of cone)

% vertical height of cone from geometry
h_cone = ((D_base - D_shaft)/2) * tand(theta_cone);

z_shaft_base = z_plinth_top + h_cone;         % level where diameter becomes 15 m

% --- Build diameter profile D(z) along the wave grid (z_CD_wave) ---
D_profile = zeros(size(z_CD_wave));

for i = 1:length(z_CD_wave)
    z = z_CD_wave(i);

    if z <= z_plinth_top
        % region of constant base diameter
        D_profile(i) = D_base;

    elseif z <= z_shaft_base
        % conical transition: linear interpolation from D_base to D_shaft
        z_rel = z - z_plinth_top;  % vertical distance measured from top of plinth
        D_profile(i) = D_base - (D_base - D_shaft) * (z_rel / h_cone);

    else
        % vertical shaft of constant diameter 15 m
        D_profile(i) = D_shaft;
    end
end

%% --- Morison equation components (per unit length) ---

% Relative horizontal velocity = (wave orbital + mean current)
u_rel = u_wave + Total_current;         % [m/s], same vertical grid as z_CD_wave

% Inertia term (uses acceleration a_wave)
F_I = rho * (pi/4) * CM .* (D_profile.^2) .* a_wave;              % [N/m]

% Drag term (uses relative velocity)
F_D = 0.5 * rho * CD .* D_profile .* u_rel .* abs(u_rel);         % [N/m]

% Total horizontal force per unit length
Fz_total = F_I + F_D;                                             % [N/m]

%% --- Resultant force and overturning moment ---

% Net horizontal force (integrated from seabed to still water level)
F_total = trapz(z_CD_wave, Fz_total);          % [N]

% Moments about seabed (z = z_seabed) and about CD (z = 0)
lever_seabed = z_CD_wave - z_seabed;           % [m]
lever_CD      = z_CD_wave;                     % [m]

M_seabed = trapz(z_CD_wave, Fz_total .* lever_seabed);   % [N·m]
M_CD     = trapz(z_CD_wave, Fz_total .* lever_CD);       % [N·m]

fprintf('Resultant horizontal Morison force (seabed–SWL): %.1f kN\n', F_total/1e3);
fprintf('Overturning moment about seabed:                 %.1f MN·m\n', M_seabed/1e6);
fprintf('Overturning moment about Chart Datum (z=0):      %.1f MN·m\n', M_CD/1e6);

%% --- Plot force distribution vs depth (optional) ---
figure;
plot(F_I, z_CD_wave, 'LineWidth', 1.5); hold on;
plot(F_D, z_CD_wave, 'LineWidth', 1.5);
plot(Fz_total, z_CD_wave, 'LineWidth', 2);
yline(0,'k--');      % CD
yline(eta,'k:');     % still water level
xlabel('Force per unit length (N/m)');
ylabel('z (m) relative to Chart Datum');
title('Morison force distribution along structure');
legend('Inertia term','Drag term','Total','Location','best');
grid on; grid minor;

%% --------------------------------------------------------------
%  STEP 9 – Resultant force & moment for discrete phases (every 30°)
% --------------------------------------------------------------

% We use linear wave theory:
% u(z,phi) = u_wave(z) * cos(phi)
% a(z,phi) = -a_wave(z) * sin(phi)
% where u_wave and a_wave are the amplitudes you already computed.

phases_deg = 0:30:330;               % 0, 30, ..., 330 degrees
phases_rad = deg2rad(phases_deg);
Nphi       = numel(phases_rad);

F_total_phase   = zeros(1, Nphi);    % resultant horizontal force [N]
M_seabed_phase  = zeros(1, Nphi);    % overturning moment about seabed [N·m]
M_CD_phase      = zeros(1, Nphi);    % moment about Chart Datum z=0 [N·m]

lever_seabed = z_CD_wave - z_seabed; % lever arm to seabed
lever_CD     = z_CD_wave;            % lever arm to CD

for iphi = 1:Nphi
    phi = phases_rad(iphi);

    % wave orbital velocity & acceleration at this phase
    u_wave_phi =  u_wave .* cos(phi);
    a_wave_phi = -a_wave .* sin(phi);  % sign consistent with du/dt

    % add steady current (wind + tide)
    u_rel = u_wave_phi + Total_current;

    % Morison terms per unit length
    F_I = rho * (pi/4) * CM .* (D_profile.^2) .* a_wave_phi;           % inertia
    F_D = 0.5 * rho * CD .* D_profile .* u_rel .* abs(u_rel);          % drag
    Fz  = F_I + F_D;                                                   % total

    % Integrate over depth
    F_total_phase(iphi)  = trapz(z_CD_wave, Fz);                       % [N]
    M_seabed_phase(iphi) = trapz(z_CD_wave, Fz .* lever_seabed);       % [N·m]
    M_CD_phase(iphi)     = trapz(z_CD_wave, Fz .* lever_CD);           % [N·m]
end

% Print a small table (force in kN, moment in MNm)
fprintf('\nPhase (deg)   F_total (kN)   M_seabed (MN·m)   M_CD (MN·m)\n');
for iphi = 1:Nphi
    fprintf('%4.0f         %10.2f      %10.3f        %10.3f\n', ...
        phases_deg(iphi), ...
        F_total_phase(iphi)/1e3, ...
        M_seabed_phase(iphi)/1e6, ...
        M_CD_phase(iphi)/1e6);
end

% Optional plots
figure;
subplot(2,1,1);
plot(phases_deg, F_total_phase/1e3, '-o', 'LineWidth', 1.5);
xlabel('Wave phase (deg)');
ylabel('Resultant force (kN)');
title('Total horizontal Morison force vs phase');
grid on;

subplot(2,1,2);
plot(phases_deg, M_seabed_phase/1e6, '-o', 'LineWidth', 1.5);
xlabel('Wave phase (deg)');
ylabel('Overturning moment about seabed (MN·m)');
title('Overturning moment vs phase');
grid on;

%% --------------------------------------------------------------
%  STEP 10 – Extract maximum total force and moment over all phases
% --------------------------------------------------------------

% Use absolute values (worst-case magnitude)
[ F_max_abs, idx_Fmax ]   = max(abs(F_total_phase));
[ M_max_abs, idx_Mmax ]   = max(abs(M_seabed_phase));
[ MCD_max_abs, idx_MCDmax ] = max(abs(M_CD_phase));

phase_Fmax   = phases_deg(idx_Fmax);
phase_Mmax   = phases_deg(idx_Mmax);
phase_MCDmax = phases_deg(idx_MCDmax);

fprintf('\n=== PEAK VALUES OVER ALL PHASES (30° steps) ===\n');
fprintf('Max |F_total|       = %.2f kN   at phase = %3.0f°\n', F_max_abs/1e3,  phase_Fmax);
fprintf('Max |M_seabed|      = %.3f MN·m at phase = %3.0f°\n', M_max_abs/1e6,  phase_Mmax);
fprintf('Max |M_CD|          = %.3f MN·m at phase = %3.0f°\n\n', MCD_max_abs/1e6, phase_MCDmax);












%% ==============================================================
% wind_load_from_profile.m  (UNIFIED VARIABLE NAMES)
%
% Requires existing workspace variables from your main script:
%   z_CD_wind   U_wind   eta   h_water
%
% Outputs:
%   F_wind_total  [N]
%   M_wind_total  [N·m] about seabed
% ==============================================================
%% -------- CHECK INPUTS ---------------------------------------
req = {'z_CD_wind','U_wind','eta','h_water'};
for k = 1:numel(req)
    if ~exist(req{k},'var')
        error('%s is missing. Run wind/wave profile script first.', req{k});
    end
end
if ~exist('g','var');  g = 9.81; end

rho_air = 1.226;
q_profile = 0.5 * rho_air .* U_wind.^2;

%% -------- STRUCTURE GEOMETRY (CD reference) ------------------
h        = 15;      % underside of deck above CD (temporary)
H_deck   = 22;      % deck vertical height
B_deck   = 50;      % deck width in wind direction

D_shaft  = 15;      % shaft diameter
n_shaft  = 2;

z_seabed = -h_water;
z_water  = eta;

z_deck_bottom = h;
z_deck_top    = h + H_deck;

%% -------- COEFFICIENTS ---------------------------------------
alpha = deg2rad(90);       % wind normal
C_rect = 2.0;              % Cs1 per your spec
C_cyl  = 1.0;              % assumption given by you

%% =============================================================
%   WIND ON RECTANGULAR DECK (integrated with profile)
% =============================================================
mask_rect = (z_CD_wind >= z_deck_bottom) & (z_CD_wind <= z_deck_top);

z_rect = z_CD_wind(mask_rect);
q_rect = q_profile(mask_rect);

width_rect = B_deck;

dF_rect = C_rect .* q_rect .* width_rect;      % N/m
F_wind_rect = trapz(z_rect, dF_rect);          % N

lever_rect = trapz(z_rect, dF_rect .* (z_rect - z_seabed)) / F_wind_rect;
M_wind_rect = F_wind_rect * lever_rect;

%% =============================================================
%   WIND ON CYLINDRICAL SHAFT
%   (between still water level eta and underside of deck h)
% =============================================================
z_exposed_low  = z_water;        % = eta
z_exposed_high = z_deck_bottom;  % = h

% logical mask selecting the heights where wind hits the shaft
mask_cyl = (z_CD_wind >= z_exposed_low) & (z_CD_wind <= z_exposed_high);

if any(mask_cyl)
    z_cyl = z_CD_wind(mask_cyl);
    q_cyl = q_profile(mask_cyl);

    width_cyl = D_shaft * n_shaft;       % 2 shafts each D_shaft wide
    dF_cyl    = C_cyl .* q_cyl .* width_cyl;   % N/m

    % Resultant force on cylinder (integral of distributed load)
    F_wind_cyl = trapz(z_cyl, dF_cyl);                % N

    % Overturning moment about seabed (no division → no Inf)
    M_wind_cyl = trapz(z_cyl, dF_cyl .* (z_cyl - z_seabed));  % N·m
else
    % No exposed shaft in this water level range
    F_wind_cyl = 0;
    M_wind_cyl = 0;
end

%% =============================================================
%   TOTAL WIND
% =============================================================
F_wind_total = F_wind_rect + F_wind_cyl;
M_wind_total = M_wind_rect + M_wind_cyl;

fprintf("\n===== WIND LOAD (Profile Integrated) =====\n");
fprintf("Rectangular block force  = %.1f kN\n", F_wind_rect/1e3);
fprintf("Cylindrical shaft force  = %.1f kN\n", F_wind_cyl/1e3);
fprintf("TOTAL wind force         = %.1f kN\n\n", F_wind_total/1e3);

fprintf("Rectangular moment       = %.2f MNm\n", M_wind_rect/1e6);
fprintf("Cylindrical moment       = %.2f MNm\n", M_wind_cyl/1e6);
fprintf("TOTAL wind moment        = %.2f MNm\n\n", M_wind_total/1e6);

%% -------- OPTIONAL PLOT ---------------------------------------
figure;
plot(q_profile, z_CD_wind, 'LineWidth', 2);
set(gca,'YDir','normal');
xlabel('Dynamic pressure q(z) [N/m^2]');
ylabel('Height above Chart Datum (m)');
title('Wind pressure profile at HAT');
grid on;

M_Total = M_max_abs/1e6 + M_wind_total/1e6
fprintf('Maximum Total Moment for wind wave combined   = %.3f MN·m at phase = %3.0f°\n\n', M_Total/1e6, phase_MCDmax);







