% ==============================================================
% wind_load_from_profile.m  (standalone)
% ==============================================================

clear; clc; close all;

% 1) Load environment data saved from main script
projectFile = 'env_data.mat';
if ~isfile(projectFile)
    error('File %s not found. Run the main wind–wave script first.', projectFile);
end

load(projectFile, 'z_CD_wind', 'U_wind', 'eta', 'h_water', 'g');



% ==============================================================
% wind_load_from_profile.m  (UNIFIED VARIABLE NAMES)
%
% Requires existing workspace variables from your main script:
%   z_CD_wind   U_wind   eta   h_water
%
% Outputs:
%   F_wind_total  [N]
%   M_wind_total  [N·m] about seabed
% ==============================================================

clearvars -except z_CD_wind U_wind eta h_water g
clc;

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

M_Total = M_max_abs/1e6 + M_wind_total
fprintf('Maximum Total Moment for wind wave combined   = %.3f MN·m at phase = %3.0f°\n\n', MM_Total/1e6, phase_MCDmax);


