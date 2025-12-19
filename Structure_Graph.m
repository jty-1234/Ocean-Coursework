% ==============================================================
% plot_structure_CD.m
%
% Plot the whole structure geometry with:
%   - CD at z = 0
%   - H (seabed to CD) and h (CD to underside of deck) shown
% ==============================================================

clearvars -except h_water eta
clc; close all;

% ---- check needed inputs from previous script ----
if ~exist('h_water','var') || ~exist('eta','var')
    error('h_water and/or eta not found. Run the wind/wave script first.');
end

%% --- Basic geometric parameters (edit if needed) ---
H          = h_water;  % in the drawing: seabed -> CD
h          = 15;       % temporary: CD -> underside of deck

deck_width = 50;       % plan width (wind direction) [m]
deck_height = 22;      % vertical height of deck [m]

D_base     = 40;       % trial base diameter at seabed [m]  <-- CHANGE IF NEEDED
D_shaft    = 15;       % circular shaft diameter [m]
h_plinth   = 5;        % height of cylindrical plinth on seabed [m]
theta_cone = 60;       % cone angle to horizontal [deg]

%% --- Key z-levels in CD coordinates ---
z_CD      = 0;             % chart datum
z_seabed  = -H;            % seabed
z_water   = eta;           % still water level (HAT)
z_plinth_top = z_seabed + h_plinth;

% cone height from geometry
h_cone = ((D_base - D_shaft)/2) * tand(theta_cone);
z_shaft_base = z_plinth_top + h_cone;      % where shaft becomes 15 m
z_deck_bottom = h;                         % underside of deck
z_deck_top    = h + deck_height;          % top of deck

%% --- Horizontal dimensions (radii / half-widths) ---
R_base  = D_base/2;
R_shaft = D_shaft/2;
R_deck  = deck_width/2;

%% --- Build outline coordinates ---

% foundation (plinth + cone + shaft up to deck bottom)
z_outline = [ ...
    z_seabed;          % 1 bottom of base
    z_seabed;          % 2
    z_plinth_top;      % 3
    z_shaft_base;      % 4
    z_deck_bottom;     % 5
    z_deck_bottom;     % 6
    z_deck_top;        % 7
    z_deck_top];       % 8

x_outline = [ ...
   -R_base;            % left at seabed
    R_base;            % right at seabed
    R_base;            % right at plinth top
    R_shaft;           % right at shaft base
    R_shaft;           % right at deck bottom
   -R_shaft;           % left at deck bottom
   -R_deck;            % left at deck top
    R_deck];           % right at deck top

% To draw closed outline we’ll plot segments individually.

%% --- Plotting ---
figure; hold on; grid on;

% seabed line
plot([-R_deck-20, R_deck+20], [z_seabed, z_seabed], 'k--', 'LineWidth', 1.2);
text(0, z_seabed-2, 'Seabed', 'HorizontalAlignment','center');

% CD line
plot([-R_deck-20, R_deck+20], [z_CD, z_CD], 'k:', 'LineWidth', 1.0);
text(0, 0.8, 'Chart Datum (CD)', 'HorizontalAlignment','center');

% still water level
plot([-R_deck-20, R_deck+20], [z_water, z_water], 'c--', 'LineWidth', 1.2);
text(0, z_water+1, 'Still water level', 'HorizontalAlignment','center');

% --- draw base plinth ---
plot([-R_base, R_base], [z_seabed, z_seabed], 'b', 'LineWidth', 2);              % bottom
plot([-R_base, -R_base], [z_seabed, z_plinth_top], 'b', 'LineWidth', 2);        % left
plot([ R_base,  R_base], [z_seabed, z_plinth_top], 'b', 'LineWidth', 2);        % right

% --- draw conical part (straight lines from base radius to shaft radius) ---
plot([-R_base, -R_shaft], [z_plinth_top, z_shaft_base], 'b', 'LineWidth', 2);
plot([ R_base,  R_shaft], [z_plinth_top, z_shaft_base], 'b', 'LineWidth', 2);

% --- draw vertical shaft up to deck bottom ---
plot([-R_shaft, -R_shaft], [z_shaft_base, z_deck_bottom], 'b', 'LineWidth', 2);
plot([ R_shaft,  R_shaft], [z_shaft_base, z_deck_bottom], 'b', 'LineWidth', 2);

% --- draw deck (rectangle) ---
plot([-R_deck,  R_deck], [z_deck_bottom, z_deck_bottom], 'b', 'LineWidth', 2);   % bottom
plot([-R_deck,  R_deck], [z_deck_top,    z_deck_top],    'b', 'LineWidth', 2);   % top
plot([-R_deck, -R_shaft], [z_deck_bottom, z_deck_top],   'b', 'LineWidth', 2);   % left outer
plot([ R_deck,  R_shaft], [z_deck_bottom, z_deck_top],   'b', 'LineWidth', 2);   % right outer

%% --- Dimension lines for H and h ---
x_dim = R_deck + 10;   % x-position for dimension lines

% H: seabed to CD
plot([x_dim, x_dim], [z_seabed, z_CD], 'k-', 'LineWidth', 1);
text(x_dim+2, (z_seabed+z_CD)/2, 'H', 'HorizontalAlignment','left');

% h: CD to underside of deck
plot([x_dim, x_dim], [z_CD, z_deck_bottom], 'k-', 'LineWidth', 1);
text(x_dim+2, (z_CD+z_deck_bottom)/2, 'h', 'HorizontalAlignment','left');

xlabel('Horizontal coordinate (m)');
ylabel('z relative to Chart Datum (m)');
title('Structure Geometry (CD = 0)');

axis equal;
xlim([-R_deck-25, R_deck+25]);
ylim([z_seabed-5, z_deck_top+5]);
set(gca,'YDir','normal');
