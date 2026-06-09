% Clear workspace and close figures
clear; clc; close all;

% Given Data for Power Scan
erosion_rate_power = [9.31345e17, 1.1035e18, 1.537e18, 1.63755e18, 1.6217e18]; % Erosion rate (atoms/s)
power_levels = [36, 72, 360, 1800, 3600]; % Power levels (kW)

% Density Scan (Fixed 360kW power)
erosion_rate_density = [1.9e17, 1.537e18]; % Erosion rates for 0cm (Ref Case) and 5cm Far
density_labels = {'2.82e15', '2.38e16'}; % Labels for density scan locations
density_positions = [2.82e15, 2.38e16]; % Place 0cm Far at the top (2), 5cm Far at the bottom (1)

% Create Figure
figure;

% Left Y-Axis: Power vs. Erosion Rate
yyaxis left;
semilogx(erosion_rate_power, power_levels, 'bo-', 'LineWidth', 2, 'MarkerSize', 10);
ylabel('Power (kW)', 'FontSize', 18);
xlabel('Erosion Rate (atoms/s)', 'FontSize', 18);
ylim([0, 4200]); % Extra space for better labeling
grid on;
hold on;

% Annotate Power Scan Points with Extra Space for Labels
annotations_power = {'Low 36kW', 'Low 72kW', 'Base Case 360kW', 'Low 1800kW', 'High 3600kW'};
for i = 1:length(erosion_rate_power)
    text(erosion_rate_power(i), power_levels(i) * 1.08, annotations_power{i}, ...
        'FontSize', 14, 'HorizontalAlignment', 'center');
end

% Right Y-Axis: Density Scan vs. Erosion Rate (0cm Far on Top, 5cm Far on Bottom)
yyaxis right;
semilogx(erosion_rate_density, density_positions, 'rs--', 'LineWidth', 2, 'MarkerSize', 10);
yticks([2.82E15, 2.38e16]);
yticklabels(density_labels);
ylabel('Density Scan Locations (at 360kW)', 'FontSize', 18);

% Annotate Density Scan Points with Extra Space
for i = 1:length(erosion_rate_density)
    text(erosion_rate_density(i), density_positions(i), sprintf('%.2e', erosion_rate_density(i)), ...
        'FontSize', 14, 'Color', 'r', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
end

% Adjust x-axis limits to include both scans with extra space
xlim([0.5*min([erosion_rate_power, erosion_rate_density]) * 0.8, max([erosion_rate_power, erosion_rate_density]) * 1.2]);

% Set x-axis to log scale
set(gca, 'XScale', 'log');

% Legend
legend({'Power vs. Erosion Rate', 'Density vs. Erosion Rate (at 360kW)'}, 'Location', 'northwest');

% Title
title('Power and Density Scans vs. Erosion Rate (Log Scale)', 'FontSize', 20);

hold off;