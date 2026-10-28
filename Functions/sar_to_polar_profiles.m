function [V_all, r_grid, stats, coverage_all] = ...
    sar_to_polar_profiles(sar_lat, sar_lon, wsp, center_lat, center_lon, ...
                          dr_km, r_max_km, daz_deg, plot_type)
% SAR_TO_POLAR_PROFILES  Polar conversion + sector profiles + quadrant stats + coverage
%
% Outputs:
%   V_all        (9 × n_r) : Row 1 = azimuthal mean
%                           Rows 2–9 = 8 sectors (0–45, 45–90, ..., 315–360)
%   r_grid       (1 × n_r) : Radial axis (km)
%   stats                  : Vmax, Rmax, R34, R50, R64 per QUADRANT
%   coverage_all (9 × n_r) : Coverage (%) — same layout as V_all
%                            Row 1     = azimuthal mean coverage
%                            Rows 2–9  = 8 sectors coverage

    %% --- Defaults ---
    if nargin < 6  || isempty(dr_km),     dr_km     = 0.5;    end
    if nargin < 7  || isempty(r_max_km),  r_max_km  = 300;    end
    if nargin < 8  || isempty(daz_deg),   daz_deg   = 1.0;    end
    if nargin < 9  || isempty(plot_type), plot_type = 'none'; end

    %% --- Convert to double ---
    sar_lat = double(sar_lat); sar_lon = double(sar_lon);
    wsp = double(wsp);
    center_lat = double(center_lat); center_lon = double(center_lon);

    %% ============================================================
    %  STEP 1: Distance & azimuth (0° = North)
    %% ============================================================
    R_earth = 6371;
    lat0 = deg2rad(center_lat); lon0 = deg2rad(center_lon);
    latr = deg2rad(sar_lat);    lonr = deg2rad(sar_lon);
    dlon = lonr - lon0;

    cos_c = sin(lat0).*sin(latr) + cos(lat0).*cos(latr).*cos(dlon);
    c = acos(min(max(cos_c, -1), 1));
    dist_km = R_earth * c;

    az_rad = atan2(cos(latr).*sin(dlon), ...
                   cos(lat0).*sin(latr) - sin(lat0).*cos(latr).*cos(dlon));
    az_deg = mod(rad2deg(az_rad), 360);

    %% ============================================================
    %  STEP 2: Build polar mesh
    %% ============================================================
    r_grid = 0:dr_km:r_max_km;
    n_r = length(r_grid);

    az_grid = linspace(-pi, pi, round(360/daz_deg)+1);
    az_grid(end) = [];
    n_az = length(az_grid);

    [R_polar, AZ_polar] = meshgrid(r_grid, az_grid);
    R_polar  = R_polar'; AZ_polar = AZ_polar';

    fprintf('Polar mesh: %d × %d\n', n_r, n_az);

    %% ============================================================
    %  STEP 3: Polar mesh → Cartesian
    %% ============================================================
    X_polar = R_polar .* sin(AZ_polar);
    Y_polar = R_polar .* cos(AZ_polar);

    %% ============================================================
    %  STEP 4: SAR pixels → Cartesian
    %% ============================================================
    X_sar = dist_km(:) .* sin(az_rad(:));
    Y_sar = dist_km(:) .* cos(az_rad(:));
    V_sar = wsp(:);

    valid = ~isnan(X_sar) & ~isnan(Y_sar) & ~isnan(V_sar);
    X_sar = X_sar(valid); Y_sar = Y_sar(valid); V_sar = V_sar(valid);

    fprintf('Valid SAR pixels: %d\n', length(V_sar));

    %% ============================================================
    %  STEP 5: Interpolate
    %% ============================================================
    F = scatteredInterpolant(X_sar, Y_sar, V_sar, 'natural', 'none');
    V_polar = F(X_polar, Y_polar);

    %% ============================================================
    %  STEP 5b: Distance mask
    %% ============================================================
    kdt = KDTreeSearcher([X_sar, Y_sar]);
    [~, dist_to_nearest] = knnsearch(kdt, [X_polar(:), Y_polar(:)]);
    dist_to_nearest = reshape(dist_to_nearest, size(V_polar));
    V_polar(dist_to_nearest > 2.5) = NaN;

    fprintf('Polar points filled: %.1f%%\n', ...
            100 * sum(~isnan(V_polar(:))) / numel(V_polar));

    %% ============================================================
    %  STEP 5c: Build binary mask (1 = valid, 0 = NaN)
    %% ============================================================
    valid_mask = ~isnan(V_polar);   % n_r × n_az logical

    %% ============================================================
    %  STEP 6: Build V_all (9 × n_r) + coverage_all (9 × n_r)
    %% ============================================================
    V_all        = nan(9, n_r);
    coverage_all = nan(9, n_r);
    az_grid_deg  = mod(rad2deg(az_grid), 360);

    % ---- Row 1: Azimuthal mean ----
    V_all(1, :)        = mean(V_polar, 2, 'omitnan');
    coverage_all(1, :) = 100 * mean(valid_mask, 2);

    % ---- Rows 2–9: 8 sectors (45° each) ----
    sector_edges_deg = [0, 45, 90, 135, 180, 225, 270, 315, 360];
    for k = 1:8
        in_sector = (az_grid_deg >= sector_edges_deg(k)) & ...
                    (az_grid_deg <  sector_edges_deg(k+1));
        idx_k = find(in_sector);
        if ~isempty(idx_k)
            V_all(1+k, :)        = mean(V_polar(:, idx_k), 2, 'omitnan');
            coverage_all(1+k, :) = 100 * mean(valid_mask(:, idx_k), 2);
        end
    end

    %% ============================================================
    %  STEP 7: Build separate QUADRANT profiles for stats
    %% ============================================================
    V_quad = nan(4, n_r);   % NE, SE, SW, NW
    quad_edges_deg = [0, 90, 180, 270, 360];

    for q = 1:4
        in_quad = (az_grid_deg >= quad_edges_deg(q)) & ...
                  (az_grid_deg <  quad_edges_deg(q+1));
        idx_q = find(in_quad);
        if ~isempty(idx_q)
            V_quad(q, :) = mean(V_polar(:, idx_q), 2, 'omitnan');
        end
    end

    %% ============================================================
    %  STEP 8: Vmax/Rmax for quadrants
    %% ============================================================
    quad_labels = {'NE', 'SE', 'SW', 'NW'};

    stats.labels = quad_labels;
    stats.Vmax   = nan(4, 1);
    stats.Rmax   = nan(4, 1);

    for q = 1:4
        [Vmax_q, idx_q] = max(V_quad(q, :));
        if ~isnan(Vmax_q)
            stats.Vmax(q) = Vmax_q;
            stats.Rmax(q) = r_grid(idx_q);
        end
    end

    %% ============================================================
    %  STEP 9: R34/R50/R64 for quadrants
    %% ============================================================
    kt_to_ms = 0.5144;
    V34 = 34 * kt_to_ms;
    V50 = 50 * kt_to_ms;
    V64 = 64 * kt_to_ms;

    stats.R34 = nan(4, 1);
    stats.R50 = nan(4, 1);
    stats.R64 = nan(4, 1);

    for q = 1:4
        V_prof = V_quad(q, :);
        Rmax_q = stats.Rmax(q);

        if isnan(Rmax_q) || all(isnan(V_prof))
            continue;
        end

        stats.R34(q) = find_first_crossing_after_rmax(r_grid, V_prof, V34, Rmax_q);
        stats.R50(q) = find_first_crossing_after_rmax(r_grid, V_prof, V50, Rmax_q);
        stats.R64(q) = find_first_crossing_after_rmax(r_grid, V_prof, V64, Rmax_q);
    end

    %% ============================================================
    %  STEP 10: Report
    %% ============================================================
    fprintf('\n');
    fprintf('=====================================================================\n');
    fprintf('  SECTOR PROFILES (V_all: 9 rows)\n');
    fprintf('=====================================================================\n');
    fprintf('  Row 1     = Azimuthal mean\n');
    fprintf('  Rows 2–9  = Sectors (0–45, 45–90, ..., 315–360)\n');

    fprintf('\n');
    fprintf('=====================================================================\n');
    fprintf('  QUADRANT STATS (stats: 4 rows)\n');
    fprintf('=====================================================================\n');
    fprintf('%-8s | %-8s | %-8s | %-8s | %-8s | %-8s\n', ...
            'Quad', 'Vmax', 'Rmax', 'R34', 'R50', 'R64');
    fprintf('%-8s-|-%8s-|-%8s-|-%8s-|-%8s-|-%8s\n', ...
            '--------', '--------', '--------', '--------', '--------', '--------');
    for q = 1:4
        fprintf('%-8s | %8s | %8s | %8s | %8s | %8s\n', ...
                quad_labels{q}, ...
                fmt_num(stats.Vmax(q)), ...
                fmt_num(stats.Rmax(q)), ...
                fmt_num(stats.R34(q)), ...
                fmt_num(stats.R50(q)), ...
                fmt_num(stats.R64(q)));
    end
    fprintf('=====================================================================\n');

    % ---- Coverage summary ----
    fprintf('\n');
    fprintf('=====================================================================\n');
    fprintf('  COVERAGE (coverage_all: 9 rows)\n');
    fprintf('=====================================================================\n');
    fprintf('  Row 1     = Azimuthal mean coverage (%%)\n');
    fprintf('  Rows 2–9  = 8 sector coverage (%%)\n');
    fprintf('  Mean coverage (row 1): %.1f%%\n', ...
            mean(coverage_all(1, :), 'omitnan'));

    %% ============================================================
    %  STEP 11: PLOT
    %% ============================================================
    if strcmpi(plot_type, 'none'), return; end

    clim_lo = 0;
    clim_hi = min(80, max(wsp(:)) + 5);

    % ---- Figure 1: Original SAR wind ----
    figure('Position', [50 50 800 700]);
    h = pcolor(sar_lon, sar_lat, wsp);
    set(h, 'EdgeColor', 'none'); shading flat;
    set(gca, 'Color', [1 1 1]); axis equal tight;
    xlabel('Longitude (°E)'); ylabel('Latitude (°N)');
    title(sprintf('Original SAR Wind\nCenter: %.2f°N, %.2f°E', ...
          center_lat, center_lon));
    colorbar; colormap('jet'); caxis([clim_lo clim_hi]);
    hold on;
    plot(center_lon, center_lat, 'k+', 'MarkerSize', 15, 'LineWidth', 2);
    hold off;

    % ---- Figure 2: Polar wind field ----
    az_wrap = [az_grid, az_grid(1) + 2*pi];
    V_wrap  = [V_polar, V_polar(:,1)];
    [AZ_w, R_w] = meshgrid(az_wrap, r_grid);
    X_w = R_w .* sin(AZ_w);
    Y_w = R_w .* cos(AZ_w);

    figure('Position', [100 100 800 750]);
    h = pcolor(X_w, Y_w, V_wrap);
    set(h, 'EdgeColor', 'none'); shading flat;
    set(gca, 'Color', [1 1 1]); axis equal;
    xlabel('East (km)'); ylabel('North (km)');
    title(sprintf('SAR Wind (Polar)\nCenter: %.2f°N, %.2f°E', ...
          center_lat, center_lon));
    colorbar; colormap('jet'); caxis([clim_lo clim_hi]);
    hold on;
    th = linspace(0, 2*pi, 300);
    for rr = 50:50:r_max_km
        plot(rr*sin(th), rr*cos(th), 'k:', 'LineWidth', 0.5);
    end
    text(0,  r_max_km*1.05, 'N', 'HorizontalAlignment','center', ...
         'FontWeight','bold', 'FontSize', 14);
    text(r_max_km*1.05, 0, 'E', 'HorizontalAlignment','center', ...
         'FontWeight','bold', 'FontSize', 14);
    text(0, -r_max_km*1.05, 'S', 'HorizontalAlignment','center', ...
         'FontWeight','bold', 'FontSize', 14);
    text(-r_max_km*1.05, 0, 'W', 'HorizontalAlignment','center', ...
         'FontWeight','bold', 'FontSize', 14);
    plot(0, 0, 'k+', 'MarkerSize', 15, 'LineWidth', 2);
    hold off;

    % ---- Figure 3: Sector profiles ----
    figure('Position', [50 50 1200 700]);
    plot(r_grid, V_all(1, :), 'k-', 'LineWidth', 3, ...
         'DisplayName', 'Azimuthal mean');
    hold on;
    sector_colors = lines(8);
    sector_labels = {'0–45°', '45–90°', '90–135°', '135–180°', ...
                     '180–225°', '225–270°', '270–315°', '315–360°'};
    for k = 1:8
        plot(r_grid, V_all(1+k, :), ...
            'Color', sector_colors(k,:), 'LineWidth', 1.5, ...
            'DisplayName', sector_labels{k});
    end
    xlabel('Radius (km)'); ylabel('Wind Speed (m/s)');
    title(sprintf('Sector Profiles (8 × 45°)\nCenter: %.2f°N, %.2f°E', ...
          center_lat, center_lon));
    legend('Location', 'best', 'NumColumns', 3, 'FontSize', 8);
    grid on; xlim([0 r_max_km]);
    hold off;

    % ---- Figure 4: Quadrant stats bar chart ----
    figure('Position', [100 100 1000 700]);
    bar_data = [stats.R34, stats.R50, stats.R64];
    b = bar(1:4, bar_data);
    b(1).FaceColor = 'g';
    b(2).FaceColor = 'm';
    b(3).FaceColor = 'r';
    set(gca, 'XTick', 1:4, 'XTickLabel', quad_labels);
    ylabel('Radius (km)');
    title('R34 / R50 / R64 per Quadrant');
    legend({'R34', 'R50', 'R64'}, 'Location', 'best');
    grid on;

    % ---- Figure 5: Coverage per sector ----
    figure('Position', [100 100 1200 700]);
    plot(r_grid, coverage_all(1, :), 'k-', 'LineWidth', 3, ...
         'DisplayName', 'Azimuthal mean');
    hold on;
    for k = 1:8
        plot(r_grid, coverage_all(1+k, :), ...
            'Color', sector_colors(k,:), 'LineWidth', 1.5, ...
            'DisplayName', sector_labels{k});
    end
    xlabel('Radius (km)'); ylabel('Coverage (%)');
    title(sprintf('Coverage per Sector (%%)\nCenter: %.2f°N, %.2f°E', ...
          center_lat, center_lon));
    legend('Location', 'best', 'NumColumns', 3, 'FontSize', 8);
    grid on; xlim([0 r_max_km]); ylim([0 105]);
    hold off;
end

%% ================================================================
%  Helper: first crossing after Rmax
%% ================================================================
function R = find_first_crossing_after_rmax(r_grid, V_profile, V_threshold, Rmax)
    R = NaN;
    if isnan(Rmax), return; end

    [~, idx_rmax] = min(abs(r_grid - Rmax));

    for i = idx_rmax:length(r_grid)
        if isnan(V_profile(i)), continue; end
        if V_profile(i) < V_threshold
            for j = i:-1:idx_rmax
                if ~isnan(V_profile(j)) && V_profile(j) >= V_threshold
                    R = r_grid(j);
                    return;
                end
            end
            return;
        end
    end
end

%% ================================================================
%  Helper: format number or NaN
%% ================================================================
function s = fmt_num(x)
    if isnan(x)
        s = 'NaN';
    else
        s = sprintf('%.1f', x);
    end
end