function [wsp_QC, ring_id, n_removed, outlier_mask, ring_stats] = ...
    sar_qc_polar_grid(lat, lon, wsp, center_lat, center_lon, ...
                      dr_km, method, k)
% SAR_QC_POLAR_GRID  Outlier removal on polar grid (ring-based QC)
%
% Creates a polar grid with dr_km radial resolution and removes outliers
% based on statistics per radial ring.
%
% Usage:
%   [wsp_QC, ~, n_removed, ~, ~] = sar_qc_polar_grid(...
%       sar_lat, sar_lon, wsp_sar_masked, ...
%       usa_lat(i,1), usa_lon(i,1), ...
%       1, 'sigma', 3);
%
% Inputs:
%   lat, lon      - 2D SAR grids (degrees)
%   wsp           - 2D wind speed (m/s)
%   center_lat    - TC center latitude
%   center_lon    - TC center longitude
%   dr_km         - Ring width (km), default 1
%   method        - 'sigma' | 'mad'   (default: 'sigma')
%   k             - threshold multiplier (default: 3)
%
% Outputs:
%   wsp_QC        - QC'd wind speed
%   ring_id       - Radial ring index for each pixel
%   n_removed     - Number of outliers removed
%   outlier_mask  - Logical mask (true = outlier)
%   ring_stats    - Struct with per-ring statistics

    %% --- Defaults ---
    if nargin < 6 || isempty(dr_km),  dr_km  = 1;       end
    if nargin < 7 || isempty(method), method = 'sigma'; end
    if nargin < 8 || isempty(k),      k      = 3;       end

    %% --- Convert to double ---
    lat = double(lat); lon = double(lon); wsp = double(wsp);
    center_lat = double(center_lat); center_lon = double(center_lon);

    %% ============================================================
    %  STEP 1: Distance & azimuth from center
    %% ============================================================
    R_earth = 6371;
    lat0 = deg2rad(center_lat); lon0 = deg2rad(center_lon);
    latr = deg2rad(lat);        lonr = deg2rad(lon);
    dlon = lonr - lon0;

    cos_c = sin(lat0).*sin(latr) + cos(lat0).*cos(latr).*cos(dlon);
    c = acos(min(max(cos_c, -1), 1));
    dist_km = R_earth * c;

    %% ============================================================
    %  STEP 2: Assign each pixel to a radial ring
    %% ============================================================
    r_max = max(dist_km(:), [], 'omitnan');
    r_max = ceil(r_max / dr_km) * dr_km;
    r_edges = 0:dr_km:r_max;
    n_r = length(r_edges) - 1;

    ring_id = discretize(dist_km(:), r_edges);
    ring_id = reshape(ring_id, size(wsp));

    fprintf('Polar grid: %d rings (%.1f km each)\n', n_r, dr_km);

    %% ============================================================
    %  STEP 3: Per-ring outlier detection
    %% ============================================================
    wsp_QC = wsp;
    outlier_mask = false(size(wsp));

    ring_stats = struct('mean',   nan(1, n_r), ...
                        'std',    nan(1, n_r), ...
                        'median', nan(1, n_r), ...
                        'mad',    nan(1, n_r), ...
                        'n',      zeros(1, n_r));

    for kk = 1:n_r
        in_ring = (ring_id == kk) & ~isnan(wsp);
        V_ring = wsp(in_ring);

        if numel(V_ring) < 5
            continue;
        end

        % Ring statistics
        ring_stats.mean(kk)   = mean(V_ring);
        ring_stats.std(kk)    = std(V_ring);
        ring_stats.median(kk) = median(V_ring);
        ring_stats.mad(kk)    = mad(V_ring, 1);
        ring_stats.n(kk)      = numel(V_ring);

        % Outlier detection
        switch lower(method)
            case 'sigma'
                if ring_stats.std(kk) == 0, continue; end
                outliers = abs(V_ring - ring_stats.mean(kk)) > ...
                           k * ring_stats.std(kk);

            case 'mad'
                if ring_stats.mad(kk) == 0, continue; end
                mod_z = 0.6745 * abs(V_ring - ring_stats.median(kk)) / ...
                        ring_stats.mad(kk);
                outliers = mod_z > k;

            otherwise
                error('Unknown method: %s', method);
        end

        % Apply outliers
        idx_in_ring = find(in_ring);
        outlier_mask(idx_in_ring(outliers)) = true;
    end

    % Apply mask
    wsp_QC(outlier_mask) = NaN;
    n_removed = sum(outlier_mask(:));

    fprintf('Polar QC (%s, k=%.1f, dr=%.1f km): removed %d outliers\n', ...
            method, k, dr_km, n_removed);

    % %% ============================================================
    % %  PLOT 1: Before vs After
    % %% ============================================================
    % figure('Position', [50 50 1400 600]);
    % clim_lo = 0;
    % clim_hi = min(80, max(wsp(:)) + 5);
    % levels = linspace(clim_lo, clim_hi, 20);
    % 
    % subplot(1,2,1);
    % contourf(lon, lat, wsp, levels, 'LineStyle', 'none');
    % axis equal tight;
    % xlabel('Longitude'); ylabel('Latitude');
    % title(sprintf('(a) Before QC — %d valid', sum(~isnan(wsp(:)))));
    % colorbar; colormap('jet'); caxis([clim_lo clim_hi]);
    % 
    % subplot(1,2,2);
    % contourf(lon, lat, wsp_QC, levels, 'LineStyle', 'none');
    % axis equal tight;
    % xlabel('Longitude'); ylabel('Latitude');
    % title(sprintf('(b) After Polar QC — %d removed', n_removed));
    % colorbar; colormap('jet'); caxis([clim_lo clim_hi]);
    % 
    % sgtitle(sprintf('Polar Grid QC — %.1f km rings', dr_km), ...
    %         'FontSize', 14, 'FontWeight', 'bold');
    % 
    % %% ============================================================
    % %  PLOT 2: Outlier map
    % %% ============================================================
    % figure('Position', [100 100 750 600]);
    % contourf(lon, lat, double(outlier_mask), [0 0.5 1], 'LineStyle', 'none');
    % axis equal tight;
    % xlabel('Longitude'); ylabel('Latitude');
    % title(sprintf('Outliers Removed (n = %d)', n_removed));
    % colorbar;
    % colormap([1 1 1; 1 0 0]);   % white = kept, red = removed
    % caxis([0 1]);

end