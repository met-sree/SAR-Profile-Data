function [Vmax, Rmax, Vmax_p99, Rmax_p99] = ...
         extract_vmax_rmax(lat, lon, wsp, center_lat, center_lon, pctl)
% EXTRACT_VMAX_RMAX  Global Vmax/Rmax plus 99th-percentile Vmax/Rmax
%
% OUTPUTS:
%   Vmax      - global maximum wind speed (m/s)
%   Rmax      - distance from center to that max pixel (km)
%   Vmax_p99  - pctl-th percentile wind speed (m/s), default pctl=99
%   Rmax_p99  - MEAN distance from center to pixels with wsp >= Vmax_p99 (km)

    if nargin < 6 || isempty(pctl), pctl = 99; end

    %% --- Convert to double ---
    lat = double(lat); lon = double(lon); wsp = double(wsp);
    center_lat = double(center_lat); center_lon = double(center_lon);

    %% ============================================================
    %  STEP 1: Distance from center
    %% ============================================================
    R_earth = 6371;
    lat0 = deg2rad(center_lat); lon0 = deg2rad(center_lon);
    latr = deg2rad(lat);        lonr = deg2rad(lon);
    dlon = lonr - lon0;

    cos_c = sin(lat0).*sin(latr) + cos(lat0).*cos(latr).*cos(dlon);
    c = acos(min(max(cos_c, -1), 1));
    dist_km = R_earth * c;

    %% ============================================================
    %  STEP 2: Global Vmax / Rmax (unchanged)
    %% ============================================================
    [Vmax, idx_max] = max(wsp(:));
    Rmax = dist_km(idx_max);

    %% ============================================================
    %  STEP 3: 99th-percentile Vmax + mean radius of that wind class
    %% ============================================================
    Vmax_p99 = NaN; Rmax_p99 = NaN;
    valid = ~isnan(wsp);
    if any(valid(:))
        Vmax_p99 = prctile(wsp(valid), pctl);
        in_class = valid & (wsp >= Vmax_p99);
        if any(in_class(:))
            Rmax_p99 = mean(dist_km(in_class));
        end
    end

    %% ============================================================
    %  STEP 4: Report
    %% ============================================================
    fprintf('\n=== TC Wind Structure ===\n');
    fprintf('  Vmax      = %6.1f m/s  (%.1f kt)\n', Vmax, Vmax/0.5144);
    fprintf('  Rmax      = %6.1f km\n', Rmax);
    fprintf('  Vmax(p%d)  = %6.1f m/s  (%.1f kt)\n', pctl, Vmax_p99, Vmax_p99/0.5144);
    fprintf('  Rmax(p%d)  = %6.1f km  (mean radius of >=p%d wind)\n', ...
            pctl, Rmax_p99, pctl);
end
