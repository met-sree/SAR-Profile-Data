function [wsp_etopo, wsp_buffered, n_removed_etopo, n_removed_buffered] = ...
    apply_etopo_landmask(z, etopo_lat, etopo_lon, ...
                         lat, lon, wsp, margin_deg)
% APPLY_ETOPO_LANDMASK  ETOPO land mask + 2.5 km dilation + 3-panel plot
%
% Outputs:
%   wsp_etopo           - Wind after ETOPO mask (land → NaN)
%   wsp_buffered        - Wind after ETOPO + 2.5 km dilation (final)
%   n_removed_etopo     - Pixels removed by ETOPO mask
%   n_removed_buffered  - Pixels removed by ETOPO + dilation

    %% --- Defaults ---
    if nargin < 7 || isempty(margin_deg), margin_deg = 1; end

    %% --- Convert to double ---
    lat = double(lat); lon = double(lon); wsp = double(wsp);
    z = double(z);
    etopo_lat = double(etopo_lat);
    etopo_lon = double(etopo_lon);

    %% --- Ensure lat/lon are 2D matrices ---
    if isvector(lat) && isvector(lon)
        [lon, lat] = meshgrid(lon, lat);
    end

    %% --- Size check ---
    assert(isequal(size(lat), size(lon), size(wsp)), ...
           'lat, lon, wsp must be same size');

    %% ============================================================
    %  STEP 0: Auto-detect z orientation
    %% ============================================================
    n_lat = length(etopo_lat);
    n_lon = length(etopo_lon);

    if size(z, 1) == n_lat && size(z, 2) == n_lon
        fprintf('ETOPO z is (lat, lon): %d × %d\n', n_lat, n_lon);
    elseif size(z, 1) == n_lon && size(z, 2) == n_lat
        z = z';
        fprintf('ETOPO z transposed to (lat, lon): %d × %d\n', n_lat, n_lon);
    else
        error('ETOPO z size mismatch');
    end

    %% ============================================================
    %  STEP 1: Normalize longitude
    %% ============================================================
    if max(etopo_lon) > 180
        etopo_lon = mod(etopo_lon + 180, 360) - 180;
        [etopo_lon, sort_idx] = sort(etopo_lon);
        z = z(:, sort_idx);
    end

    %% ============================================================
    %  STEP 2: Crop ETOPO to SAR bbox
    %% ============================================================
    lat_min = min(lat(:)) - margin_deg;
    lat_max = max(lat(:)) + margin_deg;
    lon_min = min(lon(:)) - margin_deg;
    lon_max = max(lon(:)) + margin_deg;

    lat_idx = find(etopo_lat >= lat_min & etopo_lat <= lat_max);
    lon_idx = find(etopo_lon >= lon_min & etopo_lon <= lon_max);

    z_crop   = z(lat_idx, lon_idx);
    lat_crop = etopo_lat(lat_idx);
    lon_crop = etopo_lon(lon_idx);

    fprintf('ETOPO crop: %d × %d\n', size(z_crop,1), size(z_crop,2));

    %% ============================================================
    %  STEP 3: Interpolate ETOPO to SAR grid
    %% ============================================================
    [LON_crop, LAT_crop] = meshgrid(lon_crop, lat_crop);
    F = scatteredInterpolant(LON_crop(:), LAT_crop(:), z_crop(:), ...
                             'linear', 'none');

    etopo_at_sar = nan(size(wsp));
    valid_pts = ~isnan(lon) & ~isnan(lat);

    if sum(valid_pts(:)) > 0
        etopo_at_sar(valid_pts) = F(lon(valid_pts), lat(valid_pts));
    end

    %% ============================================================
    %  STEP 4: ETOPO land mask (elevation > -2 m)
    %% ============================================================
    land_mask_etopo = etopo_at_sar > -2;
    land_mask_etopo(isnan(etopo_at_sar)) = true;

    n_removed_etopo = sum(land_mask_etopo(:) & ~isnan(wsp(:)));

    fprintf('ETOPO only:      %d pixels removed\n', n_removed_etopo);

    %% ============================================================
    %  STEP 5: Dilation — 5 iterations (~2.5 km buffer)
    %% ============================================================
    land_mask_buffered = land_mask_etopo;
    for iter = 1:5
        m = land_mask_buffered;
        m = m | [m(2:end,:);   false(1,size(m,2))];     % N
        m = m | [false(1,size(m,2)); m(1:end-1,:)];     % S
        m = m | [m(:,2:end),   false(size(m,1),1)];     % E
        m = m | [false(size(m,1),1), m(:,1:end-1)];     % W
        land_mask_buffered = m;
    end

    n_removed_buffered = sum(land_mask_buffered(:) & ~isnan(wsp(:)));

    fprintf('ETOPO + 2.5 km:  %d pixels removed (+%d)\n', ...
            n_removed_buffered, n_removed_buffered - n_removed_etopo);

    %% ============================================================
    %  STEP 6: Apply masks
    %% ============================================================
    % Wind with ETOPO-only mask
    wsp_etopo = wsp;
    wsp_etopo(land_mask_etopo) = NaN;

    % Wind with ETOPO + dilation mask (final)
    wsp_buffered = wsp;
    wsp_buffered(land_mask_buffered) = NaN;

    fprintf('Final removed:   %d pixels (%.2f%%)\n', ...
            n_removed_buffered, ...
            100*n_removed_buffered/sum(~isnan(wsp(:))));

    % %% ============================================================
    % %  PLOT: 3 panels — Before | ETOPO only | ETOPO + dilation
    % %% ============================================================
    % figure('Position', [50 50 1600 550]);
    % 
    % clim_lo = 0;
    % clim_hi = min(80, max(wsp(:)) + 5);
    % levels  = linspace(clim_lo, clim_hi, 20);
    % 
    % % --- Panel (a): Before ---
    % subplot(1,3,1);
    % contourf(lon, lat, wsp, levels, 'LineStyle', 'none');
    % axis equal tight;
    % xlabel('Longitude', 'FontSize', 10);
    % ylabel('Latitude', 'FontSize', 10);
    % title(sprintf('(a) Before — %d valid', sum(~isnan(wsp(:)))), ...
    %       'FontSize', 12);
    % colorbar; colormap('jet'); caxis([clim_lo clim_hi]);
    % 
    % % --- Panel (b): ETOPO only ---
    % subplot(1,3,2);
    % contourf(lon, lat, wsp_etopo, levels, 'LineStyle', 'none');
    % axis equal tight;
    % xlabel('Longitude', 'FontSize', 10);
    % ylabel('Latitude', 'FontSize', 10);
    % title(sprintf('(b) ETOPO only — %d removed', n_removed_etopo), ...
    %       'FontSize', 12);
    % colorbar; colormap('jet'); caxis([clim_lo clim_hi]);
    % 
    % % --- Panel (c): ETOPO + 2.5 km dilation ---
    % subplot(1,3,3);
    % contourf(lon, lat, wsp_buffered, levels, 'LineStyle', 'none');
    % axis equal tight;
    % xlabel('Longitude', 'FontSize', 10);
    % ylabel('Latitude', 'FontSize', 10);
    % title(sprintf('(c) ETOPO + 2.5 km — %d removed', n_removed_buffered), ...
    %       'FontSize', 12);
    % colorbar; colormap('jet'); caxis([clim_lo clim_hi]);
    % 
    % sgtitle('ETOPO Land Mask: Before | ETOPO only | ETOPO + 2.5 km Buffer', ...
    %         'FontSize', 14, 'FontWeight', 'bold');

end