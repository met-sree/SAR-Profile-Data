function [wsp_QC, outlier_mask, n_removed] = ...
    sar_qc_patch_mask(lat, lon, wsp, threshold, neighbor_type, R)
% SAR_QC_PATCH_MASK  Flexible patch mask with plots
%
% Removes pixels where the wind speed differs from neighbors by more
% than `threshold` m/s.
%
% Usage:
%   [wsp_QC, outlier_mask, n_removed] = ...
%       sar_qc_patch_mask(sar_lat, sar_lon, wsp_sar_masked, 50, '8');
%
% Inputs:
%   lat, lon       - 2D SAR grids (degrees)
%   wsp            - 2D wind speed (m/s)
%   threshold      - Max allowed difference to neighbors (m/s), default 50
%   neighbor_type  - '4' | '8' | '24' | 'circ'  (default '8')
%   R              - Radius for 'circ' (pixels), default 3
%
% Outputs:
%   wsp_QC         - QC'd wind speed (outliers → NaN)
%   outlier_mask   - Logical mask (true = removed)
%   n_removed      - Number of pixels removed

    %% --- Defaults ---
    if nargin < 4 || isempty(threshold),     threshold     = 50;   end
    if nargin < 5 || isempty(neighbor_type), neighbor_type = '8';  end
    if nargin < 6 || isempty(R),             R             = 3;    end

    %% --- Convert to double ---
    lat = double(lat); lon = double(lon); wsp = double(wsp);

    %% ============================================================
    %  PATCH MASK
    %% ============================================================
    outlier_mask = false(size(wsp));

    switch neighbor_type
        case '4'
            wsp_pad = padarray(wsp, [1 1], NaN);
            for di = -1:1
                for dj = -1:1
                    if di == 0 && dj == 0, continue; end
                    if abs(di) + abs(dj) == 2, continue; end
                    neighbor = wsp_pad(2+di:end-1+di, 2+dj:end-1+dj);
                    outlier_mask = outlier_mask | (abs(wsp - neighbor) > threshold);
                end
            end

        case '8'
            wsp_pad = padarray(wsp, [1 1], NaN);
            for di = -1:1
                for dj = -1:1
                    if di == 0 && dj == 0, continue; end
                    neighbor = wsp_pad(2+di:end-1+di, 2+dj:end-1+dj);
                    outlier_mask = outlier_mask | (abs(wsp - neighbor) > threshold);
                end
            end

        case '24'
            wsp_pad = padarray(wsp, [2 2], NaN);
            for di = -2:2
                for dj = -2:2
                    if di == 0 && dj == 0, continue; end
                    neighbor = wsp_pad(3+di:end-2+di, 3+dj:end-2+dj);
                    outlier_mask = outlier_mask | (abs(wsp - neighbor) > threshold);
                end
            end

        case 'circ'
            wsp_pad = padarray(wsp, [R R], NaN);
            for di = -R:R
                for dj = -R:R
                    if di == 0 && dj == 0, continue; end
                    if sqrt(di^2 + dj^2) > R, continue; end
                    neighbor = wsp_pad(R+1+di:end-R+di, R+1+dj:end-R+dj);
                    outlier_mask = outlier_mask | (abs(wsp - neighbor) > threshold);
                end
            end

        otherwise
            error('Unknown neighbor_type: %s', neighbor_type);
    end

    %% --- Apply mask ---
    wsp_QC = wsp;
    wsp_QC(outlier_mask) = NaN;
    n_removed = sum(outlier_mask(:));

    %% --- Report ---
    n_total_valid = sum(~isnan(wsp(:)));
    fprintf('Patch mask (%s, threshold = %.1f m/s): removed %d pixels (%.2f%%)\n', ...
            neighbor_type, threshold, n_removed, 100*n_removed/n_total_valid);

    %% ============================================================
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
    % xlabel('Longitude', 'FontSize', 11);
    % ylabel('Latitude', 'FontSize', 11);
    % title(sprintf('(a) Before Patch Mask — %d valid pixels', n_total_valid), ...
    %       'FontSize', 12);
    % colorbar; colormap('jet'); caxis([clim_lo clim_hi]);
    % 
    % subplot(1,2,2);
    % contourf(lon, lat, wsp_QC, levels, 'LineStyle', 'none');
    % axis equal tight;
    % xlabel('Longitude', 'FontSize', 11);
    % ylabel('Latitude', 'FontSize', 11);
    % title(sprintf('(b) After Patch Mask — %d removed', n_removed), ...
    %       'FontSize', 12);
    % colorbar; colormap('jet'); caxis([clim_lo clim_hi]);
    % 
    % sgtitle(sprintf('Patch Mask QC (%s-neighbor, threshold = %.1f m/s)', ...
    %         neighbor_type, threshold), 'FontSize', 14, 'FontWeight', 'bold');
    % 
    % %% ============================================================
    % %  PLOT 2: Outlier map
    % %% ============================================================
    % figure('Position', [100 100 750 600]);
    % contourf(lon, lat, double(outlier_mask), [0 0.5 1], 'LineStyle', 'none');
    % axis equal tight;
    % xlabel('Longitude', 'FontSize', 11);
    % ylabel('Latitude', 'FontSize', 11);
    % title(sprintf('Patch Mask Outliers (%s-neighbor, n = %d)', ...
    %       neighbor_type, n_removed), 'FontSize', 12);
    % colorbar;
    % colormap([1 1 1; 1 0 0]);   % white = kept, red = removed
    % caxis([0 1]);

end