function process_tc_case(case_name)
% PROCESS_TC_CASE  Process all SAR snapshots for a TC case
%   process_tc_case('WP242025_RAGASA');
%
% 500 m version — also computes 99th-percentile Vmax/Rmax per branch
% (SAR.FieldStatics_p99). Requires extract_vmax_rmax.m returning 4 outputs.

    clc

    %% ============================================================
    %  USER CONFIGURATION — CHANGE VALUES HERE
    %% ============================================================
    CFG.addpath_path = '/scratch/ynekkali/SAR_DATA_PreProcessing/500m/Script/Functions/';
    CFG.folder_base  = '/scratch/ynekkali/NOAA_STAR_500m/Named_TC/';
    CFG.ibfilePath   = '/scratch/ynekkali/IBTrACS.ALL.v04r01.nc';
    CFG.etopo_file   = '/scratch/ynekkali/ETOPO_2022_v1_60s_N90W180_bed.nc';
    CFG.mat_folder   = '/scratch/ynekkali/SAR_DATA_PreProcessing/500m/Processed_Data/';

    CFG.margin       = 1;
    CFG.dr_qc        = 1;
    CFG.method_qc    = 'sigma';
    CFG.k_qc         = 3;
    CFG.threshold_patch = 50;
    CFG.neighbor     = '8';
    CFG.dr_polar     = 0.5;
    CFG.rmax_polar   = 300;
    CFG.daz_polar    = 0.5;
    CFG.plot_type    = 'none';
    CFG.pctl         = 99;      % >>> ADDED: percentile for Vmax/Rmax
    % ============================================================
    %  END OF USER CONFIGURATION
    % ============================================================

    %% ============================================================
    %  SETUP — Apply configuration
    %% ============================================================
    fprintf('\n=== PROCESS TC CASE: %s ===\n', case_name);
    addpath(CFG.addpath_path);

    folder_path    = strcat(CFG.folder_base, case_name, '/');
    ibfilePath     = CFG.ibfilePath;
    etopo_file     = CFG.etopo_file;
    mat_folderPath = CFG.mat_folder;

    if ~exist(mat_folderPath, 'dir'), mkdir(mat_folderPath); end
    tc_name = case_name(5:end);

    %% ============================================================
    %  READ ETOPO
    %% ============================================================
    fprintf('\n[1/5] Reading ETOPO...\n');
    z         = ncread(etopo_file, 'z');
    etopo_lat = ncread(etopo_file, 'lat');
    etopo_lon = ncread(etopo_file, 'lon');
    fprintf('      ETOPO size: %d × %d\n', size(z, 1), size(z, 2));

    %% ============================================================
    %  READ IBTrACS
    %% ============================================================
    fprintf('\n[2/5] Reading IBTrACS...\n');
    name_ibtracs = ncread(ibfilePath, 'name');
    sid_ibtracs  = ncread(ibfilePath, 'sid');
    sid_str      = cellstr(sid_ibtracs');
    name_str     = cellstr(name_ibtracs');
    ib_year      = cellfun(@(x) x(1:min(4, length(x))), sid_str, 'UniformOutput', false);
    tc_year_name = cellfun(@(x,y) [x '_' y], ib_year, name_str, 'UniformOutput', false);
    [is_member, ~] = ismember(tc_year_name, tc_name);
    matching_tc_indices = find(is_member);

    if isempty(matching_tc_indices)
        error('TC %s not found in IBTrACS', tc_name);
    end

    IBTrACS   = create_IBTrACKS(ibfilePath, matching_tc_indices(1));
    valid_idx = ~isnan(IBTrACS.usa_lon) & ~isnan(IBTrACS.usa_lat);
    fprintf('      IBTrACS timesteps: %d\n', length(IBTrACS.datetimes));

    % ============================================================
    %  normalise USA radii to n_time x 4 and build masks/axes
    % ============================================================
    n_time_ib = numel(IBTrACS.datetimes);
    usa_r34_q = normalise_quad(IBTrACS.usa_r34, n_time_ib);
    usa_r50_q = normalise_quad(IBTrACS.usa_r50, n_time_ib);
    usa_r64_q = normalise_quad(IBTrACS.usa_r64, n_time_ib);

    idx_usa_r34 = valid_idx & any(~isnan(usa_r34_q), 2);
    idx_usa_r50 = valid_idx & any(~isnan(usa_r50_q), 2);
    idx_usa_r64 = valid_idx & any(~isnan(usa_r64_q), 2);

    if any(idx_usa_r34), t_usa_r34 = datenum(IBTrACS.datetimes(idx_usa_r34)); else, t_usa_r34 = []; end
    if any(idx_usa_r50), t_usa_r50 = datenum(IBTrACS.datetimes(idx_usa_r50)); else, t_usa_r50 = []; end
    if any(idx_usa_r64), t_usa_r64 = datenum(IBTrACS.datetimes(idx_usa_r64)); else, t_usa_r64 = []; end

    fprintf('      USA radii rows -> R34:%d  R50:%d  R64:%d\n', ...
            nnz(idx_usa_r34), nnz(idx_usa_r50), nnz(idx_usa_r64));
    % ============================================================

    %% ============================================================
    %  LIST TEXT FILES
    %% ============================================================
    fprintf('\n[3/5] Listing SAR snapshot files...\n');
    if ~exist(folder_path, 'dir')
        error('Folder not found: %s', folder_path);
    end

    txt_files = dir(fullfile(folder_path, '*.txt'));
    txt_names = {txt_files.name};

    if isempty(txt_names)
        warning('No text files found');
        return;
    end

    num_files = length(txt_names);
    fprintf('      Found %d text files\n', num_files);

    %% ============================================================
    %  PREALLOCATE — IBTrACS_Interpolated
    %% ============================================================
    IBTrACS_Interpolated.snap_time    = NaT(num_files, 1);
    IBTrACS_Interpolated.snap_cen_lat = nan(num_files, 1);
    IBTrACS_Interpolated.snap_cen_lon = nan(num_files, 1);
    IBTrACS_Interpolated.usa_lat      = nan(num_files, 1);
    IBTrACS_Interpolated.usa_lon      = nan(num_files, 1);
    IBTrACS_Interpolated.usa_rmax     = nan(num_files, 1);
    IBTrACS_Interpolated.usa_vmax     = nan(num_files, 1);
    IBTrACS_Interpolated.usa_sshs     = nan(num_files, 1);
    IBTrACS_Interpolated.usa_status   = cell(num_files, 1);
    IBTrACS_Interpolated.cma_lat      = nan(num_files, 1);
    IBTrACS_Interpolated.cma_lon      = nan(num_files, 1);
    IBTrACS_Interpolated.cma_vmax     = nan(num_files, 1);
    IBTrACS_Interpolated.tokyo_lat    = nan(num_files, 1);
    IBTrACS_Interpolated.tokyo_lon    = nan(num_files, 1);
    IBTrACS_Interpolated.tokyo_vmax   = nan(num_files, 1);

    IBTrACS_Interpolated.usa_r34_mean = nan(num_files, 1);
    IBTrACS_Interpolated.usa_r50_mean = nan(num_files, 1);
    IBTrACS_Interpolated.usa_r64_mean = nan(num_files, 1);

    %% ============================================================
    %  PREALLOCATE — SAR
    %% ============================================================
    SAR.snap_time    = NaT(num_files, 1);
    SAR.snap_cen_lat = nan(num_files, 1);
    SAR.snap_cen_lon = nan(num_files, 1);
    SAR.base_name    = cell(num_files, 1);

    SAR.RadialProfile_raw      = [];
    SAR.RadialProfile_etopo    = [];
    SAR.RadialProfile_buffered = [];
    SAR.RadialProfile_3sigma   = [];
    SAR.RadialProfile_patch    = [];

    SAR.CoverageProfile_raw      = [];
    SAR.CoverageProfile_etopo    = [];
    SAR.CoverageProfile_buffered = [];
    SAR.CoverageProfile_3sigma   = [];
    SAR.CoverageProfile_patch    = [];

    SAR.ProfileStatics_raw      = nan(4, 5, num_files);
    SAR.ProfileStatics_etopo    = nan(4, 5, num_files);
    SAR.ProfileStatics_buffered = nan(4, 5, num_files);
    SAR.ProfileStatics_3sigma   = nan(4, 5, num_files);
    SAR.ProfileStatics_patch    = nan(4, 5, num_files);

    SAR.wsp_3sigma = cell(num_files, 1);
    SAR.wsp_patch  = cell(num_files, 1);

    SAR.sar_lat = cell(num_files, 1);
    SAR.sar_lon = cell(num_files, 1);

    SAR.FieldStatics = nan(2, num_files, 5);

    % >>> ADDED: 99th-percentile Vmax/Rmax per branch
    SAR.FieldStatics_p99 = nan(2, num_files, 5);   % 1=Vmax_p99, 2=Rmax_p99
    % >>> END ADDED

    SAR.r_grid = [];

    if any(valid_idx)
        valid_times = datenum(IBTrACS.datetimes(valid_idx));
    else
        valid_times = [];
    end

    %% ============================================================
    %  LOOP OVER ALL SNAPSHOTS
    %% ============================================================
    fprintf('\n[4/5] Processing %d snapshots...\n', num_files);

    for i = 1:num_files
        fprintf('\n--- Snapshot %d/%d: %s ---\n', i, num_files, txt_names{i});

        [~, base_name, ~] = fileparts(txt_names{i});
        dt = datetime(regexp(base_name, '\d+_\d+_\d+_\d+_\d+_\d+', 'match'), ...
                      'InputFormat', 'yyyy_MM_dd_HH_mm_ss');

        [snap_cen_lat, snap_cen_lon, snap_date] = ...
            ExtractStromInfoFromTXTfile(folder_path, base_name);

        base_name_nc = strrep(base_name, '_NOAAinfo_summary', '');
        nc_file = strcat(folder_path, base_name_nc, '_wind_level2.nc');

        if ~exist(nc_file, 'file')
            warning('NetCDF not found: %s', nc_file);
            continue;
        end

        sar_lat = ncread(nc_file, 'latitude');
        sar_lon = ncread(nc_file, 'longitude');
        wsp_sar = ncread(nc_file, 'sar_wind');
        mask    = ncread(nc_file, 'mask');

        wsp_sar_masked = wsp_sar;
        wsp_sar_masked(mask == 1) = NaN;

        fprintf('  SAR size: %d × %d\n', size(wsp_sar_masked, 1), size(wsp_sar_masked, 2));

        center_lat = snap_cen_lat;
        center_lon = snap_cen_lon;
        fprintf('  Center: %.3f°N, %.3f°E\n', center_lat, center_lon);

        % ---- QC: ETOPO land mask ----
        fprintf('  ETOPO land mask (margin = %.2f°)...\n', CFG.margin);
        [wsp_etopo, wsp_buffered, ~, ~] = ...
            apply_etopo_landmask(z, etopo_lat, etopo_lon, ...
                                 sar_lat, sar_lon, wsp_sar_masked, CFG.margin);

        % ---- QC: 3-sigma polar grid ----
        fprintf('  3-sigma QC (dr = %.2f km, method = %s, k = %.1f)...\n', ...
                CFG.dr_qc, CFG.method_qc, CFG.k_qc);
        [wsp_3sigma, ~, ~, ~, ~] = ...
            sar_qc_polar_grid(sar_lat, sar_lon, wsp_buffered, ...
                              center_lat, center_lon, ...
                              CFG.dr_qc, CFG.method_qc, CFG.k_qc);

        % ---- QC: patch mask ----
        fprintf('  Patch mask QC (threshold = %.1f m/s, %s-neighbor)...\n', ...
                CFG.threshold_patch, CFG.neighbor);
        [wsp_patch, ~, ~] = ...
            sar_qc_patch_mask(sar_lat, sar_lon, wsp_buffered, ...
                              CFG.threshold_patch, CFG.neighbor);

        % ---- Polar conversion for all 5 wind fields ----
        fprintf('  Polar conversion (dr = %.2f km, rmax = %.0f km, daz = %.2f°)...\n', ...
                CFG.dr_polar, CFG.rmax_polar, CFG.daz_polar);

        fprintf('    - RAW\n');
        [RadialProfile_raw, r_grid, ProfileStatics_raw, CoverageProfile_raw] = ...
            sar_to_polar_profiles(sar_lat, sar_lon, wsp_sar_masked, ...
                                  center_lat, center_lon, ...
                                  CFG.dr_polar, CFG.rmax_polar, CFG.daz_polar, ...
                                  CFG.plot_type);

        fprintf('    - ETOPO\n');
        [RadialProfile_etopo, ~, ProfileStatics_etopo, CoverageProfile_etopo] = ...
            sar_to_polar_profiles(sar_lat, sar_lon, wsp_etopo, ...
                                  center_lat, center_lon, ...
                                  CFG.dr_polar, CFG.rmax_polar, CFG.daz_polar, ...
                                  CFG.plot_type);

        fprintf('    - BUFFERED\n');
        [RadialProfile_buffered, ~, ProfileStatics_buffered, CoverageProfile_buffered] = ...
            sar_to_polar_profiles(sar_lat, sar_lon, wsp_buffered, ...
                                  center_lat, center_lon, ...
                                  CFG.dr_polar, CFG.rmax_polar, CFG.daz_polar, ...
                                  CFG.plot_type);

        fprintf('    - 3-SIGMA\n');
        [RadialProfile_3sigma, ~, ProfileStatics_3sigma, CoverageProfile_3sigma] = ...
            sar_to_polar_profiles(sar_lat, sar_lon, wsp_3sigma, ...
                                  center_lat, center_lon, ...
                                  CFG.dr_polar, CFG.rmax_polar, CFG.daz_polar, ...
                                  CFG.plot_type);

        fprintf('    - PATCH\n');
        [RadialProfile_patch, ~, ProfileStatics_patch, CoverageProfile_patch] = ...
            sar_to_polar_profiles(sar_lat, sar_lon, wsp_patch, ...
                                  center_lat, center_lon, ...
                                  CFG.dr_polar, CFG.rmax_polar, CFG.daz_polar, ...
                                  CFG.plot_type);

        % ---- Allocate on first iteration ----
        if i == 1
            [n_prof, n_r] = size(RadialProfile_etopo);
            fprintf('  RadialProfile size: %d profiles × %d radii\n', n_prof, n_r);

            SAR.RadialProfile_raw      = nan(n_prof, n_r, num_files);
            SAR.RadialProfile_etopo    = nan(n_prof, n_r, num_files);
            SAR.RadialProfile_buffered = nan(n_prof, n_r, num_files);
            SAR.RadialProfile_3sigma   = nan(n_prof, n_r, num_files);
            SAR.RadialProfile_patch    = nan(n_prof, n_r, num_files);

            SAR.CoverageProfile_raw      = nan(n_prof, n_r, num_files);
            SAR.CoverageProfile_etopo    = nan(n_prof, n_r, num_files);
            SAR.CoverageProfile_buffered = nan(n_prof, n_r, num_files);
            SAR.CoverageProfile_3sigma   = nan(n_prof, n_r, num_files);
            SAR.CoverageProfile_patch    = nan(n_prof, n_r, num_files);

            SAR.r_grid         = r_grid;

            SAR.profile_labels = ProfileStatics_etopo.labels;
            SAR.metric_labels  = {'Vmax', 'Rmax', 'R34', 'R50', 'R64'};
        end

        % ---- Store profiles ----
        SAR.RadialProfile_raw(:,:,i)      = RadialProfile_raw;
        SAR.RadialProfile_etopo(:,:,i)    = RadialProfile_etopo;
        SAR.RadialProfile_buffered(:,:,i) = RadialProfile_buffered;
        SAR.RadialProfile_3sigma(:,:,i)   = RadialProfile_3sigma;
        SAR.RadialProfile_patch(:,:,i)    = RadialProfile_patch;

        % ---- Store coverage ----
        SAR.CoverageProfile_raw(:,:,i)      = CoverageProfile_raw;
        SAR.CoverageProfile_etopo(:,:,i)    = CoverageProfile_etopo;
        SAR.CoverageProfile_buffered(:,:,i) = CoverageProfile_buffered;
        SAR.CoverageProfile_3sigma(:,:,i)   = CoverageProfile_3sigma;
        SAR.CoverageProfile_patch(:,:,i)    = CoverageProfile_patch;

        % ---- Build 4 × 5 matrices (quadrants) ----
        M_raw      = [ProfileStatics_raw.Vmax(:),      ProfileStatics_raw.Rmax(:),      ...
                      ProfileStatics_raw.R34(:),       ProfileStatics_raw.R50(:),       ...
                      ProfileStatics_raw.R64(:)];
        M_etopo    = [ProfileStatics_etopo.Vmax(:),    ProfileStatics_etopo.Rmax(:),    ...
                      ProfileStatics_etopo.R34(:),     ProfileStatics_etopo.R50(:),     ...
                      ProfileStatics_etopo.R64(:)];
        M_buffered = [ProfileStatics_buffered.Vmax(:), ProfileStatics_buffered.Rmax(:), ...
                      ProfileStatics_buffered.R34(:),  ProfileStatics_buffered.R50(:),  ...
                      ProfileStatics_buffered.R64(:)];
        M_3sigma   = [ProfileStatics_3sigma.Vmax(:),   ProfileStatics_3sigma.Rmax(:),   ...
                      ProfileStatics_3sigma.R34(:),    ProfileStatics_3sigma.R50(:),    ...
                      ProfileStatics_3sigma.R64(:)];
        M_patch    = [ProfileStatics_patch.Vmax(:),    ProfileStatics_patch.Rmax(:),    ...
                      ProfileStatics_patch.R34(:),     ProfileStatics_patch.R50(:),     ...
                      ProfileStatics_patch.R64(:)];

        SAR.ProfileStatics_raw(:,:,i)      = M_raw;
        SAR.ProfileStatics_etopo(:,:,i)    = M_etopo;
        SAR.ProfileStatics_buffered(:,:,i) = M_buffered;
        SAR.ProfileStatics_3sigma(:,:,i)   = M_3sigma;
        SAR.ProfileStatics_patch(:,:,i)    = M_patch;

        % ---- Global Vmax/Rmax (max) + 99th-percentile Vmax/Rmax ----
        fprintf('  Global Vmax/Rmax + p99...\n');

        [Vmax_raw_g,   Rmax_raw_g,   Vp99_raw,   Rp99_raw] = ...
            extract_vmax_rmax(sar_lat, sar_lon, wsp_sar_masked, center_lat, center_lon, CFG.pctl);
        [Vmax_etopo_g, Rmax_etopo_g, Vp99_etopo, Rp99_etopo] = ...
            extract_vmax_rmax(sar_lat, sar_lon, wsp_etopo,      center_lat, center_lon, CFG.pctl);
        [Vmax_buf_g,   Rmax_buf_g,   Vp99_buf,   Rp99_buf] = ...
            extract_vmax_rmax(sar_lat, sar_lon, wsp_buffered,   center_lat, center_lon, CFG.pctl);
        [Vmax_3sig_g,  Rmax_3sig_g,  Vp99_3sig,  Rp99_3sig] = ...
            extract_vmax_rmax(sar_lat, sar_lon, wsp_3sigma,     center_lat, center_lon, CFG.pctl);
        [Vmax_patch_g, Rmax_patch_g, Vp99_patch, Rp99_patch] = ...
            extract_vmax_rmax(sar_lat, sar_lon, wsp_patch,      center_lat, center_lon, CFG.pctl);

        SAR.FieldStatics(1, i, 1) = Vmax_raw_g;   SAR.FieldStatics(2, i, 1) = Rmax_raw_g;
        SAR.FieldStatics(1, i, 2) = Vmax_etopo_g; SAR.FieldStatics(2, i, 2) = Rmax_etopo_g;
        SAR.FieldStatics(1, i, 3) = Vmax_buf_g;   SAR.FieldStatics(2, i, 3) = Rmax_buf_g;
        SAR.FieldStatics(1, i, 4) = Vmax_3sig_g;  SAR.FieldStatics(2, i, 4) = Rmax_3sig_g;
        SAR.FieldStatics(1, i, 5) = Vmax_patch_g; SAR.FieldStatics(2, i, 5) = Rmax_patch_g;

        SAR.FieldStatics_p99(1, i, 1) = Vp99_raw;   SAR.FieldStatics_p99(2, i, 1) = Rp99_raw;
        SAR.FieldStatics_p99(1, i, 2) = Vp99_etopo; SAR.FieldStatics_p99(2, i, 2) = Rp99_etopo;
        SAR.FieldStatics_p99(1, i, 3) = Vp99_buf;   SAR.FieldStatics_p99(2, i, 3) = Rp99_buf;
        SAR.FieldStatics_p99(1, i, 4) = Vp99_3sig;  SAR.FieldStatics_p99(2, i, 4) = Rp99_3sig;
        SAR.FieldStatics_p99(1, i, 5) = Vp99_patch; SAR.FieldStatics_p99(2, i, 5) = Rp99_patch;

        % ---- Store 2D QC wind fields ----
        SAR.wsp_3sigma{i} = single(wsp_3sigma);
        SAR.wsp_patch{i}  = single(wsp_patch);

        % ---- Store native lat/lon ----
        SAR.sar_lat{i} = sar_lat;
        SAR.sar_lon{i} = sar_lon;

        SAR.snap_time(i)    = snap_date;
        SAR.snap_cen_lat(i) = center_lat;
        SAR.snap_cen_lon(i) = center_lon;
        SAR.base_name{i}    = base_name_nc;

        % ---- IBTrACS_Interpolated ----
        IBTrACS_Interpolated.snap_time(i)    = snap_date;
        IBTrACS_Interpolated.snap_cen_lat(i) = center_lat;
        IBTrACS_Interpolated.snap_cen_lon(i) = center_lon;

        if any(valid_idx)
            IBTrACS_Interpolated.usa_lat(i)    = interp1(valid_times, IBTrACS.usa_lat(valid_idx),  datenum(dt), 'linear');
            IBTrACS_Interpolated.usa_lon(i)    = interp1(valid_times, IBTrACS.usa_lon(valid_idx),  datenum(dt), 'linear');
            IBTrACS_Interpolated.usa_rmax(i)   = interp1(valid_times, IBTrACS.usa_rmw(valid_idx),  datenum(dt), 'linear');
            IBTrACS_Interpolated.usa_vmax(i)   = interp1(valid_times, IBTrACS.usa_wind(valid_idx), datenum(dt), 'linear');
            IBTrACS_Interpolated.usa_sshs(i)   = interp1(valid_times, IBTrACS.sshs(valid_idx),     datenum(dt), 'nearest');

            usa_status1 = IBTrACS.status(valid_idx);
            [~, nearest_idx] = min(abs(valid_times - datenum(dt)));
            IBTrACS_Interpolated.usa_status{i} = usa_status1{nearest_idx};

            IBTrACS_Interpolated.cma_lat(i)    = interp1(valid_times, IBTrACS.cma_lat(valid_idx),   datenum(dt), 'linear');
            IBTrACS_Interpolated.cma_lon(i)    = interp1(valid_times, IBTrACS.cma_lon(valid_idx),   datenum(dt), 'linear');
            IBTrACS_Interpolated.cma_vmax(i)   = interp1(valid_times, IBTrACS.cma_wind(valid_idx),  datenum(dt), 'linear');

            IBTrACS_Interpolated.tokyo_lat(i)  = interp1(valid_times, IBTrACS.tokyo_lat(valid_idx), datenum(dt), 'linear');
            IBTrACS_Interpolated.tokyo_lon(i)  = interp1(valid_times, IBTrACS.tokyo_lon(valid_idx), datenum(dt), 'linear');
            IBTrACS_Interpolated.tokyo_vmax(i) = interp1(valid_times, IBTrACS.tokyo_wind(valid_idx),datenum(dt), 'linear');
        end

        % ============================================================
        %  mean USA R34/R50/R64 at SAR time (mean of 4 quadrants)
        % ============================================================
        if any(idx_usa_r34)
            q  = usa_r34_q(idx_usa_r34, :);
            qm = mean(q, 2, 'omitnan');
            IBTrACS_Interpolated.usa_r34_mean(i) = ...
                safe_interp(t_usa_r34, qm, datenum(dt), 'linear');
        end
        if any(idx_usa_r50)
            q  = usa_r50_q(idx_usa_r50, :);
            qm = mean(q, 2, 'omitnan');
            IBTrACS_Interpolated.usa_r50_mean(i) = ...
                safe_interp(t_usa_r50, qm, datenum(dt), 'linear');
        end
        if any(idx_usa_r64)
            q  = usa_r64_q(idx_usa_r64, :);
            qm = mean(q, 2, 'omitnan');
            IBTrACS_Interpolated.usa_r64_mean(i) = ...
                safe_interp(t_usa_r64, qm, datenum(dt), 'linear');
        end
        % ============================================================

        fprintf('  ✅ Done\n');
    end

    %% ============================================================
    %  METADATA
    %% ============================================================
    IBTrACS_Interpolated.tc_name   = tc_name;
    IBTrACS_Interpolated.case_name = case_name;
    IBTrACS_Interpolated.source    = 'IBTrACS v04r01 (interpolated at SAR times)';
    if any(valid_idx)
        IBTrACS_Interpolated.sshs_lmi = max(IBTrACS.sshs(valid_idx));
    end

    IBTrACS_Interpolated.usa_r34_mean_units = 'same as IBTrACS usa_r34 (nm by default)';
    IBTrACS_Interpolated.usa_r50_mean_units = 'same as IBTrACS usa_r50 (nm by default)';
    IBTrACS_Interpolated.usa_r64_mean_units = 'same as IBTrACS usa_r64 (nm by default)';
    IBTrACS_Interpolated.usa_radii_agg      = 'mean over 4 quadrants (NE,SE,SW,NW), then time-interpolated';

    IBTrACS.tc_name   = tc_name;
    IBTrACS.case_name = case_name;
    IBTrACS.source    = 'IBTrACS v04r01 (raw — no interpolation)';

    SAR.tc_name    = tc_name;
    SAR.case_name  = case_name;
    SAR.num_files  = num_files;
    SAR.r_grid     = r_grid;
    SAR.dr_km      = CFG.dr_polar;
    SAR.r_max_km   = CFG.rmax_polar;
    SAR.daz_deg    = CFG.daz_polar;
    SAR.RadialProfile_rows   = 'Row 1 = mean, Rows 2-9 = 8 sectors (0-45, 45-90, ..., 315-360)';
    SAR.CoverageProfile_rows = 'Row 1 = mean, Rows 2-9 = 8 sectors (coverage %)';
    SAR.center_source = 'SAR snapshot center';

    SAR.ProfileStatics_dims    = 'ProfileStatics_* are 4 × 5 × num_files (quadrants only)';
    SAR.ProfileStatics_rows    = {'NE', 'SE', 'SW', 'NW'};
    SAR.ProfileStatics_columns = {'Vmax', 'Rmax', 'R34', 'R50', 'R64'};

    SAR.FieldStatics_dims   = 'FieldStatics is 2 × num_files × 5';
    SAR.FieldStatics_rows   = {'Vmax', 'Rmax'};
    SAR.FieldStatics_fields = {'wsp_sar_masked (raw)', 'wsp_etopo', ...
                               'wsp_buffered', 'wsp_3sigma', 'wsp_patch'};

    % >>> ADDED: metadata for p99 array
    SAR.FieldStatics_p99_dims   = 'FieldStatics_p99 is 2 × num_files × 5';
    SAR.FieldStatics_p99_rows   = {'Vmax_p99', 'Rmax_p99'};
    SAR.FieldStatics_p99_fields = SAR.FieldStatics_fields;
    SAR.FieldStatics_p99_note   = ['Vmax_p99 = 99th-percentile wind over valid ', ...
        '(non-NaN) pixels; Rmax_p99 = mean distance from center to pixels with ', ...
        'wsp >= Vmax_p99'];
    % >>> END ADDED

    SAR.wsp_3sigma_dims = 'cell {num_files,1}; native 2-D single grid';
    SAR.wsp_patch_dims  = 'cell {num_files,1}; native 2-D single grid';

    SAR.CFG_used = CFG;

    %% ============================================================
    %  SAVE — 3 branches
    %% ============================================================
    fprintf('\n[5/5] Saving MAT file...\n');
    mat_file_path = fullfile(mat_folderPath, [tc_name '.mat']);
    fprintf('      Path: %s\n', mat_file_path);

    save(mat_file_path, 'IBTrACS', 'IBTrACS_Interpolated', 'SAR', '-v7.3');

    fprintf('\n✅ SAVED: %s\n', mat_file_path);
    fprintf('   Branches: IBTrACS | IBTrACS_Interpolated | SAR\n');
    fprintf('   Snapshots: %d / %d processed\n\n', num_files, num_files);
end


% ================================================================
%  LOCAL HELPERS
% ================================================================
function q = normalise_quad(raw, n_time)
% Accept any of: n_time x 4, 4 x n_time, 1 x 4 x n_time, n_time x 1 x 4, []
    q = nan(n_time, 4);
    if isempty(raw), return; end

    if numel(size(raw)) == 3
        raw = squeeze(raw);                 % collapse singleton dims
    end

    if isequal(size(raw), [n_time 4])
        q = raw;
    elseif isequal(size(raw), [4 n_time])
        q = raw.';
    elseif isequal(size(raw), [n_time 1])
        q(:,1) = raw;                       % only one quadrant available
    elseif isvector(raw) && numel(raw) == n_time
        q(:,1) = raw(:);
    else
        error('normalise_quad: unexpected size %s for n_time=%d', ...
              mat2str(size(raw)), n_time);
    end
end

function v = safe_interp(t, y, tq, method)
    if nargin < 4, method = 'linear'; end
    v = NaN;
    good = ~isnan(y);
    if nnz(good) < 2, return; end
    t = t(good); y = y(good);
    if tq < min(t) || tq > max(t), return; end
    v = interp1(t, y, tq, method);
end
