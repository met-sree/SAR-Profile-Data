function IBTrACS = create_IBTrACKS(filePath, matching_tc_indices)

sids = cellstr(ncread(filePath, 'sid')');
sid = sids{matching_tc_indices};

sshs1 = ncread(filePath, 'usa_sshs');
sshs = sshs1(:, matching_tc_indices);

basin1 = ncread(filePath, 'basin');
basin = string(char(basin1(1,1,matching_tc_indices))) + string(char(basin1(2,1,matching_tc_indices)));

times = ncread(filePath, 'time');
times = times(:, matching_tc_indices);
reference_date = datetime(1858, 11, 17, 0, 0, 0);
datetimes = reference_date + days(times);
status1 = ncread(filePath, 'usa_status');
status = string(squeeze(status1(:, :, matching_tc_indices))');
% USA variables
usa_lat = ncread(filePath, 'usa_lat');
usa_lat = usa_lat(:, matching_tc_indices);

usa_lon = ncread(filePath, 'usa_lon');
usa_lon = usa_lon(:, matching_tc_indices);

usa_wind = ncread(filePath, 'usa_wind');
usa_wind = usa_wind(:, matching_tc_indices);

usa_pres = ncread(filePath, 'usa_pres');
usa_pres = usa_pres(:, matching_tc_indices);

usa_sshs = ncread(filePath, 'usa_sshs');
usa_sshs = usa_sshs(:, matching_tc_indices);

usa_r34 = ncread(filePath, 'usa_r34');
usa_r34 = usa_r34(:, :, matching_tc_indices)';

usa_r50 = ncread(filePath, 'usa_r50');
usa_r50 = usa_r50(:, :, matching_tc_indices)';

usa_r64 = ncread(filePath, 'usa_r64');
usa_r64 = usa_r64(:, :, matching_tc_indices)';

usa_rmw = ncread(filePath, 'usa_rmw');
usa_rmw = usa_rmw(:, matching_tc_indices) * 1.852;

usa_eye = ncread(filePath, 'usa_eye');
usa_eye = usa_eye(:, matching_tc_indices);

% TOKYO (JMA) variables
tokyo_lat = ncread(filePath, 'tokyo_lat');
tokyo_lat = tokyo_lat(:, matching_tc_indices);

tokyo_lon = ncread(filePath, 'tokyo_lon');
tokyo_lon = tokyo_lon(:, matching_tc_indices);

tokyo_wind = ncread(filePath, 'tokyo_wind');
tokyo_wind = tokyo_wind(:, matching_tc_indices);

% CMA variables
cma_lat = ncread(filePath, 'cma_lat');
cma_lat = cma_lat(:, matching_tc_indices);

cma_lon = ncread(filePath, 'cma_lon');
cma_lon = cma_lon(:, matching_tc_indices);

cma_wind = ncread(filePath, 'cma_wind');
cma_wind = cma_wind(:, matching_tc_indices);

% Define variable names
var_names = {'datetimes', 'sid', 'sshs', 'basin', ...
    'usa_lat', 'usa_lon', 'usa_wind', 'usa_pres', ...
    'usa_r34', 'usa_r50', 'usa_r64', 'usa_rmw', 'usa_eye', 'status',...
    'tokyo_lat', 'tokyo_lon', 'tokyo_wind', ...
    'cma_lat', 'cma_lon', 'cma_wind', ...
    'matching_tc_indices'};

% Create empty structure
IBTrACS = struct();

% Assign each variable to the structure
for i = 1:length(var_names)
    IBTrACS.(var_names{i}) = eval(var_names{i});
end

end
