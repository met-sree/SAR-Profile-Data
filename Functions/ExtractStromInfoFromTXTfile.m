function [cen_lat, cen_lon, snap_date] = ExtractStromInfoFromTXTfile(folder_path, base_name)
% EXTRACT_STORM_INFO Extract storm information from text file
%
% Syntax:
%   [cen_lat, cen_lon, Vmax, Rmax, snap_date] = ExtractStromInfoFromTXTfile(folder_path, base_name)
%
% Inputs:
%   folder_path - Path to folder containing the text file
%   base_name   - Base name of the text file (without .txt extension)
%
% Outputs:
%   cen_lat    - Storm center latitude (degrees)
%   cen_lon    - Storm center longitude (degrees)
%   Vmax       - Maximum wind speed from all quadrants (kts)
%   Rmax       - First value of RMax range (nmi)
%   snap_date  - Acquisition date as datetime object
%
% Example:
%   [lat, lon, vmax, rmax, date] = ExtractStromInfoFromTXTfile('C:/data/', 'storm_info_2023');

    % Initialize outputs
    cen_lat = NaN;
    cen_lon = NaN;
    Vmax = NaN;
    Rmax = NaN;
    snap_date = NaT; % Not-a-Time for datetime
    
    % Construct full file path
    file_path = fullfile(folder_path, [base_name, '.txt']);
    
    % Check if file exists
    if ~isfile(file_path)
        error('File not found: %s', file_path);
    end
    
    try
        % Open and read the file
        fid = fopen(file_path, 'r');
        if fid == -1
            error('Cannot open file: %s', file_path);
        end
        
        text_content = textscan(fid, '%s', 'Delimiter', '\n');
        fclose(fid);
        
        lines = text_content{1};
        
        % Initialize quadrant wind speeds for Vmax calculation
        quadrant_winds = [];
        
        % Parse each line
        for j = 1:length(lines)
            line = lines{j};
            
            % Extract Storm Center Latitude
            if contains(line, 'Storm Center Latitude:')
                tokens = strsplit(line, ':');
                if length(tokens) >= 2
                    cen_lat = str2double(strtrim(tokens{2}));
                end
                
            % Extract Storm Center Longitude
            elseif contains(line, 'Storm Center Longitude:')
                tokens = strsplit(line, ':');
                if length(tokens) >= 2
                    cen_lon = str2double(strtrim(tokens{2}));
                end
                
            % Extract RMax - first value from range
            elseif contains(line, 'RMax (nmi):')
                tokens = strsplit(line, ':');
                if length(tokens) >= 2
                    range_str = strtrim(tokens{2});
                    % Split the range by '-' and take first value
                    range_parts = strsplit(range_str, '-');
                    if ~isempty(range_parts)
                        Rmax = str2double(strtrim(range_parts{1}));
                    end
                end
                
            % Extract quadrant wind speeds for Vmax
            elseif contains(line, 'Quadrant') && contains(line, 'VMax (kts):')
                tokens = strsplit(line, ':');
                if length(tokens) >= 2
                    wind_speed = str2double(strtrim(tokens{2}));
                    if ~isnan(wind_speed)
                        quadrant_winds = [quadrant_winds, wind_speed];
                    end
                end
                
            % Extract Acquisition Date
            elseif contains(line, 'Acquisition Date:')
                tokens = strsplit(line, ':');
                if length(tokens) >= 2
                    % Get everything after "Acquisition Date:"
                    date_str = strtrim(strjoin(tokens(2:end), ':'));
                    % Remove "UTC" if present
                    date_str = strrep(date_str, ' UTC', '');
                    try
                        snap_date = datetime(date_str, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
                    catch
                        % Try alternative parsing if format doesn't match
                        snap_date = datetime(date_str);
                    end
                end
            end
        end
        
        % Calculate Vmax as the highest value from all quadrants
        if ~isempty(quadrant_winds)
            Vmax = max(quadrant_winds);
        end
        
        % Display warnings for missing values
        if isnan(cen_lat)
            warning('Storm center latitude not found in file: %s', file_path);
        end
        if isnan(cen_lon)
            warning('Storm center longitude not found in file: %s', file_path);
        end
        if isnan(Vmax)
            warning('Maximum wind speed (Vmax) not found in file: %s', file_path);
        end
        if isnan(Rmax)
            warning('RMax not found in file: %s', file_path);
        end
        if isnat(snap_date)
            warning('Acquisition date not found in file: %s', file_path);
        end
        
    catch ME
        % Close file if still open (in case of error)
        if fid ~= -1
            fclose(fid);
        end
        rethrow(ME);
    end
    
    % Display results
    fprintf('Extracted from %s.txt:\n', base_name);
    fprintf('  Latitude:  %.4f°\n', cen_lat);
    fprintf('  Longitude: %.4f°\n', cen_lon);
    fprintf('  Vmax:      %.2f kts\n', Vmax);
    fprintf('  Rmax:      %.2f nmi\n', Rmax);
    
    if ~isnat(snap_date)
        fprintf('  Date:      %s\n', datestr(snap_date, 'yyyy-mm-dd HH:MM:SS'));
    else
        fprintf('  Date:      Not found\n');
    end
end
