function [Te, dens, RateCoeff, State] = ADF11_CX(file)
% Reads charge exchange (CX) rate coefficients from an ADF11-style file.
% Example: [Te, dens, RateCoeff, State] = ADF11_CX('ccd89_ne.dat');

    % Open file
    fileID = fopen(file, 'r');
    if fileID == -1
        error('Could not find charge exchange file: %s', file);
    end

    % Read first line (header)
    tline = fgetl(fileID);
    fprintf('DEBUG: Raw Header Line: %s\n', tline);  % Print the full header line for inspection

    % ✅ Extract numeric values only
    numValues = regexp(tline, '[-+]?\d*\.?\d+', 'match');  % Find all numbers in header

    % ✅ Ensure at least 6 numeric values are present
    if length(numValues) < 6
        error('Malformed header in %s: Expected at least 6 numeric values, found %d', file, length(numValues));
    end

    % Convert extracted numeric strings to doubles
    IZMAX = str2double(numValues{1});
    IDMAXD = str2double(numValues{2});
    ITMAXD = str2double(numValues{3});
    IZ1MIN = str2double(numValues{4});
    IZ1MAX = str2double(numValues{5});  % ✅ Now correctly extracted!

    % ✅ Validate extracted values
    if any(isnan([IZMAX, IDMAXD, ITMAXD, IZ1MIN, IZ1MAX]))
        error('Failed to extract valid numbers from header: %s', tline);
    end

    % Debug print final values
    fprintf('DEBUG: Parsed Header Values: IZMAX=%d, IDMAXD=%d, ITMAXD=%d, IZ1MIN=%d, IZ1MAX=%d\n', ...
        IZMAX, IDMAXD, ITMAXD, IZ1MIN, IZ1MAX);

    % Initialize output arrays
    RateCoeff = zeros(IDMAXD, ITMAXD, IZ1MAX);
    State = zeros(IZ1MAX, 2);

    % Skip separator line
    fgetl(fileID);

    % Read density grid
    dens = fscanf(fileID, '%f', IDMAXD);
    
    % Read temperature grid
    Te = fscanf(fileID, '%f', ITMAXD);

    % Skip metadata separator
    fgetl(fileID);

    % Read rate coefficients for each charge state
    for i = 1:IZ1MAX
        tline = fgetl(fileID);  % Read charge state header
        C = strsplit(strtrim(tline), '/');
        if length(C) < 2
            error('Malformed charge state header in file: %s', file);
        end
        D = strsplit(C{2}, '=');
        
        if length(D) < 2
            error('Could not extract charge state index from line: %s', tline);
        end

        IZ1 = str2double(D{2});  % Extract charge state
        IZ = IZ1 - 1;

        % ✅ Ensure charge state increments properly
        State(i, :) = [IZ, IZ1];

        % ✅ Debug print to verify charge states
        fprintf('DEBUG: Processing Charge State IZ1=%d at index %d\n', IZ1, i);

        % Validate charge state index
        if isnan(IZ1) || IZ1 < 0
            error('Invalid charge state index IZ1=%d in file: %s', IZ1, file);
        end

        % Read rate coefficients (correct dimensions)
        rawRates = fscanf(fileID, '%f', [ITMAXD, IDMAXD]);

        % ✅ Debugging output
        fprintf('DEBUG: Read %d elements for charge state IZ1=%d (Expected: %d x %d = %d)\n', ...
            numel(rawRates), IZ1, IDMAXD, ITMAXD, IDMAXD * ITMAXD);

        if numel(rawRates) ~= (IDMAXD * ITMAXD)
            error('Unexpected rate coefficient matrix size: Expected [%d x %d], got %d elements', ...
                IDMAXD, ITMAXD, numel(rawRates));
        end

        % ✅ Correctly assign data
        RateCoeff(:, :, i) = reshape(rawRates, [ITMAXD, IDMAXD])';

        % Skip separator line
        fgetl(fileID);
    end

    % Close file
    fclose(fileID);

end