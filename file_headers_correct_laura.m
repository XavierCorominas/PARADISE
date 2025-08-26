
%% CODE to CORRECT HEADERS OF EEG FILES that its names have been manually modified:

%specify your data folder root.
baseDir = '/mnt/projects/PARADISE/PARADISE_1/XEB/EEG_clean/';


%% Run this section on block.No need to modify anything.

% List all .VMK files and change the name for the actual file names
S = dir(fullfile(baseDir, '*.vmrk'));
if isempty(S)
    warning('No .vmrk files found in %s', baseDir);
end

for k = 1:numel(S)
    vmrkFile = fullfile(S(k).folder, S(k).name);
    [~, baseName] = fileparts(vmrkFile);
    newDataFile = baseName + ".eeg";  % target value

    try
        % Keep empty lines (use 'read' instead of 'preserve')
        lines = readlines(vmrkFile, "EmptyLineRule","read");

        % Find DataFile= line (case-insensitive); fallback to line 5 if needed
        idx = find(startsWith(strtrim(lines), "DataFile=", 'IgnoreCase', true), 1, 'first');
        if isempty(idx)
            if numel(lines) >= 5
                idx = 5;
            else
                warning('Skipping (no DataFile line and <5 lines): %s', vmrkFile);
                continue;
            end
        end

        % Optional: capture old value for logging
        oldLine = strtrim(lines(idx));
        if startsWith(oldLine, "DataFile=", 'IgnoreCase', true)
            oldVal = strtrim(extractAfter(oldLine, '='));
        else
            oldVal = '(unknown)';
        end

        % Replace and write
        lines(idx) = "DataFile=" + newDataFile;
        writelines(lines, vmrkFile);

        fprintf('✓ %s : %s  →  %s\n', S(k).name, oldVal, newDataFile);
    catch ME
        warning('Failed to update %s: %s', vmrkFile, ME.message);
    end
end

%% List all .VHDR files and change the name for the actual file names

S = dir(fullfile(baseDir, '*.vhdr'));

if isempty(S)
    warning('No .vhdr files found in %s', baseDir);
end

for k = 1:numel(S)
    vhdrFile = fullfile(S(k).folder, S(k).name);
    [~, baseName] = fileparts(vhdrFile);
    newEEG  = baseName + ".eeg";
    newVMRK = baseName + ".vmrk";

    try
        % Read lines (keep empties)
        lines = readlines(vhdrFile, "EmptyLineRule","read");

        % Find the lines (case-insensitive); fallback to 6/7 if needed
        idxData = find(startsWith(strtrim(lines), "DataFile=",  'IgnoreCase', true), 1, 'first');
        idxMark = find(startsWith(strtrim(lines), "MarkerFile=", 'IgnoreCase', true), 1, 'first');

        if isempty(idxData)
            if numel(lines) >= 6, idxData = 6; else
                warning('Skipping (no DataFile line and <6 lines): %s', vhdrFile);
                continue
            end
        end
        if isempty(idxMark)
            if numel(lines) >= 7, idxMark = 7; else
                warning('Skipping (no MarkerFile line and <7 lines): %s', vhdrFile);
                continue
            end
        end

        % (Optional) log old values
        oldData = strtrim(extractAfter(strtrim(lines(idxData)), '='));
        oldMark = strtrim(extractAfter(strtrim(lines(idxMark)), '='));

        % Replace lines
        lines(idxData) = "DataFile="  + newEEG;
        lines(idxMark) = "MarkerFile="+ newVMRK;

        % Write back
        writelines(lines, vhdrFile);

        fprintf('✓ %s : DataFile %s → %s | MarkerFile %s → %s\n', ...
            S(k).name, oldData, newEEG, oldMark, newVMRK);

    catch ME
        warning('Failed to update %s: %s', vhdrFile, ME.message);
    end
end