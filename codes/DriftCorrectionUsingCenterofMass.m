% 2D Drift Correction for Telomere Data Using Nucleus COM Drift
% GG 20250711 
%Telomere csv in nm and COM csv in um

clear; clc;

% Load nucleus center of mass drift data (CSV)
[nucleus_file, nucleus_path] = uigetfile('*.csv', 'Select nucleus COM drift CSV file');
nucleus_data = readtable(fullfile(nucleus_path, nucleus_file));

% Load telomere localization data
[telomere_file, telomere_path] = uigetfile('*.csv', 'Select telomere localization CSV file');
telomere_data = readtable(fullfile(telomere_path, telomere_file), 'VariableNamingRule', 'preserve');

% Rename columns for consistency
nucleus_data.Properties.VariableNames{'frames'} = 'frame';

% Convert COM coordinates from µm → nm
nucleus_data.XM_nm = nucleus_data.XM * 1000;  % 1 µm = 1000 nm
nucleus_data.YM_nm = nucleus_data.YM * 1000;

% Calculate relative drift (subtract first frame COM) in nm
firstFrame = min(nucleus_data.frame);
idxFirst   = nucleus_data.frame == firstFrame;

driftX_nm = nucleus_data.XM_nm - nucleus_data.XM_nm(idxFirst);
driftY_nm = nucleus_data.YM_nm - nucleus_data.YM_nm(idxFirst);

% Smooth COM drift trajectories to reduce high-frequency jitter
windowFrames  = 3;          % <-- adjust based on imaging frame-rate & noise
smoothMethod  = 'movmean';

driftX_nm = smoothdata(driftX_nm, smoothMethod, windowFrames);
driftY_nm = smoothdata(driftY_nm, smoothMethod, windowFrames);

% Plot raw COM trajectory to visually inspect drift (smoothed)
figure;
plot(nucleus_data.frame, nucleus_data.XM_nm, '-o', 'DisplayName','X (raw)', 'Color',[0.3 0.7 1]); hold on;
plot(nucleus_data.frame, nucleus_data.YM_nm, '-o', 'DisplayName','Y (raw)', 'Color',[0.6 0.3 0.8]);
plot(nucleus_data.frame, nucleus_data.XM_nm(idxFirst) + driftX_nm, '-', 'DisplayName','X (smoothed)', 'LineWidth',1.6, 'Color',[0.3 0.7 1]);
plot(nucleus_data.frame, nucleus_data.YM_nm(idxFirst) + driftY_nm, '-', 'DisplayName','Y (smoothed)', 'LineWidth',1.6, 'Color',[0.6 0.3 0.8]);

xlabel('Frame', 'FontWeight','bold');
ylabel('Position (nm)', 'FontWeight','bold');
title('Nucleus Center of Mass Trajectory (Raw vs Smoothed)', 'FontWeight','bold');
legend('Location','best'); grid on;

% Save COM plot as .PNG and .fig
[~, nucleus_base, ~] = fileparts(nucleus_file);
save_dir = nucleus_path;

png_name = fullfile(save_dir, [nucleus_base '_COM_trajectory_smoothed.png']);
fig_name = fullfile(save_dir, [nucleus_base '_COM_trajectory_smoothed.fig']);

saveas(gcf, png_name);
savefig(gcf, fig_name);

fprintf('Saved COM plot (raw + smoothed) to:\n  %s\n  %s\n', png_name, fig_name);

% Store drift in table
nucleus_data.driftX_nm = driftX_nm;
nucleus_data.driftY_nm = driftY_nm;

% Join drift with telomere localizations
merged_data = outerjoin(telomere_data, ...
                        nucleus_data(:, {'frame','driftX_nm','driftY_nm'}), ...
                        'Keys','frame','MergeKeys',true);

% Replace any missing drift with zero (for frames without drift data)
merged_data.driftX_nm(isnan(merged_data.driftX_nm)) = 0;
merged_data.driftY_nm(isnan(merged_data.driftY_nm)) = 0;

% Apply drift correction (subtract drift from telomere positions)
merged_data.x_drift_corrected_nm = merged_data.("x [nm]") - merged_data.driftX_nm;
merged_data.y_drift_corrected_nm = merged_data.("y [nm]") - merged_data.driftY_nm;

% Remove rows with any NaNs before saving
nanRows = any(ismissing(merged_data), 2);
if any(nanRows)
    fprintf('Found and removed %d row(s) with NaN values.\n', sum(nanRows));
    merged_data(nanRows, :) = [];
else
    disp('No NaN values found in corrected localization data.');
end

% Save the drift-corrected data to a new CSV
output_filename = fullfile(telomere_path, 'telomere_localizations_drift_corrected.csv');
writetable(merged_data, output_filename);
fprintf('Drift correction complete! File saved as:\n%s\n', output_filename);

% Create new table for ChromTrack with required format
chromtrack_table                       = table();
chromtrack_table.frame                 = merged_data.frame;
chromtrack_table.("x [nm]")            = merged_data.x_drift_corrected_nm;
chromtrack_table.("y [nm]")            = merged_data.y_drift_corrected_nm;
chromtrack_table.("uncertainty [nm]")  = merged_data.("uncertainty [nm]");

chromtrack_filename = fullfile(telomere_path, 'drift_corrected_telomere_localizations_format_for_chromtrack.csv');
writetable(chromtrack_table, chromtrack_filename);
fprintf('Drift correction tables saved to:\n  %s\n  %s\n', output_filename, chromtrack_filename);

% Apply drift correction to raw image stack for visual check (optional)
applyStack = questdlg('Apply drift correction to raw image stack for visual inspection?', ...
                      'Drift Correction to Image Stack', 'Yes','No','Yes');
if strcmp(applyStack,'Yes')
    % Select the raw image stack (multi-page TIFF)
    [stack_file, stack_path] = uigetfile({'*.tif;*.tiff','TIFF Stacks (*.tif, *.tiff)'}, ...
                                         'Select raw image stack');
    stack_fullpath = fullfile(stack_path, stack_file);
    info    = imfinfo(stack_fullpath);
    nFrames = numel(info);

    % Build per-frame inverse shifts (pixels) — convert nm → pixels
    pixel_size_nm = 130.7;  % pixelsize for image stack
    shiftX = -nucleus_data.driftX_nm / pixel_size_nm;
    shiftY = -nucleus_data.driftY_nm / pixel_size_nm;

    % Forward-fill any zeros after the first non-zero (optional)
    for f = 2:nFrames
        if shiftX(f)==0 && shiftX(f-1)~=0
            shiftX(f) = shiftX(f-1);
            shiftY(f) = shiftY(f-1);
        end
    end

    % --- Write drift-corrected stack safely (BigTIFF) ---
    corrected_stack_file = fullfile(stack_path, 'drift_corrected_stack.tif');
    t = Tiff(corrected_stack_file, 'w8');   % 'w8' = BigTIFF, uncompressed

    % Template image to define tags
    I_raw0 = imread(stack_fullpath, 1, 'Info', info);
    bitDepth = info(1).BitDepth;

    tagstruct.ImageLength = size(I_raw0,1);
    tagstruct.ImageWidth  = size(I_raw0,2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = bitDepth;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Compression = Tiff.Compression.None;

    for f = 1:nFrames
        I_raw  = imread(stack_fullpath, f, 'Info', info);
        I_corr = imtranslate(I_raw, [shiftX(f) shiftY(f)], ...
                             'OutputView','same','FillValues',0);

        t.setTag(tagstruct);
        t.write(I_corr);

        if f < nFrames
            t.writeDirectory(); % create new page for next frame
        end

        if mod(f,50)==0 || f==nFrames
            fprintf('Processed frame %d / %d\n', f, nFrames);
        end
    end

    t.close();
    fprintf('Drift-corrected image stack saved as:\n  %s\n', corrected_stack_file);
end
