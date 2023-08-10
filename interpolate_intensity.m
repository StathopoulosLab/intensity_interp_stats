function [A_curve, norm_A, p, tbl, stats, m] = interpolate_intensity(tf_output_file)
%INTERPOLATE_INTENSITY Interpolates the values in a measurement.
%
%   Inputs
%       filename: the filename of an excel file with the data, either in
%       the current directory, or the full path to the file
%
%   Reads in data from an excel file. Data is expected in the form of the
%   first three columns beong replicate one, the next three columns being
%   replicate two etc. For each group of three columns in a replicate, the
%   first column is expected to be the x or position, the second column is
%   the intensity of channel one at those x positions, and the third column
%   is the intensity of channel two at those x positions. The data is then
%   averaged across replicates.

% Define step size of interpolation
x_intrp = (0:1/999:1)';

% Get filename to file
filename = get_filename;

% Get names of sheets in excel file
sheets = sheetnames(filename);

% Initilize variables
data = cell(size(sheets, 1), 1);
intrp_data = cell(size(sheets, 1), 1);
mean_data = cell(size(sheets, 1), 1);
A_curve = cell(size(sheets, 1), 1);
norm_A = cell(size(sheets, 1), 2);

% For each sheet
for i = 1:size(sheets, 1)
    % Read in values from sheet i
    data{i} = readtable(filename, 'Sheet', i);
    
    % Initialize variables
    intrp_data{i} = zeros(size(x_intrp,1), size(data{i},2));
    mean_data{i} = zeros(size(x_intrp,1), 3);
    A_curve{i} = zeros(size(data{i},2)/3, 2);
    norm_A{i,1} = zeros(size(data{i},2)/3, 1);
    norm_A{i,2} = zeros(size(data{i},2)/3, 1);
    
    % For each data set (assuming data is grouped in columns where each
    % replicate data set is three columns)
    for j = 1:size(data{i},2)/3
        % Save the values where interpolation is to occur in the first
        % column, and then every third column after that (1, 4, 7, ... etc)
        % for each data set
        intrp_data{i}(:,3*j-2) = x_intrp;
        
        % Use the trapezoid rule to find the area under the curve, or total
        % intensity, using the first column as x and the second column as y
        % selecting only data that is not NaN
        A_curve{i}(j,1) = trapz(data{i}{:,3*j-2}(~isnan(data{i}{:,3*j-2})),...
            data{i}{:,3*j-1}(~isnan(data{i}{:,3*j-2})));
    
        % Use the trapezoid rule to find the area under the curve, or total
        % intensity, using the first column as x and the third column as y
        % selecting only data that is not NaN
        A_curve{i}(j,2) = trapz(data{i}{:,3*j-2}(~isnan(data{i}{:,3*j-2})),...
            data{i}{:,3*j}(~isnan(data{i}{:,3*j-2})));
        
        % Normalize total intensity by dividing the total intensity in
        % from column 3 by column 2
        norm_A{i,1}(j,1) = A_curve{i}(j,2) ./ A_curve{i}(j,1);
        
        % Interpolate the data, but only for data that is not NaN. This is
        % done for the first column and second column per replicate. First
        % column is the x or position and second column is the intenisty
        % for channel one.
        intrp_data{i}(:,3*j-1) = interp1(...
            data{i}{:,3*j-2}(~isnan(data{i}{:,3*j-2})),...
            data{i}{:,3*j-1}(~isnan(data{i}{:,3*j-2})),...
            x_intrp,...
            'makima', NaN);

        % Interpolate the data, but only for data that is not NaN. This is
        % done for the first column and third column per replicate. First
        % column is the x or position and third column is the intenisty
        % for channel two.
        intrp_data{i}(:,3*j) = interp1(...
            data{i}{:,3*j-2}(~isnan(data{i}{:,3*j-2})),...
            data{i}{:,3*j}(~isnan(data{i}{:,3*j-2})),...
            x_intrp,...
            'makima', NaN);
    end
    
    % Create a group name or label for each normalized total intensity for
    % statistical analysis
    norm_A{i,2} = repmat(sheets(i,1), [size(norm_A{i,1},1), 1]);
    
    % Save the x interpolated values in the mean_data. This is the same for
    % all samples
    mean_data{i}(:,1) = x_intrp;

    % Take the mean of the replicates. Each replicate is concatenated at
    % the end of the second dimension of the intrp_data array, so the
    % second column and every column +3 after that are the first channel
    % measurement, and the third column and every column +3 after that are
    % the second channel
    mean_data{i}(:,2) = mean(intrp_data{i}(:,2:3:size(data{i},2)),2);
    mean_data{i}(:,3) = mean(intrp_data{i}(:,3:3:size(data{i},2)),2);

    if tf_output_file
        % Open dialog box for saving excel file output, but only if true
        % was inputed
        [file, path] = uiputfile('*.xlsx');
        filename_out = fullfile(path,file);
    
        % Make an excel file with the outputs. Variable names are saved in
        % the first row. Data is added to the row below. Mean data is added
        % many columns down from the interpolated data.
        writematrix(string(data{i}.Properties.VariableNames),...
            filename_out, 'Sheet', sheets(i,1), 'Range', 'A1');
        writematrix(intrp_data{i}, filename_out, 'Sheet', sheets(i,1),...
            'Range', 'A2');
        writematrix(mean_data{i}, filename_out, 'Sheet', sheets(i,1),...
            'Range', 'Z2');
    end
end

% Concatenate the normnalized total intensities and their group labels
A = cat(1, norm_A{:,1});
g = cat(1, norm_A{:,2});

% Perform one-way ANOVA on the normalized total intensities
[p,tbl,stats] = anova1(A, g, 'off');

% Perform a multiple comparision test, currently using Tukey's method, need
% to update to dunnett's test
m = multcompare(stats, 'display', 'off');

% Code for plotting

% figure; h1 = plot(mean_data{1}(:,1),mean_data{1}(:,2),...
% mean_data{2}(:,1),mean_data{2}(:,2),...
% mean_data{3}(:,1),mean_data{3}(:,2));
% set(h1, 'linewidth', 4, 'markersize', 15);
% set(gca, 'Fontsize', 24);
% legend({'\itwt', '\itzfh1', '\itsna'}, 'Fontsize', 24);
% legend(gca, 'boxoff');
% title('HLH54F', 'Fontsize', 36);
% 
% figure; h2 = plot(mean_data{1}(:,1),mean_data{1}(:,3),...
% mean_data{2}(:,1),mean_data{2}(:,3),...
% mean_data{3}(:,1),mean_data{3}(:,3));
% set(h2, 'linewidth', 4, 'markersize', 15);
% set(gca, 'Fontsize', 24);
% legend({'\itwt', '\itzfh1', '\itsna'}, 'Location', 'northwest', 'Fontsize', 24);
% legend(gca, 'boxoff');
% title('Doc2');
end

function path = get_filename
%GET_FILENAME Gets filename using a UI.
% 
%   Outputs
%       path: path to file selected
% 
%   Overview
%       Saves the filename of the selected excel file.
    
    % Use menu to select file
    [name,folder] = uigetfile({'*.xlsx', 'Excel files (*.xlsx)'},...
        'Select one file', 'MultiSelect', 'off');
    
    % Construct full path
    path = fullfile(folder,name);
end