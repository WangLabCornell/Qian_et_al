clear all;
%close all;

%% list all files
traces = traces_from_excel("240318_Shuming_Pol2 stall torque.xlsx", "traceFolder");
trace_title = ' Pol II transcription';

file_index = 1:length(traces);
file_index = 1;

% Starting sign of torque from force
% +1 for (+) torque downstream
% -1 for (-) torque upstream
torque_sign = 1;
trace_data.torque_sign = torque_sign;

%% Define which tasks to use
% Define some named taskes to use later in the code:
% These define the task list index for specific tasks.
% All the task list index numbers have a +1 since MATLAB starts counting at 1 instead of 0.
task_plot_list = [6:11] + 1; % List of tasks to plot versus time.
task_torque_list = [7,8,9,11] + 1; %Tasks to plot the torque
task_height_finder_initial = 3+1; % Task for the initial height finder.
task_height_finder_final = 12+1; % Task for the final height finder.
task_torque_end = 11+1; % Task for the final torque level check.
task_initial_activity = 7+1; % Task to look for the replisome intitial activity where the DNA changes chirality from (-) to (+) supercoiled.
task_stall = 9+1; % Task for stalling of polII

% Plotting options
linewidth = 1;
movemean_step = 500; % size of the moving window (for showing data)

%% For loop over traces to process
for i = file_index
    
    % Import data from all tasks in a processed trace:
    % The resulting variable 'tasks' is a cell arrary of structs
    % tasks{i}.name File name of the task
    % tasks{i}.data Table of the data
    trace_data.name = traces(i).name;
    trace_data.folder = traces(i).folder;
    trace_data.index = i;
    trace_data.notes = "";
    tasks = import_trace(strcat(trace_data.folder,'\',trace_data.name));

    %Parse file and folder names to extract date, trace number, and user.
    %For this data set user name and date are in the name
    folder_parts = strsplit(trace_data.folder,"\");
    folder_parts = strsplit(folder_parts{end},"_");
    trace_name_parts = strsplit(trace_data.name,"_");
    trace_data.user = folder_parts{2};
    trace_data.date = folder_parts{1};
    trace_data.number = trace_name_parts{2}; 
    %Create a unique output filename for this trace.
    trace_data.output_file_name = sprintf("%03d_%s_%s_%s",trace_data.index,trace_data.date,trace_data.number,trace_data.user);

    time0 = tasks{task_plot_list(1)}.data.times(1);

    % Loop over all tasks to add new data columns
    % Adding:
    %   - torque from force
    for task = 1:length(tasks)
        tasks{task}.data.TorqueFromForce = torque_from_force_transcription(movmean(tasks{task}.data.FzpN, movemean_step), torque_sign);
    end

    % Get the locations of the torque sign flips:
    % Make sure we are at least 10 turns in.
    index_crop = find(tasks{task_initial_activity}.data.Turn <= -10, 1, "first");
    [M, I] = max(movmean(tasks{task_initial_activity}.data.Zextnm(index_crop:end),100));
    trace_data.fork_position.flip_times = tasks{task_initial_activity}.data.times(I+index_crop -1);
    trace_data.fork_position.extension_max = M;
   
    %% Start the main plot
    %Get the x limits for the time figures.
    xlimits = [0, tasks{task_plot_list(end)}.data.times(end) - time0];
    
    f_main = figure(1);
    clf("reset");
    f_main.Position = [10 10 1800 1500];

    % Define subplots and link x-axis:
    ax_turns = subplot(6,4, [1,2,3]);
    ax_extension = subplot (6,4, [5,6,7]);
    ax_force = subplot (6,4, [9,10,11]);
    ax_replicated_bp = subplot (6,4, [13,14,15]);
    ax_torque = subplot (6,4, [17,18,19,21,22,23]);
 
    linkaxes([ax_turns, ax_extension, ax_force, ax_replicated_bp, ax_torque],'x');

    %% plot turns
    subplot(ax_turns);
    %title(strcat(trace_data.output_file_name," (",trace_title,")"), 'Interpreter', 'none','FontSize', 18);
    title(trace_data.output_file_name, 'Interpreter', 'none','FontSize', 18);
    hold on
    for task = task_plot_list
        plot (tasks{task}.data.times - time0, tasks{task}.data.Turn, 'Linewidth', linewidth);
    end
    plot_xlines(trace_data.fork_position.flip_times - time0);
    ylabel('Turns');
    xlim(xlimits);

    %% plot extension
    subplot(ax_extension);
    hold on
    for task = task_plot_list
        plot (tasks{task}.data.times - time0, tasks{task}.data.Zextnm, 'Linewidth', linewidth);
    end
    plot_xlines(trace_data.fork_position.flip_times - time0);
    ylabel('Extension (nm)');
    xlim(xlimits);
    
    %% plot force
    subplot(ax_force);
    hold on
    for task = task_plot_list
        plot (tasks{task}.data.times -time0, tasks{task}.data.FzpN, 'Linewidth', linewidth);
    end
    
	plot_xlines(trace_data.fork_position.flip_times - time0);
    ylabel('Force (pN)');
    xlim(xlimits);
    
    %% plot movemean force
    subplot(ax_force);
    hold on
    for task = task_plot_list
        plot (tasks{task}.data.times -time0, movmean(tasks{task}.data.FzpN, movemean_step),'color',[0.7 0.7 0.7]);
    end
    [stall_force_for_display, index] = max(movmean(tasks{task_stall}.data.FzpN, movemean_step));
    time_force_max = tasks{task_stall}.data.times(index) ;
    plot(time_force_max- time0,stall_force_for_display,'.',"MarkerSize",20, "Color","Red")
    text(25,1.4,['stall force: ' num2str(stall_force_for_display)  'pN'], 'FontSize', 20);
    
    %% plot torque 
    subplot (ax_torque)
    hold on
    for task = task_plot_list
        plot(tasks{task}.data.times - time0,tasks{task}.data.TorqueFromForce, 'Linewidth', linewidth);
    end
    % Find the max torque from the stalling step.
    [max_torque, index] = max(torque_sign*tasks{task_stall}.data.TorqueFromForce);%determine the stall torque value
    torque_for_display = max_torque * torque_sign;
    time_torque_max = tasks{task_stall}.data.times(index) ;
    plot(time_torque_max- time0,torque_for_display,'.',"MarkerSize",20, "Color","Red")

	plot_xlines(trace_data.fork_position.flip_times - time0);
    grid on
    xlim(xlimits);
    xlabel('Time (s)');
    ylabel('Torque (pNnm)');
    text(25,40,['Stall Torque: ' num2str(torque_for_display)  'pN.nm'], 'FontSize', 20);
    
    %% plot height finder
    subplot (6,4,8)
    plot(tasks{task_height_finder_initial}.data.TrapZnm, tasks{task_height_finder_initial}.data.Zsignal);
    hold on
    plot(tasks{task_height_finder_final}.data.TrapZnm, tasks{task_height_finder_final}.data.Zsignal);
    %Align the height finder before and after experiment.
    [height_fit_initial, zsignal_fit_initial] = fit_height_finder(tasks{task_height_finder_initial}.data.TrapZnm, tasks{task_height_finder_initial}.data.Zsignal);
    [height_fit_end, zsignal_fit_end] = fit_height_finder(tasks{task_height_finder_final}.data.TrapZnm, tasks{task_height_finder_final}.data.Zsignal);
    surface_drift = height_fit_end - height_fit_initial;
    plot(tasks{task_height_finder_final}.data.TrapZnm - surface_drift, tasks{task_height_finder_final}.data.Zsignal);
    title('height finder');
    ylim([0 0.25]);
    xlabel('Trap height (nm)');
    legend('Initial','End',['End shifted' char newline 'Drift:' num2str(surface_drift,'%0.0f') 'nm']);
    ylabel('Zsignal');

   
    task = task_torque_end;
    j = 1;
    if(find(task_plot_list == task))
        [~,j] = find(task_plot_list == task);
    end
      [time_trimmed,toque_trimmed] = trim_torque(tasks{task}.data.times, tasks{task}.data.Turn, tasks{task}.data.TorqueFromForce);  
    plot (time_trimmed - time0, movmean(toque_trimmed, movemean_step), 'Linewidth', linewidth, 'Color', default_MATLAB_colors(j));
    xlabel('Time(s)');
    ylabel('Torque (pNnm)');
    title('Denature torque at end');
    ylim([-30 30]);
    
    saveas(f_main, strcat(trace_data.output_file_name,".png") );
    saveas(f_main, strcat(trace_data.output_file_name,".fig") );
    save(strcat(trace_data.output_file_name,"_trace_data.mat"),"trace_data");
    
    %% Special plots for further ananlysis
    %polII_stall_torque(tasks,trace_data,task_plot_list, task_stall, movemean_step)
    polII_stall(tasks,trace_data)
end


%% Functions
% ####################################################################
% ####################################################################
% ####################################################################
% ####################################################################
% ####################################################################
function polII_stall(tasks,trace_data)

    % The core of this function uses Yifeng's LUT DNA HAT curve calculator
    % to convert force and extension to PolII position on the template.
    movmean_step = 500;
    movmean_step_torque = 500;
    decimate_step = 50;
    %task = 9+1;
    
    % Collect all the data from the steps to plot.
    task_list = [6,7,8,9] + 1;
    task_index= [];
    time_raw = [];
    force_raw = [];
    ext_raw = [];
    torque_raw = [];
    turns_raw = [];
    for k = 1:length(task_list)
        task = task_list(k);
        task_index = [task_index; length(time_raw)+1];
        time_raw = [time_raw; tasks{task}.data.times];
        force_raw = [force_raw; tasks{task}.data.FzpN];
        ext_raw = [ext_raw; tasks{task}.data.Zextnm];
        torque_raw = [torque_raw; tasks{task}.data.TorqueFromForce];
        turns_raw = [turns_raw; tasks{task}.data.Turn];
    end
    

    time0 = time_raw(task_index(4)); %Shift the time = 0 to the stalling step.
    time_raw = time_raw - time0;

    %Smoothe the variables
    force_smooth = movmean(force_raw,movmean_step);
    ext_smooth = movmean(ext_raw,movmean_step);
    torque_smooth = movmean(torque_raw,movmean_step_torque);


    time = decimate_with_average(time_raw, decimate_step);
    force =  decimate_with_average(force_smooth, decimate_step);
    extension = decimate_with_average(ext_smooth, decimate_step);
    torque = decimate_with_average(torque_smooth, decimate_step);
    turns = decimate_with_average(turns_raw, decimate_step);
    
    % Find the stall torque as the max torque during the stalling step.
    [trace_data.stall_torque.torque, trace_data.stall_torque.index] = max(torque_smooth(task_index(4):end));
    trace_data.stall_torque.time = time_raw(trace_data.stall_torque.index+task_index(4)-1);

    %Find the time of the peak in extension in the sign flip step

    [ext_max, index_max] = max(ext_smooth(task_index(3):(task_index(4)-1)));
    index_max = index_max+task_index(3)-1;
    time_max = time_raw(index_max);

    %initialize the RNAP position array
    position = zeros(size(force));
    torque_lookup = zeros(size(force));

    initial_turns = tasks{task}.data.Turn(1)
    a = 0.338; % [nm/bp] the counter length per bp of dsDNA.
    N_initial = 5000; % number of basepairs in template.

    % Fit the peak of the extension when the PolII goes over the hat curve
    % to get the initial turns.  This accounts for PolII activity prior to
    % grabbing the tether.
    index_max = find(time>=time_max,1,"first")
    force_max = force_smooth(index_max)
    %ext_max = extension(index_max)
    [sigma, normalized_extension, ~] = normalized_extension_LUT(force_max);
    n_max = ext_max/(a*normalized_extension(1))

    turn_offset = turns(index_max);

     figure(3)
        n = (n_max - (turns(index_max) - turn_offset)*10.5 )./(sigma+1);
        extension_LUT = a * n .* normalized_extension;
        position_max = interp1(extension_LUT', n, ext_max);
        plot(n,extension_LUT)

        hold on
        sigma = -flip(sigma);
        normalized_extension = flip(normalized_extension);
        n = (n_max - (turns(index_max) - turn_offset)*10.5)./(sigma+1);
        extension_LUT = a * n .* normalized_extension;
        plot(n,extension_LUT, "Color", default_MATLAB_colors(1));
        yline(ext_max, '--')
         
        plot(position_max,ext_max,'.',"MarkerSize",30)
        hold off
      
   size(force)
    parfor i = 1:length(force)
        %Make the nomalized extension vs sigma LUT for this force data
        %point
        [sigma, normalized_extension, torque_LUT] = normalized_extension_LUT(force(i));

        %If on negative side flip the lookup table to the negative side
        if(time(i) < time_max)
            sigma = -flip(sigma);
            normalized_extension = flip(normalized_extension);
            torque_LUT = -flip(torque_LUT);
        end
        
        %find number of bp transcribed from LUT
        % interp1(x,y,X) given y = y(x) find y(X)
        n = (n_max + (turns(i) - turn_offset)*10.5)./(sigma+1);
        extension_LUT = a * n .* normalized_extension;
        
        try
            [max_ext_LUT, max_index] = max(extension_LUT);
            if extension(i) > max_ext_LUT
                out = [n(max_index), torque_LUT(max_index)];
            else
                out = interp1(extension_LUT', [n, torque_LUT], extension(i));
            end
            
            position(i) = out(1);
            torque_lookup(i) = out(2);
            exit = false;
        catch
            position(i) = NaN;
            torque_lookup(i) = NaN;
            exit = true
        end
        if isnan(position(i))
            exit = true
        end
        %Debugging plot.
        % figure(3)
        % plot(n,extension_LUT)
        % hold on
        % sigma = -flip(sigma);
        % normalized_extension = flip(normalized_extension);
        % n = (n_max + (turns(i) - turn_offset)*10.5)./(sigma+1);
        % extension_LUT = a * n .* normalized_extension;
        % plot(n,extension_LUT, "Color", default_MATLAB_colors(1));
        % yline(ext_max, '--')
        % yline(extension(i))
        % plot(position(i),extension(i),'.',"MarkerSize",30)
        % hold off
        % xlim([4000,5400])
        % ylim([0, 1400])
        % pause(0.01)
        % %End debugging plot
        % if exit
        %     answer = inputdlg("Continue?");
        % end
    end
    
    position = N_initial - position;
    index = find(time >= trace_data.stall_torque.time, 1, "first");
    trace_data.stall_torque.position = position(index);

    sigma = (n_max - position)/10.5 - (turns-turn_offset);

    fig = figure(2);
    clf
    fig.Position = [1813 198 1400 1200];
    ax_position = subplot(4,1,1:2);
    ax_ext = subplot(4,1,3);
    ax_torque = subplot(4,1,4);
    linkaxes([ax_torque,ax_position, ax_ext],'x')


    subplot(ax_ext)
    plot(time,extension, "LineWidth",2);
    hold on 
    plot(time_raw,ext_smooth);
    plot(time_max,ext_max,'.',"MarkerSize",30)
    yline(ext_max,'--', "LineWidth",2);
    xline(time_max,'--', "LineWidth",2);
    ylabel("Extension (nm)");

    
    subplot(ax_position)
    plot(time,position, "LineWidth",2)
   
    title(trace_data.output_file_name, 'Interpreter', 'none','FontSize', 14);
    ylim([-300, 1200]);
    xlim([-120,120])
    text(trace_data.stall_torque.time+2,trace_data.stall_torque.position + 100, sprintf("Stall position = %0.0f nt", trace_data.stall_torque.position));
    xline(trace_data.stall_torque.time,'--', "LineWidth",2);
    xline(time_max,'--', "LineWidth",2);
    yline(trace_data.stall_torque.position,'--', "LineWidth",2);
    ylabel("RNAP position (nt)")
    

    subplot(ax_torque)
    hold on
    plot(time_raw,torque_raw,"Color", lighten_plot_color(default_MATLAB_colors(1),0.8));
    plot(time,torque,"LineWidth",2, "Color",default_MATLAB_colors(1) );
    plot(time,torque_lookup,"LineWidth",2, "Color","Black" );
    
    xline(trace_data.stall_torque.time,'--', "LineWidth",2);
    yline(trace_data.stall_torque.torque,'--', "LineWidth",2);
    text(trace_data.stall_torque.time+2,trace_data.stall_torque.torque + 5, sprintf("Stall time = %0.1f s\nStall torque = %0.1f pNnm", trace_data.stall_torque.time,trace_data.stall_torque.torque));
    ylim([-8, 18])
    ylabel("Torque (pNnm)")
    xlabel("Time from start of stalling step (s)")

    set(findall(gcf,'-property','FontSize'),'FontSize',18)

    % Save data
    saveas(gcf, strcat(trace_data.output_file_name,"_position.png") );
    time = time';
    position = position';
    torque = torque';
    torque_lookup = torque_lookup';
    T = table(time, position, torque, torque_lookup);
    writetable(T,strcat(trace_data.output_file_name,"_stall.txt"));
    save(strcat(trace_data.output_file_name,"_trace_data.mat"),"trace_data");
    
end

function polII_stall_torque(tasks,trace_data,task_plot_list, task_stall, movmean_step)
    
    task = task_plot_list(4);
    time = tasks{task}.data.times;
    time = time - time(1);
    force_raw = tasks{task}.data.FzpN;
    force_50 = movmean(tasks{task}.data.FzpN,50);
    force = movmean(tasks{task}.data.FzpN,movmean_step);
   
    torque_raw = torque_from_force_transcription(tasks{task_stall}.data.FzpN, trace_data.torque_sign);
    torque_50 = torque_from_force_transcription(movmean(tasks{task_stall}.data.FzpN, 50), trace_data.torque_sign);
    torque = torque_from_force_transcription(movmean(tasks{task_stall}.data.FzpN, movmean_step), trace_data.torque_sign);
    
    [trace_data.stall_torque.torque, trace_data.stall_torque.index] = max(trace_data.torque_sign*torque);
    trace_data.stall_torque.torque = trace_data.stall_torque.torque*trace_data.torque_sign;
    trace_data.stall_torque.time = time(trace_data.stall_torque.index);
    trace_data.stall_torque.force = force(trace_data.stall_torque.index);
    trace_data.stall_torque.total_time = time(end);

    
    fig = figure(2);
    clf
    fig.Position = [10 10 1400 1200];
    
    ax_force = subplot(2,1,1);
    ax_torque = subplot(2,1,2);
    linkaxes([ax_torque,ax_force],'x')

    subplot(ax_force)
    hold on
    plot(time,force_raw,"Color", lighten_plot_color(default_MATLAB_colors(1),0.8));
    plot(time,force_50,"Color", lighten_plot_color(default_MATLAB_colors(1),0.3));
    plot(time,force,"LineWidth",2, "Color", default_MATLAB_colors(1));
    xline(trace_data.stall_torque.time,'--', "LineWidth",2);
    yline(trace_data.stall_torque.force,'--', "LineWidth",2);
    text(trace_data.stall_torque.time+2,trace_data.stall_torque.force + 0.2, sprintf("Stall force = %0.2f pN", trace_data.stall_torque.force));
    ylim([0,1.5])
    ylabel("Force (pN)")
    title(trace_data.output_file_name, "Interpreter","none");
    legend(["Raw data (100 Hz)", "Smoothed 0.5 second window", sprintf("Smoothed %0.1f second window", movmean_step/100)])

    subplot(ax_torque)
    hold on
    plot(time,torque_raw,"Color", lighten_plot_color(default_MATLAB_colors(1),0.8));
    plot(time,torque_50,"Color", lighten_plot_color(default_MATLAB_colors(1),0.3));
    plot(time,torque,"LineWidth",2, "Color",default_MATLAB_colors(1) );
    xline(trace_data.stall_torque.time,'--', "LineWidth",2);
    yline(trace_data.stall_torque.torque,'--', "LineWidth",2);
    text(trace_data.stall_torque.time+2,trace_data.stall_torque.torque + 2, sprintf("Stall time = %0.2f s\nTime after stall = %0.2f s\nStall torque = %0.2f pNnm", trace_data.stall_torque.time,trace_data.stall_torque.total_time - trace_data.stall_torque.time, trace_data.stall_torque.torque));
    ylim([-18, 0])
    xlim([0,200])
    ylabel("Torque (pNnm)")

    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    saveas(fig, strcat(trace_data.output_file_name,"_polII_stall.png") );
    save(strcat(trace_data.output_file_name,"_trace_data.mat"),"trace_data");

end

function data = importfile(filename, dataLines)
%IMPORTFILE Import data from a text file

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 69);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["V_AV", "V_DV", "V_PV", "V_SV", "V_LV", "V_RV", "V_sumV", "V_xV", "V_yV", "V_pxV", "V_pyV", "V_pzV", "dts", "Vpx_AOV", "Vpy_AOV", "Vpz_AOV", "scanbacklog", "LockInPhase", "TrapXnm", "TrapYnm", "TrapZnm", "V_AV1", "V_DV1", "V_PV1", "V_SV1", "V_LV1", "V_RV1", "V_sumV1", "V_xV1", "V_yV1", "PpowermW", "SpowermW", "anglecompensatedPpowerW", "anglecompensatedSpowerW", "anglecompensatedV_L", "anglecompensatedV_R", "anglecompensatedV_sumV", "anglecompensatedV_xV", "anglecompensatedV_yV", "Psirad", "LaserpowerbeforeobjectiveW", "E_in", "E_out", "DeltaE", "Xsignal", "Ysignal", "Zsignal", "TrapXAOnm", "TrapYAOnm", "TrapZAOnm", "FzpN", "Zextnm", "TaupNnm", "Turn", "FxpN", "Xextnm", "FypN", "Yextnm", "times", "DeltaZsignal", "DeltaZbeadnm", "Zbeadnm", "DeltaTrapZnm", "DeltaXsignal", "Xbeadnm", "DeltaYsignal", "Ybeadnm", "StopCondition", "VarName69"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "VarName69", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "VarName69", "EmptyFieldRule", "auto");

% Import the data
data = readtable(filename, opts);
end

function tasks = import_trace(folder)
    files = dir(folder);
    files = files(3:end);
    tasks = {};
    for j = 1:length(files)
        tasks{j} = struct('name',files(j).name,'data',importfile(strcat(files(j).folder,'\',files(j).name)));
    end
end

function [time_trimmed, torque_trimmed] = trim_torque(time, turns, torque)
    % Function to trim the beginning of the torque data
    % The angle is changing just a bit at the beginning
    bool =(mod(turns,1) == 0);
    idx = find(bool);
    if ~isempty(idx)
        starting_index = idx(1);
    else
        starting_index = 1;
    end
    torque_trimmed = torque(starting_index:end);
    time_trimmed = time(starting_index:end);
end

function torque_from_force = torque_from_force_transcription(force, torque_sign)
    %Get the torque from the force using a pre-calibrated conversion
    %of naked DNA in transcription buffer.  See:
    % 240430_ji33_Torque from Force in transcription buffer.pptx
    A = 12.07;
    B = 0.6409;
    F_melt = 0.8202;
    
    if torque_sign > 0
        torque_from_force = A*(abs(force).^B);
    else
        torque_from_force = - heaviside(F_melt - force) .* A.*(abs(force).^B) - heaviside(force - F_melt)* A*(abs(F_melt).^B);
    end
end

function [x_fit, y_fit] = fit_height_finder(x,y)
    zsignal = 0.15; %Z signal to align height finder
    zrange = 0.01; %range for linear fit to smooth noise.
    slice_index = find( y < (zsignal+zrange) & y > (zsignal-zrange));
    y_slice = y(slice_index);
    x_slice = x(slice_index);
    fit = polyfit(y_slice,x_slice,1);
    x_fit = polyval(fit,zsignal);
    y_fit = zsignal;
end

function plot_xlines(xlines)
    for i= 1:length(xlines)
        xline(xlines(i),'--');
    end
end

function plot_color_light = lighten_plot_color(plot_color,percent)
    %Lighten plot color by mixing with white.
    plot_color_light = plot_color + percent*(1-plot_color);
end

function color = default_MATLAB_colors(index)
    if index < 1
        index = 1;
    end
    index = mod(index-1,7)+1;
    plot_colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
    color = plot_colors{index};
end

function out = decimate_with_average(y,num_points)
    out = mean(reshape([y(:); nan(mod(-numel(y),num_points),1)],num_points,[]),"omitnan");
end

function [sigma, normalized_extension, torque] = normalized_extension_LUT(applied_force)

% Parameters used for the MMS model.
L_p = 43;
K_0 = 1200;
kBT = 4.09;

%%
%caluclate buckling transition
force = [0.02215 0.04122 0.10551 0.23913 0.25 0.5 1 2 5];%pN
Ceff = [20.69 28.07 42.08 65.33 68 79.81 92.39 94.73 100.7];%Xiang PRL, nm
% figure(1)
% semilogx(force,Ceff)
% hold on

[fitresult1, gof1] = Ceff_force_LUT(force, Ceff);


%torque_slope = fitresult1(applied_force) * kBT * 2 * pi / L0;%pNnm/turn
torque_slope = fitresult1(applied_force) * kBT * 2 * pi / (0.338*10.5);%pNnm/turn

a = 13;%fit from using 6481 bp DNA
n = 0.6849;%fit from using 6481 bp DNA
buckling_transition_sigma = a * applied_force^n/ torque_slope;

%%
%parameters used for hat curve generation
%pre-buckling region: a0*(sigma)^2+1
force2 = [0.3 0.5 0.75 1 1.5 3];
coefficient_pre_buckling = [-576.2076 -212.4287 -105.7476 -73.9055 -29.9995 -8.9695];
[fitresult2, gof2] = coefficient_pre_buckle_LUT(force2, coefficient_pre_buckling);

%slope of post-buckling region
%post-buckling region: -a1*(𝜎−𝑎2)+𝑎0*𝑎2^2+1
ext_slope_post_buckling = [22.9814 17.2886 14.2171 12.3404 10.2941 7.6653];
[fitresult3, gof3] = slope_post_buckling_LUT(force2, ext_slope_post_buckling);


%hat curve

sigma1 = 0:0.001:buckling_transition_sigma;
factor_of_extension1 = (fitresult2(applied_force) * (sigma1).^2 + 1);

% Find sigma_max when the extension goes to zero
sigma_max = (fitresult2(applied_force) * buckling_transition_sigma^2 + 1 + buckling_transition_sigma*fitresult3(applied_force))/fitresult3(applied_force);

sigma2 = buckling_transition_sigma:0.001:sigma_max;
factor_of_extension2 = (-fitresult3(applied_force) * (sigma2 - buckling_transition_sigma) + fitresult2(applied_force) * buckling_transition_sigma^2 + 1);

factor_of_extension = [factor_of_extension1 factor_of_extension2]';
sigma = [sigma1 sigma2]';

normalized_extension = factor_of_extension * x_MMS(applied_force, kBT, L_p, K_0);
%%
%torque
%torque1 = torque_slope(k) * sigma1 * Lk0;
torque1 = torque_slope * sigma1 ;
torque2 = a * applied_force^n * ones(1,length(sigma2));
torque = [torque1 torque2]';

end

function [fitresult, gof] = Ceff_force_LUT(force, Ceff)
%CREATEFIT(FORCE,CEFF)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : force
%      Y Output: Ceff
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 22-Aug-2023 18:53:28


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( force, Ceff );

% Set up fittype and options.
ft = 'linearinterp';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'Ceff vs. force', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'force', 'Interpreter', 'none' );
% ylabel( 'Ceff', 'Interpreter', 'none' );
% grid on
end

function [fitresult, gof] = coefficient_pre_buckle_LUT(force2, coefficient_pre_buckling)
%CREATEFIT(FORCE2,COEFFICIENT_PRE_BUCKLING)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : force2
%      Y Output: coefficient_pre_buckling
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 30-Aug-2023 15:26:58


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( force2, coefficient_pre_buckling );

% Set up fittype and options.
ft = 'linearinterp';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'coefficient_pre_buckling vs. force2', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'force2', 'Interpreter', 'none' );
% ylabel( 'coefficient_pre_buckling', 'Interpreter', 'none' );
% grid on
end

function [fitresult, gof] = slope_post_buckling_LUT(force2, ext_slope_post_buckling)
%CREATEFIT(FORCE2,EXT_SLOPE_POST_BUCKLING)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : force2
%      Y Output: ext_slope_post_buckling
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 30-Aug-2023 15:31:58


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( force2, ext_slope_post_buckling );

% Set up fittype and options.
ft = 'linearinterp';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'ext_slope_post_buckling vs. force2', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'force2', 'Interpreter', 'none' );
% ylabel( 'ext_slope_post_buckling', 'Interpreter', 'none' );
% grid on
end

function x = x_MMS(F,kBT,Lp,K0)
    % Formula is from MDW, et al, Biophys. J. 72, 1335 (1997).
    % See p. 1342.

    % Returns the normalized extension x := extension / L_0
    c = K0*Lp/kBT;
    f = F/K0;
    
    % a0 + a1 * x + a2 * x**2 + x**3 = 0
    a0 = 1/4-(1/4+(c+1)*f).*((1+f).^2);
    a1 = (1+f).*((1+f)+1/2+2*(c+1)*f);
    a2 = -(2*(1+f)+1/4+(c+1)*f);

    x = zeros(size(a0));
    for i = 1:length(x)
        r = roots([1,a2(i),a1(i),a0(i)]);
        % This function has three roots for force < 0.176 pN and then three real roots after that point
        % The correct root is always the one with the smallest real part.
        x(i) = min(real(r));
    end
end


function traces = traces_from_excel(filename, column_name)
    T = readtable(filename);
    paths = T.(column_name);

    traces = [];
    for k = 1:length(paths)
        path = paths{k};
        if ~isnan(path)
            split_path = split(path,"\");
            traces(end+1).name = string(split_path(end));
            traces(end).folder = string(join(split_path(1:end-1),"\")); 
        end
    end
end