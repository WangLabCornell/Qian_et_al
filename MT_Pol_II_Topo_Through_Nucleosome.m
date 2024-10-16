clear;
close all;
clc;
%c = @cmu.colors;
% IMPORTANT! This program only works if there are two forces!

% This program is supposed to choose multiple files, but they should have
% the same protocols;

load 'Surface_Shift_Info.txt';

pickedcolor = [255, 0, 0]/255;
EndingT = 11;
tweezer2 = 0;
tweezer2 = 2;
resuming_index = 1;
smooth_time = 10;
ita = 0.3;
max_plot_height = 1600;
dw = 0.5973;
dh = 0.0414;
Li = 13017;
Lf = 12667;
L_cutoff = 12614;

flag = 0; % To indicate if we choose any files at all, even if we choose one file, this would be 1;
file_index = 0; % Number of files chosen;
key = '1 Wind';
Maximum_allowed_files = 1;
Maximum_allowed_boxes = 80;
Data_points_number = zeros(Maximum_allowed_files,1);
position_and_force = zeros(4,Maximum_allowed_boxes,Maximum_allowed_files);
target_index = [];

essential_columns = {{'DAQ time (s)'}; {'Time (ms)'}; {'tasklist (index)'};{'Magnet Angle (Turns)'};{'Magnet Height (mm)'}};% This is the number of columns essential for our hatcurve;

for ii = 1:Maximum_allowed_files
    position_and_force(2,:,ii) = (1:Maximum_allowed_boxes)-1;
end


while (1) % while not stopped, this loop keep on choosing files
    %% Access data files
    
    Address = 'Y:\Magnetic Tweezers\Jin';
    if tweezer2 == 2
        Address = 'Y:\Magnetic Tweezers\Jin\240702_Pol II Nucleosome Array scTopoII 10 mins Test\_240705_184217_PROCESSED\Data Files';
    end
    
    [filename, pathname] = uigetfile('*.*','Select the data file',Address); % Choose MT data files
    if filename == 0 % If you do not choose a file, then the program proceeds;
        break;
    else flag = 1;
    end
    new_file_name = strrep(filename,'_',' ');
    fileID = fopen ([pathname,filename],'r'); % Open the data file
    file_index = file_index + 1;
    new_folder = [pwd,'\',[filename,' plots']];
    mkdir([filename,' plots']);
    id_Date = str2double(filename(1:6));
    id_Time = str2double(filename(8:13));
    
    %% Access force calibration files
    pathname_Force = pathname(1:(length(pathname)-length('Data Files\'))); % This is to access the force calibration results
    filename_Force = 'Force Cal Results.txt'; % Find the force calibration file
    fileID_Force = fopen([pathname_Force,filename_Force],'r'); % Open the force calibration file
    
    %% Get Title lines (the first line of each data file should be word descriptions, ditch that)
    titleline = fgetl(fileID);
    essential_columns_position = find_essential_columns_by_name( titleline, essential_columns );
    Box_cell_position  = find_boxes( titleline );
    num_of_box_cells = length(Box_cell_position);
    numofboxes = num_of_box_cells/8;
    end_of_array = Box_cell_position(end);
    Box_Ext_positions = (Box_cell_position(1)+6):8:end_of_array;
    % This is the actually num of boxes, not the maximum box index. If the
    % numofboxes = 10, the maximum box index would be 9.
    
    %% Figure out which box to plot
    % The boxes being plot should have the category of "1 Wind", which is
    % the content of 'key'.
    TCindex = (1:numofboxes) - 1; % All of the box indexes
    return_matrix = find_TC_index (fileID_Force, key);
    plotting_index = return_matrix(1,:); % These are the <box indexes> that are good traces
    Start_force = return_matrix(2,:); % This is the first calibrated force
    Jump_force = return_matrix(3,:); % This is the second calibrated force
    position_and_force(1,plotting_index+1,file_index) = 1;
    position_and_force(3,plotting_index+1,file_index) = Start_force;
    position_and_force(4,plotting_index+1,file_index) = Jump_force;
    
    %% Format Spec for reading data files
    good_box_index = plotting_index+1;
    good_box_ext_position = Box_Ext_positions(good_box_index);
    reading_position = [essential_columns_position, good_box_ext_position];
    format_spec = [];
    for ii = 1:end_of_array
        if isempty(find(reading_position==ii))
            format_spec = [format_spec, '%*f'];
        else
            format_spec = [format_spec, '%f'];
        end
    end
    Data_all_cell = textscan(fileID, format_spec, 'MultipleDelimsAsOne',true, 'Delimiter','[;', 'HeaderLines',0);
    %% Close all the files
    fclose (fileID);
    fclose (fileID_Force);
    disp (['File ', num2str(file_index), ': ', filename]);
    
    %% A precaution not to exceed the maximum allowed number of files
    if file_index == Maximum_allowed_files
        break;
    end
end



if flag == 1 % This means we actually choose some files
    %% Find Time spots
    time_for_timespot = Data_all_cell{2}/1000; %This is only used for determining the steps; We use the first chosen file for that.
    task_for_timespot = Data_all_cell{3}; % The task number for the first file selected
    magnet_height_for_subplot = Data_all_cell{5}; % The magnetic height for the first file selected, used for plotting subplot;
    Turns = Data_all_cell{4}; % The magnetic turns for the first file selected, used for plotting subplot;
    Timespot_index = timespot_finding (task_for_timespot);% Time spots given by the first file selected;
    Timespot = time_for_timespot(Timespot_index);
    EndingT = length(Timespot);
    
    
    %%
    for which_file = 1:file_index
        Extension_index = (21+tweezer2):8:4000;
        plotting_index = find (position_and_force(1,:,which_file)==1);
        Time = Data_all_cell{2}/1000;
        dT = mean(diff(Time));
        output = [];
        hat_curve_record = cell(3,1);
        height_record = zeros(3,1);
        width_record = zeros(3,1);
        if isempty(target_index)
            target_index=plotting_index-1;
        end
        for boxes = resuming_index:length(plotting_index)
            %for boxes = 5:5
            if ~isempty(find(target_index==plotting_index(boxes)-1))
                clc;
                disp([num2str(boxes),'/',num2str(length(plotting_index))]);
                output_temp = [id_Date, id_Time, plotting_index(boxes)-1];
                Extension = 1000*Data_all_cell{boxes+length(essential_columns_position)};
                surface_shift = lookup_surface_shift (Surface_Shift_Info, ...
                    id_Date, id_Time, plotting_index(boxes)-1);
                Extension = Extension+surface_shift*1000;
                Extension_notfiltered = Extension;
                Extension = Smoothed_trace (Time, Extension, task_for_timespot, smooth_time);
                
                FIG_y = figure (2*boxes-1);
                pos = [100 100 1800 900];
                set(FIG_y, 'Pos', pos);
                TILE = tiledlayout(4,1, 'Padding', 'none', 'TileSpacing', 'compact');
                S0 = title (TILE, [new_file_name(1:end-4),' box ',num2str(plotting_index(boxes)-1)]);
                S0.FontSize = 20;
                S0.FontWeight = 'bold';
                nexttile
                plot (time_for_timespot, Turns,'k','LineWidth',4);
                axis([0 inf -inf inf])
                ylabel('Turns');
                xlim ([0 Timespot(EndingT)]);
                set(gca,'FontSize',15,'LineWidth',1.5);
                set(gca,'Color', 'None');
                set(gca,'TickDir','out');
                set(gca,'TickLength',[0.015, 0.01]);
                set(gca,'XColor','k', 'YColor','k');
                set(gca, 'Layer', 'top');
                ax_all = gca;
                box(ax_all,'off');
                nexttile(2, [3 1]);
                plot (Time, Extension_notfiltered, 'Color', [0.7,0.7,0.7]);
                hold on;
                plot (Time, Extension, 'Color', [0,0,1], 'LineWidth',4);
                hold on;
                for ii = 1:length(Timespot)-1
                    time_i = Timespot (ii);
                    yy = 0:10:2000;
                    xx = time_i*ones(size(yy));
                    plot(xx,yy,'k--','LineWidth',1.5);
                    hold on;
                end
                xlabel ('Time (s)');
                ylabel ('Extension (nm)');
                xlim ([0 Timespot(EndingT)]);
                ylim ([0 max_plot_height]);
                set(gca,'FontSize',15,'LineWidth',1.5);
                set(gca,'Color', 'None');
                set(gca,'TickDir','out');
                set(gca,'TickLength',[0.015, 0.01]);
                set(gca,'XColor','k', 'YColor','k');
                set(gca, 'Layer', 'top');
                ax_all = gca;
                box(ax_all,'off');
                saveas(gcf,[new_folder,['\Box ', num2str(plotting_index(boxes)-1),'_ext.fig']]);
                saveas(gcf,[new_folder,['\Box ', num2str(plotting_index(boxes)-1),'_ext.png']]);
                saveas(gcf,[new_folder,['\Box ', num2str(plotting_index(boxes)-1),'_ext.emf']]);
                close;
                
                FIG_x = figure (2*boxes);
                pos = [0 100 2250 900];
                set(FIG_x, 'Pos', pos);
                TILE = tiledlayout(2,4, 'Padding', 'none', 'TileSpacing', 'compact');
                S0 = title (TILE, [new_file_name(1:end-4),' box ',num2str(plotting_index(boxes)-1)]);
                S0.FontSize = 20;
                S0.FontWeight = 'bold';
                
                %% Plotting all the good traces
                pos_hat_range = task_for_timespot==7 | task_for_timespot==11;
                neg_hat_range = task_for_timespot==9;
                
                [turns_p, ext_p] = binning_hat(Turns(pos_hat_range), Extension_notfiltered(pos_hat_range));
                [turns_n, ext_n] = binning_hat(Turns(neg_hat_range), Extension_notfiltered(neg_hat_range));
                [turns_p, ext_p, turns_n, ext_n] = restricted_binning (turns_p, ext_p, turns_n, ext_n);
                
                nexttile(1)
                plot (turns_p, ext_p, 'Color', [0.7,0.7,0.7]);
                hold on;
                plot (turns_n, ext_n, 'Color', [0.7,0.7,0.7]);
                hold on;
                xlim ([-70 70]);
                ylim ([0 max_plot_height]);
                ylabel('Extension (nm)');
                xlabel('Turns');
                set(gca,'FontSize',15,'LineWidth',1.5);
                set(gca,'Color', 'None');
                set(gca,'TickDir','out');
                set(gca,'TickLength',[0.015, 0.01]);
                set(gca,'XColor','k', 'YColor','k');
                set(gca, 'Layer', 'top');
                ax_all = gca;
                box(ax_all,'off');
                d_ext_during_harcurve = ext_p-ext_n;
                dEXT_during = mean(d_ext_during_harcurve);
                
                if 1
                    ext_np = 0.5*(ext_n+ext_p);
                    C0_fit = fit5piece_Jin(turns_p, ext_np/1000);
                    hat_curve_record{1} = C0_fit;
                    zz = f_5piece(C0_fit,turns_p);
                    plot (turns_p, ext_np, 'Color', [0 0 0], 'LineWidth', 4);
                    hold on;
                    plot (turns_p,zz*1000,'g--','Linewidth',2);
                    hold on;
                    [nucQuality_temp,n_nuc_temp] = getNucQuality(Lf/Li*C0_fit);
                    new_C0 = Predict_Hatcurve_NonTC (C0_fit, Li, Lf, n_nuc_temp*Li/Lf, 0);
                    C0_cutoff = Predict_Hatcurve_NonTC (C0_fit, Li, L_cutoff, n_nuc_temp*Li/L_cutoff, 0);
                    [nucQuality,n_nuc] = getNucQuality(new_C0);
                    [Height, xp, wn, wp, ~, ~] = HCpara(C0_fit);
                    [Height_cutoff, ~, ~, ~, ~, ~] = HCpara(C0_cutoff);
                    height_record(1) = Height;
                    width_record(1) = wp-wn;
                    string_0 = ['Initial hat curve, N_{nuc} = ', num2str(n_nuc, '%.0f'),...
                        ' (', nucQuality, ');',...
                        newline,...
                        'Height = ', num2str(1000*Height, '%.0f'), ' nm; ',...
                        'Width = ', num2str(wp-wn, '%.1f'), ' turns.'];
                    t0 = title (string_0);
                    t0.FontSize = 12;
                    t0.FontWeight = 'normal';
                    flagg = 1;
                    
                else
                    C0_fit = [-16.0859    1.0231
                        -1.5938    1.1486
                        38.0039    0.9968
                        57.1456    0.6842];
                    zz = f_5piece(C0_fit,turns_p);
                    nexttile(1)
                    plot (turns_p,zz*1000,'c--','Linewidth',4);
                    ylim ([0 1600]);
                    hold on;
                    nexttile(2)
                    plot (turns_p,zz*1000,'c--','Linewidth',4);
                    ylim ([0 1600]);
                    hold on;
                    string_0 = ['Initial hat curve',...
                        newline,...
                        'Bad'];
                    t0 = title (string_0);
                    t0.FontSize = 12;
                    t0.FontWeight = 'normal';
                    flagg = 0;
                end
                nexttile(3, [2 2])
                if flagg == 1
                    zz_new = f_5piece(new_C0,turns_p);
                    d_turns_p = xp;
                    plot (turns_p-d_turns_p, 1000*zz_new, 'Color', [0.9290 0.6940 0.1250],'LineStyle','--', 'LineWidth', 2);
                    hold on;
                    ax_all = gca;
                    ax_all.ColorOrderIndex = 1;
                else
                    d_turns_p = 0;
                end
                centered_initial = [turns_p-d_turns_p, ext_p];
                plot (turns_p-d_turns_p, ext_p, 'LineWidth', 6);
                hold on;
                xlim ([-70 70]);
                ylim ([0 max_plot_height]);
                ylabel('Extension (nm)');
                xlabel('Turns');
                set(gca,'FontSize',15,'LineWidth',1.5);
                set(gca,'Color', 'None');
                set(gca,'TickDir','out');
                set(gca,'TickLength',[0.015, 0.01]);
                set(gca,'XColor','k', 'YColor','k');
                set(gca, 'Layer', 'top');
                ax_all = gca;
                box(ax_all,'off');
                t_all = title ('Overlaying hat curves (Centered)');
                t_all.FontSize = 12;
                t_all.FontWeight = 'normal';
                hold on;
                if flagg == 1
                    nexttile(2)
                    zz_new = f_5piece(new_C0,turns_p);
                    plot (turns_p, 1000*zz_new, 'k', 'LineWidth', 4);
                    hold on;
                    plot (turns_p, 1000*zz_new, 'g--', 'LineWidth', 2);
                    hold on;
                    xlim ([-70 70]);
                    ylim ([0 max_plot_height]);
                    ylabel('Extension (nm)');
                    xlabel('Turns');
                    set(gca,'FontSize',15,'LineWidth',1.5);
                    set(gca,'Color', 'None');
                    set(gca,'TickDir','out');
                    set(gca,'TickLength',[0.015, 0.01]);
                    set(gca,'XColor','k', 'YColor','k');
                    set(gca, 'Layer', 'top');
                    ax_all = gca;
                    box(ax_all,'off');
                    hold on;
                    [nucQuality,n_nuc] = getNucQuality(new_C0);
                    [Height_new, xp_new, wn_new, wp_new, ~, ~] = HCpara(new_C0);
                    hat_curve_record{2} = new_C0;
                    height_record(2) = Height_new;
                    width_record(2) = wp_new-wn_new;
                    string_new = ['Predicted hat curve after 350 bp, N_{nuc} = ', num2str(n_nuc, '%.0f'),...
                        ' (', nucQuality, ');',...
                        newline,...
                        'Height = ', num2str(1000*Height_new, '%.0f'), ' nm; ',...
                        'Width = ', num2str(wp_new-wn_new, '%.1f'), ' turns.'];
                    tnew = title (string_new);
                    tnew.FontSize = 12;
                    tnew.FontWeight = 'normal';
                    hold on;
                    
                    
                    range_after_flush = task_for_timespot==21;
                    [turns_temp, ext_temp] = binning_hat(Turns(range_after_flush),...
                        Extension_notfiltered(range_after_flush));
                    nexttile(5)
                    plot (turns_temp, ext_temp, 'Color', [0 0 0], 'LineWidth', 4);
                    xlim ([-70 70]);
                    ylim ([0 max_plot_height]);
                    ylabel('Extension (nm)');
                    xlabel('Turns');
                    set(gca,'FontSize',15,'LineWidth',1.5);
                    set(gca,'Color', 'None');
                    set(gca,'TickDir','out');
                    set(gca,'TickLength',[0.015, 0.01]);
                    set(gca,'XColor','k', 'YColor','k');
                    set(gca, 'Layer', 'top');
                    ax_all = gca;
                    box(ax_all,'off');
                    hold on;
                    [x_c, ~] = ginput();
                    x_c = sort(x_c);
                    range_chosen = false(size(turns_temp));
                    for xci = 1:(length(x_c)/2)
                        range_chosen_i = turns_temp>=x_c(2*xci-1) & turns_temp<=x_c(2*xci);
                        range_chosen = range_chosen | range_chosen_i;
                    end
                    
                    turns_after_flush = turns_temp(range_chosen);
                    ext_after_flush = ext_temp(range_chosen);
                    
                    Cf_fit = fit_for_shifted_hatcurve (turns_after_flush, ext_after_flush);
                    
                    
                    
                    [Heightf, ~, ~, ~, ~, ~] = HCpara(Cf_fit);
                    dH = 1000*(Height_new - Heightf);
                    [fitresult_f, gof_f] = fit_for_ita(turns_after_flush, ext_after_flush, dH, new_C0);
                    ita = fitresult_f.ita;
                    dn_ita = fitresult_f.dn;
                    Cf_fit_2 = Hat_curve_after_transcribing_through_ita_nucleosomes (ita, dn_ita, dH, new_C0);
                    hat_curve_record{3} = Cf_fit_2;
                    zf = f_5piece(Cf_fit_2,turns_after_flush);
                    [nucQualityj,n_nucj] = getNucQuality(Cf_fit_2);
                    [Heightf, xpf, wnf, wpf, ~, ~] = HCpara(Cf_fit_2);
                    height_record(3) = Heightf;
                    width_record(3) = wpf-wnf;
                    
                    plot (turns_after_flush,zf*1000,'g--','Linewidth',2);
                    hold on;
                    string_0 = ['Final hat curve',...
                        ', N_{nuc} = ', num2str(n_nucj, '%.0f'),...
                        ' (', nucQualityj, ');',...
                        newline,...
                        'Height = ', num2str(1000*Heightf, '%.0f'), ' nm; ',...
                        'Width = ', num2str(wpf-wnf, '%.1f'), ' turns.'];
                    t0 = title (string_0);
                    t0.FontSize = 12;
                    t0.FontWeight = 'normal';
                    nexttile(3, [2 2])
                    plot (turns_after_flush-xpf, ext_after_flush, 'LineWidth', 6);
                    hold on;
                    legend ('Predicted hat curve after 350 bp', 'Initial hat curve',...
                        'Final hat curve' ,'location','northeast');
                    hold on;
                    
                    centered_final = [turns_after_flush-xpf, ext_after_flush];
                    
                    nexttile(6)
                    ww = 0:10:100;
                    hh = 1000*ww*(height_record(2))/width_record(2);
                    xlim ([round(width_record(2))-15 round(width_record(2))+15]);
                    ylim ([round(height_record(2)*10)*100-150, round(height_record(2)*10)*100+150]);
                    ylabel('Hat curve height (nm)');
                    xlabel('Hat curve width');
                    set(gca,'FontSize',15,'LineWidth',1.5);
                    set(gca,'Color', 'None');
                    set(gca,'TickDir','out');
                    set(gca,'TickLength',[0.015, 0.01]);
                    set(gca,'XColor','k', 'YColor','k');
                    set(gca, 'Layer', 'top');
                    ax_all = gca;
                    box(ax_all,'off');
                    ax_all.ColorOrderIndex = 1;
                    hold on;
                    r_max = 5;
                    for rr = r_max:-1:1
                        hhp = hh*(1-rr*dh/height_record(2));
                        hhn = hh*(1+rr*dh/height_record(2));
                        wwp = ww*(1+rr*dw/width_record(2));
                        wwn = ww*(1-rr*dw/width_record(2));
                        Cr = rr/r_max*[1 1 1];
                        Pr = patch ([wwp, flip(wwn)], [hhp, flip(hhn)], Cr);
                        Pr.LineStyle = 'none';
                        hold on;
                    end
                    for kk = [1 3 2]
                        Sk = scatter(width_record(kk), 1000*height_record(kk));
                        Sk.MarkerFaceColor = 'flat';
                        Sk.SizeData = 75;
                        hold on;
                    end
                    plot (ww, 1000*Height_cutoff*ones(size(ww)), 'r--', 'LineWidth', 0.5);
                    hold on;
                    %lambda_final = lambda;
                    Ratio_vec = [width_record(3), height_record(3)]./[width_record(2), height_record(2)];
                    Ratio_final = sqrt(sum(Ratio_vec.^2))/sqrt(2);
                    delta_n = n_nuc*(Ratio_final-1);
                    t_conclusion = title(['\Delta H = ',num2str(dH, '%.1f'), ' nm;',...
                        newline,...
                        '\delta n = ', num2str(delta_n, '%.3f'),'; \eta = ', num2str(ita, '%.3f')]);
                    t_conclusion.FontSize = 12;
                    output_temp = [output_temp, n_nuc, ...
                        width_record(2), 1000*height_record(2),...
                        width_record(3), 1000*height_record(3),...
                        1000*Height_cutoff, ...
                        delta_n, reshape(C0_fit, [1 8]),  reshape(Cf_fit, [1 8]), ita, dn_ita];
                end
                
                saveas(gcf,[new_folder,['\Box ', num2str(plotting_index(boxes)-1),'_hat.fig']]);
                saveas(gcf,[new_folder,['\Box ', num2str(plotting_index(boxes)-1),'_hat.png']]);
                saveas(gcf,[new_folder,['\Box ', num2str(plotting_index(boxes)-1),'_hat.emf']]);
                close;
                output = [output; output_temp];
            end
        end
        
        
        
        
    end
    close all;
    
    %% Plotting time devide lines
    
else % There is a possibility that one doesn't choose any files
    disp('You did not choose any files!');
end


%% Functions
function return_matrix = find_TC_index (fileID_Force, key)
% This function is to deal with the force calibration file. We want to
% get which traces are good traces (in the previously catagorized good
% trace catagory), and what are the two forces (there supposed to be
% two forces);

titleline = fgetl(fileID_Force); % Get the title line, since the first line is just text

ii = 1;
box_index = 0;
while feof(fileID_Force)==0
    string = fgetl(fileID_Force);
    is_TC = strfind(string,key);
    if isempty(is_TC) == 1
    else
        TC_bead_index_matrix (ii) = box_index;
        starting = is_TC+length(key)+1;
        shorter_string = string(starting:end);
        shorter_data = str2num (shorter_string);
        Force (1,ii) = shorter_data(2);
        Force (2,ii) = shorter_data(1);
        ii = ii+1;
    end
    box_index = box_index+1;
end
return_position = TC_bead_index_matrix (1:ii-1);
return_matrix = [return_position; Force(:,1:ii-1)];
end

function index = timespot_finding (task)
% This is a function to find the time for changing to different task
% lists.
index = length(task);
task_temp = max(task)-task;
while (1)
    mask = sign(task_temp);
    if mask == 0
        break;
    end
    index = [sum(mask);index];
    task_temp = task_temp-1;
    task_temp = task_temp.*mask;
end
%timespot = time(index);
end

function [fitresult, gof] = Linear_Fit(xx, yy)

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( xx, yy );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );
end

function [ essential_columns_position ] = find_essential_columns_by_name( title_line, essential_columns )
% This function returns the column numbers of the columns we are
% caring, searching by name.
tabs_pos = find(title_line == 9); % Find were are the "tab" signs are
tabs_pos = [tabs_pos, (length(title_line)+1)];
[sizex, sizey] = size(essential_columns);
essential_columns_position = ones(1,sizex);
for ii = 1:sizex
    key = essential_columns{ii};
    find_key = strfind(title_line,key);
    essential_columns_position(ii) = find(tabs_pos>=find_key(1),1);
end
%essential_columns_position = sort(essential_columns_position);
end

function [ Box_cell_position ] = find_boxes( title_line )
% This function returns the column numbers of the columns we are
% caring, searching by name.
tabs_pos = find(title_line == 9); % Find were are the "tab" signs are
tabs_pos = [tabs_pos, (length(title_line)+1)];
key = 'Box';
find_key = strfind(title_line,key);
pp_boxes = length(find_key);
for pp = 1:pp_boxes
    Box_cell_position(pp) = find(tabs_pos>=find_key(pp),1);
end

%essential_columns_position = sort(essential_columns_position);
end

function [turns_p_output, ext_p_output, turns_n_output, ext_n_output] = restricted_binning (turns_p, ext_p, turns_n, ext_n)
[~,ip,in] = intersect(turns_p, turns_n);
turns_p_output = turns_p(ip);
turns_n_output = turns_n(in);
ext_p_output = ext_p(ip);
ext_n_output = ext_n(in);
end

function surface_shift = lookup_surface_shift (surface_shift_info, Date, Time, Box)

Data_temp_1 = surface_shift_info;
Date_index = Data_temp_1(:,1) == Date;
Data_temp_2 = Data_temp_1(Date_index, :);
Time_index = Data_temp_2(:,2) == Time;
Data_temp_3 = Data_temp_2(Time_index, :);
Box_index = Data_temp_3(:,3) == Box;
surface_shift = Data_temp_3(Box_index, 4);

end

function output = Smoothed_trace(Time, Extension, Task, smoothing_time)
dT = mean(diff(Time));
smoothing_factor = smoothing_time/dT;
reasonable_region = (Extension>=-1000 & Extension <=3000);
Time_reasonable = Time(reasonable_region);
Extension_reasonable = Extension(reasonable_region);
Extension_interp = interp1(Time_reasonable,Extension_reasonable,Time,'linear', 'extrap');
output = Smooth_data_each_task (Extension_interp, Task, smoothing_factor);
end

function output = Smooth_data_each_task (Extension, Task, smoothing_factor)
output = zeros(size(Extension));
for ii = min(Task):max(Task)
    current_region = Task==ii;
    if sum(current_region)
        output(current_region) = smoothdata(Extension(current_region), 'gaussian', smoothing_factor);
    end
end
end

function [output_turns, output_ext] = binning_hat(Turns, Extension)
min_turns = fix(min(Turns));
%min_turns = max(-40, min_turns);
max_turns = fix(max(Turns));
data = [Turns, Extension];
data_sort = sortrows(data,1);
Turns_sort = data_sort(:,1);
Extension_sort = data_sort(:,2);
output_turns = [];
output_ext = [];
for ii = min_turns:max_turns
    sub_index = find(abs(Turns_sort-ii)<=0.1);
    output_turns = [output_turns; ii];
    mean_ext = mean(Extension_sort(sub_index));
    output_ext = [output_ext; mean_ext];
end
%output = [output_turns, output_ext];
end

function R0 = fit5piece_Jin(xx, zz)
[xData, yData] = prepareCurveData( xx, zz );
ft = fittype( 'max (0, heaviside (( x2 + dx_12 ) - x)* ( ( 2*(dy_12 - k_32*dx_12 )/dx_12 + k_32 ) *(x - ( x2 + dx_12 )) + ( y2 + dy_12 )) + heaviside (x - ( x2 + dx_12 ))* heaviside (x2 - x)*  ( ( dy_12 - k_32*dx_12 )/dx_12^2 *(x - x2)^2 + k_32 *(x - x2) + y2) + heaviside (x - x2)* heaviside (( x2 + dx_32 ) - x)* ( k_32*(x-x2) + y2) + heaviside (x - ( x2 + dx_32 ))* heaviside (( x2 + dx_32 + dx_43 ) - x)*  ( ( dy_43 - k_32*dx_43 )/dx_43^2 *(x - ( x2 + dx_32 ))^2 + k_32 *(x - ( x2 + dx_32 )) + ( y2 + ( k_32*dx_32 ) )) + heaviside (x - ( x2 + dx_32 + dx_43 ))* ( ( 2*(dy_43 - k_32*dx_43 )/dx_43 + k_32 ) *(x - ( x2 + dx_32 + dx_43 )) + ( y2 + ( k_32*dx_32 ) + dy_43 )))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-20 10 2 -0.3 -0.8 -0.006 -15 0.5];
opts.StartPoint = [-15 38 20 -0.13 -0.31 -0.004 -1 1.2];
opts.Upper = [-2 50 40 -0.01 -0.01 -0.002 5 2.5];

[fitresult, gof] = fit( xData, yData, ft, opts );
x2 = fitresult.x2;
y2 = fitresult.y2;
x1 = x2 + fitresult.dx_12;
x3 = x2 + fitresult.dx_32;
x4 = x3 + fitresult.dx_43;
y1 = y2 + fitresult.dy_12;
dy_32 = fitresult.dx_32*fitresult.k_32;
y3 = y2 + dy_32;
y4 = y3 + fitresult.dy_43;
R0 = [x1, y1; x2, y2; x3, y3; x4, y4];

end

function z = f_5piece(r0,x)
% r is the Catersian coordinates of 4 points that define the 5 piece
% function
r = reshape(r0,4,2);
% top piece
y3 = r(2,2) + (r(3,2)-r(2,2))/(r(3,1)-r(2,1)).*(x-r(2,1));
%% the two curving parts
%% top left curve
M_tl = [2*r(2,1) 1 0;
    r(1,1)^2  r(1,1) 1;
    r(2,1)^2  r(2,1) 1];
C_tl = [(r(3,2)-r(2,2))/(r(3,1)-r(2,1)) ;
    r(1,2);
    r(2,2)];
para_tl= M_tl\C_tl;
y2 = para_tl(1) * x.^2 + para_tl(2) * x + para_tl(3);

%% top right curve
M_tr = [2*r(3,1) 1 0;
    r(4,1)^2  r(4,1) 1;
    r(3,1)^2  r(3,1) 1];
C_tr = [(r(3,2)-r(2,2))/(r(3,1)-r(2,1)) ;
    r(4,2);
    r(3,2)];
para_tr = M_tr\C_tr;
y4 = para_tr(1) * x.^2 + para_tr(2) * x + para_tr(3);

%% calculating alpha_n and alpha_p
% lower left piece
alpha_n = 2 * para_tl(1) * r(1,1) + para_tl(2);
y1 = r(1,2) + alpha_n * (x-r(1,1));

% lower right piece
alpha_p = 2 * para_tr(1) * r(4,1) + para_tr(2);
y5 = r(4,2) + alpha_p * (x-r(4,1));

%% whole curve
z = y1 .* (x<r(1,1)) + y2 .* (x>=r(1,1) & x<r(2,1)) + y3 .* (x>=r(2,1) & x<r(3,1)) + y4 .* (x>=r(3,1) & x<r(4,1)) + y5.* (x>=r(4,1));

end

function [Height, hatcenter, buckling_negative, buckling_positive, slope_negative, slope_positive] = HCpara(r0)
%obtains hat curve parameters from the result of a 5 piece fit
% buckling positive is buckling point relative to the hat center


r = reshape(r0,4,2);
M_tl = [2*r(2,1) 1 0;  r(1,1)^2  r(1,1) 1; r(2,1)^2  r(2,1) 1];
C_tl = [(r(3,2)-r(2,2))/(r(3,1)-r(2,1)) ;  r(1,2);  r(2,2)];
para_tl= M_tl\C_tl;

M_tr = [2*r(3,1) 1 0;  r(4,1)^2  r(4,1) 1;  r(3,1)^2  r(3,1) 1];
C_tr = [(r(3,2)-r(2,2))/(r(3,1)-r(2,1)) ;  r(4,2);  r(3,2)];
para_tr = M_tr\C_tr;
alpha_n = 2 * para_tl(1) * r(1,1) + para_tl(2);
alpha_p = 2 * para_tr(1) * r(4,1) + para_tr(2);
slope_negative = alpha_n;
slope_positive = alpha_p;
alpha_c = (r(3,2)-r(2,2))/(r(3,1)-r(2,1));

buckling_negative = (r(2,2) - r(1,2) + alpha_n * r(1,1) - alpha_c * r(2,1))/(alpha_n-alpha_c);
buckling_positive = (r(4,2) - r(3,2) - alpha_p * r(4,1) + alpha_c * r(3,1))/(alpha_c-alpha_p);

%% top left curve
M_tl = [2*r(2,1) 1 0;
    r(1,1)^2  r(1,1) 1;
    r(2,1)^2  r(2,1) 1];
C_tl = [(r(3,2)-r(2,2))/(r(3,1)-r(2,1)) ;
    r(1,2);
    r(2,2)];
para_tl= M_tl\C_tl;
Height = para_tl(3) - para_tl(2)^2/(4*para_tl(1));
hatcenter =  -para_tl(2)/2/para_tl(1);
buckling_positive = buckling_positive - hatcenter;
buckling_negative = buckling_negative - hatcenter;
end

function [nucQuality,n_nuc] = getNucQuality(r0)
%obtains goodness of a nucleosome array based on Tung's method described in
%190119_Tung_Selection criteria for single_double chromatin array data on
%AOT_MT.pptx assumes nuc array is from 197-64ex plasmid
%modified based on 190205_Seong ha_Hat curve dimensions in different buffer conditions.pptx
%load('x_f_WLC.mat');

r = reshape(r0,4,2);
%load('nVLatF.mat');
%slope = interp1(nVLatF(:,1),nVLatF(:,2),F); % N Vs Length relationship obtained on 190520
slope = -0.0414;
%intercept = interp1(nVLatF(:,1),nVLatF(:,3),F);
intercept = 3.2306;
[Height, ~, ~, buckling_positive, ~,~] = HCpara(r);
%dL = 0.200; % length measurement uncertainty on MT
n_nuc =  (Height - intercept)/slope;
n_nucRange = [n_nuc-5 n_nuc+5];
widthRange = 22.7936+ 0.5973*n_nucRange;

if buckling_positive > widthRange(1) && buckling_positive < widthRange(2) %abs(n_nuc  - n_nuc_width) < dn_nuc_width
    nucQuality = 'good';
else
    nucQuality = 'bad';
end

end

function new_C0 = Predict_Hatcurve_NonTC (C0_i, Li, Lf, p_i, dp)
L_standard = 12667;
C0_scaled = C0_i/Li*L_standard;
dp_scaled = dp/Lf*L_standard+p_i*L_standard*(1/Lf-1/Li);
new_C0_scaled = delta_hat_curve (C0_scaled, dp_scaled);
new_C0 = new_C0_scaled*Lf/L_standard;

end

function new_C0 = delta_hat_curve (C0, dp)

[~, ~, wL, wR, kL, kR] = HCpara(C0);

kM = (C0(2,2)-C0(3,2))/(C0(2,1)-C0(3,1));
Zero_turn_height = f_5piece(C0,0);

k_height = -42.6492/1000;
k_y2 = -38.4353/1000;
k_kL = -(61.76892-17.56061)/65/1000;
k_kM = -(-4.39694+4.17612)/65/1000;
k_kR = -(-57.32513+17.63203)/65/1000;
k_wL = 0.03915;
k_wR = 0.54257;
k_x2 = (1.23/4.97*25)/(13.22/6.4*70);
k_x3 = (6.39/5.28*60)/(11.38/6.96*70);

new_kL = kL+k_kL*dp;
new_kR = kR+k_kR*dp;
new_kM = kM+k_kM*dp;
new_wL = wL+k_wL*dp;
new_wR = wR+k_wR*dp;
new_x2 = C0(2,1)+k_x2*dp;
new_x3 = C0(3,1)+k_x3*dp;
new_height = Zero_turn_height+k_height*dp;

xmu = -wL;
QL = new_kL*xmu;
QM = new_kM*xmu;
ALPHAL = 0.25*(new_kL^2-new_kM^2)/(QM-QL);
%CL = QL+new_kL^2/(4*ALPHAL);
%x1_tilde = new_kL/(2*ALPHAL)+xmu;
x2_tilde = new_kM/(2*ALPHAL)+xmu;
wJL = new_x2-x2_tilde;
wJR = new_wR-new_wL+wJL;

new_C0 = construct_yshifted_hatcurve_xpara (new_x2, new_x3, new_kL, new_kM, new_kR, wJL, wJR);
new_C0(:,2) = new_C0(:,2)+new_height;

end

function C_output = construct_yshifted_hatcurve_xpara (new_x2, new_x3, new_kL, new_kM, new_kR, wJL, wJR)

YL = new_kM*wJL;
bL = YL-new_kL*wJL;
qL = bL+new_kL*new_x2;
qM1 = new_kM*new_x2;
new_x1 = new_x2 + 2*(qM1-qL)/(new_kL-new_kM);
alpha_L = 0.5*(new_kL-new_kM)/(new_x1-new_x2);

YR = new_kM*wJR;
bR = YR-new_kR*wJR;
qR = bR+new_kR*new_x3;
qM2 = new_kM*new_x3;
new_x4 = new_x3 + 2*(qM2-qR)/(new_kR-new_kM);
alpha_R = 0.5*(new_kR-new_kM)/(new_x4-new_x3);

turns = [new_x1, new_x2, new_x3, new_x4, 0];

SI = heaviside(new_x1-turns).*(new_kL*turns+bL);
SII = heaviside(turns-new_x1).*heaviside(new_x2-turns).*( alpha_L*(turns-new_x2).^2+new_kM*(turns-new_x2)+qM1 );
SIII = heaviside(turns-new_x2).*heaviside(new_x3-turns).*( new_kM*turns );
SIV = heaviside(turns-new_x3).*heaviside(new_x4-turns).*( alpha_R*(turns-new_x3).^2+new_kM*(turns-new_x3)+qM2 );
SV = heaviside(turns-new_x4).*(new_kR*turns+bR);
S_total = SI+SII+SIII+SIV+SV;

mock_z0 = S_total(5);
mock_y1 = S_total(1);
mock_y2 = S_total(2);
mock_y3 = S_total(3);
mock_y4 = S_total(4);
new_y1 = mock_y1-mock_z0;
new_y2 = mock_y2-mock_z0;
new_y3 = mock_y3-mock_z0;
new_y4 = mock_y4-mock_z0;

C_output = [new_x1, new_y1; new_x2, new_y2; new_x3, new_y3; new_x4, new_y4];

end

function C0_fit = fit_for_shifted_hatcurve (turns, ext)

maximum_height_turns = mean(turns(find(ext == max(ext))));
dn0 = maximum_height_turns+2;
C0_fit = fit5piece_Jin(turns-dn0, ext/1000);
C0_fit(:,1) = C0_fit(:,1)+dn0;

end


function [fitresult, gof] = fit_for_ita(turns_after_flush, ext_after_flush, dH, C0)
lnuc = 530.70309/64;
[xData, yData] = prepareCurveData( turns_after_flush, ext_after_flush );

C0 = reshape(C0,4,2);

string_0 = '[';
for ii = 1:4
    string_0 = [string_0,num2str(C0(ii,1)),','];
    string_0 = [string_0,num2str(C0(ii,2)),';'];
end
string_0 = string_0(1:end-1);
string_0 = [string_0, ']'];

string_1 = ['1000*Nucleosome_array_hat_curve('...
    'Hat_curve_after_transcribing_through_ita_nucleosomes (ita, dn, ', ...
    num2str(dH), ', ',...
    string_0,...
    ')',...
    ', x)'];


% Set up fittype and options.
ft = fittype( string_1, 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-100 0];
opts.StartPoint = [-20 1];
opts.Upper = [10 max(dH/lnuc,0)];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

end
