clear;
close all;
clc;
%c = @cmu.colors;
% IMPORTANT! This program only works if there are two forces!

% This program is supposed to choose multiple files, but they should have
% the same protocols;

target_box = 152;

pickedcolor = [255, 0, 0]/255;
EndingT = 11;
tweezer2 = 0;
tweezer2 = 2;
resuming_index = 1;
%ext_to_bp_factor = 4.24629351960929; % 0.3 pN
ext_to_bp_factor = 3.8372; % 0.5 pN
starting_time = 85;
shift_buffer = 0;

flag = 0; % To indicate if we choose any files at all, even if we choose one file, this would be 1;
file_index = 0; % Number of files chosen;
key = '1 Wind';
Maximum_allowed_files = 1;
Maximum_allowed_boxes = 80;
Data_points_number = zeros(Maximum_allowed_files,1);
%Data_all = zeros (32768,1024,Maximum_allowed_files);
position_and_force = zeros(4,Maximum_allowed_boxes,Maximum_allowed_files);

essential_columns = {{'DAQ time (s)'}; {'Time (ms)'}; {'tasklist (index)'};{'Magnet Angle (Turns)'};{'Magnet Height (mm)'}};% This is the number of columns essential for our hatcurve;

nucleosome_position = 403; %in bp;
RNAP_downstream_protects = 18; % # of bps protected downstream +1 site
Li = 13014+nucleosome_position-403; % Initial template length, from RNA9 to the downstream BstXI cutting site

for ii = 1:Maximum_allowed_files
    position_and_force(2,:,ii) = (1:Maximum_allowed_boxes)-1;
end


while (1) % while not stopped, this loop keep on choosing files
    %% Access data files
    
    Address = 'Y:\Magnetic Tweezers\Jin';
    if tweezer2 == 2
        Address = 'Y:\Magnetic Tweezers 2\Jin\230914_Single Nucleosome 500 nM TFIIS F3 1\_230914_185030_PROCESSED\Data Files';
    end
    
    [filename, pathname] = uigetfile('*.*','Select the data file',Address); % Choose MT data files
    if filename == 0 % If you do not choose a file, then the program proceeds;
        break;
    else flag = 1;
    end
    new_file_name = strrep(filename,'_',' ');
    fileID = fopen ([pathname,filename],'r'); % Open the data file
    file_index = file_index + 1;
    new_folder = [pwd,'/',[filename,' plots']];
    mkdir([filename,' plots']);
    
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
    magnet_turns_for_subplot = Data_all_cell{4}; % The magnetic turns for the first file selected, used for plotting subplot;
    Timespot_index = timespot_finding (task_for_timespot);% Time spots given by the first file selected;
    Timespot = time_for_timespot(Timespot_index);
    EndingT = length(Timespot);
    
    %%
    for which_file = 1:file_index
        Extension_index = (21+tweezer2):8:4000;
        plotting_index = find (position_and_force(1,:,which_file)==1);
        Time = Data_all_cell{2}/1000;
        mov = 0*plotting_index;
        mov = mov';
        mov_rate = mov;
        start_moving_time = mov;
        first_slip = mov;
        finish_flowing = mov+Timespot(8);
        finish_waiting = mov+Timespot(11);
        start_flowing = mov+Timespot(7);
        breaking_time = mov;
        start_mov_rec = mov;
        finish_mov_rec = mov;
        EndingT = length(Timespot);
        
        Time_range_0 = (Timespot_index(1):Timespot_index(6))';
        Time_range_1 = (Timespot_index(10):Timespot_index(11))';
        Time_range_2 = (Timespot_index(13):Timespot_index(14))';
        Fast_winding_range = (Timespot_index(7):Timespot_index(8))';
        
        pos_hat_1_range = (Timespot_index(1):Timespot_index(2))';
        pos_hat_2_range = (Timespot_index(5):(Timespot_index(6)+10))';
        neg_hat_range = (Timespot_index(3):Timespot_index(4))';
        pos_hat_range = [pos_hat_1_range; pos_hat_2_range];
        
        start_ext_range = (1:Timespot_index(1))';
        
        
        target_index = find(plotting_index== target_box+1);
        
        for boxes = resuming_index:length(plotting_index)
            %for boxes = [target_index]
            extension_increase_during_flow = 0;
            Nuc_info = [1];
            fig_flag = 2;
            clc;
            disp(num2str(100*boxes/length(plotting_index), '%.1f'));
            total_shift = 0;
            while (fig_flag > 0)
                FIG_x = figure (boxes);
                pos = [100 100 1500 700];
                set(FIG_x, 'Pos', pos);
                %% Plotting subplots: Magnetic height and magnetic turns
                
                subplot(6, 5, [1 3]);
                %subplot_M_turns = subplot (5,1,1);
                plot (time_for_timespot, magnet_turns_for_subplot,'k','LineWidth',1.5);
                axis([0 inf -inf inf])
                ylabel('M-turns');
                xlim ([0 Timespot(EndingT)]);
                string_3 = ['Trace been shifted down by ', num2str(total_shift, '%.0f'),' nm.'];
                title ([new_file_name(1:(length(filename)-4)),' box ',num2str(plotting_index(boxes)-1)...
                    char(10), string_3]);
                h = gca;
                h.XAxis.Visible = 'off';
                set(gca,'FontSize',13,'LineWidth',1.5);
                set(gca, 'Color', 'None');
                ax_all = gca;
                %box(ax_all,'off');
                %subplot_M_height = subplot (5,1,2);
                subplot(6, 5, [6 8]);
                plot (time_for_timespot, magnet_height_for_subplot,'k','LineWidth',1.5);
                axis([0 inf -inf inf])
                ylabel('M-height [mm]');
                xlim ([0 Timespot(EndingT)]);
                h = gca;
                h.XAxis.Visible = 'off';
                set(gca,'FontSize',13,'LineWidth',1.5);
                set(gca, 'Color', 'None');
                ax_all = gca;
                box(ax_all,'off');
                
                %% Plotting all the good traces
                %subplot_main = subplot (5,1,[3,5]);
                subplot(6, 5, [11 28]);
                
                if fig_flag == 2
                    Extension = 1000*Data_all_cell{boxes+length(essential_columns_position)};
                    Extension_notfiltered = Extension;
                    Extension = movmean(Extension,20);
                    
                else
                    Extension(Timespot_index(7):end) = Extension(Timespot_index(7):end)-extension_increase_during_flow;
                    Extension_notfiltered(Timespot_index(7):end) = Extension_notfiltered(Timespot_index(7):end)-extension_increase_during_flow;
                end
                
                Tracename = ['File',num2str(which_file),' Box',num2str(plotting_index(boxes)-1)];
                plot (Time, Extension_notfiltered, 'Color', [0.7,0.7,0.7]);
                hold on;
                plot (Time, Extension, 'Color', [0,0,1], 'DisplayName',Tracename, 'Linewidth', 2);
                xlim ([0 Timespot(EndingT)]);
                initial_range = find(Time>=0 & Time<=Timespot(1));
                initial_ext = mean(Extension(initial_range));
                initial_ext_round = round(initial_ext/100)*100;
                ylim ([0 4000]);
                h = gca;
                hold on;
                
                
                for ii = 1:length(Timespot)-1
                    time_i = Timespot (ii);
                    yy = 0:10:5000;
                    xx = time_i*ones(size(yy));
                    plot(xx,yy,'k--','LineWidth',1.5);
                    hold on;
                end
                
                h = gca;
                set(gca,'FontSize',13,'LineWidth',1.5);
                set(gca, 'Color', 'None');
                ax_all = gca;
                ylabel('Extension [nm]');
                xlabel('Time [s]');
                hold on;
                
                
                tile_y = subplot(6, 5, [4 15]);
                %string_0 = ['Box ',num2str(plotting_index(boxes)-1)];
                string_0 = [];
                [output_turns, pos_hat_ext] = binning_hat(magnet_turns_for_subplot(pos_hat_2_range),...
                    Extension_notfiltered(pos_hat_2_range));
                
                end_ext_range = (Timespot_index(6):Timespot_index(6)+10)';
                end_ext = mean(Extension_notfiltered(end_ext_range));
                start_ext = mean(Extension_notfiltered(start_ext_range));
                
                dEXT_ini_fin = end_ext - start_ext;
                
                
                [turns_p, ext_p] = binning_hat(magnet_turns_for_subplot(pos_hat_range),...
                    Extension_notfiltered(pos_hat_range));
                turns_p_index = find(turns_p>=-44 & turns_p<=69);
                turns_p = turns_p(turns_p_index);
                ext_p = ext_p(turns_p_index);
                [turns_n, ext_n] = binning_hat(magnet_turns_for_subplot(neg_hat_range),...
                    Extension_notfiltered(neg_hat_range));
                turns_n_index = find(turns_n>=-44 & turns_n<=69);
                turns_n = turns_n(turns_n_index);
                ext_n = ext_n(turns_n_index);
                
                d_ext_during_harcurve = ext_p-ext_n;
                dEXT_during = mean(d_ext_during_harcurve);
                
                string_1 = ['Ext(T) - Ext(0) = ', num2str(dEXT_ini_fin, '%.0f'), ' nm; Mean hysteresis = ',...
                    num2str(abs(dEXT_during), '%.0f'), ' nm.'];
                
                flagg = 0;
                
                if abs(dEXT_ini_fin)<=500 && dEXT_during<=50
                    ext_np = 0.5*(ext_n+ext_p);
                    C0_fit = fitSingle_m_Jin(turns_p, ext_np/1000);
                    zz = f_Single(C0_fit,turns_p);
                    plot (turns_p,zz*1000,'g--','Linewidth',1.5);
                    hold on;
                    flagg = 1;
                    pre_title = 'Good';
                    legend_script = 'Hat curve fit';
                    
                else
                    field1 = 'a';
                    field2 = 'b';
                    field3 = 'km';
                    field4 = 'kp';
                    field5 = 'n0';
                    field6 = 'nsm';
                    field7 = 'nsp';
                    value1 = -5.100481346514311e-04;
                    value2 = 3.4402;
                    value3 = 0.0479;
                    value4 = 0.0479;
                    value5 = 0;
                    value6 = 16.287056412899766;
                    value7 = 16.287056412899766;
                    
                    C0_fit = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7);
                    
                    zz = f_Single(C0_fit,turns_p);
                    
                    string_1 = [string_1, ' (Bad)'];
                    plot (turns_p,zz*1000,'c--','Linewidth',1.5);
                    hold on;
                    pre_title = 'Bad';
                    legend_script = 'Standard hat curve';
                    
                end
                
                xx = -200:0.1:200;
                dL1 = nucleosome_position+147; % # of bps from +1 site to the exit of 1st NPE
                dL1 = dL1-RNAP_downstream_protects;
                dn1 = dL1/10.5; % # of turns polymerase placed
                L1 = Li-dL1; % Remaining downstream template length
                dp1_enter = nucleosome_depletion_calculation(nucleosome_position, nucleosome_position, Nuc_info);
                dp1_exit = nucleosome_depletion_calculation(nucleosome_position+147, nucleosome_position, Nuc_info);
                new_C01 =  Predict_Hatcurve_Single (C0_fit, Li, L1, dp1_exit);
                new_C01_enter =  Predict_Hatcurve_Single (C0_fit, Li, L1+147, dp1_enter);
                
                
                yy1 = f_Single(new_C01,xx);
                yy1_enter = f_Single(new_C01_enter,xx);
                upper_curve1 = max([yy1; yy1_enter]);
                lower_curve1 = min([yy1; yy1_enter]);
                upper_curve1 = max(upper_curve1, 0);
                lower_curve1 = max(lower_curve1, 0);
                SP1 = patch([xx, flip(xx)], [upper_curve1, flip(lower_curve1)]*1000, [1 0 0]);
                SP1.FaceAlpha = 0.25;
                SP1.LineStyle = 'none';
                hold on;
                
                
                set(gca,'ColorOrderIndex',2)
                
                PP = plotting_initial_hatcurve(Extension_notfiltered, magnet_turns_for_subplot, pos_hat_range, FIG_x, tile_y);
                PP = plotting_initial_hatcurve(Extension_notfiltered, magnet_turns_for_subplot, neg_hat_range, FIG_x, tile_y);
                PP = plotting_initial_hatcurve (Extension_notfiltered, magnet_turns_for_subplot, Time_range_1 , FIG_x, tile_y);
                PP = overlaying_hatcurve (Extension_notfiltered, magnet_turns_for_subplot,Fast_winding_range, FIG_x, tile_y);
                PP.LineStyle = '--';
                PP.Color = 'b';
                set(gca,'ColorOrderIndex',4)
                
                ax = gca;
                ax.TitleFontSizeMultiplier = 0.8;
                title(string_1);
                
                ylim ([0 4000]);
                xlim ([-150 80]);
                ylabel('Extension [nm]');
                xlabel('Magnet Turns');
                %
                if fig_flag == 1
                    plot (To_be_fitted(:,1), QQ, 'c--', 'LineWidth', 1.5);
                    hold on;
                    lgd = legend (legend_script, 'Inside NPE',...
                        'Initial hat curve (- to +)', 'Initial hat curve (+ to -)',...
                        'Hat curve at 20 minutes', 'Fast Hat curve after flowing', 'Fit of final hat curve',...
                        'location','northwest');
                else
                    lgd = legend (legend_script, 'Inside NPE',...
                        'Initial hat curve (- to +)', 'Initial hat curve (+ to -)',...
                        'Hat curve at 20 minutes', 'Fast Hat curve after flowing', 'location','northwest');
                end
                
                set(gca,'FontSize',13,'LineWidth',1.5);
                set(gca, 'Color', 'None');
                ax_all = gca;
                box(ax_all,'off');
                lgd.FontSize = 8;
                set(lgd,'color','none');
                ax_all = gca;
                box(ax_all,'off');
                hold on;
                
                Waiting_range1 = find( Time>=Timespot(9) & Time<=Timespot(10) );
                Waiting_range = Waiting_range1;
                
                if fig_flag ==2
                    qstring = 'Good Trace?';
                    choice = questdlg(qstring,'Question?','Yes','Nope','Yes');
                    extension_increase_during_flow = 0;
                end
                
                
                if (strcmp(choice,'Yes'))
                    [T_pos, X_pos, T_neg, X_neg, Critical, Hat_record] = ...
                        Single_position_new_version(Time, Extension, ...
                        task_for_timespot, magnet_turns_for_subplot, nucleosome_position);
                    T_pos = T_pos(~isnan(X_pos));
                    X_pos = X_pos(~isnan(X_pos));
                    %X_pos = movmean(X_pos, 40);
                    T_neg = T_neg(~isnan(X_neg));
                    X_neg = X_neg(~isnan(X_neg));
                    %X_neg = movmean(X_neg, 40);
                    T_c = Critical(:,1);
                    
                    subplot(6, 5, [19 30]);
                    plot (T_pos, movmean(X_pos, 40), 'k', 'LineWidth', 0.5);
                    hold on;
                    plot (T_neg, movmean(X_neg, 40), 'b', 'LineWidth', 3);
                    hold on;
                    if fig_flag ==2
                        end_position_theory = X_pos(end);
                    end
                    S_end = scatter (T_pos(end),end_position_theory, 50 );
                    S_end.MarkerEdgeColor = 'm';
                    
                    hold on;
                    xlabel ('Time [s]');
                    ylabel ('Pol II Position [bp]');
                    ylim ([0 1000]);
                    xlim ([0 1300]);
                    hold on;
                    
                    ax_all = gca;
                    box(ax_all,'off');
                    unit_vec = ones(size(T_c));
                    for ii = 1:1
                        enter_ii = nucleosome_position+(ii-1)*197;
                        exit_ii = enter_ii+147;
                        plot (T_c-T_c(1), enter_ii*unit_vec, 'r--', 'LineWidth',2);
                        hold on;
                        plot (T_c-T_c(1), exit_ii*unit_vec, 'g--', 'LineWidth',2);
                        hold on;
                    end
                    hold off;
                    set(gca,'FontSize',13,'LineWidth',1.5);
                    set(gca, 'Color', 'None');
                    ax_all = gca;
                    box(ax_all,'off');
                    
                    
                    subplot(6, 5, [11 28]);
                    [~, size_cy] = size(Critical);
                    for critn = 2:2:size_cy-1
                        yy_enter = Critical(:,critn);
                        yy_exit = Critical(:,critn+1);
                        plot (T_c, yy_enter, 'r:', 'Linewidth', 1.5);
                        hold on;
                        plot (T_c, yy_exit, 'g:', 'Linewidth', 1.5);
                        hold on;
                    end
                    
                    
                    if fig_flag == 1
                        qstring_2 = 'Satisfied?';
                        choice_2 = 'Yes';
                        %choice_2 = questdlg(qstring_2,'Question?','Yes','Nope','Yes');
                        
                        if (strcmp(choice_2,'Yes'))
                            
                            fig_flag = 0;
                            % open your file for writing
                            file_name_write = [new_folder,['\',new_file_name(1:14),...
                                'box ', num2str(plotting_index(boxes)-1),' Position.txt']];
                            file_name_trace = [new_folder,['\',new_file_name(1:14),...
                                'box ', num2str(plotting_index(boxes)-1),' Re-aligned.txt']];
                            fid_write = fopen(file_name_write,'wt');
                            fid_trace = fopen(file_name_trace,'wt');
                            % write the matrix
                            
                            length_max = max(length(T_neg), length(T_pos));
                            myData = ones(length_max, 5)*-100;
                            myData(1:length(T_neg),1) = T_neg;
                            myData(1:length(T_neg),2) = X_neg;
                            myData(1:length(T_pos),3) = T_pos;
                            myData(1:length(T_pos),4) = X_pos;
                            myData(2:4, 5) = Nuc_info.';
                            myData(5:6, 5) = x_end.';
                            myData(1,5) = total_shift;
                            
                            myData =myData.';
                            myData_trace =[Time Extension_notfiltered/1000 ...
                                task_for_timespot magnet_turns_for_subplot].';
                            if fid_write > 0
                                fprintf(fid_write,'%f\t%f\t%f\t%f\t%f\n',myData);
                                fclose(fid_write);
                            end
                            if fid_trace > 0
                                fprintf(fid_trace,'%f\t%f\t%d\t%f\n',myData_trace);
                                fclose(fid_trace);
                            end
                            
                        end
                    else
                        [x_end_0, ~] = ginput(2);
                        fig_flag = 1;
                    end
                    
                    if fig_flag >0
                        [end_turns, end_ext] = binning_hat...
                            ( magnet_turns_for_subplot(Time_range_1), Extension_notfiltered(Time_range_1));
                        end_turns_0 = magnet_turns_for_subplot(Timespot_index(10));
                        %x_end = [end_turns_0, x_end_0];
                        x_end = x_end_0;
                        end_select = end_turns<=max(x_end) & end_turns>=min(x_end);
                        end_turns_select = end_turns(end_select);
                        end_ext_select = end_ext(end_select)/1000;
                        To_be_fitted = [end_turns_select, end_ext_select];
                        
                        [~, ~, ~, ~, ~, Hat_record] = ...
                            Single_position_new_version(Time, Extension, task_for_timespot,...
                            magnet_turns_for_subplot, nucleosome_position);
                        [~, ~, ~, end_position_C, ~]...
                            = Fitting_N_Correction_Single(C0_fit,To_be_fitted,Hat_record);
                        [kink_y, kink_y_theory, end_position_theory, end_position_C, QQ]...
                            = Fitting_N_Correction_Single(end_position_C,To_be_fitted,Hat_record);
                        QQ = QQ*1000;
                        
                        extension_increase_during_flow = 1000*(kink_y-kink_y_theory);
                        QQ = QQ-extension_increase_during_flow;
                        total_shift = total_shift+extension_increase_during_flow;
                        
                        clf(FIG_x);
                        
                    end
                else
                    fig_flag = 0;
                end
                
            end
            %}
            saveas(gcf,[new_folder,['\',pre_title,'_Box ', num2str(plotting_index(boxes)-1),'_H.fig']]);
            saveas(gcf,[new_folder,['\',pre_title,'_Box ', num2str(plotting_index(boxes)-1),'_H.png']]);
            saveas(gcf,[new_folder,['\',pre_title,'_Box ', num2str(plotting_index(boxes)-1),'_H.svg']]);
            hold on;
            close;
        end
        close all;
        output = [plotting_index'-1, start_moving_time, first_slip, start_flowing, finish_flowing, finish_waiting, mov, mov_rate];
        %output = [plotting_index'-1, mov, mov_rate, start_moving_time, breaking_time, start_mov_rec, finish_mov_rec];
    end
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

function PP = overlaying_hatcurve (Ext, Turns, Range, Fig_handle, tile_handle)
figure (Fig_handle);
subplot(tile_handle)
yy = Ext(Range);
xx = Turns(Range);
xx_f = movmean(xx,10);
yy_f = movmean(yy,10);
PP = plot (xx_f, yy_f, 'Linewidth', 2);
hold on;
end

function [output_turns, output_ext] = binning_hat(Turns, Extension)
min_turns = fix(min(Turns));
max_turns = fix(max(Turns));
data = [Turns, Extension];
data_sort = sortrows(data,1);
Turns_sort = data_sort(:,1);
Extension_sort = data_sort(:,2);
output_turns = [];
output_ext = [];
for ii = min_turns:max_turns
    sub_index = abs(Turns_sort-ii)<=0.1;
    output_turns = [output_turns; ii];
    mean_ext = mean(Extension_sort(sub_index));
    output_ext = [output_ext; mean_ext];
end
end

function PP = plotting_initial_hatcurve(Ext, Turns, Range, Fig_handle, tile_handle)
[output_turns, output_ext] = binning_hat(Turns(Range), Ext(Range));
data = [output_turns, output_ext];
data_sort = sortrows(data,1);
figure (Fig_handle)
subplot(tile_handle)
PP = plot (data_sort(:,1), data_sort(:,2), 'Linewidth', 1.5);
hold on;
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

function fitresult_symm = fitSingle_m_Jin(turns, z)
%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( turns, z );

% Set up fittype and options.
ft = fittype( 'heaviside((n0-ns)-x)*(a*ns^2+b+k*(x+ns-n0))+heaviside(x-(n0+ns))*(a*ns^2+b-k*(x-ns-n0))+heaviside(-(n0-ns)+x)*heaviside(-x+(n0+ns))*(a*(x-n0)^2+b)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';


opts.Lower = [-0.007 0 0 -10 0];
opts.StartPoint = [-0.001 3.5 40 0 15];
opts.Upper = [0 4.5 100 10 25];

% Fit model to data.
[fitresult, ~] = fit( xData, yData, ft, opts );
field1 = 'a';
field2 = 'b';
field3 = 'km';
field4 = 'kp';
field5 = 'n0';
field6 = 'nsm';
field7 = 'nsp';
value1 = fitresult.a;
value2 = fitresult.b;
value3 = fitresult.k;
value4 = fitresult.k;
value5 = fitresult.n0;
value6 = fitresult.ns;
value7 = fitresult.ns;
fitresult_symm = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7);


end

function z = f_Single (fitresult, x)

a = fitresult.a;
b = fitresult.b;
km = fitresult.km;
kp = fitresult.kp;
n0 = fitresult.n0;
nsm = fitresult.nsm;
nsp = fitresult.nsp;
z = heaviside((n0-nsm)-x).*(a*nsm^2+b+km*(x+nsm-n0))+heaviside(x-(n0+nsp)).*(a*nsp^2+b-kp*(x-nsp-n0))+heaviside(-(n0-nsm)+x).*heaviside(-x+(n0+nsp)).*(a*(x-n0).^2+b);

end

function dp = nucleosome_depletion_calculation(front_edge, Nuc_position_bp, Nuc_info)

kk = 1:1:10;
kk = kk.';
Nuc_info_extended = ones(size(kk));

for ll = 1:length(Nuc_info)
    if Nuc_info(ll) == 0
        Nuc_info_extended(ll) = 0;
    end
end

exiting = cumsum(Nuc_info_extended);
entering = exiting-Nuc_info_extended;

entering_position = Nuc_position_bp+(kk-1)*197;
exiting_position = entering_position+147;



total_data_y = [entering; exiting];
total_data_x = [entering_position; exiting_position];

total_data = [total_data_x total_data_y];
total_data = [[0 0]; total_data];
total_data = sortrows(total_data);

dp = -1*interp1(total_data(:,1),total_data(:,2),front_edge,'linear');

end

function new_C0 = Predict_Hatcurve_Single (fitresult, Li, Lf, dp)
dn = (Li-Lf)/10.5+dp;
scaling_factor = (Lf-(dp+1)*147)/(Li-147);

a_new = fitresult.a;
b_new = scaling_factor*fitresult.b;
km_new = fitresult.km;
kp_new = fitresult.kp;
n0_new = fitresult.n0-dn;
nsm_new = scaling_factor*fitresult.nsm;
nsp_new = scaling_factor*fitresult.nsp;

field1 = 'a';
field2 = 'b';
field3 = 'km';
field4 = 'kp';
field5 = 'n0';
field6 = 'nsm';
field7 = 'nsp';
value1 = a_new;
value2 = b_new;
value3 = km_new;
value4 = kp_new;
value5 = n0_new;
value6 = nsm_new;
value7 = nsp_new;
new_C0 = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7);


end

function [T_pos, X_pos, T_neg, X_neg, Critical, Hat_record]...
    = Single_position_new_version(Time, Extension, task, magnet_turns_for_subplot, Nuc_position_bp)

Critial_position_bp = [Nuc_position_bp; Nuc_position_bp+147];
Li = 13014+Nuc_position_bp-403;
RNAP_downstream_protects = 0;

Timespot_index = timespot_finding (task);
Timespot = Time(Timespot_index);

pos_hat_1_range = (Timespot_index(1):Timespot_index(2))';
pos_hat_2_range = (Timespot_index(5):Timespot_index(6))';
neg_hat_range = (Timespot_index(3):Timespot_index(4))';
pos_hat_range = [pos_hat_1_range; pos_hat_2_range];

start_ext_range = (1:Timespot_index(1))';

[turns_p, ext_p] = binning_hat(magnet_turns_for_subplot(pos_hat_range),...
    Extension(pos_hat_range));
turns_p_index = find(turns_p>=-44 & turns_p<=69);
turns_p = turns_p(turns_p_index);
ext_p = ext_p(turns_p_index);
[turns_n, ext_n] = binning_hat(magnet_turns_for_subplot(neg_hat_range),...
    Extension(neg_hat_range));
turns_n_index = find(turns_n>=-44 & turns_n<=69);
turns_n = turns_n(turns_n_index);
ext_n = ext_n(turns_n_index);
ext_np = 0.5*(ext_n+ext_p);
C0_fit = fitSingle_m_Jin(turns_p, ext_np/1000);

Waiting_range_continuous = (Timespot_index(9):Timespot_index(10)).';
Waiting_turns_continuous = magnet_turns_for_subplot(Waiting_range_continuous);
Initial_turn_range = (Timespot_index(9):Timespot_index(10)).';
Initial_turns = mean(magnet_turns_for_subplot(Initial_turn_range));
%Final_turn_range = (Timespot_index(9):Timespot_index(10)).';
%Final_turns = mean(magnet_turns_for_subplot(Final_turn_range));
Time_waiting = Time(Waiting_range_continuous);
Extension_waiting = Extension(Waiting_range_continuous);
%unit_vec = ones(size(Time_waiting));

Waiting_turns = Initial_turns;
X_pos = zeros(size(Waiting_turns_continuous));
T_pos = X_pos;
X_neg = X_pos;
T_neg = X_pos;
Critical = zeros(length(X_pos),length(Critial_position_bp)+1);
Size_pos = 0;
Size_neg = 0;
Size_critical = 0;

dL1 = 0:1:1500;
dL2 = [];
dL1 = dL1.';
dL2 = dL2.';
dL1 = [dL1; Critial_position_bp; dL2];
dL1 = unique(dL1);
dL1 = sort(dL1);
Hat_record = zeros(length(dL1),8);
Hat_record(:,1) = dL1;
front_edge = dL1+RNAP_downstream_protects;
dp = nucleosome_depletion_calculation(front_edge, Nuc_position_bp, 1);

Ext_th = zeros(size(dL1));
critical_height = zeros(length(Critial_position_bp),1);
for ii = 1:length(dL1)
    L1 = Li-dL1(ii);
    new_C01 =  Predict_Hatcurve_Single (C0_fit, Li, L1, dp(ii));
    [~, new_C01_array] = Single_fit_Struct2Array(new_C01, 0);
    Hat_record(ii,2:end) = new_C01_array(1:end);
    Ext_th(ii) = 1000*f_Single(new_C01, Waiting_turns);
    critial_index = find(Critial_position_bp == dL1(ii));
    if ~isempty(critial_index)
        critical_height(critial_index) = Ext_th(ii);
    end
end

dL1 = [-1000; dL1];
Ext_th = [-100; Ext_th];
max_Ext_th = max(Ext_th);
apex_index = find(Ext_th == max_Ext_th);
Ext_th_pos = Ext_th(apex_index:end);
dL_pos = dL1(apex_index:end);
Ext_th_neg = Ext_th(1:apex_index);
dL_neg = dL1(1:apex_index);


data_temp_pos = [Ext_th_pos, dL_pos];
data_temp_pos = sortrows(data_temp_pos);
Ext_th_pos = data_temp_pos(:,1);
dL_pos = data_temp_pos(:,2);

data_temp_neg = [Ext_th_neg, dL_neg];
data_temp_neg = sortrows(data_temp_neg);
Ext_th_neg = data_temp_neg(:,1);
dL_neg = data_temp_neg(:,2);

Current_range = find(abs(Waiting_turns_continuous-Waiting_turns)<0.01);
Ext_current = Extension_waiting(Current_range);
Time_current = Time_waiting(Current_range);
Ext_current = min(Ext_current, max_Ext_th);


d_size_critical = length(Current_range);
insert_index_critical = ((Size_critical+1):(Size_critical+d_size_critical)).';
unit_vector_current = ones(size(insert_index_critical));
Critical(insert_index_critical, 1) = Time_current;
for jj = 1:size(Critial_position_bp)
    Critical(insert_index_critical,  jj+1) = critical_height(jj)*unit_vector_current;
end
Size_critical = Size_critical+d_size_critical;

X_pos_temp = interp1(Ext_th_pos,dL_pos,Ext_current,'linear');
T_pos_temp = Time_current-Timespot(9);
d_size_pos = length(X_pos_temp);
insert_index_pos = ((Size_pos+1):(Size_pos+d_size_pos)).';
X_pos(insert_index_pos) = X_pos_temp;
T_pos(insert_index_pos) = T_pos_temp;
Size_pos = Size_pos+d_size_pos;

X_neg_temp = interp1(Ext_th_neg,dL_neg,Ext_current,'linear');
T_neg_temp = Time_current-Timespot(9);
d_size_neg = length(X_neg_temp);
insert_index_neg = ((Size_neg+1):(Size_neg+d_size_neg)).';
X_neg(insert_index_neg) = X_neg_temp;
T_neg(insert_index_neg) = T_neg_temp;
Size_neg = Size_neg+d_size_neg;


X_pos = X_pos(1:Size_pos);
T_pos = T_pos(1:Size_pos);

X_neg = X_neg(1:Size_pos);
T_neg = T_neg(1:Size_pos);

Critical = Critical (1:Size_critical,:);

end

function [C0_output, C0_array_output] = Single_fit_Struct2Array(C0, C0_array)

if ~isa(C0,'struct') % Convert array to structure
    C0_array_output = C0_array;
    field1 = 'a';
    field2 = 'b';
    field3 = 'km';
    field4 = 'kp';
    field5 = 'n0';
    field6 = 'nsm';
    field7 = 'nsp';
    value1 = C0_array(1);
    value2 = C0_array(2);
    value3 = C0_array(3);
    value4 = C0_array(4);
    value5 = C0_array(5);
    value6 = C0_array(6);
    value7 = C0_array(7);
    C0_output = struct(field1,value1,field2,value2,field3,value3,...
        field4,value4,field5,value5,field6,value6,field7,value7);
else % Convert structure to array
    C0_output = C0;
    C0_array_output = zeros(7,1);
    C0_array_output(1) = C0.a;
    C0_array_output(2) = C0.b;
    C0_array_output(3) = C0.km;
    C0_array_output(4) = C0.kp;
    C0_array_output(5) = C0.n0;
    C0_array_output(6) = C0.nsm;
    C0_array_output(7) = C0.nsp;
end

end

function [center_y, center_y_theory, end_position_theory, end_position_C, QQ] = Fitting_N_Correction_Single(C0_fit,To_be_fitted,Hat_record)

turns_to_be_fitted = To_be_fitted(:,1);
Ext_to_be_fitted = To_be_fitted(:,2);

[fitresult, ~, ~] = Single_restricted_fit(turns_to_be_fitted, Ext_to_be_fitted, C0_fit);
center_x = fitresult.n0;
center_y = fitresult.b;
x = turns_to_be_fitted;
QQ = f_Single(fitresult, x);

[size_x, ~] = size(Hat_record);
RL_record = zeros(size_x, 3);
for ii = 1:size_x
    C_temp = Hat_record(ii,2:end);
    dL = Hat_record(ii,1);
    RL_record(ii,1) = dL;
    RL_record(ii,2) = C_temp(5);
    RL_record(ii,3) = C_temp(2);
end
center_y_theory = interp1(RL_record(:,2),RL_record(:,3),center_x);
end_position_theory = interp1(RL_record(:,2),RL_record(:,1),center_x);
end_position_index = find(Hat_record(:,1)>=floor(end_position_theory),1);
end_position_C_array = Hat_record(end_position_index, 2:end);
[end_position_C,~] = Single_fit_Struct2Array(0, end_position_C_array);
end

function [fitresult_symm, z_alpha, zz_delta_B] = Single_restricted_fit(turns_to_be_fitted, Ext_to_be_fitted, paras)
%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( turns_to_be_fitted, Ext_to_be_fitted );

% Set up fittype and options.
ft = fittype( 'zz_delta_B+heaviside((n0-(ns*z_alpha))-x)*(a*(ns*z_alpha)^2+(b*z_alpha)+k*(x+(ns*z_alpha)-n0))+heaviside(x-(n0+(ns*z_alpha)))*(a*(ns*z_alpha)^2+(b*z_alpha)-k*(x-(ns*z_alpha)-n0))+heaviside(-(n0-(ns*z_alpha))+x)*heaviside(-x+(n0+(ns*z_alpha)))*(a*(x-n0)^2+(b*z_alpha))',...
    'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';

a_th = paras.a;
b_th = paras.b;
k_th = paras.km;
ns_th = paras.nsm;
n0_th = paras.n0;

opts.Lower = [a_th b_th k_th -200 ns_th 0.75 -0.1];
opts.StartPoint = [a_th b_th k_th -80 ns_th 0.9 0];
opts.Upper = [a_th b_th k_th n0_th ns_th 1 0.1];

% Fit model to data.
[fitresult, ~] = fit( xData, yData, ft, opts );
z_alpha = fitresult.z_alpha;
zz_delta_B = fitresult.zz_delta_B;
field1 = 'a';
field2 = 'b';
field3 = 'km';
field4 = 'kp';
field5 = 'n0';
field6 = 'nsm';
field7 = 'nsp';
value1 = fitresult.a;
value2 = fitresult.b*z_alpha+zz_delta_B;
value3 = fitresult.k;
value4 = fitresult.k;
value5 = fitresult.n0;
value6 = fitresult.ns*z_alpha;
value7 = fitresult.ns*z_alpha;

fitresult_symm = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7);


end