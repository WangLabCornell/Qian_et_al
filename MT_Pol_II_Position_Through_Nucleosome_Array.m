clear;
close all;
clc;
%c = @cmu.colors;
% IMPORTANT! This program only works if there are two forces!

% This program is supposed to choose multiple files, but they should have
% the same protocols;

target_box = 103;

pickedcolor = [255, 0, 0]/255;
EndingT = 11;
tweezer2 = 0;
tweezer2 = 2;
resuming_index = 2;
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
        Address = 'Y:\Magnetic Tweezers 2\Jin\230914_Nucleosome Array 500 nM TFIIS F3 3\_230915_161314_PROCESSED\Data Files';
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
            Nuc_info = [1,1,1];
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
                string_3 = ['Trace been shifted down by ', num2str(extension_increase_during_flow, '%.0f'), ' nm.'];
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
                ylim ([0 2000]);
                h = gca;
                hold on;
                
                
                for ii = 1:length(Timespot)-1
                    time_i = Timespot (ii);
                    yy = 0:10:2000;
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
                    C0_fit = fit5piece_Jin(turns_p, ext_np/1000);
                    zz = f_5piece(C0_fit,turns_p);
                    plot (turns_p,zz*1000,'g--','Linewidth',3);
                    hold on;
                    
                    %% Plotting w+ and w-
                    r_plot = reshape(C0_fit,4,2);
                    x_plot = (r_plot(1,1):0.1:r_plot(4,1))';
                    [RL, RR, Rp, z_plot] = Hat_Curve_Skeleton(C0_fit,x_plot);
                    plot (x_plot, 1000*z_plot, 'k--', 'LineWidth', 1, 'HandleVisibility','off');
                    hold on;
                    y_L = 1000*RL(2):1:2000;
                    y_R = 1000*RR(2):1:2000;
                    y_p = 1000*Rp(2):1:2000;
                    x_L = RL(1)*ones(size(y_L));
                    x_R = RR(1)*ones(size(y_R));
                    x_p = Rp(1)*ones(size(y_p));
                    plot (x_L, y_L, 'k--', 'LineWidth', 1, 'HandleVisibility','off');
                    hold on;
                    plot (x_R, y_R, 'k--', 'LineWidth', 1, 'HandleVisibility','off');
                    hold on;
                    plot (x_p, y_p, 'k--', 'LineWidth', 1, 'HandleVisibility','off');
                    hold on;
                    text_y = mean(y_p);
                    text_x_v = RL(1):0.1:RR(1);
                    text_y_v = text_y*ones(size(text_x_v));
                    plot (text_x_v, text_y_v, 'k--', 'LineWidth', 1, 'HandleVisibility','off');
                    hold on;
                    w_m = RL(1)-Rp(1);
                    w_p = RR(1)-Rp(1);
                    text(RL(1)-5, text_y+80, num2str(w_m, '%.1f'));
                    text(Rp(1)+15, text_y+80, num2str(w_p, '%.1f'));
                    critical_points = [RL RR Rp]';
                    S0 = scatter (critical_points(:,1), 1000*critical_points(:,2));
                    S0.MarkerFaceColor = 'r';
                    S0.MarkerEdgeColor = 'r';
                    S0.SizeData = 20;
                    S0.HandleVisibility = 'off';
                    hold on;
                    %%
                    [Nnuc_h, Nnuc_w, flagg] = Initial_hat_Nnuc_Relaxed (C0_fit, Li);
                    string_1 = [string_1, ' (Good)'];
                    string_2 = ['N_{nuc, h} = ', num2str(Nnuc_h, '%.1f'),'; N_{nuc, w+} = : ', num2str(Nnuc_w, '%.1f')];
                    
                    if flagg==1
                        string_2 = [string_2, ' (Good)'];
                        pre_title = 'good';
                    else
                        string_2 = [string_2, ' (Bad)'];
                        pre_title = 'bad';
                    end
                    flagg = 1;
                    Nnuc = min(Nnuc_h,64);
                    legend_script = '5-piece Fit';
                    
                else
                    string_1 = [string_1, ' (Bad)'];
                    string_2 = 'Nucleosome number: N/A';
                    Nnuc = 50;
                    C0_fit = [-16.0859    1.0231
                        -1.5938    1.1486
                        38.0039    0.9968
                        57.1456    0.6842];
                    pre_title = 'bad';
                    zz = f_5piece(C0_fit,turns_p);
                    plot (turns_p,zz*1000,'c--','Linewidth',3);
                    hold on;
                    legend_script = 'Standard hat curve';
                    
                end
                
                xx = -200:0.1:200;
                dL1 = nucleosome_position+147; % # of bps from +1 site to the exit of 1st NPE
                dL1 = dL1-RNAP_downstream_protects;
                dn1 = dL1/10.5; % # of turns polymerase placed
                L1 = Li-dL1; % Remaining downstream template length
                dp1_enter = nucleosome_depletion_calculation(nucleosome_position, nucleosome_position, Nuc_info);
                dp1_exit = nucleosome_depletion_calculation(nucleosome_position+147, nucleosome_position, Nuc_info);
                new_C01 = Predict_Hatcurve (C0_fit, Li, L1, Nnuc, dp1_exit);
                new_C01_enter = Predict_Hatcurve (C0_fit, Li, L1+147, Nnuc, dp1_enter);
                yy1 = f_5piece(new_C01,xx);
                yy1_enter = f_5piece(new_C01_enter,xx);
                
                upper_curve1 = max([yy1; yy1_enter]);
                lower_curve1 = min([yy1; yy1_enter]);
                upper_curve1 = max(upper_curve1, 0);
                lower_curve1 = max(lower_curve1, 0);
                SP1 = patch([xx, flip(xx)], [upper_curve1, flip(lower_curve1)]*1000, [1 0 0]);
                SP1.FaceAlpha = 0.25;
                SP1.LineStyle = 'none';
                hold on;
                
                dL2 = dL1+197;
                dn2 = dL2/10.5;
                L2 = Li-dL2;
                dp2_enter = nucleosome_depletion_calculation(nucleosome_position+197, nucleosome_position, Nuc_info);
                dp2_exit = nucleosome_depletion_calculation(nucleosome_position+147+197, nucleosome_position, Nuc_info);
                new_C02 = Predict_Hatcurve (C0_fit, Li, L2, Nnuc, dp2_exit);
                new_C02_enter = Predict_Hatcurve (C0_fit, Li, L2+147, Nnuc, dp2_enter);
                yy2 = f_5piece(new_C02,xx);
                yy2_enter = f_5piece(new_C02_enter,xx);
                
                upper_curve2 = max([yy2; yy2_enter]);
                lower_curve2 = min([yy2; yy2_enter]);
                upper_curve2 = max(upper_curve2, 0);
                lower_curve2 = max(lower_curve2, 0);
                SP2 = patch([xx, flip(xx)], [upper_curve2, flip(lower_curve2)]*1000, [0 0 1]);
                SP2.FaceAlpha = 0.25;
                SP2.LineStyle = 'none';
                hold on;
                
                dL3 = dL1+197*2;
                dn3 = dL3/10.5;
                L3 = Li-dL3;
                dp3_enter = nucleosome_depletion_calculation(nucleosome_position+197*2, nucleosome_position, Nuc_info);
                dp3_exit = nucleosome_depletion_calculation(nucleosome_position+147+197*2, nucleosome_position, Nuc_info);
                new_C03 = Predict_Hatcurve (C0_fit, Li, L3, Nnuc, dp3_exit);
                new_C03_enter = Predict_Hatcurve (C0_fit, Li, L3+147, Nnuc, dp3_enter);
                yy3 = f_5piece(new_C03,xx);
                yy3_enter = f_5piece(new_C03_enter,xx);
                
                upper_curve3 = max([yy3; yy3_enter]);
                lower_curve3 = min([yy3; yy3_enter]);
                upper_curve3 = max(upper_curve3, 0);
                lower_curve3 = max(lower_curve3, 0);
                SP3 = patch([xx, flip(xx)], [upper_curve3, flip(lower_curve3)]*1000, [0 1 0]);
                SP3.FaceAlpha = 0.25;
                SP3.LineStyle = 'none';
                hold on;
                
                set(gca,'ColorOrderIndex',2)
                
                PP = plotting_initial_hatcurve(Extension_notfiltered, magnet_turns_for_subplot, pos_hat_range, FIG_x, tile_y);
                PP = plotting_initial_hatcurve(Extension_notfiltered, magnet_turns_for_subplot, neg_hat_range, FIG_x, tile_y);
                PP = plotting_initial_hatcurve (Extension_notfiltered, magnet_turns_for_subplot, Time_range_1 , FIG_x, tile_y);
                PP = overlaying_hatcurve (Extension_notfiltered, magnet_turns_for_subplot,Fast_winding_range, FIG_x, tile_y);
                PP.LineStyle = '--';
                PP.Color = 'b';
                set(gca,'ColorOrderIndex',4)
                %PP = plotting_initial_hatcurve (Extension_notfiltered, magnet_turns_for_subplot, Time_range_2 , FIG_x, tile_y);
                %PP.LineStyle = ':';
                %PP.HandleVisibility = 'off';
                ax = gca;
                ax.TitleFontSizeMultiplier = 0.8;
                title([string_1, char(10), string_2]);
                
                ylim ([0 2000]);
                xlim ([-200 80]);
                ylabel('Extension [nm]');
                xlabel('Magnet Turns');
                %
                
                lgd = legend (legend_script, 'Inside 1st NPE', 'Inside 2nd NPE','Inside 3rd NPE',...
                    'Initial hat curve (- to +)', 'Initial hat curve (+ to -)',...
                    'Hat curve at 20 minutes', 'Fast Hat curve after flowing', 'location','northwest');
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
                else
                    
                end
                
                
                if (strcmp(choice,'Yes'))
                    [T_pos, X_pos, T_neg, X_neg, Critical, Hat_record] = ...
                        Array_Position_TC_Continuous_Single_Turn_New_Alignment(Time, Extension, ...
                        task_for_timespot, magnet_turns_for_subplot, Nuc_info, nucleosome_position);
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
                    S1 = scatter (T_pos(end),end_position_theory );
                    S1.MarkerEdgeColor = 'm';
                    hold on;
                    xlabel ('Time [s]');
                    ylabel ('Pol II Position [bp]');
                    ylim ([0 1000]);
                    xlim ([0 1200]);
                    hold on;
                    
                    ax_all = gca;
                    box(ax_all,'off');
                    unit_vec = ones(size(T_c));
                    for ii = 1:3
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
                        choice_2 = questdlg(qstring_2,'Question?','Yes','Nope','Yes');
                        if (strcmp(choice_2,'Yes'))
                            
                            fig_flag = 0;
                            % open your file for writing
                            file_name_write = [new_folder,['\',new_file_name(1:14),...
                                'box ', num2str(plotting_index(boxes)-1),' Position.txt']];
                            file_name_trace = [new_folder,['\',new_file_name(1:14),...
                                'box ', num2str(plotting_index(boxes)-1),' Extension.txt']];
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
                        fig_flag = 1;
                    end
                    
                    if fig_flag >0
                        prompt = {'First Nucleosome:','Second Nucleosome:', 'Third Nucleosome:'};
                        title_prompt = 'NPE info';
                        num_line = 1;
                        defaultvals = {num2str(Nuc_info(1), '%.0f'), ...
                            num2str(Nuc_info(2), '%.0f'), num2str(Nuc_info(3), '%.0f')};
                        answer = inputdlg(prompt,title_prompt,num_line,defaultvals);
                        Nuc_info = cellfun(@str2num, answer)';
                        [x_end_0, ~] = ginput(1);
                        [end_turns, end_ext] = binning_hat...
                            ( magnet_turns_for_subplot(Time_range_1), Extension_notfiltered(Time_range_1));
                        end_turns_0 = magnet_turns_for_subplot(Timespot_index(10));
                        x_end = [end_turns_0, x_end_0];
                        end_select = end_turns<=max(x_end) & end_turns>=min(x_end);
                        end_turns_select = end_turns(end_select);
                        end_ext_select = end_ext(end_select)/1000;
                        To_be_fitted = [end_turns_select, end_ext_select];
                        [~, ~, ~, end_position_C]...
                            = Fitting_N_Adjusting_neg(C0_fit,To_be_fitted,Hat_record);
                        [kink_y, kink_y_theory, end_position_theory, end_position_C]...
                            = Fitting_N_Adjusting_neg(end_position_C,To_be_fitted,Hat_record);
                        extension_increase_during_flow = 1000*(kink_y-kink_y_theory);
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


function [RL, RR, Rp, z] = Hat_Curve_Skeleton(r0,x)
% r is the Catersian coordinates of 4 points that define the 5 piece
% function
r = reshape(r0,4,2);
% top piece
y3 = r(2,2) + (r(3,2)-r(2,2))/(r(3,1)-r(2,1)).*(x-r(2,1));
K_m = (r(3,2)-r(2,2))/(r(3,1)-r(2,1));
B_m = r(2,2)-r(2,1)*K_m;
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
Peak_position_x = - para_tl(2)/(2*para_tl(1));
Peak_position_y = para_tl(3)-(para_tl(2))^2/(4*para_tl(1));

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
K_n = alpha_n;
B_n = r(1,2)-K_n*r(1,1);

% lower right piece
alpha_p = 2 * para_tr(1) * r(4,1) + para_tr(2);
y5 = r(4,2) + alpha_p * (x-r(4,1));
K_p = alpha_p;
B_p = r(4,2)-K_p*r(4,1);

%% Calculating the intersects
x_L = -(B_n-B_m)/(K_n-K_m);
x_R = -(B_p-B_m)/(K_p-K_m);
y_L = K_m*x_L+B_m;
y_R = K_m*x_R+B_m;

RL = [x_L; y_L];
RR = [x_R; y_R];
Rp = [Peak_position_x; Peak_position_y];

%% whole curve
z = y1 .* (x<x_L) + y3 .* (x>=x_L & x<x_R) + y5.* (x>=x_R);

end

function [Nnuc_h, Nnuc_w, status] = Initial_hat_Nnuc_Relaxed(C0, Li)
Ls = 12667;
Cs = C0*Ls/Li;
[Height_s,~, ~, buckling_positive_s, ~, ~] = HCpara(Cs);
slope_s = -0.0414; % Seong ha's fit
intercept_s = 3.2306; % Seong ha's fit
W_slope_s = 0.597; % Seong ha's fit
W_intercept_s = 22.794; % Seong ha's fit
Nnuc_h_s =  (Height_s - intercept_s)/slope_s;

Nnuc_w_s = (buckling_positive_s-W_intercept_s)/W_slope_s;
Nnuc_h = Nnuc_h_s*Li/Ls;
Nnuc_w = Nnuc_w_s*Li/Ls;

Ws_lower = 0.521 * Nnuc_h_s + 19.260; % Seong ha's lower bound for width
Ws_upper = 0.673 * Nnuc_h_s + 26.330; % Seong ha's higher bound for width

if buckling_positive_s>=Ws_lower && buckling_positive_s<= Ws_upper
    status = 1; % Good array
else
    status = 0; % Not good
end

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

function new_C0 = Predict_Hatcurve (C0_i, Li, Lf, p_i, dp)
dn = (Li-Lf)/10.5;
L_standard = 12667;
C0_scaled = C0_i/Li*L_standard;
dp_scaled = dp/Lf*L_standard+p_i*L_standard*(1/Lf-1/Li);
new_C0_scaled = delta_hat_curve (C0_scaled, dp_scaled);
new_C0 = new_C0_scaled*Lf/L_standard;
new_C0(:,1) = new_C0(:,1)-dn;

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

function [T_pos, X_pos, T_neg, X_neg, Critical, Hat_record] = ...
    Array_Position_TC_Continuous_Single_Turn_New_Alignment...
    (Time, Extension, task, magnet_turns_for_subplot, Nuc_info, Nuc_position_bp)

Critial_position_bp = [Nuc_position_bp; Nuc_position_bp+147; ...
    Nuc_position_bp+197; Nuc_position_bp+197+147;...
    Nuc_position_bp+197*2; Nuc_position_bp+197*2+147];
Li = 13014+Nuc_position_bp-403;
RNAP_downstream_protects = 0;


Timespot_index = timespot_finding (task);
Timespot = Time(Timespot_index);

pos_hat_1_range = (Timespot_index(1):Timespot_index(2))';
pos_hat_2_range = (Timespot_index(5):Timespot_index(6))';
neg_hat_range = (Timespot_index(3):Timespot_index(4))';
pos_hat_range = [pos_hat_1_range; pos_hat_2_range];


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
C0_fit = fit5piece_Jin(turns_p, ext_np/1000);
[Nnuc_h, ~, ~] = Initial_hat_Nnuc_Relaxed (C0_fit, Li);
Nnuc = min(Nnuc_h,64);

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
Hat_record = zeros(length(dL1),9);
Hat_record(:,1) = dL1;
front_edge = dL1+RNAP_downstream_protects;
dp = nucleosome_depletion_calculation(front_edge, Nuc_position_bp, Nuc_info);

Ext_th = zeros(size(dL1));
critical_height = zeros(length(Critial_position_bp),1);
for ii = 1:length(dL1)
    L1 = Li-dL1(ii);
    new_C01 =  Predict_Hatcurve (C0_fit, Li, L1, Nnuc, dp(ii));
    Hat_record(ii,2:end) = new_C01(1:end);
    Ext_th(ii) = 1000*f_5piece(new_C01, Waiting_turns);
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

function [kink_y, kink_y_theory, end_position_theory, end_position_C, QQ] = Fitting_N_Adjusting_neg(C0_fit,To_be_fitted,Hat_record)

x0 = -200:0.1:200;
turns_to_be_fitted = To_be_fitted(:,1);
Ext_to_be_fitted = To_be_fitted(:,2);

paras = prepare_paras_for_3_piece_fitting_neg(C0_fit);
[fitresult, ~] = Three_piece_restricted_fit(turns_to_be_fitted, Ext_to_be_fitted, paras, 0);
zz_B = fitresult.zz_B;
zz_x_0 = fitresult.zz_x_0;
A_original = fitresult.A_original;
Y_alpha = fitresult.Y_alpha;
Y_beta = fitresult.Y_beta;
k_I = fitresult.k_I;
k_III = fitresult.k_III;
x_alpha = fitresult.x_alpha;
x_beta = fitresult.x_beta;

x = turns_to_be_fitted;
QQ = heaviside (x_alpha+zz_x_0-x).*( k_I*(x-(zz_x_0+x_alpha))+Y_alpha ) + heaviside (x-x_beta-zz_x_0).*( k_III*(x-(zz_x_0+x_beta))+Y_beta ) + heaviside (x-x_alpha-zz_x_0).* heaviside (x_beta+zz_x_0-x).*(A_original*(x-zz_x_0).^2)+zz_B;

kink_x = k_I*(zz_x_0+x_alpha)-k_III*(zz_x_0+x_beta)-(Y_alpha-Y_beta);
kink_x = kink_x/(k_I-k_III);
kink_y = k_I*(kink_x-(zz_x_0+x_alpha))+Y_alpha+zz_B;

[size_x, ~] = size(Hat_record);
RL_record = zeros(size_x, 3);
for ii = 1:size_x
    C_temp = Hat_record(ii,2:end);
    dL = Hat_record(ii,1);
    C_temp = reshape(C_temp, 4,2);
    [RL_temp, ~, ~, ~] = Hat_Curve_Skeleton(C_temp,x0);
    RL_record(ii,1) = dL;
    RL_record(ii, 2:end) = RL_temp.';
end
kink_y_theory = interp1(RL_record(:,2),RL_record(:,3),kink_x);
end_position_theory = interp1(RL_record(:,2),RL_record(:,1),kink_x);
end_position_index = find(Hat_record(:,1)>=floor(end_position_theory),1);
end_position_C = Hat_record(end_position_index, 2:end);
end_position_C = reshape(end_position_C, 4,2);
end

function [kink_y, kink_y_theory, end_position_theory, end_position_C, QQ] = Fitting_N_Adjusting_pos(C0_fit,To_be_fitted,Hat_record)

x0 = -200:0.1:200;
turns_to_be_fitted = To_be_fitted(:,1);
Ext_to_be_fitted = To_be_fitted(:,2);

paras = prepare_paras_for_3_piece_fitting_pos(C0_fit);
[fitresult, ~] = Three_piece_restricted_fit(turns_to_be_fitted, Ext_to_be_fitted, paras, -50);
zz_B = fitresult.zz_B;
zz_x_0 = fitresult.zz_x_0;
A_original = fitresult.A_original;
Y_alpha = fitresult.Y_alpha;
Y_beta = fitresult.Y_beta;
k_I = fitresult.k_I;
k_III = fitresult.k_III;
x_alpha = fitresult.x_alpha;
x_beta = fitresult.x_beta;

x = turns_to_be_fitted;
QQ = heaviside (x_alpha+zz_x_0-x).*( k_I*(x-(zz_x_0+x_alpha))+Y_alpha ) + heaviside (x-x_beta-zz_x_0).*( k_III*(x-(zz_x_0+x_beta))+Y_beta ) + heaviside (x-x_alpha-zz_x_0).* heaviside (x_beta+zz_x_0-x).*(A_original*(x-zz_x_0).^2)+zz_B;

kink_x = k_I*(zz_x_0+x_alpha)-k_III*(zz_x_0+x_beta)-(Y_alpha-Y_beta);
kink_x = kink_x/(k_I-k_III);
kink_y = k_I*(kink_x-(zz_x_0+x_alpha))+Y_alpha+zz_B;

[size_x, ~] = size(Hat_record);
RR_record = zeros(size_x, 3);
for ii = 1:size_x
    C_temp = Hat_record(ii,2:end);
    dL = Hat_record(ii,1);
    C_temp = reshape(C_temp, 4,2);
    [~, RR_temp, ~, ~] = Hat_Curve_Skeleton(C_temp,x0);
    RR_record(ii,1) = dL;
    RR_record(ii, 2:end) = RR_temp.';
end
kink_y_theory = interp1(RR_record(:,2),RR_record(:,3),kink_x);
end_position_theory = interp1(RR_record(:,2),RR_record(:,1),kink_x);
end_position_index = find(Hat_record(:,1)>=floor(end_position_theory),1);
end_position_C = Hat_record(end_position_index, 2:end);
end_position_C = reshape(end_position_C, 4,2);
end

function [paras] = prepare_paras_for_3_piece_fitting_neg(C0_fit)
x1 = C0_fit(1,1);
y1 = C0_fit(1,2);
x2 = C0_fit(2,1);
y2 = C0_fit(2,2);
x3 = C0_fit(3,1);
y3 = C0_fit(3,2);

dy_12 = y1-y2;
dx_12 = x1-x2;
k_32 = (y3-y2)/(x3-x2);

A_original = ( dy_12 - k_32*dx_12 )/dx_12^2;
B_original = -2*A_original*x2+k_32;
para_center_original = -B_original/(2*A_original);
dx_1p = x1-para_center_original;
dx_2p = x2-para_center_original;
k_I = ( 2*(dy_12 - k_32*dx_12 )/dx_12 + k_32 );
k_III = k_32;
x_alpha = dx_1p;
x_beta = dx_2p;
Y_alpha = A_original*x_alpha^2;
Y_beta = A_original*x_beta^2;
paras = [A_original Y_alpha Y_beta  k_I  k_III  x_alpha  x_beta];
end

function [paras] = prepare_paras_for_3_piece_fitting_pos(C0_fit)
x2 = C0_fit(2,1);
y2 = C0_fit(2,2);
x3 = C0_fit(3,1);
y3 = C0_fit(3,2);
x4 = C0_fit(4,1);
y4 = C0_fit(4,2);

dy_43 = y4-y3;
dx_43 = x4-x3;
k_32 = (y3-y2)/(x3-x2);

A_original = ( dy_43 - k_32*dx_43 )/dx_43^2 ;
B_original = -2*A_original*x3+k_32;
para_center_original = -B_original/(2*A_original);
dx_3p = x3-para_center_original;
dx_4p = x4-para_center_original;
k_I = k_32;
k_III = ( 2*(dy_43 - k_32*dx_43 )/dx_43 + k_32 );
x_alpha = dx_3p;
x_beta = dx_4p;
Y_alpha = A_original*x_alpha^2;
Y_beta = A_original*x_beta^2;
paras = [A_original Y_alpha Y_beta  k_I  k_III  x_alpha  x_beta];
end

function [fitresult, gof] = Three_piece_restricted_fit(turns_to_be_fitted, Ext_to_be_fitted, paras, dx_0)

[xData, yData] = prepareCurveData( turns_to_be_fitted, Ext_to_be_fitted );

% Set up fittype and options.
ft = fittype( 'heaviside(x_alpha+zz_x_0-x)*(k_I*(x-(zz_x_0+x_alpha))+Y_alpha)+heaviside(x-x_beta-zz_x_0)*(k_III*(x-(zz_x_0+x_beta))+Y_beta)+heaviside(x-x_alpha-zz_x_0)*heaviside(x_beta+zz_x_0-x)*(A_original*(x-zz_x_0)^2)+zz_B', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [paras 0 -150];
opts.StartPoint = [paras 1 -10+dx_0];
opts.Upper = [paras 2 10];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

end
