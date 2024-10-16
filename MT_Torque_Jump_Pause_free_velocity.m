clear;
close all;
clc;
%c = @cmu.colors;
% IMPORTANT! This program only works if there are two forces!

% This program is supposed to choose multiple files, but they should have
% the same protocols;

target_box = [
   11
78
84
];

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
        Address = 'Y:\Magnetic Tweezers 2\Jin';
    end
    
    [filename, pathname] = uigetfile('*.*','Select the data file',Address); % Choose MT data files
    if filename == 0 % If you do not choose a file, then the program proceeds;
        break;
    else flag = 1;
    end
    new_file_name = strrep(filename,'_',' ');
    fileID = fopen ([pathname,filename],'r'); % Open the data file
    file_index = file_index + 1;
    new_folder = [pwd,'/',[filename,' plots_No Shift']];
    mkdir([filename,' plots_No Shift']);
    
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
    Start_force = return_matrix(3,:); % This is the first calibrated force
    Jump_force = return_matrix(2,:); % This is the second calibrated force
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
    id_A = filename(1:6);
    id_B = filename(8:13);
    id_A_num = round(str2double(id_A));
    id_B_num = round(str2double(id_B));
    
    
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
        date_rec = id_A_num+mov;
        DAQ_time_rec = id_B_num+mov;
        force_rec = mov;
        output = zeros(length(mov), 16)-100;
        
        EndingT = length(Timespot);
        
        Time_range_0 = (Timespot_index(1):Timespot_index(6))';
        Time_range_1 = (Timespot_index(10):Timespot_index(11))';
        Time_range_2 = (Timespot_index(13):Timespot_index(14))';
        Fast_winding_range = (Timespot_index(14):Timespot_index(15))';
        pos_hat_range = (Timespot_index(3):Timespot_index(4))';
        neg_hat_range = (Timespot_index(5):Timespot_index(6))';
                
        start_ext_range = (1:Timespot_index(1))';
        
        
        target_index_tf = ismember(plotting_index,target_box+1);
        target_index = 1:length(plotting_index);
        target_index = target_index(target_index_tf);

        for boxes = resuming_index:length(plotting_index)
        %for boxes = 1:10
            
            Nuc_info = [0];
            fig_flag = 2;
            clc;
            disp(num2str(100*boxes/length(plotting_index), '%.1f'));
            total_shift = 0;
            
            FIG_x = figure (2*boxes-1);
            pos = [100 100 1500 700];
            set(FIG_x, 'Pos', pos);
            %% Plotting subplots: Magnetic height and magnetic turns
            
            subplot(6, 5, [1 3]);
            %subplot_M_turns = subplot (5,1,1);
            plot (time_for_timespot, magnet_turns_for_subplot,'k','LineWidth',1.5);
            axis([0 inf -inf inf])
            ylabel('M-turns');
            xlim ([0 Timespot(EndingT)]);
            current_force = Jump_force(boxes);
            force_rec(boxes) = current_force;
            string_3 = ['Force: ', num2str(current_force, '%.2f'), ' pN'];
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
            
            Extension = 1000*Data_all_cell{boxes+length(essential_columns_position)};
            Extension_notfiltered = Extension;
            Extension = movmean(Extension,20);
            
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
            
            string_1 = ['Mean hysteresis = ', num2str(abs(dEXT_during), '%.0f'), ' nm.'];
            flagg = 0;
            
            if dEXT_during<=50
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
                value6 = 8;
                value7 = 8;
                
                C0_fit = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7);
                
                zz = f_Single(C0_fit,turns_p);
                
                string_1 = [string_1, ' (Bad)'];
                plot (turns_p,zz*1000,'c--','Linewidth',1.5);
                hold on;
                pre_title = 'Bad';
                legend_script = 'Standard hat curve';
                
            end
                        
            xx = -100:0.1:50;
            dL1 = nucleosome_position+147; % # of bps from +1 site to the exit of 1st NPE
            dL1 = dL1-RNAP_downstream_protects;
            dn1 = dL1/10.5; % # of turns polymerase placed
            L1 = Li-dL1; % Remaining downstream template length
            dp1_enter = nucleosome_depletion_calculation(nucleosome_position, nucleosome_position, Nuc_info);
            dp1_exit = nucleosome_depletion_calculation(nucleosome_position+147, nucleosome_position, Nuc_info);
            new_C01 =  Predict_Hatcurve_Single (C0_fit, Li, L1, dp1_exit);
            new_C01_enter =  Predict_Hatcurve_Single (C0_fit, Li, L1+147, dp1_enter);
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
            
            ylim ([0 2000]);
            xlim ([-25 35]);
            ylabel('Extension [nm]');
            xlabel('Magnet Turns');
            lgd = legend (legend_script, 'Inside NPE',...
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
            
            
            qstring = 'Good Trace?';
            choice = questdlg(qstring,'Question?','Yes','Nope','Yes');
            extension_increase_during_flow = 0;
            
            
            if (strcmp(choice,'Yes'))
                [T_pos, X_pos, T_neg, X_neg, Critical, Hat_record] = ...
                Stall_torque_new_alignment(Time, Extension, ...
                task_for_timespot, magnet_turns_for_subplot, nucleosome_position);
                T_pos = T_pos(~isnan(X_pos));
                X_pos = X_pos(~isnan(X_pos));
                %X_pos = movmean(X_pos, 40);
                T_neg = T_neg(~isnan(X_neg));
                X_neg = X_neg(~isnan(X_neg));
                %X_neg = movmean(X_neg, 40);
                T_c = Critical(:,1);
                
                subplot(6, 5, [19 30]);
                plot (T_pos, movmean(X_pos, 40), 'k', 'LineWidth', 3);
                hold on;
                plot (T_neg, movmean(X_neg, 40), 'b', 'LineWidth', 3);
                hold on;
                
                hold on;
                xlabel ('Time [s]');
                ylabel ('Pol II Position [bp]');
                ylim ([0 600]);
                xlim ([0 1300]);
                hold on;
                
                ax_all = gca;
                box(ax_all,'off');
                unit_vec = ones(size(T_c));
                
                hold off;
                set(gca,'FontSize',13,'LineWidth',1.5);
                set(gca, 'Color', 'None');
                ax_all = gca;
                box(ax_all,'off');
                
                saveas(gcf,[new_folder,['\',pre_title,'_Box ', num2str(plotting_index(boxes)-1),'_H.fig']]);
                saveas(gcf,[new_folder,['\',pre_title,'_Box ', num2str(plotting_index(boxes)-1),'_H.png']]);
                saveas(gcf,[new_folder,['\',pre_title,'_Box ', num2str(plotting_index(boxes)-1),'_H.svg']]);
                hold on;
                close;
                
                
                parameters_temp = replisome_position_vs_time_adapted_to_Jin_func_v2...
                    ([T_pos, X_pos], 0.3, [new_file_name(1:14),...
                    'box ', num2str(plotting_index(boxes)-1), ' Pause analysis']);
                
                file_name_write = [new_folder,['\',new_file_name(1:14),...
                    'box ', num2str(plotting_index(boxes)-1),' Position.txt']];

                % write the matrix

                myData = zeros(length(T_pos), 3)-100;
                myData(:,1) = T_pos;
                myData(:,2) = X_pos;
                myData(1,3) = id_A_num;
                myData(2,3) = id_B_num;
                myData(3,3) = current_force;
                myData(4:(length(parameters_temp)+3),3) = parameters_temp.';
                output(boxes, 5:end) = parameters_temp;

                myData =myData.';
                myData_trace =[Time Extension_notfiltered/1000 ...
                    task_for_timespot magnet_turns_for_subplot].';
                fid_write = fopen(file_name_write,'wt');
                
                if fid_write > 0
                    fprintf(fid_write,'%f\t%f\t%f\n',myData);
                    fclose(fid_write);
                end
            else
                saveas(gcf,[new_folder,['\',pre_title,'_Box ', num2str(plotting_index(boxes)-1),'_H.fig']]);
                saveas(gcf,[new_folder,['\',pre_title,'_Box ', num2str(plotting_index(boxes)-1),'_H.png']]);
                saveas(gcf,[new_folder,['\',pre_title,'_Box ', num2str(plotting_index(boxes)-1),'_H.svg']]);
                hold on;
                close;
            end
            
            output(boxes, 1:4) = [id_A_num, id_B_num, plotting_index(boxes)-1, current_force];
            
            %}
            
        end
        close all;
        %output = [date_rec, DAQ_time_rec, plotting_index'-1, ];
       
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

function fitresult_symm = fitSingle_m_Jin(turns, z)
%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( turns, z );

% Set up fittype and options.
ft = fittype( 'heaviside((n0-ns)-x)*(a*ns^2+b+k*(x+ns-n0))+heaviside(x-(n0+ns))*(a*ns^2+b-k*(x-ns-n0))+heaviside(-(n0-ns)+x)*heaviside(-x+(n0+ns))*(a*(x-n0)^2+b)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';


opts.Lower = [-0.007 0 0 -2 3];
opts.StartPoint = [-0.001 1.3 40 0 7];
opts.Upper = [0 2.5 100 2 12];

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
    = Stall_torque_new_alignment(Time, Extension, task, magnet_turns_for_subplot, Nuc_position_bp)

Critial_position_bp = [Nuc_position_bp; Nuc_position_bp+147];
Li = 5000;
RNAP_downstream_protects = 0;

Timespot_index = timespot_finding (task);
Timespot = Time(Timespot_index);           

neg_hat_range = (Timespot_index(5):Timespot_index(6))';
pos_hat_range = (Timespot_index(3):Timespot_index(4))';

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

Waiting_range_continuous = (Timespot_index(17):Timespot_index(18)).';
Waiting_turns_continuous = magnet_turns_for_subplot(Waiting_range_continuous);
Initial_turn_range = (Timespot_index(15):Timespot_index(16)).';
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

dL1 = 0:1:800;
dL2 = [];
dL1 = dL1.';
dL2 = dL2.';
dL1 = [dL1; Critial_position_bp; dL2];
dL1 = unique(dL1);
dL1 = sort(dL1);
Hat_record = zeros(length(dL1),8);
Hat_record(:,1) = dL1;
front_edge = dL1+RNAP_downstream_protects;


Ext_th = zeros(size(dL1));
critical_height = zeros(length(Critial_position_bp),1);
for ii = 1:length(dL1)
    L1 = Li-dL1(ii);
    new_C01 =  Predict_Hatcurve_Single (C0_fit, Li, L1, 0);
    [~, new_C01_array] = Single_fit_Struct2Array(new_C01, 0);
    Hat_record(ii,2:end) = new_C01_array(1:end);
    Ext_th(ii) = 1000*f_Single(new_C01, Waiting_turns);
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
    T_pos_temp = Time_current-Timespot(17);
    d_size_pos = length(X_pos_temp);
    insert_index_pos = ((Size_pos+1):(Size_pos+d_size_pos)).';
    X_pos(insert_index_pos) = X_pos_temp;
    T_pos(insert_index_pos) = T_pos_temp;
    Size_pos = Size_pos+d_size_pos;
    
    X_neg_temp = interp1(Ext_th_neg,dL_neg,Ext_current,'linear');
    T_neg_temp = Time_current-Timespot(17);
    d_size_neg = length(X_neg_temp);
    insert_index_neg = ((Size_neg+1):(Size_neg+d_size_neg)).';
    X_neg(insert_index_neg) = X_neg_temp;
    T_neg(insert_index_neg) = T_neg_temp;
    Size_neg = Size_neg+d_size_neg;
    

X_pos = X_pos(1:Size_pos);
T_pos = T_pos(1:Size_pos)+Timespot(17);

X_neg = X_neg(1:Size_pos);
T_neg = T_neg(1:Size_pos)+Timespot(17);

Critical = Critical (1:Size_critical,:);

end

function output = replisome_position_vs_time_adapted_to_Jin_func_v2(Test_Data, dwell_thresh, file_name)


x_raw_temp = Test_Data(:,1)-Test_Data(1,1);
y_raw_temp = Test_Data(:,2);

fig_0 = figure (1000);
fig_0.Position = [10 100 1800 800];
set(gca,'FontSize',13,'LineWidth',1.5);
set(gca, 'Color', 'None');
ax_all = gca;
box(ax_all,'off');
hold on
plot(x_raw_temp,y_raw_temp, 'b', 'LineWidth',2);
ylabel('Pol II Position (bp)')
xlabel('Time (s)')
ylim([0 600]);
xlim([0 800]);
hold on;
[x_c, ~] = ginput(2);
close(fig_0);
x_c = sort(x_c);
in_range = x_raw_temp>=x_c(1) & x_raw_temp<=x_c(2);
x_raw = x_raw_temp(in_range);
y_raw = y_raw_temp(in_range);
x_raw = x_raw-x_raw(1);
time = y_raw;

indirect_torque_raw = zeros(size(time));
dt = (x_raw(101)-x_raw(1))/100;

%Trim down the data to the max position
[y_max,I_max] = max(y_raw);
I_end = length(y_raw);

%May need to tweak the ending parameter

%Detect fork regression as a maximum followed by regression of at least 100
%bp.  If found, cut the trace at the maximum.
[y_min,~] = min(y_raw(I_max:end));
if(y_max-y_min > 50)
    I_end = I_max;
end

I_start = find(indirect_torque_raw<10,1,"first");
if(I_end <= I_start)
    I_end = I_start+1;
end
%I_start
%I_end
y = y_raw(I_start:end);
x = x_raw(I_start:end);
indirect_torque = indirect_torque_raw(I_start:end);

% Tweak parameters

time_cutoff = x(end);
time_cutoff2 = x(end);


%[y_smoothed,v] = calculate_velocity_position(x,y,0.5);



% Make figure and define subplots
fig = figure(1001);
fig.Position = [10 100 1800 800];

% Plot position vs time
ax_position = subplot(1,5,[1 3]);
title(file_name);
set(gca,'FontSize',13,'LineWidth',1.5);
set(gca, 'Color', 'None');
ax_all = gca;
box(ax_all,'off');
hold on
P_a = plot(x_raw,y_raw, 'b', 'LineWidth',2);
ylabel('Pol II Position (bp)')
xlabel('Time (s)')
xlim([0 150]);
ylim([0 600]);
hold on;
%title(trace_data.output_file_name, "Interpreter","none");

y_smoothed = smoothdata(y,'gaussian',80);
y_smoothed_original = y_smoothed;
y_smoothed = Transcript_length_no_TFIIS(y_smoothed);
plot(x,y_smoothed, 'g', 'LineWidth',3);
hold on;

% Dwell time histogram
ax_dwell = subplot(1,5,4);
set(gca,'FontSize',13,'LineWidth',1.5);
set(gca, 'Color', 'None');
ax_all = gca;
box(ax_all,'off');
hold on;
bins = 0:600;
h = histogram(y_smoothed,bins,'Orientation', 'horizontal');
h.BinCounts = h.BinCounts*dt;
h.LineStyle = 'none';
h.FaceColor = 'b';
xline(dwell_thresh,'--');
xlabel('Dwell Time (s/bp)')
title('Dwell Histogram');
xlim([0 3]);
ylim ([0 600]);



pauses = get_pauses(h.BinCounts,bins,dwell_thresh);


% Create and index array to select all the regions with pausing.
index = false(size(y_smoothed));

for i = 1:length(pauses)
    min_index = find((y_smoothed > pauses(i).start),1);
    max_index = find((y_smoothed <= pauses(i).end),1,'last');
    index = index | ((x > x(min_index)) & (x < x(max_index)));
end
subplot(ax_position)
hold on;
plot(x(index),y_smoothed(index),'.','MarkerSize',10,'Color','Red')
%tracedata.paused = index;
xline(time_cutoff);


% Find pauses less than cutoff:

index_cutoff = find(x < time_cutoff,1,'last');
[pause_regions, numRegions] = bwlabel(index(1:index_cutoff));
pause_durations = zeros(numRegions,1);
pause_start_index = zeros(numRegions,1);
for i = 1:numRegions
    x_region = x(pause_regions == i);
    pause_durations(i) = x_region(end)-x_region(1);
    pause_start_index(i) = find(pause_regions == i,1,'first');
end

pause_duration = mean(pause_durations);

index_cutoff = find(x < time_cutoff2,1,'last');
[active_regions, numRegions] = bwlabel(~index(1:index_cutoff));
active_distances = zeros(numRegions,1);

for i = 1:numRegions
    
    y_region = y_smoothed(active_regions == i);
    active_distances(i) = y_region(end)-y_region(1);

end

active_distance = mean(active_distances);


%Find the procesivity
processivity = max(y)-y(1);

% Find the overall velocity from a linear fit to the position vs time.

durations = [20, 30, 40, 60];

subplot(ax_position)
hold on


overall_velocities = zeros(size(durations));
overall_velocities_direct = overall_velocities;
overall_torques = zeros(size(durations));
for i = 1:length(durations)
    time_end = x(1)+durations(i);
    i_end = find(x>time_end,1,"first");
    if isempty(i_end) 
        i_end = length(x);
    end
    overall_x = x(1:i_end);
    overall_y = y_smoothed_original(1:i_end);
    fit = polyfit(overall_x,overall_y,1);
    fit_y = polyval(fit,overall_x);
    plot(overall_x,fit_y, '--', 'LineWidth',2, 'Color', 'Black');
    overall_velocities(i) = fit(1);
    overall_velocities_direct(i) = (overall_y(end) - overall_y(1))/(overall_x(end) - overall_x(1));
    %overall_torques(i) = mean(indirect_torque(1:i_end));
end


plot_text = [...
    sprintf('Processivity = %0.1f bp', processivity)...
    newline...
    sprintf('Overall Velocity (%d s) = %0.1f bp/s', durations(1), overall_velocities(1))...
    newline...
    sprintf('Overall Velocity (%d s) = %0.1f bp/s', durations(2), overall_velocities(2))...
    newline...
    sprintf('Overall Velocity (%d s) = %0.1f bp/s', durations(3), overall_velocities(3))...
    newline...
    sprintf('Overall Velocity (%d s) = %0.1f bp/s', durations(4), overall_velocities(4))...
    ];
text(50,300,plot_text);



% What is this?
[active_regions, numRegions] = bwlabel(~index);

tstart = x_raw(I_start);

% Add the vertical lines to all the plots.
subplot(ax_position)
yline(y(1));
yline(y(1)+processivity);
xline(tstart)
for i = 1:length(durations)
    xline(x(1)+durations(i),'--');
end



ax_velocity_regions = subplot(1,5,5);


velocity_regions = [];
for i = 1:numRegions
    region = active_regions == i;
    x_cropped = x(region);
    y_cropped = y_smoothed(region);

    fit = polyfit(x_cropped,y_cropped,1);
    velocity_regions(i).distance = polyval(fit,x_cropped(end)) - polyval(fit,x_cropped(1));
    velocity_regions(i).velocity = fit(1);
    velocity_regions(i).indirect_torque = mean(indirect_torque(region));
end
hold off;
plot([velocity_regions.distance],[velocity_regions.velocity],'b.','MarkerSize',30);
average_velocity = sum([velocity_regions.distance].*[velocity_regions.velocity])/sum([velocity_regions.distance]);
average_torque = sum([velocity_regions.distance].*[velocity_regions.indirect_torque])/sum([velocity_regions.distance]);
title(['Weighted Mean = ', num2str(average_velocity, '%.1f') ,' bp/s']);
ylabel('Velocity of region (bp/s)');
xlabel('Distance of active region (bp)');
set(gca,'FontSize',13,'LineWidth',1.5);
set(gca, 'Color', 'None');
ax_all = gca;
box(ax_all,'off');
%legend(sprintf("Weighted Mean = %0.1f bp/s" + newline + "Weighted torque = %0.1f pNnm",average_velocity, average_torque),'Location','Northwest');

trace_data.velocity.pause_free_velocity = average_velocity;
trace_data.velocity.pause_free_torque = average_torque;

set(findall(gcf,'-property','FontSize'),'FontSize',18);
saveas(gcf,[file_name, '.fig']);
saveas(gcf,[file_name, '.png']);
hold off;
close all;
output = [overall_velocities, overall_velocities_direct, average_velocity, processivity, ...
    x_c(1)+Test_Data(1,1), x_c(2)+Test_Data(1,1)];

end

function output = Transcript_length_no_TFIIS(x)
    NN = length(x);
    output = x;
    for jj = 2:NN
        output(jj) = max(output(1:jj));
    end
end

function pauses = get_pauses(dwell_hist,bins,dwell_thresh)
    
    ind = dwell_hist > dwell_thresh;
    
    % Label regions with unique value
    [labeledVector, numRegions] = bwlabel(ind);
    
    pauses = [];
    for i = 1:numRegions
        bins_selected = bins(labeledVector == i);
        pauses(i).start = min(bins_selected);
        pauses(i).end = max(bins_selected)+1;
    end
   
end