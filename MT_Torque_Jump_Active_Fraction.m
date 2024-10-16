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
smoothing_time = 10;
L0 = 5000;
Minimal_moving_distance = 10.5;
Cutoff_time = [30; 60; 120; 180];

flag = 0; % To indicate if we choose any files at all, even if we choose one file, this would be 1;
file_index = 0; % Number of files chosen;
key = '1 Wind';
Maximum_allowed_files = 1;
Maximum_allowed_boxes = 80;
Data_points_number = zeros(Maximum_allowed_files,1);
%Data_all = zeros (32768,1024,Maximum_allowed_files);
position_and_force = zeros(4,Maximum_allowed_boxes,Maximum_allowed_files);

essential_columns = {{'DAQ time (s)'}; {'Time (ms)'}; {'tasklist (index)'};{'Magnet Angle (Turns)'};{'Magnet Height (mm)'}};% This is the number of columns essential for our hatcurve;


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
    %filename = '240229_180649_Extended force jump 6 pNnm no TFIIS.txt';
    %pathname = 'Y:\Magnetic Tweezers 2\Jin\240229_Extended force jump 6 pNnm no TFIIS\_240301_162002_PROCESSED\Data Files\';
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
    Turns = Data_all_cell{4}; % The magnetic turns for the first file selected, used for plotting subplot;
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
        torque_rec = mov;

        pos_hat_range = (Timespot_index(3):Timespot_index(4))';
        neg_hat_range = (Timespot_index(5):Timespot_index(6))';
                
        start_ext_range = (1:Timespot_index(1))';
        
        

        for boxes = resuming_index:length(plotting_index)
        %for boxes = 1:10
            
            
            clc;
            disp(num2str(100*boxes/length(plotting_index), '%.1f'));
            
            FIG_x = figure (boxes);
            pos = [100 100 1300 700];
            set(FIG_x, 'Pos', pos);
            TILE = tiledlayout(5,3, 'Padding', 'none', 'TileSpacing', 'compact');
            current_force = Jump_force(boxes);
            force_rec(boxes) = current_force;
            current_torque = 11.9946213491175*(current_force).^0.656556817954428;
            torque_rec(boxes) = current_torque;
            string_force = ['Force: ', num2str(current_force, '%.2f'), ' pN, Torque: ',...
                num2str(current_torque, '%.1f'), ' pN{\cdot}nm'];
            title (TILE,[new_file_name(1:(length(filename)-4)),' box ',num2str(plotting_index(boxes)-1)...
                char(10), string_force]);
            %% Plotting subplots: Magnetic turns
            nexttile(1, [1, 2]);
            
            P_magnet_turns = plot (time_for_timespot, Turns,'k','LineWidth',4);
            hold on;
            yLimits = get(gca,'YLim');
            xlim ([0 max(Time)]);
            ylabel('M-turns');
            for ii = 1:length(Timespot)-1
                time_i = Timespot (ii);
                yy = yLimits(1):0.1:yLimits(2);
                xx = time_i*ones(size(yy));
                plot(xx,yy,'k--','LineWidth',0.5);
                hold on;
            end
            
            ax_all = gca;
            ax_all.XAxis.Visible = 'off';
            set(gca,'FontSize',15,'LineWidth',1.5);
            set(gca,'Color', 'None');
            set(gca,'TickDir','out');
            set(gca,'TickLength',[0.015, 0.01]);
            set(gca,'XColor','k', 'YColor','k');
            set(gca, 'Layer', 'top');
            box(ax_all,'off');
            %% Plotting subplots: Magnet heights
            nexttile(4, [1,2]);
            P_magnet_heights = plot (time_for_timespot, magnet_height_for_subplot,'k','LineWidth',4);
            hold on;
            ylabel('M-height (mm)');
            xlim ([0 max(Time)]);
            ylim ([0 7]);
            for ii = 1:length(Timespot)-1
                time_i = Timespot (ii);
                yy = yLimits(1):0.1:yLimits(2);
                xx = time_i*ones(size(yy));
                plot(xx,yy,'k--','LineWidth',0.5);
                hold on;
            end
            
            ax_all = gca;
            ax_all.XAxis.Visible = 'off';
            set(gca,'FontSize',15,'LineWidth',1.5);
            set(gca,'Color', 'None');
            set(gca,'TickDir','out');
            set(gca,'TickLength',[0.015, 0.01]);
            set(gca,'XColor','k', 'YColor','k');
            set(gca, 'Layer', 'top');
            box(ax_all,'off');
    
            %% Plotting all the good traces
            
            nexttile(7,[3 2]);
            
            Extension = 1000*Data_all_cell{boxes+length(essential_columns_position)};
            Extension_notfiltered = Extension;
            Extension = Smoothed_trace(Time, Extension, task_for_timespot, smoothing_time);
           
            
            plot (Time, Extension_notfiltered, 'Color', [0.7,0.7,0.7]);
            hold on;
            plot (Time, Extension, 'Color', [0,0,1], 'Linewidth', 4);
            xlim ([0 max(Time)]);
            ylim ([0 1500]);
            hold on;
            
            for ii = 1:length(Timespot)-1
                time_i = Timespot (ii);
                yy = 0:0.01:4;
                yy = yy*1000;
                xx = time_i*ones(size(yy));
                plot(xx,yy,'k--','LineWidth',0.5);
                hold on;
            end
            
            ylabel('Extension (nm)');
            xlabel('Time (s)');
            set(gca,'FontSize',15,'LineWidth',1.5);
            set(gca, 'color', 'none');
            set(gca,'TickDir','out');
            set(gca,'TickLength',[0.015, 0.01]);
            set(gca,'XColor','k', 'YColor','k');
            set(gca, 'Layer', 'top');
            ax_all = gca;
            box(ax_all,'off');

            %% Plot initial hat curve and fit
            nexttile(9,[3 1]);
            [turns_p, ext_p] = binning_hat(Turns(pos_hat_range),...
                Extension_notfiltered(pos_hat_range));
            [turns_n, ext_n] = binning_hat(Turns(neg_hat_range),...
                Extension_notfiltered(neg_hat_range));
            plot (turns_p, ext_p, 'LineWidth', 4);
            hold on;
            plot (turns_n, ext_n, 'LineWidth', 4);
            hold on;
            
            dEXT_during = mean(abs(ext_p-ext_n));
            string_1 = ['Mean hysteresis = ', num2str(abs(dEXT_during), '%.0f'), ' nm.'];
            flagg = 0;
            
            if dEXT_during<=50
                ext_np = 0.5*(ext_n+ext_p);
                C0_fit = fitSingle_m_Jin(turns_p, ext_np/1000);
                zz = f_Single(C0_fit,turns_p);
                plot (turns_p,zz*1000,'g','Linewidth',1.5);
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
                plot (turns_p,zz*1000,'c','Linewidth',1.5);
                hold on;
                pre_title = 'Bad';
                legend_script = 'Standard hat curve';
                
            end
                        
            
            ylim ([0 1500]);
            xlim ([0 turns_p(end)]);
            ylabel('Extension [nm]');
            xlabel('Magnet Turns');
            lgd = legend ('Initial hat curve (- to +)', 'Initial hat curve (+ to -)', ...
                legend_script, 'location','southwest');
            lgd.FontSize = 12;
            set(lgd,'color','none');
            
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
            
            %% Calculating where it reaches buckling
            Waiting_range = (Timespot_index(17):Timespot_index(18)).';
            Wait_turns = mean(Turns(Waiting_range));
            n0 = C0_fit.n0;
            nsp = C0_fit.nsm;
            kp = C0_fit.kp*1000;
            N_0 = L0/10.5;
            sigma_s = nsp/N_0;
            Z_s = f_Single(C0_fit,nsp+n0)*1000;
            npp = Z_s/kp+nsp;
            sigma_p = npp/N_0;
            kp_corrected = kp*(1+sigma_p);
            Z_buckling = (1+Wait_turns/N_0)*Z_s/(sigma_s+1);
            Ext_after_jump = mean(Extension(Waiting_range(1:10)));
            Ext_high = min(Z_buckling,Ext_after_jump);
            Ext_low = Ext_high-kp_corrected*(Minimal_moving_distance/10.5);
            
            nexttile(7,[3 2]);
            unit_vec = ones(size(Waiting_range));
            plot (Time(Waiting_range), Ext_high*unit_vec, 'r', 'LineWidth', 1.5);
            hold on;
            plot (Time(Waiting_range), Ext_low*unit_vec, 'g', 'LineWidth', 1.5);
            hold on;
            yy = (0:100:5000).';
            unit_vec_y = ones(size(yy));
            for qq = 1:length(Cutoff_time)
                Cutoff_t = Cutoff_time(qq)+Timespot(17);
                plot (Cutoff_t*unit_vec_y, yy, 'm--', 'LineWidth', 1.5);
                hold on;
            end
            nexttile(3,[2 1]);
            plot (Time, Extension_notfiltered, 'Color', [0.7,0.7,0.7]);
            hold on;
            plot (Time, Extension, 'Color', [0,0,1], 'Linewidth', 4);
            hold on;
            plot (Time(Waiting_range), Ext_high*unit_vec, 'r', 'LineWidth', 1.5);
            hold on;
            plot (Time(Waiting_range), Ext_low*unit_vec, 'g', 'LineWidth', 1.5);
            hold on;
            for ii = 1:length(Timespot)-1
                time_i = Timespot (ii);
                yy = 0:0.01:4;
                yy = yy*1000;
                xx = time_i*ones(size(yy));
                plot(xx,yy,'k--','LineWidth',0.5);
                hold on;
            end
            yy = (0:100:5000).';
            unit_vec_y = ones(size(yy));
            for qq = 1:length(Cutoff_time)
                Cutoff_t = Cutoff_time(qq)+Timespot(17);
                plot (Cutoff_t*unit_vec_y, yy, 'm--', 'LineWidth', 1.5);
                hold on;
            end
            xlim ([Timespot(17)-20 Timespot(17)+60]);
            ylim ([Ext_high-200 Ext_high+200]);
            ylabel('Extension (nm)');
            xlabel('Time (s)');
            set(gca,'FontSize',15,'LineWidth',1.5);
            set(gca,'Color', 'None');
            set(gca,'TickDir','out');
            set(gca,'TickLength',[0.015, 0.01]);
            set(gca,'XColor','k', 'YColor','k');
            set(gca, 'Layer', 'top');
            ax_all = gca;
            box(ax_all,'off');
            
            saveas(gcf,[new_folder,['\',pre_title,'_Box ', num2str(plotting_index(boxes)-1),'.fig']]);
            saveas(gcf,[new_folder,['\',pre_title,'_Box ', num2str(plotting_index(boxes)-1),'.png']]);
            saveas(gcf,[new_folder,['\',pre_title,'_Box ', num2str(plotting_index(boxes)-1),'.svg']]);
            close;
        end
        close all;
        output = [date_rec, DAQ_time_rec, plotting_index.'-1, force_rec, torque_rec ];
       
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
opts.StartPoint = [-0.001 1.3 40/1000 0 7];
opts.Upper = [0 2.5 100/1000 2 12];

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

function new_C0 = Predict_Hatcurve_Single (fitresult, Li, Lf, dp)
dn = (Li-Lf)/10.5+dp;
scaling_factor = Lf/Li;

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

function output = Smooth_data_each_task (Extension, Task, smoothing_factor)
    output = zeros(size(Extension));
    for ii = min(Task):max(Task)
        current_region = Task==ii;
        if sum(current_region)
            output(current_region) = smoothdata(Extension(current_region), 'gaussian', smoothing_factor);
        end
    end
end

function output = Smoothed_trace(Time, Extension, Task, smoothing_time)
    dT = mean(diff(Time));
    smoothing_factor = smoothing_time/dT;
    reasonable_region = (Extension>=-100 & Extension <=2000);
    Time_reasonable = Time(reasonable_region);
    Extension_reasonable = Extension(reasonable_region);
    Extension_interp = interp1(Time_reasonable,Extension_reasonable,Time,'linear', 'extrap');
    output = Smooth_data_each_task (Extension_interp, Task, smoothing_factor);
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