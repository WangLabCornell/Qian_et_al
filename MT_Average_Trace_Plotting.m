clear;
close all;
clc;

Nuc_pos_chart = [403, 403+197, 403+197+197];
INDEX = 1;

Fit_min = 0;
Fit_max = 950;
Nuc_pos = Nuc_pos_chart(INDEX);
dx_fit = 50;
ita = 0.3;
%% Selected files for analysis
[selectedFiles, path] = uigetfile('.txt',...
    'Select One or More Files','Multiselect', 'on');
if isequal(selectedFiles,0)
   disp('User selected Cancel');
end

max_t_plotting = 1000;
dT = 1;
N_trace = length(selectedFiles);
t_plotting = (0:dT:max_t_plotting).';
N_plotting = zeros(size(t_plotting));
Average_x = N_plotting;
Std_x = N_plotting;
x_plotting = zeros(length(t_plotting), N_trace)-1;
B_record = zeros(length(selectedFiles),3);

Hist_Time = [250 500 750];
Hist_record = zeros(N_trace, length(Hist_Time))-10000;


%% Show each trace
for ii = 1:N_trace
    percentage = ii/N_trace;
    disp_text = ['Percentage: ', num2str(100*percentage, '%.0f'), '%;'];
    clc;
    disp(disp_text);
    id_i = [path selectedFiles{1,ii}];
    trace_i = load(id_i);
    tt_i = trace_i(:,1);
    xx_i = trace_i(:,2);
    %xx_i = smoothdata(xx_i,'gaussian',180);
    start_t = trace_i(2*INDEX,3);
    Nuc_info = trace_i(1,3);
    b = de2bi( Nuc_info,3 );
    B_record(ii,:) = b;
        
    if start_t>=0 && b(1) == 1
        
        consider_index = find(tt_i>=start_t);
        tt_i_consider = tt_i(consider_index);
        xx_i_consider = xx_i(consider_index);
        tt_i_consider = tt_i_consider-tt_i_consider(1);
        max_t_i = max(tt_i_consider);
        max_t_i = min(max_t_i, max_t_plotting);
        t_plotting_i = (0:dT:max_t_i).';
        x_plotting_i = interp1(tt_i_consider,xx_i_consider,t_plotting_i);
        %%
        x_plotting_i = Transcript_length_no_TFIIS(x_plotting_i);
        %%
        N_i = length(t_plotting_i);
        N_plotting(1:N_i) = N_plotting(1:N_i)+1;
        x_plotting(1:N_i, ii) = x_plotting_i;
        
        for ww = 1:length(Hist_Time)
            current_hist_index = find(t_plotting_i==Hist_Time(ww));
            if ~isempty(current_hist_index)
                Hist_record(ii, ww) = x_plotting_i(current_hist_index);
            end
        end
        
    end
end

for jj = 1:length(t_plotting)
    x_jj = x_plotting(jj, :);
    x_jj_pos_index = find(x_jj>=0);
    x_jj_pos = x_jj(x_jj_pos_index);
    Average_x(jj) = mean(x_jj_pos);
    Std_x(jj) = std(x_jj_pos);
end

Average_x = Average_x-18;
Average_x = max(Average_x,0);
error_x = Std_x;
SEM_x = Std_x./sqrt(N_plotting);
lower_x = Average_x-error_x*ita;
upper_x = Average_x+error_x*ita;
lower_x_SEM = Average_x-SEM_x;
upper_x_SEM = Average_x+SEM_x;

Output_Xmas = [t_plotting, Average_x, upper_x, lower_x, upper_x_SEM, lower_x_SEM, N_plotting];


patch_t = [t_plotting; flip(t_plotting)];
patch_x = [lower_x_SEM; flip(upper_x_SEM)];



figure (1)
P1 = patch(patch_t,patch_x,'b');
P1.FaceAlpha = 0.3;
P1.LineStyle = 'none';
hold on;
plot (t_plotting, Average_x, 'r', 'LineWidth', 3);
hold on;
enter_i = Nuc_pos*ones(size(t_plotting));
plot (t_plotting, enter_i, 'k--', 'LineWidth', 3);
hold on;
plot (t_plotting, enter_i+147, 'k--', 'LineWidth', 3);
hold on;
xlim([0 1000]);
ylim([0 1000]);
xlabel ('Time [s]');
ylabel ('Pol II Position [bp]');
set(gca,'FontSize',13,'LineWidth',1.5);
set(gca, 'Color', 'None');
ax_all = gca;
box(ax_all,'off');
hold off;

zero_x = Average_x-Nuc_pos;

function output = Transcript_length_no_TFIIS(x)
    NN = length(x);
    output = x;
    for jj = 2:NN
        output(jj) = max(output(1:jj));
    end
end

