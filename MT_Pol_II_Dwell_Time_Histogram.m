clear;
close all;
clc;

Nuc_pos_chart = [403, 403+197, 403+197+197];
INDEX = 1;
JMAX = 10;
Nuc_pos = Nuc_pos_chart(INDEX);
dyad_pos = Nuc_pos+147/2;
RNAP_protects = 20;

%% Selected files for analysis
[selectedFiles, path] = uigetfile('.txt',...
    'Select One or More Files','Multiselect', 'on');
if isequal(selectedFiles,0)
   disp('User selected Cancel');
end

max_t_plotting = 1000;
max_x_counting = 1200;
dx = 3;
dT = 1;
N_trace = length(selectedFiles);
single_file_flag = isa(selectedFiles,'char');
if single_file_flag
    N_trace = 1;
end
Time_edge = (-dT/2:dT:max_t_plotting+dT/2).';
t_plotting = (Time_edge(1:end-1)+Time_edge(2:end))/2;
x_counting_edge = 0:dx:200;
x_counting_edge = [-x_counting_edge, x_counting_edge];
x_counting_edge = unique(x_counting_edge);
x_counting_edge = sort(x_counting_edge);
x_counting_edge = x_counting_edge+dyad_pos;
x_counting_edge = x_counting_edge.';
x_center = ( x_counting_edge(1:end-1)+x_counting_edge(2:end) )*0.5;
T_counting = zeros(size(x_center));
T_counting_survived = T_counting;
dyad_index = find(x_center>=dyad_pos-RNAP_protects,1);

pp = 0;
qq = 0;

%% Show each trace
for ii = 1:N_trace
    percentage = ii/N_trace;
    disp_text = ['Percentage: ', num2str(100*percentage, '%.0f'), '%;'];
    clc;
    disp(disp_text);
    if single_file_flag
        filename = selectedFiles;
    else
        filename = selectedFiles{1,ii};
    end
    id_i = [path filename];
    trace_i = load(id_i);
    tt_i = trace_i(:,1);
    xx_i = trace_i(:,2);
    Nuc_occupancy = trace_i(1,3);
    Nuc_info = flip(dec2bin(Nuc_occupancy,3));
    
    if Nuc_info(1) =='1' 
        [Nt,~,bint] = histcounts(tt_i, Time_edge);
        x_plotting_i = zeros(size(t_plotting));
        for dti = 1:length(t_plotting)
            x_plotting_i(dti) = mean(xx_i(bint==dti));
        end
        
        if JMAX >0
            
        x_plotting_i = smoothdata(x_plotting_i,'gaussian',180);
        
        end
        x_plotting_i = Transcript_length_no_TFIIS(x_plotting_i);
        [N_i, ~, bin_i] = histcounts(x_plotting_i, x_counting_edge);
        N_i = N_i.';
        T_i = N_i*dT;
        T_counting = T_counting+T_i;
        pp = pp+1;
        passed_dyad = find(bin_i>=dyad_index);
        if ~isempty(passed_dyad)
            T_counting_survived = T_counting_survived+T_i;
            qq=qq+1;
        end
    end
end

t_per_x = T_counting/dx/pp;
t_per_x_survived = T_counting_survived/dx/qq;
x_center_origin = x_center-dyad_pos;


figure (1)
%plot (x_center_origin , t_per_x, 'b','LineWidth', 3);
%hold on;
plot (x_center_origin , t_per_x_survived, 'r','LineWidth', 3);
hold on;
set(gca,'FontSize',13,'LineWidth',1.5);
set(gca,'TickDir','out');
set(gca, 'Color', 'none');
ax_all = gca;
Y_Lim = ax_all.YLim;
box(ax_all,'off');
plot ([0 0], [0 1000], 'k--','LineWidth', 3);
hold on;
xlabel ('Active site position [bp]');
ylabel ('Dwell Time per bp [s/bp]');
xlim ([-100 100]);
ylim ([0 max(Y_Lim)]);
hold off;

output = [x_center_origin, t_per_x_survived];

function output = Transcript_length_no_TFIIS(x)
    NN = length(x);
    output = x;
    for jj = 2:NN
        output(jj) = max(output(1:jj));
    end
end

