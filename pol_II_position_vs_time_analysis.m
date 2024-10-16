clear
traces = dir("../*force_clamp.txt");


for k = 28
filename = traces(k).name;

trace_name = regexp(string(filename),".*(?=_force_clamp.txt)","match");
tracedata = readtable("../"+filename);
load("../"+trace_name + "_trace_data.mat");

time = tracedata.time;
x_raw = tracedata.time;
y_raw = -tracedata.turns;
torque_from_force_raw = tracedata.torque;
dt = x_raw(2)-x_raw(1);



%Trim down the data to the max position
[y_max,I_max] = max(y_raw);

manual_trace_time_end = [...
    15.49,...
    52.87,...
    time(end),...
    12.49,...
    70,...
    62,... %trace 6
    28.9,...
    time(end),...
    85.8,...
    58.47
    ];


I_end = find(time >= time(end), 1, "first");





I_start = 1;
if(I_end <= I_start)
    I_end = I_start+1;
end

I_start
I_end
y = y_raw(I_start:I_end);
x = x_raw(I_start:I_end);
torque_from_force = torque_from_force_raw(I_start:I_end);

% Tweak parameters
dwell_thresh = 2; %Threshold for pauses in the dwell histogram in s/turna
time_cutoff = x(end);
time_cutoff2 = x(end);


[y_smoothed,v] = calculate_velocity_position(x,y,0.5);


% Make figure and define subplots
figure(1);
clf
rows = 3;
cols = 4;
ax_position = subplot(rows,cols,[1 2]);
ax_velocity = subplot(rows,cols,[9 10]);
ax_torque = subplot(rows,cols,[5 6]);
ax_dwell = subplot(rows,cols,[3 4]);
%ax_velocity_binned = subplot(rows,cols,4);
%ax_velocity_hist = subplot(rows,cols,[7,8]);
ax_velocity_regions = subplot(rows,cols,11:12);
linkaxes([ax_position ax_dwell],'y');
linkaxes([ax_position ax_velocity, ax_torque],'x');

% Plot position vs time
subplot(ax_position);
hold on
%plot(x,y, 'LineWidth',2)
plot(x_raw,y_raw, 'LineWidth',2)
%plot(x,y_smoothed, 'LineWidth',2)
ylabel('Number of turns (turn)')
xlabel('Time (s)')
title(trace_data.output_file_name, "Interpreter","none");

% Plot Velocity vs time
subplot(ax_velocity);
plot(x,v, 'LineWidth',2);
ylabel('Velocity (turn/s)');
xlabel('Time (s)');

title('Velicity vs Time');

% Plot torque vs time
subplot(ax_torque);

plot(x_raw, torque_from_force_raw, 'LineWidth',2);
ylim([0,22]);
ylabel("Torque (pN nm)");
xlabel("Time (s)");
title("Torque vs Time");
xline(time_cutoff, '--');


% Dwell time histogram
subplot(ax_dwell);
bin_width = 1;
bins = 1:bin_width:120;


h = histogram(y_smoothed,bins,'Orientation', 'horizontal');
h.BinCounts = h.BinCounts*dt;
xline(dwell_thresh,'--');
xlabel('Dwell Time (s/turn)')
title('Dwell Histogram');


% %  Plot the velocity binned by the postion data.
% ind = round(y_smoothed);
% binned = accumarray(ind,v)./accumarray(ind,ones(size(v)));
% subplot(ax_velocity_binned);
% plot(binned,1:length(binned), 'LineWidth',2);
% xline(0,'--');
% xlabel('Velocity (bp/s)')
% title('Velicity binned by position');
% 
% % Plot the histogram of velocity binned by position
% subplot(ax_velocity_hist);
% histogram(binned);
% title('Binned velocity histogram')
% xlabel('Velocity (bp/s)');
% ylabel('Count')
% average_velocity_binned = mean(binned(binned>10));
% legend(sprintf("Mean (v>10) = %0.1f",average_velocity_binned))
% % subplot(ax_position);
% % hold on
% % position_pauses = y_smoothed()




pauses = get_pauses(h.BinCounts,bins,dwell_thresh)

fig = figure(1);
fig.Position = [10 10 1800 1500];
subplot(ax_position);
hold on

% Create and index array to select all the regions with pausing.
index = false(size(y_smoothed));

for i = 1:length(pauses)
    min_index = find((y_smoothed > pauses(i).start),1);
    max_index = find((y_smoothed <= pauses(i).end + bin_width),1,'last');
    index = index | ((x > x(min_index)) & (x <= x(max_index)));
end
plot(x(index),y_smoothed(index),'.','MarkerSize',10,'Color','Red')


xline(time_cutoff, '--');


tracedata.paused = index;
tracedata.cylinder_turns_pause = y-y(1);
tracedata.cylinder_turns_active = tracedata.cylinder_turns_pause;
tracedata.cylinder_turns_pause(index) = NaN;
tracedata.cylinder_turns_active(~index) = NaN;

writetable(tracedata, trace_data.output_file_name + "pause.txt");

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
trace_data.velocity.processivity = processivity;

% Find the overall velocity from a linear fit to the position vs time.



subplot(ax_position)
hold on



fit = polyfit(x,y,1);
fit_y = polyval(fit,x);
plot(x,fit_y, '--', 'LineWidth',2, 'Color', 'Black');
overall_velocity = fit(1);
overall_torque = mean(torque_from_force);

trace_data.velocity.duration = duration;
trace_data.velocity.overall_velocity = overall_velocity;
trace_data.velocity.overall_torque = overall_torque;


plot_text = [...
    sprintf('Processivity = %0.1f turn', processivity)...
    newline...
    sprintf('Overall Velocity = %0.1f turn/s', overall_velocity)...
    ];
text(10,110,plot_text);



% What is this?
[active_regions, numRegions] = bwlabel(~index);

tstart = x_raw(I_start);

% Add the vertical lines to all the plots.
subplot(ax_position)
yline(y(1));
yline(y(1)+processivity);
xline(tstart)

subplot(ax_torque)
xline(tstart)

subplot(ax_velocity)
xline(tstart)



subplot(ax_torque)
hold on
average_torque = overall_torque;
plot(x,ones(size(x))*average_torque, '--', 'LineWidth',2, 'Color', 'Black')

plot_text = [sprintf('Average torque = %0.1f pNnm',overall_torque)];
text(10,15,plot_text);
xlim([0 120])

subplot(ax_velocity_regions)
hold on

velocity_regions = [];
for i = 1:numRegions
    region = active_regions == i;
    x_cropped = x(region);
    y_cropped = y_smoothed(region);

    fit = polyfit(x_cropped,y_cropped,1);
    velocity_regions(i).distance = polyval(fit,x_cropped(end)) - polyval(fit,x_cropped(1));
    velocity_regions(i).velocity = fit(1);
    velocity_regions(i).torque_from_force = mean(torque_from_force(region));

    
end
hold off
plot([velocity_regions.distance],[velocity_regions.velocity],'.','MarkerSize',30);
average_velocity = sum([velocity_regions.distance].*[velocity_regions.velocity])/sum([velocity_regions.distance]);
average_torque = sum([velocity_regions.distance].*[velocity_regions.torque_from_force])/sum([velocity_regions.distance]);
title('Velocity of active regions')
ylabel('Velocity of region (turn/s)');
xlabel('Distance of active region (turn)')
legend(sprintf("Weighted Mean = %0.1f turn/s" + newline + "Weighted torque = %0.1f pNnm",average_velocity, average_torque),'Location','SouthEast');

trace_data.velocity.pause_free_velocity = average_velocity;
trace_data.velocity.pause_free_torque = average_torque;



set(findall(gcf,'-property','FontSize'),'FontSize',18)
saveas(fig, strcat(trace_data.output_file_name,"_velocity.png") );

save(trace_data.output_file_name+"_trace_data.mat","trace_data");



end


function pauses = get_pauses(dwell_hist,bins,dwell_thresh)
    
    ind = dwell_hist > dwell_thresh;
    
    % Label regions with unique value
    [labeledVector, numRegions] = bwlabel(ind);
    
    pauses = [];
    for i = 1:numRegions
        bins_selected = bins(labeledVector == i);
        pauses(i).start = min(bins_selected);
        pauses(i).end = max(bins_selected);
    end
   
    % figure(2)
    % subplot(2,1,1)
    % plot(labeledVector)

end