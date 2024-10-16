data = readtable("028_240320_0038_Shumingpause.txt");

data.transcribed_turns = -data.turns; 

% Filter cylinder turns to make it more smooth
data.cylinder_turns = movmean(data.turns,100);

% Found the source of the phase delay
% figure(3)
% hold on
% plot(data.time,data.cylinder_turns)
% temp = data.replicated_bp/10.5;
% plot(data.time,temp - temp(1))

%data.cylinder_turns = data.cylinder_turns - data.cylinder_turns(1);

data.transcribed_turns_moving = data.transcribed_turns;
data.transcribed_turns_moving(logical(data.paused)) = NaN;
data.transcribed_turns_not_moving = data.transcribed_turns;
data.transcribed_turns_not_moving(~logical(data.paused)) = NaN;


fig = figure(1);
set(gcf,'Position',[100 100 1280 720])
clf;
set(gcf,'color','w');

%draw_frame(data,4000,5000);

vid = VideoWriter("animation_16x9.avi");
vid.Quality = 100;
frameRate = 30.0;
vid.FrameRate = frameRate;
open(vid);



t_start = 25;
t_end = 40;
duration = t_end - t_start;

dt = data.time(2)-data.time(1);
n_frames = int16(duration*frameRate);
i_start = int16(t_start/dt)

for i_frame = 1:n_frames
    i_data = int16(double(i_frame)/frameRate/dt) + i_start;
    clf;
    draw_frame(data, i_start, i_data);
    writeVideo(vid,getframe(gcf));
end

vid.close();


function draw_frame(data, i_start, i_data)

    color_moving = [1 0 0]; % Red
    color_not_moving = [0 0 0]; %Black

    if ~data.paused(i_data)
        color = color_moving;
        linewidth = 4;
    else
        color = color_not_moving;
        linewidth = 2;
    end

    


    
    % Plot the trace
    hold on
    
    
    % Plot entire trace with lighter color
    percent_lighter = 0.5;
    plot(data.time,data.transcribed_turns_not_moving,'.',"MarkerSize",6,"Color",lighten_plot_color(color_not_moving,percent_lighter));
    plot(data.time,data.transcribed_turns_moving,'.',"MarkerSize",6,"Color",lighten_plot_color(color_moving,percent_lighter));
    
    % Plot the animated region in bolder color.  Also larger font?
    plot(data.time(i_start:i_data),data.transcribed_turns_not_moving(i_start:i_data),'.',"MarkerSize",10,"Color",color_not_moving);
    plot(data.time(i_start:i_data),data.transcribed_turns_moving(i_start:i_data),'.',"MarkerSize",10,"Color",color_moving);

    xline(data.time(i_start),'--',"Color", [0.1 0.1 0.1], "LineWidth",1.5);
    xline(data.time(i_data),'--',"Color", [0.1 0.1 0.1], "LineWidth",1.5);
    plot(data.time(i_data),data.transcribed_turns(i_data),'.',"MarkerSize",30,"Color",color)
    

    % 240401 - Add scale bar
    position = [43,15];
    rectangle_position = [position, 0.2, 100/10.5*1.041]; % [x y w h] 1.041 is the "Jin factor" at 0.1 pN
    rectangle("Position",rectangle_position, "FaceColor", "Black");
    text(position(1) + 0.3, position(2) + 5, " 100 bp", "FontSize", 18)

    xlim([0,47]);
    ylim([0,75]);
    %yticks([80:20:160])
    xlabel("Time (s)");
    ylabel("Pol II Rotation of DNA (turns)");
    set(gca,'linewidth',2)
    set(gca,"FontSize",18)
    box off

    %Plot the cylinder
    inset = axes("Position", [0.11 0.50 0.4 0.4]);
    axis off;
    pos = [-1 -1 2 2];
    rectangle('Position',pos,'Curvature',[1 1], ...
        'LineWidth',2);
    hold on
    p1 = [0 0];
    theta = data.cylinder_turns(i_data)*2*pi;
    p2 = [cos(theta) sin(theta)];
    dp = p2-p1;

    
    quiver(p1(1),p1(2),dp(1),dp(2),0,"LineWidth",linewidth,"Color",color, ...
        "MaxHeadSize",0.8);

    %Draw static arrow over cylinder to indicate overall cylinder rotation
   
    arc_arrow();

    ylim([-1,1.5])
    axis equal

   
end

function plot_color_light = lighten_plot_color(plot_color,percent)
    %Lighten plot color by mixing with white.
    plot_color_light = plot_color + percent*(1-plot_color);
end