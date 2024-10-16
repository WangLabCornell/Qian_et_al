clear;
close all;
clc;

file_name = 'Example Trace';
file_name = [file_name, '.txt'];
Data = load (file_name);

Li = 5000;
smoothing_length = 30;
N_fit = 80;

Time = Data(:,1);
Extension = Data(:,2);
Task = Data(:,3);
Turns = Data(:,4);

Timespot_index = timespot_finding (Task);
Timespot = Time(Timespot_index);

neg_hat_range_high = (Timespot_index(5):Timespot_index(6))';
pos_hat_range_high = (Timespot_index(3):Timespot_index(4))';
neg_hat_range_low = (Timespot_index(10):Timespot_index(11))';
pos_hat_range_low_1 = (Timespot_index(8):Timespot_index(9))';
pos_hat_range_low_2 = (Timespot_index(12):Timespot_index(13))';
pos_hat_range_low = [pos_hat_range_low_1; pos_hat_range_low_2];

[turns_p_high, ext_p_high] = binning_hat(Turns(pos_hat_range_high),...
    Extension(pos_hat_range_high));
[turns_n_high, ext_n_high] = binning_hat(Turns(neg_hat_range_high),...
    Extension(neg_hat_range_high));
ext_np_high = 0.5*(ext_n_high+ext_p_high);
C0_fit_high = fitSingle_m_Jin(turns_p_high, ext_np_high/1000);

[turns_p_low, ext_p_low] = binning_hat(Turns(pos_hat_range_low),...
    Extension(pos_hat_range_low));
[turns_n_low, ext_n_low] = binning_hat(Turns(neg_hat_range_low),...
    Extension(neg_hat_range_low));
ext_np_low = 0.5*(ext_n_low+ext_p_low);
C0_fit_low = fitSingle_m_Jin(turns_p_low, ext_np_low/1000);

Range_low = (Timespot_index(15):Timespot_index(16)).';
Range_high = (Timespot_index(17):Timespot_index(18)).';
Turns_low = mean(Turns(Range_low));
Turns_high = mean(Turns(Range_high));
Time_high = Time(Range_high);
Time_low = Time(Range_low);
Extension_high = Extension(Range_high);
Extension_low = Extension(Range_low);
Extension_high = smoothdata(Extension_high, "gaussian", smoothing_length);
Extension_low = smoothdata(Extension_low, "gaussian", smoothing_length);


X_pos_high = zeros(size(Range_high));
T_pos_high = X_pos_high;
dL_high = (0:0.5:800).';
Hat_record_high = zeros(length(dL_high),8);
Hat_record_high(:,1) = dL_high;
Ext_th_high = zeros(size(dL_high));

X_pos_low = zeros(size(Range_low));
T_pos_low = X_pos_low;
dL_low = dL_high;
Hat_record_low = zeros(length(dL_low),8);
Hat_record_low(:,1) = dL_low;
Ext_th_low = zeros(size(dL_low));

for ii = 1:length(dL_high)
    L_high = Li-dL_high(ii);
    new_C01_high =  Predict_Hatcurve_Single (C0_fit_high, Li, L_high, 0);
    [~, new_C01_array_high] = Single_fit_Struct2Array(new_C01_high, 0);
    Hat_record_high(ii,2:end) = new_C01_array_high(1:end);
    Ext_th_high(ii) = 1000*f_Single(new_C01_high, Turns_high);
end

for ii = 1:length(dL_low)
    L_low = Li-dL_low(ii);
    new_C01_low =  Predict_Hatcurve_Single (C0_fit_low, Li, L_low, 0);
    [~, new_C01_array_low] = Single_fit_Struct2Array(new_C01_low, 0);
    Hat_record_low(ii,2:end) = new_C01_array_low(1:end);
    Ext_th_low(ii) = 1000*f_Single(new_C01_low, Turns_low);
end


jump_time = 0.5*(Time_low(end)+Time_high(1));

max_Ext_th_high = max(Ext_th_high);
apex_index_high = find(Ext_th_high == max_Ext_th_high);
Ext_th_pos_high = Ext_th_high(apex_index_high:end);
dL_pos_high = dL_high(apex_index_high:end);


data_temp_pos_high = [Ext_th_pos_high, dL_pos_high];
data_temp_pos_high = sortrows(data_temp_pos_high);
Ext_th_pos_high = data_temp_pos_high(:,1);
dL_pos_high = data_temp_pos_high(:,2);

X_int_high = interp1(Ext_th_pos_high,dL_pos_high,Extension_high,'linear');

max_Ext_th_low = max(Ext_th_low);
apex_index_low = find(Ext_th_low == max_Ext_th_low);
Ext_th_pos_low = Ext_th_low(apex_index_low:end);
dL_pos_low = dL_low(apex_index_low:end);
Ext_th_neg_low = [-1000; Ext_th_low(1:apex_index_low)];
dL_neg_low = [-1000; dL_low(1:apex_index_low)];

data_temp_pos_low = [Ext_th_pos_low, dL_pos_low];
data_temp_pos_low = sortrows(data_temp_pos_low);
Ext_th_pos_low = data_temp_pos_low(:,1);
dL_pos_low = data_temp_pos_low(:,2);

data_temp_neg_low = [Ext_th_neg_low, dL_neg_low];
data_temp_neg_low = sortrows(data_temp_neg_low);
Ext_th_neg_low = data_temp_neg_low(:,1);
dL_neg_low = data_temp_neg_low(:,2);

X_int_pos_low = interp1(Ext_th_pos_low,dL_pos_low,Extension_low,'linear');
X_int_neg_low = interp1(Ext_th_neg_low,dL_neg_low,Extension_low,'linear');
dX_low = X_int_pos_low-X_int_neg_low;
exchange_index = (dX_low==min(dX_low));
exchange_time = mean(Time_low(exchange_index));
X_int_low = [X_int_neg_low(Time_low<=exchange_time); X_int_pos_low(Time_low>exchange_time)];

X_all = [X_int_low; X_int_high];
X_all = smoothdata(X_all, "gaussian",40);
Time_all = [Time_low; Time_high]-jump_time;
output = [Time_all, X_all];
output_low = [Time_all(Time_all<0), X_all(Time_all<0)];
output_high = [Time_all(Time_all>=0), X_all(Time_all>=0)];

figure (1)
plot (Time_all(Time_all<0), X_all(Time_all<0), 'b', 'LineWidth', 3);
hold on;
plot (Time_all(Time_all>=0), X_all(Time_all>=0), 'r', 'LineWidth', 3);
hold on;
xlim ([-50 100]);
ylim ([0 600]);
hold on;

figure (2)
plot (Time_low-jump_time, Extension_low, 'b', 'LineWidth', 3);
hold on;
plot (Time_high-jump_time, Extension_high, 'r', 'LineWidth', 3);
hold on;
xlim([-50 20])
hold on;

output_ext_high = [Time_high-jump_time, Extension_high];
output_ext_low = [Time_low-jump_time, Extension_low];

function index = timespot_finding (task)
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
    sub_index = find(abs(Turns_sort-ii)<=0.1);
    output_turns = [output_turns; ii];
    mean_ext = mean(Extension_sort(sub_index));
    output_ext = [output_ext; mean_ext];
end

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