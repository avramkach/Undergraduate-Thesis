clear
load('lecm_data/emf.mat'); %corresponding VOC-SOC

%Load actual submodule test voltages and currents
%submodule 1
% load('lecm_data/drive_profile_lecm_submodule_1');
% load('lecm_data/exp_vterm_lecm_submodule_1.mat'); %actual voltages
% load('lecm_data/drive_profile_baseline_submodule_1');
% load('lecm_data/exp_vterm_baseline_submodule_1.mat');

%submodule 5
%load('lecm_data/drive_profile_lecm_submodule_5');
%load('lecm_data/exp_vterm_lecm_submodule_5.mat');
load('lecm_data/drive_profile_baseline_submodule_5'); %better fit without charge curve initial
load('lecm_data/exp_vterm_baseline_submodule_5.mat');

%setup ecm parameters
%model from EIS data
%load('ecm_params/ecmparams_eis_ext_v2_submodule_1.mat'); 
%load('ecm_params/ecmparams_eis_ext_c_v1_submodule_1.mat');
%load('ecm_params/ecmparams_eis_ext_v2_submodule_1_HhighR.mat');
%load('ecm_params/ecmparams_eis_ext_v2_submodule_5.mat');
%load('ecm_params/ecmparams_eis_ext_c_v1_submodule_5.mat');
%load('ecm_params/ecmparams_eis_ext_c_v2_submodule_5.mat');
load('ecm_params/ecmparams_eis_ext_c_v2_extrap_submodule_5.mat');
%load('ecm_params/ecmparams_eisanalyser_submodule_5.mat'); %RC-RC

temp_array=[15 25 35];
% voc_array=[3.2 3.6 4];
voc_array=[3 3.6 4.2];
%temp left/right, soc up/down top, left corner smallest
Ri_array=[ecmparams(1:3,1)'; ecmparams(4:6,1)'; ecmparams(7:9,1)'];
C1_array=[ecmparams(1:3,2)'; ecmparams(4:6,2)'; ecmparams(7:9,2)'];
R1_array=[ecmparams(1:3,3)'; ecmparams(4:6,3)'; ecmparams(7:9,3)'];
C2_array=[ecmparams(1:3,4)'; ecmparams(4:6,4)'; ecmparams(7:9,4)'];
R2_array=[ecmparams(1:3,5)'; ecmparams(4:6,5)'; ecmparams(7:9,5)'];

%vq1 = interp1(x,v,xq); %default method is 'linear'
%x: known x values (voltages)
%v: known y vales y = v(x) (SOCs)
%xq: query points (interpolated x point with data)

%make initial SOC the same
%Look up table emf.mat structures: voc_soc_lut (voltage) 2.62-4.20V,
%voc_soc_lut_index (SOC) 0-100%


%find max index of vterm, compare discharge and charge

%Peaks around 2415 index at 1208
[maxv, maxi] = max(exp_vterm); %find max voltage index

%Average current to this index
avgChargeI = mean(drive_profile(1:maxi));

%Find voltage increase due to circuit after time
circuitChargeV = avgChargeI*(ecmparams(1,1)+ecmparams(1,3)+ecmparams(1,5));

%Set initital vterm as the first experimental voltage - circuit
%voltage(t=inf)
initial_vterm = exp_vterm(1)-circuitChargeV;

%initial_vterm = exp_vterm(1)-0.05; %makes it more accurate
%initial_vterm = exp_vterm(1)-0.1; %for submodule 5
%initial_vterm = exp_vterm(1)-0.01;
initial_soc = interp1(voc_soc_lut, voc_soc_lut_index, initial_vterm);

%reset to 0
drive_profile_index = drive_profile_index - 1;
exp_vterm_index = exp_vterm_index - 1;


%setup other sim params
%Model S submodule capacity 3.2Ah per cell * 74 cells in parallel = 236.8Ah
%3.4Ah*74 = 251.6
%Model S full module: 236.8Ah per module * 3.6V (nominal) = 852.48Wh per module
%852.48Wh per module * 6 modules (each 74 cells) = 5114.88Wh

%charging slope, Q_total lower = steeper slope, higher = flatter slope
%Q_total = 350; %350Ah and -0.1 on vterm initial works well for lecm submodule 5
Q_total = 350;
Q_cell = 3.2;
initial_temp=25;
T_est=0.02;
%T_est=0.002;

%setup initial model values
initial_voc=interp1(voc_soc_lut_index, voc_soc_lut, initial_soc);
Ri=interp2(temp_array, voc_array, Ri_array, initial_temp, initial_voc);
C1=interp2(temp_array, voc_array, C1_array, initial_temp, initial_voc);
R1=interp2(temp_array, voc_array, R1_array, initial_temp, initial_voc);
C2=interp2(temp_array, voc_array, C2_array, initial_temp, initial_voc);
R2=interp2(temp_array, voc_array, R2_array, initial_temp, initial_voc);

%set time of simulation to be last time of drive_profile_index
sim_time = drive_profile_index(end);

%plot voc vs soc
figure;
plot(voc_soc_lut, voc_soc_lut_index);
ylabel('SOC (%)')
xlabel('Voltage (V)')
title('SOC VOC Curve')
