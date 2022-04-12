clear
%load('drive_profile_ftp.mat');
%load('emf_curve_single.mat');
%load('exp_vterm.mat');
load('lecm_data/drive_profile_lecm_submodule_1');
load('lecm_data/emf.mat'); %corresponding VOC-SOC
load('lecm_data/exp_vterm_lecm_submodule_1.mat'); %actual voltages

%setup ecm parameters
%model from EIS data
%load('ecmparams_eis.mat');
load('ecm_params/ecmparams_eis_ext_v2_submodule_1.mat');
%load('ecmparams_eis_ext_v2_submodule_1_lowR.mat'); %R0 = 1e-4
%load('ecmparams_eis_ext_v2_submodule_1_highR.mat'); %R0 = 4e-4
%load('ecmparams_eis_ext_v2_submodule_1_HhighR.mat'); %R0 = 6e-4
temp_array=[15 25 35];
% voc_array=[3.2 3.6 4];
voc_array=[3 3.6 4.2];
%temp left/right, soc up/down top, left corner smallest
Ri_array=[ecmparams(1:3,1)'; ecmparams(4:6,1)'; ecmparams(7:9,1)'];
C1_array=[ecmparams(1:3,2)'; ecmparams(4:6,2)'; ecmparams(7:9,2)'];
R1_array=[ecmparams(1:3,3)'; ecmparams(4:6,3)'; ecmparams(7:9,3)'];
C2_array=[ecmparams(1:3,4)'; ecmparams(4:6,4)'; ecmparams(7:9,4)'];
R2_array=[ecmparams(1:3,5)'; ecmparams(4:6,5)'; ecmparams(7:9,5)'];
%model from SOC-VOC test data
% load('ecmparams_vsc.mat');
% temp_array=[-25 0 25 50];
% voc_array=[3.074209315 3.442787708 3.510374873 3.585462578 3.637769453 ...
%     3.680661885 3.739318176 3.841407965 3.943322913 4.053565395 4.183012835];


%vq1 = interp1(x,v,xq); %default method is 'linear'
%x: known x values (voltages)
%v: known y vales y = v(x) (SOCs)
%xq: query points (interpolated x point with data)

%make initial SOC the same
%Look up table emf.mat structures: voc_soc_lut (voltage) 2.62-4.20V,
%voc_soc_lut_index (SOC) 0-100%

%initial_vterm = exp_vterm(1)-0.05; %makes it more accurate
initial_vterm = exp_vterm(1);
initial_soc = interp1(voc_soc_lut, voc_soc_lut_index, initial_vterm);

%setup other sim params
%Model S submodule capacity 3.2Ah per cell * 74 cells in parallel = 236.8Ah
%Model S full module: 236.8Ah per module * 3.6V (nominal) = 852.48Wh per module
%852.48Wh per module * 6 modules (each 74 cells) = 5114.88Wh

%charging slope, Q_total lower = steeper slope, higher = flatter slope
Q_total = 236.8;
Q_cell = 3.2;
initial_temp=25;
%T_est=0.02;
T_est=0.02;
%setup initial model values
initial_voc=interp1(voc_soc_lut_index, voc_soc_lut, initial_soc);
Ri=interp2(temp_array, voc_array, Ri_array, initial_temp, initial_voc);
C1=interp2(temp_array, voc_array, C1_array, initial_temp, initial_voc);
R1=interp2(temp_array, voc_array, R1_array, initial_temp, initial_voc);
C2=interp2(temp_array, voc_array, C2_array, initial_temp, initial_voc);
R2=interp2(temp_array, voc_array, R2_array, initial_temp, initial_voc);

% Needs to run the entire drive profile, NO
% If set T_est = 2s use gain 0.5 on clock, NO

%arbitrary topology settings
ns=1;
nl=1;
nc=1;
nm=1;
cell_temp=25*ones(ns,nl,nc);
cell_soc_init=initial_soc*ones(ns,nl,nc);

%set time of simulation to be last time of drive_profile_index
sim_time = drive_profile_index(end);

%% Calculate ECM parameters
for s = 1:ns
    for l = 1:nl
        for c = 1:nc
            cell_voc_init(s,l,c)=interp1(voc_soc_lut_index, voc_soc_lut, cell_soc_init(s,l,c));
            cell_r0_init(s,l,c)=interp2(temp_array, voc_array, Ri_array, cell_temp(s,l,c), cell_voc_init(s,l,c));
            cell_c1_init(s,l,c)=interp2(temp_array, voc_array, C1_array, cell_temp(s,l,c), cell_voc_init(s,l,c));
            cell_r1_init(s,l,c)=interp2(temp_array, voc_array, R1_array, cell_temp(s,l,c), cell_voc_init(s,l,c));
            cell_c2_init(s,l,c)=interp2(temp_array, voc_array, C2_array, cell_temp(s,l,c), cell_voc_init(s,l,c));
            cell_r2_init(s,l,c)=interp2(temp_array, voc_array, R2_array, cell_temp(s,l,c), cell_voc_init(s,l,c));
        end
    end
end