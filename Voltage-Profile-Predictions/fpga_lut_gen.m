load('emf_curve_single.mat');

lut_size=11;

voc_soc_lut_index_interp=floor(interp(1000000*voc_soc_lut_index,2));
voc_soc_lut_index_interp_hex=dec2hex(voc_soc_lut_index_interp',8);

voc_soc_lut_interp=floor(1000000*interp(voc_soc_lut,2));
voc_soc_lut_interp_hex=dec2hex(voc_soc_lut_interp',8);

plot(voc_soc_lut_index_interp,voc_soc_lut_interp)

% load('ecmparams_eis.mat');
% temp_array=[15 25 35];
% % voc_array=[3.2 3.6 4];
% voc_array=[3 3.6 4.2];
% %temp left/right, soc up/down top, left corner smallest
% Ri_array=[ecmparams(1:3,1)'; ecmparams(4:6,1)'; ecmparams(7:9,1)'];
% C1_array=[ecmparams(1:3,2)'; ecmparams(4:6,2)'; ecmparams(7:9,2)'];
% R1_array=[ecmparams(1:3,3)'; ecmparams(4:6,3)'; ecmparams(7:9,3)'];
% C2_array=[ecmparams(1:3,4)'; ecmparams(4:6,4)'; ecmparams(7:9,4)'];
% R2_array=[ecmparams(1:3,5)'; ecmparams(4:6,5)'; ecmparams(7:9,5)'];
% 
% temp_min=min(temp_array);
% temp_max=max(temp_array);
% temp_range=temp_max-temp_min;
% voc_min=min(voc_array);
% voc_max=max(voc_array);
% voc_range=voc_max-voc_min;
% lut_size=8;
% 
% [temp_array_interp,voc_array_interp]=meshgrid(temp_min:temp_range/(lut_size-1):temp_max,voc_min:voc_range/(lut_size-1):voc_max);
% 
% ri_interp=round(1000000*interp2(temp_array,voc_array,Ri_array,temp_array_interp,voc_array_interp));
% r1_interp=round(1000000*interp2(temp_array,voc_array,R1_array,temp_array_interp,voc_array_interp));
% r2_interp=round(1000000*interp2(temp_array,voc_array,R2_array,temp_array_interp,voc_array_interp));
% c1_interp=round(interp2(temp_array,voc_array,C1_array,temp_array_interp,voc_array_interp));
% c2_interp=round(interp2(temp_array,voc_array,C2_array,temp_array_interp,voc_array_interp));
% temp_index_interp=round(1000*(temp_min:temp_range/(lut_size-1):temp_max))';
% voc_index_interp=round(1000000*(voc_min:voc_range/(lut_size-1):voc_max))';
% 
% ri_interp_hex=dec2hex(ri_interp,8);
% r1_interp_hex=dec2hex(r1_interp,8);
% r2_interp_hex=dec2hex(r2_interp,8);
% c1_interp_hex=dec2hex(c1_interp,8);
% c2_interp_hex=dec2hex(c2_interp,8);
% temp_index_interp_hex=dec2hex(temp_index_interp,8)
% voc_index_interp_hex=dec2hex(voc_index_interp,8)