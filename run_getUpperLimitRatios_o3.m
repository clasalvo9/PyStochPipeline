%addpath(genpath('C:\Users\welco\stochmat-radiometer'))
%Universal Parameters


startHLa = 1238166018;
endHLa = 1253977218;
startHLb = 1256663956;
endHLb = 1269363617;
timesHLa = 1238166018 : 7200 : 1253977218 ;
timesHLb = 1256663956 : 7200 : 1269363617 ;
timesHL = cat(2, timesHLa, timesHLb);
ndaysHL = (endHLa-startHLa + endHLb - startHLb)/(60*60*24);

snrs = -8 : 0.5 : 8;

calibration_error_H = 0.0696;
calibration_error_L = 0.0637;
calibration_error_HL = sqrt(calibration_error_H^2 + calibration_error_L^2 + (calibration_error_H*calibration_error_L)^2);

conf = 0.95;

%SN
%fprintf('SN 1987a\n');
%ra = 5.591119444444444;
%dec = -69.26994166666667;
%ul_r = get_upper_limit_ratios(snrs, ra, dec, timesHL, calibration_error_HL, conf);
%!move ul_ratios.mat ul_ratios_sn1987a_HL_O3_mat.mat

%AntiVela
fprintf('Anti-Vela Jr\n');
ra = 2.3213897160;
dec = -46.298055555;
ul_r = get_upper_limit_ratios(snrs, ra, dec, timesHL, calibration_error_HL, conf);
!move ul_ratios.mat ul_ratios_AntiVela_HL_O3.mat

%CA
fprintf('Anti-Cas A\n');
ra = 6.1237704239609;
dec = 58.8117;
ul_r = get_upper_limit_ratios(snrs, ra, dec, timesHL, calibration_error_HL, conf);
!move ul_ratios.mat ul_ratios_AntiCasA_HL_O3.mat
