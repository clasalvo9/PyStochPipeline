%Universal Parameters


startHL = 800000000;
endHL =  800000000+14*30*60*60*24;
timesHL = startHL : 7200 : endHL ;
ndaysHL = (endHL-startHL)/(60*60*24)

snrs = -8 : 0.5 : 8;

calibration_error_H = 0.0696;
calibration_error_L = 0.0637;
calibration_error_HL = sqrt(calibration_error_H^2 + calibration_error_L^2 + (calibration_error_H*calibration_error_L)^2)

conf = 0.95;

%SN
fprintf('SN 1987a\n');
ra = 5.591119444444444;
dec = -69.26994166666667;
ul_r = get_upper_limit_ratios(snrs, ra, dec, timesHL, calibration_error_HL, conf);
!move ul_ratios.mat ul_ratios_sn1987a_HL_TEST14_mat.mat
