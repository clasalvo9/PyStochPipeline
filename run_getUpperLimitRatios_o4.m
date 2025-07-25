%addpath(genpath('C:\Users\welco\stochmat-radiometer'))
%Universal Parameters

%O4a Time array for the three baselines
startHLa = 1368980142;
endHLa = 1389427879;
%startHLb = 1256663956;
%endHLb = 1269363617;
timesHLa = startHLa : 7200 : endHLa ;
%timesHLb = 1256663956 : 7200 : 1269363617 ;
%NOTE H AND L HARDCODED!!! 
timesHL = timesHLa; %cat(2, timesHLa, timesHLb);
ndaysHL = (endHLa-startHLa)/(60*60*24);

snrs = -8 : 0.5 : 8;

calibration_error_H = 0.0693;
calibration_error_L = 0.0410;
calibration_error_HL = sqrt(calibration_error_H^2 + calibration_error_L^2 + (calibration_error_H*calibration_error_L)^2);

conf = 0.95;

%AntiVela
fprintf('Anti-Vela Jr\n');
ra = 2.3213897160;
dec = -46.298055555;
ul_r = get_upper_limit_ratios(snrs, ra, dec, timesHL, calibration_error_HL, conf);
!move ul_ratios.mat ul_ratios_AntiVela_HL_O4a.mat

%CA
fprintf('Anti-Cas A\n');
ra = 6.1237704239609;
dec = 58.8117;
ul_r = get_upper_limit_ratios(snrs, ra, dec, timesHL, calibration_error_HL, conf);
!move ul_ratios.mat ul_ratios_AntiCasA_HL_O4a.mat
