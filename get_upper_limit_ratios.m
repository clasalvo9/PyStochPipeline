function [ul_ratio] = get_upper_limit_ratios(snrs, ra, dec, times, calibration_error, conf)
% get_upper_limit_ratios(snrs, ra, dec, times, calibration_error, conf)

    ul_ratio = zeros(length(snrs),1);
    ul_ratio_one_sigma = zeros(length(snrs),1);
    for ii = 1:length(snrs)
        snr = snrs(ii);
        [palpha, palpha_circ, h0] = skyfactor(snr, 1, ra, dec, times, calibration_error);
        ul_ratio(ii) = get_ul_from_posterior(palpha, h0, conf) ...
            / get_ul_from_posterior(palpha_circ, h0,conf);
        ul_ratio_one_sigma(ii) = get_ul_from_posterior(palpha, h0, 0.68) ...
            / get_ul_from_posterior(palpha_circ, h0, 0.68);
        if mod(ii, 10) == 0
            fprintf('done with snr %4.2f\n', snr);
        end
    end
save('ul_ratios','snrs','ul_ratio','ra','dec','times','calibration_error','conf', 'ul_ratio_one_sigma')
end


function ul = get_ul_from_posterior(p, h0, conf)
    % get cdf, normalize to integrate to 1
    c = cumsum(p);
    c = c / c(end);
    [c,I] = unique(c);
    h0 = h0(I);
    ul = interp1(c, h0, conf);
end
    
