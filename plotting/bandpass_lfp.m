function y = bandpass_lfp(x, fs, band)
% BANDPASS_LFP  Zero-phase FIR band-pass, same design as the detectors (fir1, order 6*fs/bandwidth).
order = round(6 * fs / diff(band));
b = fir1(order, band / (fs / 2), 'bandpass');
y = filtfilt(b, 1, double(x));
end
