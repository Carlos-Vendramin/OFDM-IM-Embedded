%Frequencia e largura de banda

function [Fc, Bw, fs, F_corte, a_filt,b_filt, H_freq, w, H_eff] = freq_bw(N,Fc, Bw, fs, F_corte)


[b_filt,a_filt]=butter(4,F_corte/(fs/2));
[H_freq, w] = freqz(b_filt, a_filt, N, 'whole');
H_eff = H_freq .^ 2; 

end