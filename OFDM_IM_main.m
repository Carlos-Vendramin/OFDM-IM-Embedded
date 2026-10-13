clear all; clc; close all;
%% Configuraçao do transmissor OFDM-IM
N = 64;              % total de subportadoras
G = 16;               
n = N / G;           
k = 3;               
M = 16;             % Ordem da modulação (M-QAM)
num_blocos = 5000;  % Numero de amostras a serem transmitidas no loop
nCicP = 64;          %tamanho do prefixo ciclico

EbNo_vec = 0:1:20;
BER = zeros(1, length(EbNo_vec));

%----------------------------------%
%       Correçao da potencia       %
%----------------------------------%
Pot_N = n/k;
Fator_Potencia = sqrt(Pot_N);

%-------------------------------------------%
%       Frequencia e largura de banda       %
%-------------------------------------------%
[Fc, Bw, fs, F_Corte, a_filt,b_filt, H_freq, w, H_eff] = freq_bw(N, 2.4e9, 20e6, 20e6, 9.9e6);

%----------------------------------

%-------------------------------------------%
%              Modelagem CFO                %
%-------------------------------------------%
%Carrier frequency offset
%[cfo_hz] = cfo(Crystal_ppm, Fc)
[cfo_hz] = cfo(10e-6, Fc);
%----------------------------------

%----------------------------------------------------%
%            Declaraçao de vetores RX/TX             %
%----------------------------------------------------%
[x_total, y_total, y_eq] = rxtx_Vec(N, num_blocos);

%----------------------------------

%-------------------------------------------%
%            Lookup table (LUT)             %
%-------------------------------------------%
[p1, p2, p_grupo, combinacoes_LUT, LUT_var] = LUT(n, k, M);
p_total = G * p_grupo; 


%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------%
%           Catalogo de verossimilhança (Maximum Likelihood - ML)             %
%-----------------------------------------------------------------------------%
combinacoes = 2^p_grupo;
Catalogo_X = zeros(n, combinacoes);
Catalogo_Bits = zeros(combinacoes, p_grupo);
for c = 1:combinacoes
    [Catalogo_X, qam_sum, Catalogo_Bits, b_indice, b_sim] = ML(p_grupo, c, Catalogo_X, Catalogo_Bits, LUT_var, p1, M, k, Fator_Potencia, n);
end


%-----------------------------------------------------------------------------

%-------------------------------------------%
%            Transmissao                    %
%-------------------------------------------%
disp('--- TRANSMISSOR ---');
disp('Fazendo um loop de transmissao, pra poder gerar uma grande amostra')
disp('Logica monte carlo para a curva de BER');

max_tx_hardware = 4.0;
max_rx_hardware = 4.0;
bits_adc = 12;
bits_dac = 8;
Delta_tx = 2 * max_tx_hardware / (2^bits_dac);
Delta_rx = 2 * max_rx_hardware / (2^bits_adc);
 
%painel de controle
Habilita_CFO         = true;   
Compensa_CFO         = true;
Habilita_Quantizacao = true;
Habilita_Filtro      = true;
Habilita_Rayleigh    = true; %Canal Rayleigh seletivo na frequencia


if ~Habilita_Filtro
    H_eff = ones(N, 1);
end

for snr_indice = 1:length(EbNo_vec)
    
    EbNo_dB = EbNo_vec(snr_indice);
    SNR_dB = EbNo_dB + 10*log10(p_total/(N+nCicP));
    
    erros = 0; 
    tempo = 0;
    txstate = zeros(max(length(a_filt)-1, length(b_filt)-1), 1);
    rxstate = zeros(max(length(a_filt)-1, length(b_filt)-1), 1);
    
%-----------------------------------------------------------%
%               variaveis de sensoriamento                  %
%-----------------------------------------------------------%
%Resoluçao dos bins de sensoriamento
% \Delta R = \frac{c}{2B}, Sendo B a largura de banda do sinal. No caso estudado, Delta R = 7,5
    distancia=1;  %distancia do alvo em bins. [1 bin = 1 Delta R] => 1 bin = 7,5 metros
    velocidade=2; %velocidade do alvo em bins. [1 bin = 1 Delta R] => 1 bin = 7,5 metros

    Radar_tx_Matriz = zeros(N, num_blocos);
    Radar_rx_Matriz = zeros(N, num_blocos);
    tempo_radar = 0;
%----------------------------------------------%

    

    for b = 1:num_blocos
        if Habilita_Rayleigh
   
            n_taps = 10; %numero de taps
            
            h_T = (randn(n_taps, 1) + 1i * randn(n_taps, 1)) / sqrt(2 * n_taps);
            
            % Resposta em frequência do canal
            H_channel = fft(h_T, N);
        else
            % Se Rayleigh = false, o impulso do canal é ideal e plano
            h_T = 1; 
            H_channel = ones(N, 1);
        end
        H_eff_total = H_eff .* H_channel;
                
        [x_bloco, y_bloco, bits_tx, txstate, rxstate, tempo, x_Cic_Final] = TxRx_Cascata(fs, Fator_Potencia, p_grupo, p_total, SNR_dB, LUT_var, N, n, G, nCicP, p1, M, k, a_filt, b_filt, txstate, rxstate, tempo, cfo_hz, Delta_tx, Delta_rx, Habilita_CFO, Compensa_CFO, Habilita_Quantizacao, Habilita_Filtro, h_T);
        
        [bits_rx, erros, x_total, y_total, y_eq] = eqlt(b, N, n, G, p_grupo, x_bloco, y_bloco, H_eff_total, Catalogo_X, Catalogo_Bits, bits_tx, erros, x_total, y_total, y_eq);
        
    %-----------------------------------------%
    %              Sensoriamento              %
    %-----------------------------------------%
    
        % fase_range = exp(-1i * 2 * pi * (0:N-1).' * distancia / N);

        % Atraso circular funciona pois a distância (1 bin) está dentro do CP de 64 amostras
        alvo_tempo_atrasado = circshift(x_Cic_Final, distancia);
        fase_doppler = exp(1i * 2 * pi * velocidade * b / num_blocos);
        % alvo_refletido_freq = (x_bloco .* fase_range) * fase_doppler;
        % 
        alvo_tempo = alvo_tempo_atrasado.*fase_doppler;
        
        Vetor_tempo_radar = tempo_radar + (0:length(x_Cic_Final)-1).'/fs;
        tempo_radar = Vetor_tempo_radar(end) + (1/fs);
        
        if Habilita_CFO && Compensa_CFO
            Cfo_radar = exp(-1i * 2 * pi * cfo_hz * Vetor_tempo_radar);
            alvo_tempo = alvo_tempo .* Cfo_radar;
           
        end

        if Habilita_Filtro
            alvo_tempo = filter(b_filt, a_filt, alvo_tempo);
        end

        if Habilita_Quantizacao
            rx_re_radar = round(real(alvo_tempo)/Delta_rx)*Delta_rx;
            rx_im_radar = round(imag(alvo_tempo)/Delta_rx)*Delta_rx;
            alvo_tempo = rx_re_radar + (1i*rx_im_radar);
        end
        %Remoçao do prefixo ciclico no receptor do radar. Antes da FFT
        alvo_tempo = awgn(alvo_tempo, SNR_dB, 'measured');
        alvo_tempo_sem_cp = alvo_tempo(nCicP+1:end);

        y_radar_bloco = fft(alvo_tempo_sem_cp, N);
        
        Radar_tx_Matriz(:, b) = x_bloco;
        Radar_rx_Matriz(:, b) = y_radar_bloco;
        
    %-----------------------------------------%

    end
    filtro_casado = zeros(N, num_blocos);
    mascara_ativas = abs(Radar_tx_Matriz) > 0;
    filtro_casado(mascara_ativas) = Radar_rx_Matriz(mascara_ativas) ./ Radar_tx_Matriz(mascara_ativas);
    
    
    Range = ifft(filtro_casado, N, 1);
    Doppler = fftshift(fft(Range, num_blocos, 2), 2);


    %-----------------------------------------%
    %      Iniciando a calculo da BER         %
    %-----------------------------------------%
    total_bits_TX = num_blocos * p_total;
    BER(snr_indice) = erros / total_bits_TX;
    
    disp(['SNR = ', num2str(SNR_dB), ' dB | BER = ', num2str(BER(snr_indice))]);
end


%-----------------------------------------%
%               Plotagem                  %
%-----------------------------------------%

timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));

fig_constelacao=figure;
%pontos recebidos com ruido (vermelho)
plot(real(y_eq), imag(y_eq), 'r.', 'MarkerSize', 8);
hold on;
%transmitidos e ideiais (circulos pretos)
plot(real(x_total), imag(x_total), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');

%formataçao do grafico
grid on;
xlabel('Fase (I) - Parte Real');
ylabel('Quadratura (Q) - Parte Imaginária');
%PS: o K e o numero de subportadoras ativas, esqueci disso.
title(['Constelação OFDM-IM (N=', num2str(N), ', K=', num2str(G*k), ', ', num2str(M), '-QAM) - SNR = ', num2str(SNR_dB), ' dB']);
legend('Sinal Recebido (com ruído)', 'Sinal Transmitido (Ideal)', 'Location', 'best');

% Fixa os limites dos eixos X e Y para -2 a +2 para visualização simétrica
axis([-4 4 -4 4]);
disp('--- RECEPTOR ---');
disp('Sinal recebido com ruído na Frequência:');
disp(y_bloco);
png_constelacao = ['constelacao_ofdm_im_HardwareLimitado_', timestamp, '.png'];
exportgraphics(fig_constelacao, png_constelacao, 'Resolution', 300);


% 4. Curva BER
fig_BER=figure;
% semilogy é usado porque a BER decai exponencialmente (eixo Y logarítmico)
semilogy(EbNo_vec, BER, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
grid on;
xlabel('E_B/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title(['Desempenho de BER do OFDM-IM (N=', num2str(N), ', K=', num2str(G*k), ', ', num2str(M), '-QAM)'])
axis([min(EbNo_vec) max(EbNo_vec) 1e-5 1]); % Fixa os limites do gráfico
png_BER = ['BER_ofdm_im_HardwareLimitado_', timestamp, '.png'];
exportgraphics(fig_BER, png_BER, 'Resolution', 300);

%-----------------------------------------------------------------------------

% 5. Mapa range doppler
fig_RangeDoppler=figure;
doppler_axis = linspace(-num_blocos/2, num_blocos/2 - 1, num_blocos);
range_axis = 0:(N-1);

mesh(doppler_axis, range_axis, 10*log10(abs(Doppler).^2));
view(2);
xlabel('Frequência Doppler (Velocidade)');
ylabel('Range (Distância)');
title(['Mapa Range-Doppler OFDM-IM | SNR = ', num2str(SNR_dB), ' dB']);
colorbar;

axis([-30 40 0 30]); 

png_RD = ['RangeDoppler_ofdm_im_HardwareLimitado_', timestamp, '.png'];
exportgraphics(fig_RangeDoppler, png_RD, 'Resolution', 300);