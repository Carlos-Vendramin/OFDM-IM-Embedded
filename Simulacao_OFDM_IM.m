
clear all; clc; close all;

N = 4;              % total de subportadoras
K = 2;              % Número de subportadoras ativas
M = 16;              % Ordem da modulação (4-QAM)
num_blocos = 100;  % Numero de amostras a serem transmitidas no loop
%EbNo_vec= 0:1:2;

%----------------------------------
%Passei a utilizar energia/bit (correçao matematica)

EbNo_vec = 0:2:10               
%----------------------------------
%Esse vetor BER vai ser o responsavel por guardar os valores da BER pra
%cada valor de SNR
BER = zeros(1, length(EbNo_vec));




%vetores pra receber os dados dos blocos gerados
%x_total e o total de amostras transmitidas
x_total = zeros(N * num_blocos, 1); % Sinal transmitido perfeito
%y_total e o total de amostras recebidas com ruido. Sao vetores incialmente
%sem conteudo para receberem posterirmente valores em suas celulas
y_total = zeros(N * num_blocos, 1); % Sinal recebido com ruído

%-----------------------------------------------------------------------------
% Quantidade de bits, automatizada para se adaptar a valores variados
p1 = floor(log2(nchoosek(N,K)));             %indice
%combinaçeos possiveis(2^K)
p2 = K*log2(M);             %qam
%cada subportadora vai enviar 4 simbolos, que tambem carrega 2 bits, logo
%2bits de
%portadoras ligadas * 2 bits de simbolos

p = p1 + p2;        %total de bits
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% Look up table
combinacoes_LUT = nchoosek(1:N, K);
LUT = combinacoes_LUT(1:2^p1,:);
%-----------------------------------------------------------------------------


%Catalogo de verossimilhança (Maximum Likelihood - ML)
%logica monte carlo
combinacoes=2^p;

Catalogo_X = zeros(N, combinacoes);
Catalogo_Bits = zeros(combinacoes, p);


for c = 1:combinacoes
%gerando os bits a partir dos indices

%-------------------------------------------
%Conversao decimal->binario
   bits_teste = int2bit(c-1, p).';       
%------------------------------------------%


%-------------------------------------------------------------------------
%esse pega cada valor de C, e substitui toda a coluna pelo sua respectiva
%linha de bits
Catalogo_Bits(c, :) = bits_teste;
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
%Separando os bits em bits de indice e bits de transmissao (genericamente
%falando, sao bits 'uteis')
b_indice = bits_teste(1:p1);
b_sim = bits_teste(p1+1:end);
%sim de simbolo
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
%separando por +V e -V. Caso o bit seja 1, a saida e +V e caso seja 0, a
%saida e -V
indice_dec = bit2int(b_indice.', p1);
indices_ativos = LUT(indice_dec+1, :);
%a linha de indices ativos busca na LUT os correspondente
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
%vetor de indices, nesse caso, 2 bits
%aqui estou pegando o vetor de simbolos b_sim e transformando em uma matriz
%b_sim_matriz = reshape(b_sim, log2(M), K).';
b_sim_matriz = reshape(b_sim, log2(M), K);
simbolos_dec = bit2int(b_sim_matriz, log2(M));
qam_sum = qammod(simbolos_dec, M, 'UnitAveragePower', true);
%tambem passei a utilizar o qammod para otimizar o codigo. Fazer tudo
%'barebone' acabou ficando mais lendo do que o esperado.
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
%Aqui era a separaçao dos bits e a geraçao do QAM, passei a usar o qammod
%for i = 1:K
 %   
  %      bit_atual = b_sim((i-1)*2 +1: i*2);
   %     parte_real = 2 * bit_atual(1) - 1;
    %    parte_imag = 2 * bit_atual(2) - 1;
     %   qam_sum(i) = (parte_real + 1i * parte_imag) / sqrt(2);
%
%
%end
%deixei tudo comentado pra futuras consultas.
%-------------------------------------------------------------------------

vec_test_BER = zeros(N,1);
vec_test_BER(indices_ativos) = qam_sum;
Catalogo_X(:, c) = vec_test_BER; %Salva o bloco de bits ideiais


end
%% Transmissao
disp('--- TRANSMISSOR ---');
disp('Fazendo um loop de transmissao, pra poder gerar uma grande amostra')
disp('Logica monte carlo para a curva de BER');

for snr_indice = 1:length(EbNo_vec)
 
    
 %------------------------------------------------   
 %Passei a utilizar energia/bit
 EbNo_dB = EbNo_vec(snr_indice)                 %
 SNR_dB = EbNo_dB + 10*log10(p/N);              %
    %snr_linear=10^(SNR_dB / 10);               %
%------------------------------------------------


    erros=0;
%a ideia aqui e fazer com que essa minha transmissao e recepçao ocorram em
%loops, pra pode gerar um erro mais verossimil.
for b = 1:num_blocos
bits_tx = randi([0 1], 1, p);

%separando os bits, entre o indice e simbolos
bits_indice  = bits_tx(1:p1);%aqui a gente pega os dois primeiros bits e usa pra indice
%os 2 primeiros sao de indice e o restante de simbolos
bits_simbolo = bits_tx(p1+1:end); %aqui pega os restantes  e usa pra simbolos (informaçao)


%convetendo os bits de indice que anteriormente estavam em binario p/
%Decimal

dec_indice = bit2int(bits_indice.', p1);
indices_ativos = LUT(dec_indice + 1, :); % +1 pois o MATLAB começa no índice 1

%------------------------------------------------------------------------------------------------------------------------
% Mapeamento 4-QAM (Manual) -> passei a utilizar com toolbox, deixei aqui
% pra consultas futuras.
% Para 4-QAM, pegamos de 2 e    m 2 bits. 
% O primeiro bit define o eixo Real (I) e o segundo o eixo Imaginário (Q)
% Regra: bit '0' vira amplitude -1, bit '1' vira amplitude +1
% simbolos_qam = zeros(1, K);
% 
% for i = 1:K
%     % Pega os 2 bits do símbolo atual
%     b_atual = bits_simbolo( (i-1)*2 + 1 : i*2 ); 
% 
% 
%     parte_real = 2 * b_atual(1) - 1; %nesse passo, se o bit e 0, a tensao e -1, se e 1, a tensao e +1
%     parte_imag = 2 * b_atual(2) - 1;
% 
%     % Monta o número complexo e divide por sqrt(2) para normalizar a energia
%     simbolos_qam(i) = (parte_real + 1i * parte_imag) / sqrt(2);
% end
%------------------------------------------------------------------------------------------------------------------------

%------------------------------------------------------------------------------------------------------------------------

%Mapeamento QAM, utilizando toolbox


    b_sim_matriz=reshape(bits_simbolo, log2(M), K);
    simbolos_dec=bit2int(b_sim_matriz, log2(M));
    simbolos_qam=qammod(simbolos_dec, M, 'UnitAveragePower', true);


%------------------------------------------------------------------------------------------------------------------------


%------------------------------------------------------------------------------------------------------------------------

% --- Formação do Bloco OFDM-IM ---
x_bloco = zeros(N, 1);
x_bloco(indices_ativos) = simbolos_qam; 

%------------------------------------------------------------------------------------------------------------------------



%% Canal e IFFT
x_tempo = ifft(x_bloco, N);
%------------------------------------------------------------------------------------------------------------------------

% 
% % --- Adição de Ruído AWGN (Manual) ---
% %SNR_dB = 15;
% % Converte SNR de dB para escala linear
% snr_linear = 10^(SNR_dB / 10); 

%adiçao do ruido
% % Calcula a potência média do sinal transmitido
% potencia_sinal = mean(abs(x_tempo).^2); 
% 
% % Calcula a variância (potência) necessária para o ruído
% potencia_ruido = potencia_sinal / snr_linear;
% 
% % Gera ruído complexo Gaussiano (randn gera distribuição normal)
% ruido_real = sqrt(potencia_ruido / 2) * randn(N, 1);
% ruido_imag = sqrt(potencia_ruido / 2) * randn(N, 1);
% ruido_complexo = ruido_real + 1i * ruido_imag;
% 
% % Soma o ruído ao sinal
% y_tempo = x_tempo + ruido_complexo;


%todo esse bloco foi substituido somente pela funçao awgn


y_tempo = awgn(x_tempo, SNR_dB, 'measured');
%------------------------------------------------------------------------------------------------------------------------


%% 4. RECEPTOR (RX)
y_bloco = fft(y_tempo, N);


%------------------------------------------------------------------------------------------------------------------------
%salvando as iformaçoes geradas num vetor

VetorBeg = (b-1) * N+1;
VetorEnd = b*N;

x_total(VetorBeg:VetorEnd) = x_bloco;
y_total(VetorBeg:VetorEnd) = y_bloco;
%------------------------------------------------------------------------------------------------------------------------



%detector de ML
%distancias = zeros(1, combinacoes);
% Calcula a distância euclidiana ao quadrado entre o y_bloco
% recebido e TODOS os 64 blocos possíveis do Dicionário
% 
% for c=1:combinacoes
% 
%     distancias(c)=sum(abs(y_bloco-Catalogo_X(:,c)).^2);
% 
% 
% end
%-------------------------------------------------------------------------------


distancias = sum(abs(y_bloco - Catalogo_X).^2,1);
%-------------------------------------------------------------------------------

%encontrando o indice de menor distancia

[~, indice_menor]=min(distancias);

bits_rx=Catalogo_Bits(indice_menor, :);


%-----------------------------------------%
%      Iniciando a contagem de erros      %
%-----------------------------------------%


erros_bloco=sum(bits_tx ~= bits_rx);
erros=erros + erros_bloco;

end
%iniciando calculo da BER
total_bits_TX=num_blocos*p;
BER(snr_indice)=erros/total_bits_TX;

disp(['SNR = ', num2str(SNR_dB), ' dB | BER = ', num2str(BER(snr_indice))]);
end



figure;
%pontos recebidos com ruido (vermelho)
plot(real(y_total), imag(y_total), 'r.', 'MarkerSize', 8);
hold on;
%transmitidos e ideiais (circulos pretos)
plot(real(x_total), imag(x_total), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');

%formataçao do grafico
grid on;
xlabel('Fase (I) - Parte Real');
ylabel('Quadratura (Q) - Parte Imaginária');
title(['Constelação OFDM-IM (N=4, K=2, 4-QAM) - SNR = ', num2str(SNR_dB), ' dB']);
legend('Sinal Recebido (com ruído)', 'Sinal Transmitido (Ideal)', 'Location', 'best');

% Fixa os limites dos eixos X e Y para -2 a +2 para visualização simétrica
axis([-2 2 -2 2]);
disp('--- RECEPTOR ---');
disp('Sinal recebido com ruído na Frequência:');
disp(y_bloco);

%% 4. PLOTAR A CURVA DE BER
figure;
% semilogy é usado porque a BER decai exponencialmente (eixo Y logarítmico)
semilogy(EbNo_vec, BER, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('Desempenho de BER do OFDM-IM (N=4, K=2, 4-QAM)');
axis([min(EbNo_vec) max(EbNo_vec) 1e-5 1]); % Fixa os limites do gráfico