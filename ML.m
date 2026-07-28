%-----------------------------------------------------------------------------%
%           Catalogo de verossimilhança (Maximum Likelihood - ML)             %
%-----------------------------------------------------------------------------%


function[Catalogo_X, qam_sum, Catalogo_Bits, b_indice, b_sim] = ML(p, c, Catalogo_X, Catalogo_Bits, LUT, p1,M,K,Fator_Potencia,N)

 bits_teste = int2bit(c-1, p).';       
   Catalogo_Bits(c, :) = bits_teste;
   b_indice = bits_teste(1:p1);
   b_sim = bits_teste(p1+1:end);
%"sim" de simbolo

%separando por +V e -V. Caso o bit seja 1, a saida e +V e caso seja 0, a
%saida e -V
indice_dec = bit2int(b_indice.', p1);
indices_ativos = LUT(indice_dec+1, :);
%a linha de indices ativos busca na LUT os correspondente
b_sim_matriz = reshape(b_sim, log2(M), K);
simbolos_dec = bit2int(b_sim_matriz, log2(M));
qam_sum = qammod(simbolos_dec, M, 'UnitAveragePower', true)*Fator_Potencia;

vec_test_BER = zeros(N,1);
vec_test_BER(indices_ativos) = qam_sum;
Catalogo_X(:, c) = vec_test_BER; %Salva o bloco de bits ideiais


end