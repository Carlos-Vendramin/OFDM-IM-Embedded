%-------------------------------------------%
%            Lookup table (LUT)             %
%-------------------------------------------%
function [p1, p2, p, combinacoes_LUT, LUT] = LUT(N,K,M)

p1 = floor(log2(nchoosek(N,K))); %indice
p2 = K*log2(M);                  %qam/simbolos
p = p1 + p2;                     %total de bits

combinacoes_LUT = nchoosek(1:N, K);
LUT = combinacoes_LUT(1:2^p1,:);


end