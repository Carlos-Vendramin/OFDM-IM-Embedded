%----------------------------------------------------%
%            Declaraçao de vetores RX/TX             %
%----------------------------------------------------%


function [x_total, y_total, y_eq] = rxtx_Vec(N, num_blocos)

%vetores pra receber os dados dos blocos gerados
%x_total e o total de amostras transmitidas
x_total = zeros(N * num_blocos, 1); % Sinal transmitido perfeito
%y_total e o total de amostras recebidas com ruido. Sao vetores incialmente
y_total = zeros(N * num_blocos, 1); % Sinal recebido com ruído
%sem conteudo para receberem posterirmente valores em suas celulas
%-----------------------------------------------------------------------------
y_eq = zeros(N * num_blocos, 1); % Para visualização da constelação


end