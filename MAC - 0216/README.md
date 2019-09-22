Projeto MAC-0216 Técnicas de Programação 1 - Projeto 1

Alunos:
Vitor Daisuke Tamae                NUSP 10705620
Jeniffer Florinda Martins da Silva NUSP 10377966
Giulia NUSP ??

O projeto é composto por um arquivo C que recebe como argumento na linha de comando o TimeStep da simulação e na entrada padrão os dados listados no enunciado do ep. Os testes fornecidos estão na pasta input.

A unidade temporal utilizada foi o "segundo" para familiarização das unidades com o SI.

O seguinte comando foi executado para compilação e execução:

gcc sw1.c -o sw1
./sw1 "timestep" < inputs/"arquivo".txt

Um exemplo de teste é dado no makefile, onde é executada a leitura de 7.txt com um timestep de 0.01 segundos.

Foi utilizado no makefile o programa Valgrind para auxílio no controle da memória alocada e desalocada.

Os resultados da simulação são impressos na saída padrão, mostrando o tempo, as posições e velocidades das naves e dos projéteis.


