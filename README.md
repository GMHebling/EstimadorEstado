# EstimadorEstado


## Branches ativos:
   1. master: Estimador de Estado Trifásico WLS com rotinas de esparsidade e Multifrontal QR
   2. Referencia: Branch com códigos complementares ao master, para servir como base para novas implementações e consultas
   3. AnaliseResidual:
   4. BayesianFusion:
   5. pythonMC:
   6. lagrangeIgualdades: Estimador de Estado Trifásico utilizando restrições de igualdade, rotinas de esparsidade e Multifrontal QR
   7. SparseHatchel: Estimador de Estado Trifásico utilizando a matriz aumentada de Hatchel, rotinas de esparsidade e Multifrontal QR 
   8. BranchCurrent: Estimador baseado no Branch Current em coordenadas retangulares

## Guia de Implementação
  * Faça o pull do master, este é o código mais recente e com maior número de melhorias de desempenho. O branch possui também os sistemas de teste já validados.
  * Crie um novo branch, específico para a funcionalidade que será implementada
  * Se necessário, o branch Referencia pode ser consultado para buscar exemplos. 
  
 ### Guia de instalação do SuiteSparse
 Para utilizar o código no branch master, é necessário fazer a instalação do pacote SuiteSparse. Até o momento, só foi possível fazer a instalação em máquinas Unix/Linux. O código abaixo não é garantia que a instalação será concluída 100%. É necessário verificar se a compilação do código do estimador é feita corretamente.
 
 ```
sudo apt install g++ wget m4 cmake libopenblas-dev liblapack-dev libcholmod3

wget https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/v5.6.0.zip

unzip v5.6.0.zip

cd SuiteSparse

make

sudo make install INSTALL=/usr/local

LD_LIBRARY_PATH=“/lib:/usr/lib:/usr/local/lib:/usr/local/lib/lib:/usr/local/lib/include"

 ```
