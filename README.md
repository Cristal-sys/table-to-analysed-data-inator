# table-to-analysed-data-inator
It takes data from lab essay and analyses it for scientific use. the through radme.txt file in portuguese is given below:

Olá! Seja bem-vind@ ao table-to-analysed-data-inator! :3

Trata-se de um programa que lê seus dados brutos de um .txt, e produz uma tabela .txt com os dados analisados! Para entender como usá-lo, precisamos passar por esses pontos:

I) como colocar os dados em um arquivo txt
II) como escrever a função final no arquivo txt
III)como utilizar o .exe


************************************** I : como colocar os dados brutos em um arquivo .txt? ***************************


cada variável é composta de:

a) um nome, simbolizado por um caractere ascii. prefira letras simples (a, b, c...), evite simbolos estranhos (@, #, $, ...) e jamais utilise caracteres não-ascii (à, â, ã,...)

b) uma incerteza instrumental; um número real, a incerteza do equipamento utilizado para medir.

c) medidas; aquelas que o seu equipamento registra nos experimentos, e que pode ser mais de uma.

dessa maneira, colocaremos, no arquivo, a) o nome da variável, b) sua incerteza instrumental, separada por enter ou espaço, c) suas medidas, separadas por enter ou espaço entre si, e d) um 0, para sinalizar que as medidas de a acabaram. Segue um exemplo:
-------------------------------------------
a
0.05
1.33 1.32 1.45
0
b1
0.001
0.067 0.045 0.039
0*
a*b1
---------------------------------------------

IMPORTANTE:
* separe as casas decimais com pontos ( . ), nunca com vírgula.
** nunca coloque uma medição nula; isso não é físico, e o numero 0 serve somente para dizer que uma tabela acabou.
*** escolha unidades convenientes para escrevê-las na tabela. prefira 1.3 um, ao invés de 0.00013cm. 
****(ESSENCIAL) o último 0, que sinaliza o fim de todas as tabelas, deve ser seguido de um *, para indicar o início da expressão da 
função.


******************************* II : como escrever a função final no arquivo txt? *********************************************


Agora é a hora de escrever uma função, isto é, a expressão final, que usa as outras variáveis e suas incertezas. É essencial seguir alguns passos simples para respeitar a formatação da função, e para que nenhum erro ocorra. Por favor, memorize as seguintes regras:

i) uma expressão é composta por um operando, um operador, e um operando.
ii) só cubra com parênteses uma expressão se ela for um operando de uma outra expressão.
iii) jamais cubra um só número ou uma só variável com parênteses.
vi) podemos usar, até a presente versão, os operadores + - * / ^ (mais, menos, vezes, dividido, elevado a)
v) Sen(a), Cos(a) e Tan(a) são escritos, até a presente versão, como 1Sa, 1Ca, 1Ta, e são consideradas operações, assim como arcsen(a), arccos(a) e arctan(a), que são escritos como 0Sa, 0Ca, 0Ta, respectivamente.

seguem alguns exemplos de expressões corretas e erradas:

a+b                     correto!
a/b                     correto!
(a)+(b)                 ERRADO! 
(a+b)*c                 correto!
a+(b/c)                 correto!
a+b+c                   ERRADO!
a+(b+c)                 correto!
a*b*b                   ERRADO!
a*(b*c)                 correto!
a^(1/2)                 correto! (não trabalhamos com sqrt, é assim que tomamos a raiz de um número aqui!)
1Sa                     correto!
1Sa+1Cb                 ERRADO!
(1Sa)+(1Cb)             correto!
2*1Ta                   ERRADO!
2*(1Ta)                 correto!

IMPORTANTE: 
*o número de euler e pi já estão implementados aqui, com um bom número de casas decimais. basta representá-los por "e" e "p", caso já não existam variáveis com esses nomes.

******************************* II : como utilizar o .exe? *********************************************************

Agora que toda a tabela foi escrita, é simples! basta abrir o arquivo .exe, e digitar o nome do arquivo .txt com os dados brutos. Isso vai gerar um arquivo .txt com os seus dados já tratados, com média e incerteza final de cada variável! A incerteza é calculada considerando a incerteza instrumental que você forneceu, e o desvio-padrão das várias medida que vpcê registrou. segue a tabela gerada pelo arquivo do exemplo (I):


-------------------------------------------------------------------------------------------------
a 1.366666666666666666782314898398453806294 ± 0.06514940095230687765978953504553139453037
b 0.05033333333333333333304421275400386548426 ± 0.008569973421454960642025359962570263405723

a*b

0.06878901111100000000136327552024972931122 ± 0.01216271434964022345783787523720853585019
-------------------------------------------------------------------------------------------------

colocamos várias casas decimais para a sua comodidade, basta arredondar para a precisão que quiser!



Parabéns! Você chegou ao fim desse tutorial! Agora, você está pronto para utilizar o table-to-analysed-data-inator! Se estiver interessado em novas atualizações do programa, tiver dúvidas ou sugestões, contate o desenvolvedor desse programa:


Bons exercícios de laboratório! :3

