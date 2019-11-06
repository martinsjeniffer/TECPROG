// Este arquivo: spwar.h
// Breve descrição da biblioteca: ...
/////////////////////////////////////

#ifndef _SPWAR_H
#define _SPWAR_H

// Inclusão de outras interfaces
// necessárias para entender esta
// interface.


// Definições das macros que
// a biblioteca bib usa.

#define G 6.673e-11 //constante gravitacional
#define WW 1500 //altura da janela
#define WH 800  //largura da janela
#define PI 3.14159265359    //pi
#define NumIma 4 //numero de imagens por quarto de rotação
#define RB 10   // raio de alcance do projétil antes de explodir
#define BOMBRB 50   //raio do alcance da bomba, ou no caso o projétil que explodirá
#define BOMINTERVAL 2 //duração do tempo final, onde o disparo se tornará uma bomba
#define DISPLIM 10 //tempo limite de sobrevivência do disparo


// Protótipos das funções da
// biblioteca.
/////////////////////////////////////

/**************************************/
/* Aritmética com ponteiros double*   */
/**************************************/


/* Quadrado do módulo do vetor c */
double AbsSqrd(double* c);

/* Multiplica o vetor c por uma constante k */
double* Multk(double* c, double k);

/* Adiciona o conteúdo de c2 em c */
double* Add(double* c,double* c2);

/* Subtrai o conteúdo de c2 em c */
double* Sub(double* c,double* c2);

/* Retorna um vetor cópia de c */
double* Copia(double* c);

//retorna modulo de x
double AbsVal(double x);

//retorna o maximo entre x e y
double max(double x, double y);

//retorna o minimo entre x e y
double min(double x, double y);

//escolhe qual imagem do vetor tem a rotação mais apropriada e retorna sua posição no vetor
int Orientacao(double *a);

//checa se dois objetos em coordenadas c1 e c2 se colidiram ou não, bomb é um inteiro que indica se o objeto é uma bomba
//disptime é o tempo limite para o projétil não explodir. Desde tempo até seu fim este estará com um alcance maior de colisão
int CheckCollision(double* c1, double* c2,int bomb,int disptime);

/*Se a coordenada estiver fora dos limites do mapa, colocar o objeto dentro dos limites*/
double* CheckLimits(double*c,double bounds[]);




// A função xxx recebe ... e
// devolve ... tal que ...
//
int xxxx (int a, char c, int v[]);

// A função yyyy recebe ... e
// devolve ... tal que ...
//
void yyyy (double a, int b);

#endif
