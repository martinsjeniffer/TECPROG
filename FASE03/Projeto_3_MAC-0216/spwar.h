// Este arquivo: spwar.h
// Breve descrição da biblioteca:
/*
Biblioteca com comandos usados no arquivo spacewar.c. As funções da biblioteca são esenciais para o funcionamento do jogo, mas não para a compreensão dele, assim, a biblioteca contém apenas as funções auxiliares assim como definições e etc.

*/
/////////////////////////////////////

#ifndef _SPWAR_H
#define _SPWAR_H

// Inclusão de outras interfaces
// necessárias para entender esta
// interface.

#include "xwc.h"
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
#define AC 10   // aceleração gerada pelos controles manuais na nave
#define INITACB 600 //velocidade inicial do projétil lançado


/*********************************/
/*   Imprementação dos Structs   */
/*********************************/

typedef struct bullet *Bullet; /* struct do projetil */

typedef struct nave *Nave; /* struct da nave */

typedef struct planet *Planet; /* struct do Planeta */

typedef struct sprite *Sprite;    /* struct do sprite, que armazenará a imagem do objeto */

struct bullet{
    double m;
    double *c,*v,*a;
    int o;
    struct sprite* spt;
};

struct nave{
    char* nome;
    double m;
    double *c,*v,*a;
    int o;
    struct sprite* spt;
};

struct planet{
    double r,m,t;
};

struct sprite {
  PIC *P;
  PIC *Aux;
  MASK *Msk;
};

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


/**************************************/
/*  Definição das fronteiras do mapa  */
/**************************************/

void DefineBoundaries(Nave* Ns,double bounds[]);
/*Se a coordenada estiver fora dos limites do mapa, colocar o objeto dentro dos limites*/
double* CheckLimits(double*c,double bounds[]);


/*********************************/
/*    Funções Construtoras       */
/*********************************/

void PrepareMask(WINDOW* w1, PIC P, MASK Msk);
/*cria vetor de struct contendo a picture e a mascara correspondente para cada uma das possiveis orientações*/
Sprite CriaSprite(int norb, WINDOW* w1);
Nave CriaNave(double m, double x, double y, double vx, double vy, char* nome,WINDOW* w1,int numnave);
/*Criação de apenas 2 naves porém pode ser facilmente generalizada para n naves*/
Nave* CriaNaves(double**naves, char**nomes,WINDOW* w1);
Planet CriaPlaneta(double r,double m,double t);
Bullet CriaBullet(double m,double x, double y,double vx,double vy,WINDOW* w1);
/* Criação dos projéteis fornecidos na entrada padrão */
Bullet* CriaBullets(double**bullets, int nb,WINDOW* w1);


/**************************************/
/*      Desenho dos Objetos           */
/**************************************/

PIC DrawShip(WINDOW* w1,PIC Paux,PIC Pbkg,PIC Pplnt, PIC P,MASK Msk,int i,int j,double bounds[]);
void DrawWindow(WINDOW* w1, PIC Pf, PIC Pbkg,PIC Pplnt, Nave* Ns,double bounds[],Bullet* Bs,int nb,double t,double tb,int* dispclk);
/*imprime dados no prompt*/
int ImprimeDados(int bo,double t, double tb,int nb, Nave* Ns,Bullet* Bs);


/**************************************/
/*      Forças gravitacionais         */
/**************************************/

/* Força gravicional entre Planeta e Projétil */
double* gravityPB(Planet P, Bullet B);
/* Interação entre o Planeta e Nave */
double* gravityPN(Planet P, Nave N);
/* Interação entre as naves */
double* gravityNN(Nave N1,Nave N2);
/* Interação entre os Projétil e Nave */
double* gravityBN(Bullet B, Nave N);
/* Interação entre Projéteis */
double* gravityBB(Bullet B1,Bullet B2);


/*********************************/
/*    Liberação de memória       */
/*********************************/

/*funções para desalocar a memória no final do programa*/
void FreeNave(Nave N);

void FreeBullet(Bullet B);

void FreeSprite(Sprite S, int n);


#endif
