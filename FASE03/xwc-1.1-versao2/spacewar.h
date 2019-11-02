//BIBLIOTECA SPACEWAR (HEADER)
// Este arquivo: spacewar.h
// Breve descrição da biblioteca: ...
/////////////////////////////////////

#ifndef _SPACEWAR_H
#define _SPACEWAR_H

// Inclusão de outras interfaces
// necessárias para entender esta
// interface.

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/select.h>
#include <termios.h>

// Definições das macros que
// a biblioteca bib usa.

#define G 6.673e-11

#define RaioN 2 //raio da nave
#define RaioB 1 //raio das bullets

#define NumIma 4 //numero de imagens para cada lado rotacao
#define NumBul 5 //numero máx de bullets ativos para cada nave
#define TBul 1 //tempo de vida das bullets



struct termios orig_termios;

typedef struct bullet *Bullet; /* struct do projetil */

typedef struct nave *Nave; /* struct da nave */

typedef struct planet *Planet; /* struct do Planeta */

struct nave{
    char* nome;
    double m;
    double *c,*v,*a;
    int o;
    struct bullet *lb; //lista de bullets
};


struct bullet{
    double m;
    double *c,*v,*a;
    int o;
};

struct planet{
    double r,m,t;
};

// Protótipos das funções da
// biblioteca.
/////////////////////////////////////

//configurações para nao esperar input
void reset_terminal_mode();
void set_conio_terminal_mode();
int kbhit();
int getch();

void freeNave(Nave N);

void FreeBullet(Bullet B);

Nave CriaNave(double m, double x, double y, double vx, double vy, char* nome);

/*Criação de apenas 2 naves porém pode ser facilmente generalizada para n naves*/
Nave* CriaNaves(double**naves, char**nomes);

Planet CriaPlaneta(double r,double m,double t);

Bullet CriaBullet(double m,double x, double y,double vx,double vy);

/* Criação dos projéteis fornecidos na entrada padrão */
Bullet* CriaBullets(double**bullets, int nb);

void DefineBoundaries(Nave* Ns,double bounds[]);

int ImprimeDados(int bo,double t, double tb,int nb, Nave* Ns,Bullet* Bs);

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

double AbsVal(double x);
/*Se a coordenada estiver fora dos limites do mapa, colocar o objeto dentro dos limites*/
double* CheckLimits(double*c,double bounds[]);

double max(double x, double y);

double min(double x, double y);

//0 aponta para cima, 4 para direita, 8 para baixo, 12 para esquerda etc
//só quem controla é o jogador
int Orientacao(int o, int n);

double Dist(double *d1, double *d2);

int Colisao(double *d1, double *d2, double r1, double r2);

double* Zera(double* f);

//checa se c pertence a s e devolve i+1; 0 se não pertence
int Valid(int c, char s[10]);

#endif
