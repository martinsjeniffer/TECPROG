//    BIBLIOTECA DO JOGO SPACEWAR
//    spwar.c
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "spwar.h"

/**************************************/
/* Aritmética com ponteiros double*   */
/**************************************/


/* Quadrado do módulo do vetor c */
double AbsSqrd(double* c){
    return (c[0])*(c[0])+(c[1])*(c[1]);
}

/* Multiplica o vetor c por uma constante k */
double* Multk(double* c, double k){
    double newx = c[0]*k, newy=c[1]*k;
    c[0]=newx;
    c[1]=newy;
    return c;
}

/* Adiciona o conteúdo de c2 em c */
double* Add(double* c,double* c2){
    c[0]=c[0]+c2[0];
    c[1]=c[1]+c2[1];
    return c;
}

/* Subtrai o conteúdo de c2 em c */
double* Sub(double* c,double* c2){
    (c[0])-=(c2[0]);
    (c[1])-=(c2[1]);
    return c;
}

/* Retorna um vetor cópia de c */
double* Copia(double* c){
    double* c2=malloc(sizeof(double)*2);
    c2[0]=c[0];
    c2[1]=c[1];
    return c2;
}

//retorna modulo de x
double AbsVal(double x){
    if(x>=0)return x;
    return -x;
}

//retorna o maximo entre x e y
double max(double x, double y){
    if(x>y)return x;
    return y;
}

//retorna o minimo entre x e y
double min(double x, double y){
    if(x<y)return x;
    return y;
}

//escolhe qual imagem do vetor tem a rotação mais apropriada e retorna sua posição no vetor
int Orientacao(double *a)
{
    /*
    Vetor de imagens com a seguinte organização (de 0 a 2n angulos a direita e 2n a 4n angulos a esquerda):
    0... |n  ...  |2n... |3n  .....|4n
    cima |direita |baixo |esquerda |cima
    Devolve uma int orientacao
    */
    double segm;
    double alpha; //angulo entre vetor a e (0,1)
    int r;
    if(a[0]==0 && a[1]<=0)alpha=0;
    else if(a[0]==0 && a[1]>0)alpha=PI;
    else alpha = atan((double)(-a[1]/a[0]));
    segm = 2*PI/(4*NumIma);
    r=a[0]>0 ? -(alpha-(double)(PI/2))/segm : -(alpha-(double)(PI/2))/segm+8;
    return r;
}

//checa se dois objetos em coordenadas c1 e c2 se colidiram ou não, bomb é um inteiro que indica se o objeto é uma bomba
//disptime é o tempo limite para o projétil não explodir. Desde tempo até seu fim este estará com um alcance maior de colisão
int CheckCollision(double* c1, double* c2,int bomb,int disptime){
    double* c=Copia(c1);
    Sub(c,c2);
    if(AbsSqrd(c)<RB+(bomb==1)*(disptime>DISPLIM-BOMINTERVAL)*BOMBRB){
        free(c);
        return 1;
    }
    return 0;
}

/*Se a coordenada estiver fora dos limites do mapa, colocar o objeto dentro dos limites*/
double* CheckLimits(double*c,double bounds[]){
    while(c[1]>bounds[0])
        (c[1])-=(bounds[0]-bounds[1]);
    while(c[1]<bounds[1])
        (c[1])+=(bounds[0]-bounds[1]);
    while(c[0]>bounds[2])
        (c[0])-=(bounds[2]-bounds[3]);
    while(c[0]<bounds[3])
        (c[0])+=(bounds[2]-bounds[3]);
    return c;
}
