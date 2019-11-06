//    BIBLIOTECA DO JOGO SPACEWAR
//    spwar.c
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "xwc.h"
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

/**************************************/
/*  Definição das fronteiras do mapa  */
/**************************************/

void DefineBoundaries(Nave* Ns,double bounds[]){
    /*Os limites são definidos de acordo com as posições das naves*/
    /*O limite superior, por exemplo, é definido como o máximo entre 0 e o dobro da maior das coordenadas y das naves */
    /*Analogamente, o mínimo é definido o mínimo entre 0 e o dobro da menor das coordenadas y das naves  */
    /*O mesmo é feito para os limites direito e esquerdo*/
    bounds[0]= max(0,2*max(Ns[0]->c[1],Ns[1]->c[1])); /*UpperBound*/
    bounds[1]= min(0,2*min(Ns[0]->c[1],Ns[1]->c[1])); /*LowerBound*/
    bounds[2]= max(0,2*max(Ns[0]->c[0],Ns[1]->c[0])); /*RightBound*/
    bounds[3]= min(0,2*min(Ns[0]->c[0],Ns[1]->c[0])); /*LeftBound */
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

/*********************************/
/*    Funções Construtoras       */
/*********************************/

void PrepareMask(WINDOW* w1, PIC P, MASK Msk){
  /* recebe
    window w1
    pic auxiliar para aplicar mascara P
    mascara a ser aplicada Msk

    aplica Msk em P
  */

  PutPic(w1, P,  0, 0, WW, WH, 0, 0);
  SetMask(P,Msk);
}


/*cria vetor de struct contendo a picture e a mascara correspondente para cada uma das possiveis orientações*/
Sprite CriaSprite(int norb, WINDOW* w1){
    Sprite S = malloc(sizeof(struct sprite));
    int i;
    char str[50],strint[10],strnave[2];
    S->P=malloc(4*NumIma*(sizeof(PIC)));
    S->Aux=malloc(4*NumIma*(sizeof(PIC)));
    S->Msk=malloc(4*NumIma*(sizeof(MASK)));
    if(norb==0){
        for(i=0;i<4*NumIma;i++){
            S->Msk[i] = NewMask(w1,WW,WH);
            strcpy(str,"images/bullet/tiro");
            sprintf(strint,"%d",i);
            strcat(str,strint);
            strcat(str,".xpm");
            S->P[i] = ReadPic(w1, str, S->Msk[i]);
            S->Aux[i] = NewPic(w1, WW, WH);
            PrepareMask(w1, S->Aux[i], S->Msk[i]);
        }
        return S;
    }
    for(i=0;i<4*NumIma;i++){
        S->Msk[i] = NewMask(w1,WW,WH);
        strcpy(str,"images/nave/nave");
        sprintf(strnave,"%d",norb);
        sprintf(strint,"%d",i);
        strcat(str,strnave);
        strcat(str,"-");
        strcat(str,strint);
        strcat(str,".xpm");
        S->P[i] = ReadPic(w1, str, S->Msk[i]);
        S->Aux[i] = NewPic(w1, WW, WH);
        PrepareMask(w1, S->Aux[i], S->Msk[i]);
    }
    return S;
}

Nave CriaNave(double m, double x, double y, double vx, double vy, char* nome,WINDOW* w1,int numnave){
    Nave N = malloc(sizeof(struct nave));
    double *c,*v,*a;

    c = malloc(sizeof(double)*2); /* vetor posição    */
    v = malloc(sizeof(double)*2); /* vetor velocidade */
    a = malloc(sizeof(double)*2); /* vetor aceleração */
    N->spt= CriaSprite(numnave,w1); /*vetor de sprite*/
    N->nome=nome;
    N->m=m;
    N->c=c; /*coordenadas*/
    N->v=v; /*velocidade*/
    N->a=a; /*aceleração*/
    N->o=0; //começa apontando para cima
    N->c[0]=x;
    N->c[1]=y;

    N->v[0]=vx;
    N->v[1]=vy;

    N->a[0]=0;
    N->a[1]=0;

    return N;
}

/*Criação de apenas 2 naves porém pode ser facilmente generalizada para n naves*/
Nave* CriaNaves(double**naves, char**nomes,WINDOW* w1){
    int i = 0;
    Nave* NaveList = malloc(2*sizeof(Nave));
    for(i = 0; i < 2; i++){
        NaveList[i] = CriaNave(naves[i][0],naves[i][1],naves[i][2],naves[i][3],naves[i][4],nomes[i],w1,i+1);
    }
    return NaveList;
}

Planet CriaPlaneta(double r,double m,double t){
    Planet P = malloc(sizeof(struct planet));
    P->m=m; /* massa */
    P->t=t; /* tempo  de simulação*/
    P->r=r; /* raio  */
    return P;
}

Bullet CriaBullet(double m,double x, double y,double vx,double vy,WINDOW* w1){
    Bullet B = malloc(sizeof(struct bullet));
    double *c,*v,*a;
    int i=0; //contador
    struct dirent *de;  // Pointer for directory entry
    c = malloc(sizeof(double)*2); /* vetor posição    */
    v = malloc(sizeof(double)*2); /* vetor velocidade */
    a = malloc(sizeof(double)*2); /* vetor aceleração */
    B->spt= CriaSprite(0,w1); //sprite das imagens
    //identico a nave
    B->m=m;
    B->c=c;
    B->v=v;
    B->a=a;
    B->o=0;

    B->c[0]=x;
    B->c[1]=y;

    B->v[0]=vx;
    B->v[1]=vy;

    B->a[0]=0;
    B->a[1]=0;

    B->o=0;

    return B;
}


/* Criação dos projéteis fornecidos na entrada padrão */
Bullet* CriaBullets(double**bullets, int nb,WINDOW* w1){
    int i = 0;
    Bullet* BulletList = malloc((nb+2)*sizeof(Bullet));
    for(i = 0; i < nb; i++){
        BulletList[i] = CriaBullet(bullets[i][0], bullets[i][1], bullets[i][2], bullets[i][3], bullets[i][4],w1);
    }
    return BulletList;
}

/**************************************/
/*      Desenho dos Objetos           */
/**************************************/

PIC DrawShip(WINDOW* w1,PIC Paux,PIC Pbkg,PIC Pplnt, PIC P,MASK Msk,int i,int j,double bounds[]){
    /* recebe
    window w1
    pic auxiliar Paux
    pic do background Pbkg
    pic do objeto P
    mask do objeto Msk

    desenha o objeto em Paux na posição (i,j) e retorna Paux
    */
    PutPic(Pbkg, Pplnt,  0, 0,400,400,(-bounds[3])*WW/(bounds[2]-bounds[3])-200,
    (-bounds[1])*WH/(bounds[0]-bounds[1])-200);

    UnSetMask(Paux);
    PutPic(Paux, Pbkg,  0, 0, WW, WH, 0, 0);
    SetMask(Paux,Msk);
    PutPic(Paux, P,  0, 0, 100, 100, i, j);
    return Paux;
}

void DrawWindow(WINDOW* w1, PIC Pf, PIC Pbkg,PIC Pplnt, Nave* Ns,double bounds[],Bullet* Bs,int nb,double t,double tb,int* dispclk){
    /* recebe
    window w1
    pic transitorio para sobrepor todos os objetos Pf
    pic da nave1 P1
    pic da nave2 P2
    pic auxiliar da nave1 Paux1
    pic do background Pbkg
    pic auxiliar da nave2 Paux2
    mask da nave1 Msk1
    mask da nave2 Msk2
    Vetor contendo as duas naves Ns
    fronteiras do mapa bounds
    vetor de pics auxiliares dos projeteis Paux
    vetor de pics dos projeteis Pproj
    vetor de masks dos projeteis MaskProj
    Vetor contendo projeteis Bs
    numero de projeteis nb
    tempo atual t
    tempo limite dos projeteis tb

    desenha todos os objetos em w1

    as posições calculadas são recalculadas para as dimensões da tela
    */

    int i;

    Pf=DrawShip(w1,Ns[0]->spt->Aux[Ns[0]->o],Pbkg,Pplnt,Ns[0]->spt->P[Ns[0]->o],Ns[0]->spt->Msk[Ns[0]->o],(-bounds[3]+Ns[0]->c[0])*WW/(bounds[2]-bounds[3]),
    (-bounds[1]+Ns[0]->c[1])*WH/(bounds[0]-bounds[1]),bounds);

    Pf=DrawShip(w1,Ns[1]->spt->Aux[Ns[1]->o],Pf,Pplnt,Ns[1]->spt->P[Ns[1]->o],Ns[1]->spt->Msk[Ns[1]->o],(-bounds[3]+Ns[1]->c[0])*WW/(bounds[2]-bounds[3]),
    (-bounds[1]+Ns[1]->c[1])*WH/(bounds[0]-bounds[1]),bounds);
    for(i=0;i<nb+2;i++){
        if((t<tb && i<nb) || (i>=nb && dispclk[i-nb])){
            Pf=DrawShip(w1,Bs[i]->spt->Aux[Bs[i]->o],Pf,Pplnt,Bs[i]->spt->P[Bs[i]->o],Bs[i]->spt->Msk[Bs[i]->o],(-bounds[3]+Bs[i]->c[0])*WW/(bounds[2]-bounds[3]),
            (-bounds[1]+Bs[i]->c[1])*WH/(bounds[0]-bounds[1]),bounds);
        }
    }
    PutPic(w1, Pf,  0, 0, WW, WH, 0, 0);
}


/*imprime dados no prompt*/
int ImprimeDados(int bo,double t, double tb,int nb, Nave* Ns,Bullet* Bs){
    int i;
    if(bo && t>tb){
            bo=0;
            printf("Tempo Limite de simulação dos projéteis alcançado!\n\n");
        }
    /* resultados impressos na saída padrão */
    printf("Tempo %.3lf:\n",t);
    for(i=0;i<2;i++){
            printf("%5s:      x: %12.4e | y: %12.4e |",Ns[i]->nome,Ns[i]->c[0],Ns[i]->c[1]);
            printf(" vx: %12.4e | vy: %12.4e\n",Ns[i]->v[0],Ns[i]->v[1]);
    }
    if(bo)
            for(i=0;i<nb;i++){
                printf("projétil%d:  x: %12.4e | y: %12.4e |",i+1,Bs[i]->c[0],Bs[i]->c[1]);
                printf(" vx: %12.4e | vy: %12.4e \n",Bs[i]->v[0],Bs[i]->v[1]);
            }
    printf("\n");
    return bo;
}



/**************************************/
/*      Forças gravitacionais         */
/**************************************/


/* Força gravicional entre Planeta e Projétil */
double* gravityPB(Planet P, Bullet B){
    double* g=Copia(B->c);
    double dist= (double)(-G*P->m*B->m/(pow(AbsSqrd(B->c),3/2)));
    g=Multk(g,dist);
    return g;
}

/* Interação entre o Planeta e Nave */
double* gravityPN(Planet P, Nave N){
    double* g = Copia(N->c);

    /*         -G * m * M         */
    /*      (x**2 + y**2)^3/2    */
    double dist = (double)(-G* P->m * N->m/(pow(AbsSqrd(N->c),3/2)));

    g = Multk(g, dist);
    return g;
}

/* Interação entre as naves */
double* gravityNN(Nave N1,Nave N2){
    double* g = Sub(Copia(N2->c),N1->c),
            dist = (double)(-G*N2->m*N1->m/(pow(AbsSqrd(g),3/2)));

    g = Multk(g,dist);
    return g;
}

/* Interação entre os Projétil e Nave */
double* gravityBN(Bullet B, Nave N){
    double* g = Sub(Copia(N->c),B->c),
            dist = (double)(-G*B->m*N->m/(pow(AbsSqrd(g),3/2)));

    g = Multk(g,dist);
    return g;
}

/* Interação entre Projéteis */
double* gravityBB(Bullet B1,Bullet B2){
    double* g = Sub(Copia(B2->c),B1->c),
            dist = (double)(-G*B2->m*B1->m/(pow(AbsSqrd(g),3/2)));

    g = Multk(g,dist);
    return g;
}





/*********************************/
/*    Liberação de memória       */
/*********************************/
/*funções para desalocar a memória no final do programa*/
void FreeNave(Nave N){
    free(N->c);
    free(N->v);
    free(N->a);
    free(N->nome);
    FreeSprite(N->spt);
    free(N);
}
void FreeBullet(Bullet B){
    free(B->c);
    free(B->v);
    free(B->a);
    FreeSprite(B->spt);
    free(B);
}
void FreeSprite(Sprite S){
    int i;
    for(i=0;i<16;i++){
        free(S->P[i]);
        free(S->Aux[i]);
    }
    free(S);
}
