#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <dirent.h>
#include "xwc.h"

#define G 6.673e-11
#define WW 1500
#define WH 800
#define PI 3.14159265359

#define NumIma 4 //numero de imagens para cada lado rotacao


typedef struct bullet *Bullet; /* struct do projetil */

typedef struct nave *Nave; /* struct da nave */

typedef struct planet *Planet; /* struct do Planeta */

typedef struct sprite *Sprite;    /* struct do sprite */

struct sprite {
  PIC *P;
  PIC *Aux;
  MASK *Msk;
};

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

/*funções para desalocar a memória no final do programa*/
void freeNave(Nave N){
    free(N->c);
    free(N->v);
    free(N->a);
    free(N->nome);
    free(N);
}

void FreeBullet(Bullet B){
    free(B->c);
    free(B->v);
    free(B->a);
    free(B);
}

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
    N->spt= CriaSprite(numnave,w1);
    N->nome=nome;
    N->m=m;
    N->c=c;
    N->v=v;
    N->a=a;
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
    Bullet B = malloc(sizeof(struct nave));
    double *c,*v,*a;
    int i=0; //contador
    struct dirent *de;  // Pointer for directory entry
    c = malloc(sizeof(double)*2); /* vetor posição    */
    v = malloc(sizeof(double)*2); /* vetor velocidade */
    a = malloc(sizeof(double)*2); /* vetor aceleração */
    B->spt= CriaSprite(0,w1);
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

/*void DestoiPics (PIC P[4*NumIma+1])
{
  int i;
  for (i=0; i<=4*NumIma; i++)
    FreePic(i)
}*/
/* Criação dos projéteis fornecidos na entrada padrão */
Bullet* CriaBullets(double**bullets, int nb,WINDOW* w1){
    int i = 0;
    Bullet* BulletList = malloc(nb*sizeof(Bullet));
    for(i = 0; i < nb; i++){
        BulletList[i] = CriaBullet(bullets[i][0], bullets[i][1], bullets[i][2], bullets[i][3], bullets[i][4],w1);
    }

    return BulletList;
}

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

void Update(Bullet* Bs, Nave* Ns, Planet P, int nb, double dt, double t, double tb){
    /* recebe
    vetor com projeteis Bs
    vetor com naves Ns
    planeta P
    numero de projeteis nb
    timestep dt
    tempo atual de simulação t
    tempo maximo dos projeteis tb

    atualiza os valores das posições das naves e projeteis
    */
    int i, j;
    double* gpn1 = gravityPN(P, Ns[0]), *gpn2=gravityPN(P, Ns[1]); /* cálculo da gravidade entre o planeta e as naves 1 e 2 */
    double* gnn  = gravityNN(Ns[0], Ns[1]); /* cálculo da gravidade entre as duas naves */
    double* gbn1,*gbn2,*gpb;                /* cálculo da gravidade para os projéteis */

    for(i = 0; i < 2; i++){
        Ns[i]->a[0] = 0;
        Ns[i]->a[1] = 0;
    }

    Ns[0]->a = Add(Ns[0]->a, Multk(gnn, (double)(-1/Ns[0]->m)));

    gnn = Multk(gnn, (double)(-Ns[0]->m));

    Ns[1]->a = Add(Ns[1]->a, Multk(gnn,(double)(1/Ns[1]->m)));
    Ns[0]->a = Add(Ns[0]->a, Multk(gpn1,(double)(1/Ns[0]->m)));
    Ns[1]->a = Add(Ns[1]->a, Multk(gpn2,(double)(1/Ns[1]->m)));

    /* liberação da memória alocada */
    free(gpn1);
    free(gpn2);
    free(gnn);

    for(i = 0; i < nb; i++){
        Bs[i]->a[0] = 0;
        Bs[i]->a[1] = 0;
    }

    if(t < tb){ /* enquanto o tempo de simulação for menor que o tempo limite dos projéteis */
        for(i = 0; i < nb; i++){
            /* Cálculo da aceleração no i+1-ésimo projétil */
            for(j = i + 1; j < nb; j++){
                gbn1 = gravityBB(Bs[j], Bs[i]);
                Bs[i]->a = Add(Bs[i]->a, Multk(gbn1, (double)(1/Bs[i]->m)));

                gbn1 = Multk(gbn1, (double)(Bs[i]->m));
                Bs[j]->a = Add(Bs[j]->a, Multk(gbn1, (double)(-1/Bs[j]->m)));
                free(gbn1);
            }

            /* Projétil i+1 com Naves e Planeta */
            gbn1 = gravityBN(Bs[i], Ns[0]);
            gbn2 = gravityBN(Bs[i], Ns[1]);
            gpb = gravityPB(P, Bs[i]);

            /* Atualiza aceleração do projétil e das naves */
            Ns[0]->a = Add(Ns[0]->a, Multk(gbn1,(double)(1/Ns[0]->m)));
            Multk(gbn1,Ns[0]->m);
            Ns[1]->a = Add(Ns[1]->a, Multk(gbn2,(double)(1/Ns[1]->m)));
            Multk(gbn2,Ns[0]->m);

            Bs[i]->a = Add(Bs[i]->a, Multk(gbn1,(double)(-1/Bs[i]->m)));
            Bs[i]->a = Add(Bs[i]->a, Multk(gbn2,(double)(-1/Bs[i]->m)));
            Bs[i]->a = Add(Bs[i]->a, Multk(gpb,(double)(1/Bs[i]->m)));

            /* Atualiza orientacao das naves e Projéteis */
            

            /* Atualização da velocidade e posição do i+1-ésimo projétil */

            Bs[i]->v = Add(Bs[i]->v, Multk(Bs[i]->a,dt));
            Bs[i]->c = Add(Bs[i]->c, Multk(Bs[i]->v,dt));
            Bs[i]->v = Multk(Bs[i]->v,(double)(1/dt));

            free(gbn1);
            free(gbn2);
            free(gpb);
        }
        for(i=0;i<nb;i++)
            Bs[i]->o = Orientacao(Bs[i]->a);
    }
    Ns[0]->o = Orientacao(Ns[0]->a);
    Ns[1]->o = Orientacao(Ns[1]->a);

    /* atualiza velocidade e posição das naves */
    Ns[0]->v = Add(Ns[0]->v,Multk(Ns[0]->a,dt));
    Ns[1]->v = Add(Ns[1]->v,Multk(Ns[1]->a,dt));
    Ns[0]->c = Add(Ns[0]->c,Multk(Ns[0]->v,dt));
    Ns[1]->c = Add(Ns[1]->c,Multk(Ns[1]->v,dt));
    Ns[0]->v=Multk(Ns[0]->v,(double)1/dt);
    Ns[1]->v=Multk(Ns[1]->v,(double)1/dt);
}

double AbsVal(double x){
    if(x>=0)return x;
    return -x;
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

double max(double x, double y){
    if(x>y)return x;
    return y;
}

double min(double x, double y){
    if(x<y)return x;
    return y;
}




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

void DrawWindow(WINDOW* w1, PIC Pf, PIC Pbkg,PIC Pplnt, Nave* Ns,double bounds[],Bullet* Bs,int nb,double t,double tb){
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
    if(t<tb){
        for(i=0;i<nb;i++){
            Pf=DrawShip(w1,Bs[i]->spt->Aux[Bs[i]->o],Pf,Pplnt,Bs[i]->spt->P[Bs[i]->o],Bs[i]->spt->Msk[Bs[i]->o],(-bounds[3]+Bs[i]->c[0])*WW/(bounds[2]-bounds[3]),
            (-bounds[1]+Bs[i]->c[1])*WH/(bounds[0]-bounds[1]),bounds);
        }
    }
    PutPic(w1, Pf,  0, 0, WW, WH, 0, 0);
}

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

int ImprimeDados(int bo,double t, double tb,int nb, Nave* Ns,Bullet* Bs){
    /*autoexplicativo*/
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

void Simulate(Bullet* Bs, Nave* Ns, Planet P, int nb, double dt, double tf, double tb,double bounds[],WINDOW* w1,PIC Pbkg,PIC Pplnt, PIC Pf){
    /* recebe
    Vetor com projeteis Bs
    vetor com naves Ns
    Planeta P
    numero de projeteis nb
    timestep dt
    tempo final tf
    tempo maximo dos projeteis tb
    vetor com fronteiras do mapa bounds

    simula em t no intervalo [0,tf] em stepsizes dt,
    calculando as posições dos projéteis e desenhando em uma janela
    */
    #ifdef NOXPM
    puts("Este programa só funciona com a biblioteca Xpm!");
    #else
    double t=0;
    int i, bo=1;
    printf("\nSimulação Iniciada!\n");
    
    for(t = 0; t < tf; t += dt){

        bo=ImprimeDados(bo,t,tb,nb,Ns,Bs);                      /*impressão dos dados*/

        Update(Bs, Ns, P, nb, dt, t, tb);                       /*atualiza valor das posições dos objetos*/

        for(i=0;i<nb;i++){
            Bs[i]->c=CheckLimits(Bs[i]->c,bounds);              /*se a posição nova estiver fora da fronteira*/
        }                                                       /*o objeto reaparece*/
        for(i=0;i<2;i++){                                       /*do outro lado do mapa*/
            Ns[i]->c=CheckLimits(Ns[i]->c,bounds);
        }


        DrawWindow(w1,Pf,Pbkg,Pplnt,Ns,bounds,Bs,nb,t,tb);     /*desenha objetos*/
        usleep(20000);
                                                /*intervalo entre frames*/
    }  
    getchar();
    CloseGraph();
    WDestroy(w1);

    printf("Simulação Finalizada!\n");
    #endif
}

int main(int argc, char*argv[]){
    int i,j,nb;
    double dt = atof(argv[1]) /*timestep dado por linha de comando*/, tb, bounds[4]/*limites do mapa a serem definidos*/;
    double** DadosN = malloc(sizeof(double*)*2), *DadosP=malloc(sizeof(double*)*3), **DadosB;
    char** Nomes = malloc(sizeof(char*)*2);
    Nave* Ns;
    Planet P;
    Bullet* Bs;
  	WINDOW *w1;                           /*janela a ser utilizada*/
    PIC Pf,                               /*Pic final utilizado para sobrepor os objetos 1 a 1*/
        Pplnt,                            /*pic do planeta*/
        Pbkg,                             /*pic do cenario*/
        Paux2;
    MASK
        Mskplnt;                          /*mascara do planeta*/


    w1 = InitGraph(WW,WH, "SpaceWar");
        /*criação das mascaras para tornam o fundo dos objetos transparente*/
    Mskplnt = NewMask(w1,WW,WH);

    Pbkg = ReadPic(w1, "images/background.xpm", NULL); //pic do background
    Pplnt = ReadPic(w1, "images/planeta_rosa.xpm", Mskplnt); //pic do planeta
    Pf = NewPic(w1, WW, WH);

    /*input para criação do Planeta: massa, tempo total e raio*/
    for(i = 0; i < 3; i++) {
        if(!scanf("%lf", &DadosP[i])){
            printf("A leitura dos dados do planeta falhou!");
            return 1;
        }
    }

    P = CriaPlaneta(DadosP[0], DadosP[1], DadosP[2]);
    
    PrepareMask(w1,Pbkg,Mskplnt);


    /*input para criação das Naves*/
    for(i = 0; i < 2; i++) {
        DadosN[i] = malloc(sizeof(double)*5); /* Massa, pos_x, pos_y, vel_x, vel_y */
        Nomes[i]  = malloc(sizeof(char)*20);  /* Nome das naves */
        if(!scanf("%s", Nomes[i])){
            printf("A leitura dos nomes das naves falhou!");
            return 1;
        }


        for(j = 0; j < 5; j++){
            if(!scanf("%lf", &DadosN[i][j])){
                printf("A leitura dos dados das naves falhou");
                return 1;
            }
        }
    }

    Ns = CriaNaves(DadosN, Nomes,w1);

    DefineBoundaries(Ns,bounds);

    PutPic(Pbkg, Pplnt,  0, 0,400,400,(-bounds[3])*WW/(bounds[2]-bounds[3])-200,
    (-bounds[1])*WH/(bounds[0]-bounds[1])-200);

    /*input para criação dos Projéteis*/
    if(!scanf("%d", &nb)){
        printf("A leitura do numero de projeteis falhou");
        return 1;
    }  /* número de projéteis */
    if(!scanf("%lf", &tb)){
        printf("A leitura do tempo de simulação dos projéteis falhou");
        return 1;
    }
     /* tempo de simulação dos projéteis */

    DadosB = malloc(sizeof(double*)*nb);
    for(i = 0; i < nb ; i++){
        DadosB[i] = malloc(sizeof(double)*5);
        for(j = 0; j < 5; j++){
            if(!scanf("%lf", &DadosB[i][j])){
                printf("A leitura dos dados dos projéteis falhou");
                return 1;
            }
        }
    }

    Bs = CriaBullets(DadosB, nb,w1);

    for(i=0;i<nb;i++){
        Bs[i]->c=CheckLimits(Bs[i]->c,bounds);
    }

    printf("Simulação prestes a ser iniciada...\n");
    printf("TimeStep: %.4lf\n",dt);
    printf("Número de projéteis: %d\n",nb);
    printf("Tempo de simulação: %.3lf\n",P->t);
    printf("Limite Superior do Mapa: %.3e\n",bounds[0]);
    printf("Limite Inferior do Mapa: %.3e\n",bounds[1]);
    printf("Limite Direito do Mapa: %.3e\n",bounds[2]);
    printf("Limite Esquerdo do Mapa: %.3e\n",bounds[3]);

    /*simulação do movimento*/
    Simulate(Bs, Ns, P, nb, dt, DadosP[2], tb,bounds,w1,Pbkg,Pplnt,Pf);

    /*desalocando os vetores*/
    for(i = 0; i < 2; i++){
        freeNave(Ns[i]);
        free(DadosN[i]);
    }

    free(Ns);
    free(DadosN);

    for(i = 0; i < nb; i++){
        FreeBullet(Bs[i]);
        free(DadosB[i]);
    }
    free(DadosB);
    free(Bs);
    free(DadosP);
    free(P);
    free(Nomes);

    return 0;
}
