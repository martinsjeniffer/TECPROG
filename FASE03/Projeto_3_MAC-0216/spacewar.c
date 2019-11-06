#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "spwar.h"
#include "xwc.h"

#define G 6.673e-11 //constante gravitacional
#define WW 1500 //altura da janela
#define WH 800  //largura da janela
#define PI 3.14159265359    //pi
#define AC 10   // aceleração gerada pelos controles manuais na nave
#define INITACB 600 //velocidade inicial do projétil lançado


/*********************************/
/*   Imprementação dos Structs   */
/*********************************/


typedef struct bullet *Bullet; /* struct do projetil */

typedef struct nave *Nave; /* struct da nave */

typedef struct planet *Planet; /* struct do Planeta */

typedef struct sprite *Sprite;    /* struct do sprite, que armazenará a imagem do objeto */


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


/*********************************/
/*    Liberação de memória       */
/*********************************/

/*funções para desalocar a memória no final do programa*/

void FreeSprite(Sprite S){
    int i;
    for(i=0;i<16;i++){
        free(S->P[i]);
        free(S->Aux[i]);
    }
    free(S);
}

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


/**************************************/
/*  Colisões, Orientação e Disparos   */
/**************************************/

//orientação em spwar.h

//Confere se o teclado foi apertado. Caso positivo, envia a devida orientação para prosseguir.
int CheckKB(WINDOW *w1){
    int b;
    int key;
    b=WCheckKBD(w1);

    if(b){
        key=WGetKey(w1);
        if(key==25)return 1;
        else if(key==40)return 2;
        else if(key==39)return 3;
        else if(key==38)return 4;
        else if(key==111)return 5;
        else if(key==114)return 6;
        else if(key==116)return 7;
        else if(key==113)return 8;
    }
    return 0;
}

//aceleração extra na nave, gerada por comando de teclado
void BoostN(Nave N,double bounds[]){
    double *v=malloc(sizeof(double)*2);
    double ang=(double)N->o*PI/8-PI/2;
    v[0]=AC*cos(ang);
    v[1]=AC*sin(ang);
    N->v = Add(N->v,v);
    free(v);
}


//aceleração extra do projétil disparado pela nave
void BoostB(Bullet B,double bounds[]){
    double *v=malloc(sizeof(double)*2);
    double ang=(double)B->o*PI/8-PI/2;
    v[0]=INITACB*cos(ang);
    v[1]=INITACB*sin(ang);
    B->v = Add(B->v,v);
    free(v);
}


//checa colisão com planeta
int CheckCollisionPlnt(Planet P,double* c1){
    if(AbsSqrd(c1)<P->r*P->r)
        return 1;
    return 0;
}


//cria um disparo a partir de uma nave e do comando de teclado
void CriaDisparo(WINDOW* w1,Bullet* Bs, Nave* Ns, int i,int nb,double bounds[]){
    Bs[nb+i]=CriaBullet(10,Ns[i]->c[0],Ns[i]->c[1],0,0,w1);
    Bs[nb+i]->o=Ns[i]->o;
    BoostB(Bs[nb+i],bounds);
}

//ações a serem tomadas de acordo com o input do teclado
void ActKB(int dir,int* dispclk, Nave* Ns, Bullet* Bs, int nb,double bounds[],WINDOW* w1){
    switch (dir)
    {
    case 1:                 // acelera nave 1, botão w
        BoostN(Ns[0],bounds);
        break;
    case 2:             // rotaciona nave 1 no sentido horário, botão d
        if(Ns[0]->o==15)Ns[0]->o=-1;
        (Ns[0]->o)++;
        break;
    case 3:              // nave 1 dispara projétil, botão s
        if(dispclk[0]==0)CriaDisparo(w1,Bs,Ns,0,nb,bounds);
        dispclk[0]=1;
        break;
    case 4:             //nave 1 roda no sentido anti-horário, botão a
        if(Ns[0]->o==0)Ns[0]->o=16;
        (Ns[0]->o)--;
        break;
    case 5:             // acelera nave 2, botão uparrow
        BoostN(Ns[1],bounds);
        break;
    case 6:             // rotaciona nave 2 no sentido horário, botão rightarrow
        if(Ns[1]->o==15)Ns[1]->o=-1;
        (Ns[1]->o)++;
        break;
    case 7:             // nave 2 dispara projétil, botão downarrow
        if(dispclk[1]==0)CriaDisparo(w1,Bs,Ns,1,nb,bounds);
        dispclk[1]=1;
        break;
    case 8:             // rotaciona nave 2 no sentido anti-horário, botão leftarrow
        if(Ns[1]->o==0)Ns[1]->o=16;
        (Ns[1]->o)--;
        break;
    default:    //caso contrario, ignorar ação a ser tomada.
        break;
    }
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
/*  Simulação do sistema dinâmico     */
/**************************************/

int Update(Bullet* Bs, Nave* Ns, Planet P, int nb, double dt, double t, double tb,int* dispclk,double bounds[], WINDOW *w1){
    /* recebe
    vetor com projeteis Bs
    vetor com naves Ns
    planeta P
    numero de projeteis nb
    timestep dt
    tempo atual de simulação t
    tempo maximo dos projeteis tb
    vetor cronômetro dos disparos das naves dispclk
    vetor com limites do mapa bounds
    janela w1

    atualiza os valores das posições e orientações das naves e projeteis, e checa se houve colisão no step analisado
    */
    int i, j, dir,col=0;
    double* gpn1 = gravityPN(P, Ns[0]), *gpn2=gravityPN(P, Ns[1]); /* cálculo da gravidade entre o planeta e as naves 1 e 2 */
    double* gnn  = gravityNN(Ns[0], Ns[1]); /* cálculo da gravidade entre as duas naves */
    double* gbn1,*gbn2,*gpb;                /* cálculo da gravidade para os projéteis */

    //checa se o tempo limite dos disparos passou. Caso positivo, destroi projétil.
    if(dispclk[0]>DISPLIM)dispclk[0]=0;
    else if(dispclk[0]>0) (dispclk[0])++;
    if(dispclk[1]>DISPLIM)dispclk[1]=0;
    else if(dispclk[1]>0) (dispclk[1])++;

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

     /* enquanto o tempo de simulação for menor que o tempo limite dos projéteis */
    for(i = 0; i < nb+2; i++){  //os ultimos dois projeteis correspondem aos projéteis lançados pelas naves
        if((t<tb && i<nb) || (dispclk[i-nb] && i>=nb)){
        /* Cálculo da interação entre o i-ésimo projétil e os outros projéteis*/
            for(j = i + 1; j < nb; j++){
                gbn1 = gravityBB(Bs[j], Bs[i]);
                Bs[i]->a = Add(Bs[i]->a, Multk(gbn1, (double)(1/Bs[i]->m)));

                gbn1 = Multk(gbn1, (double)(Bs[i]->m));
                Bs[j]->a = Add(Bs[j]->a, Multk(gbn1, (double)(-1/Bs[j]->m)));
                free(gbn1);

            }
            if(i<nb){   //calculo da interação entre projéteis e naves/planetas
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
                free(gbn1);
                free(gbn2);
                free(gpb);
            }
            if(i<nb || dispclk[i-nb]){  //atualiza posição e velocidade dos projéteis
                Bs[i]->v = Add(Bs[i]->v, Multk(Bs[i]->a,dt));
                Bs[i]->c = Add(Bs[i]->c, Multk(Bs[i]->v,dt));
                Bs[i]->v = Multk(Bs[i]->v,(double)(1/dt));
            }
            if(i<nb){   //checa colisões
                col=CheckCollision(Ns[0]->c,Bs[i]->c,0,0);
                if(col)return col;
                col=CheckCollision(Ns[1]->c,Bs[i]->c,0,0);
                if(col)return col;
            }
            else if(dispclk[i-nb]){ //projéteis não atingem a nave que a disparou
                if(i!=nb)col=CheckCollision(Ns[0]->c,Bs[i]->c,1,dispclk[0]);
                if(col)return col;
                if(i!=nb+1)col=CheckCollision(Ns[1]->c,Bs[i]->c,1,dispclk[1]);
                if(col)return col;
            }
        }
    }
    for(i=0;i<nb;i++)   //atualiza orientações
        Bs[i]->o = Orientacao(Bs[i]->v);


    dir=CheckKB(w1);    //gera ações captadas no teclado
    ActKB(dir,dispclk,Ns,Bs,nb,bounds,w1);

    /* atualiza velocidade e posição das naves */
    Ns[0]->v = Add(Ns[0]->v,Multk(Ns[0]->a,dt));
    Ns[1]->v = Add(Ns[1]->v,Multk(Ns[1]->a,dt));
    Ns[0]->c = Add(Ns[0]->c,Multk(Ns[0]->v,dt));
    Ns[1]->c = Add(Ns[1]->c,Multk(Ns[1]->v,dt));
    Ns[0]->v=Multk(Ns[0]->v,(double)1/dt);
    Ns[1]->v=Multk(Ns[1]->v,(double)1/dt);
    col=CheckCollision(Ns[0]->c,Ns[1]->c,0,0); //checa colisões das naves
    if(col)return col;
    col=CheckCollisionPlnt(P,Ns[0]->c);
    if(col)return col;
    col=CheckCollisionPlnt(P,Ns[1]->c);
    return col;
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
    janela w1
    Pic do plano de fundo Pbkg
    Pic do planeta Pplnt
    Pic a ser utilizado para sobreposições de imagens Pf

    simula em t no intervalo [0,tf] em stepsizes dt,
    calculando as posições dos projéteis e desenhando em uma janela
    */
    #ifdef NOXPM
    puts("Este programa só funciona com a biblioteca Xpm!");
    #else
    double t=0;
    int i, bo=1,col,*dispclk=calloc(2,sizeof(int));
    printf("\nSimulação Iniciada!\n");

    for(t = 0; t < tf; t += dt){

        bo=ImprimeDados(bo,t,tb,nb,Ns,Bs);                      /*impressão dos dados*/

        col=Update(Bs, Ns, P, nb, dt, t, tb,dispclk,bounds, w1);                       /*atualiza valor das posições dos objetos*/

        for(i=0;i<nb;i++){
            Bs[i]->c=CheckLimits(Bs[i]->c,bounds);              /*se a posição nova estiver fora da fronteira*/
        }                                                       /*o objeto reaparece*/
        for(i=0;i<2;i++){
            if(dispclk[i]){
                Bs[nb+i]->c=CheckLimits(Bs[nb+i]->c,bounds);
            }
        }
        for(i=0;i<2;i++){                                       /*do outro lado do mapa*/
            Ns[i]->c=CheckLimits(Ns[i]->c,bounds);
        }


        DrawWindow(w1,Pf,Pbkg,Pplnt,Ns,bounds,Bs,nb,t,tb,dispclk);     /*desenha objetos*/
        usleep(60000);
                                                /*intervalo entre frames*/
        if(col){
            printf("\nHouve Colisão!\n");
            break;
        }
    }
    CloseGraph();
    WDestroy(w1);

    printf("Simulação Finalizada!\n");
    #endif
}

/**************************************/
/*     Inicialização da rotina        */
/**************************************/

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
    InitKBD(w1);
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
        FreeNave(Ns[i]);
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
