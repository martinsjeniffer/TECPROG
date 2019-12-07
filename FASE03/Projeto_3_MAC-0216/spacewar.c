#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "spwar.h"


/**************************************/
/*  Colisões, Orientação e Disparos   */
/**************************************/

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
//Cria de devolve Sprite de explosão
Sprite CriaExplo (WINDOW* w1)
{
  Sprite S = malloc(sizeof(struct sprite));
  int i;
  char str[50],strint[10];
  S->P=malloc(4*NumIma*(sizeof(PIC)));
  S->Aux=malloc(4*NumIma*(sizeof(PIC)));
  S->Msk=malloc(4*NumIma*(sizeof(MASK)));
  for(i=0;i<10;i++){                  //10 foi o número escolhido de frames
      S->Msk[i] = NewMask(w1,WW,WH);
      strcpy(str,"images/explo/frame");
      sprintf(strint,"%d",i);
      strcat(str,strint);
      strcat(str,".xpm");
      S->P[i] = ReadPic(w1, str, S->Msk[i]);
      S->Aux[i] = NewPic(w1, WW, WH);
      PrepareMask(w1, S->Aux[i], S->Msk[i]);
  }
  return S;
}

//cria e devolve sprite da intro
Sprite CriaIntro (WINDOW* w1)
{
  Sprite S = malloc(sizeof(struct sprite));
  int i;
  char str[50],strint[10];
  S->P=malloc(3*(sizeof(PIC)));
  S->Aux=malloc(3*(sizeof(PIC)));
  S->Msk=malloc(3*(sizeof(MASK)));
  for(i=0;i<3;i++){                  //10 foi o número escolhido de frames
      S->Msk[i] = NewMask(w1,WW,WH);
      strcpy(str,"images/intro/intro");
      sprintf(strint,"%d",i);
      strcat(str,strint);
      strcat(str,".xpm");
      S->P[i] = ReadPic(w1, str, S->Msk[i]);
      S->Aux[i] = NewPic(w1, WW, WH);
      PrepareMask(w1, S->Aux[i], S->Msk[i]);
  }
  return S;
}

//devolve pic com P2 em (i,j) sobre P1, Msk é a mascara já aplicada em P2
PIC DrawOver (PIC P1, PIC P2, PIC Paux, MASK Msk, int i, int j)
{
  UnSetMask(Paux);
  PutPic(Paux, P1,  0, 0, WW, WH, 0, 0);
  SetMask(Paux,Msk);
  PutPic(Paux, P2,  0, 0, WW, WH, i, j);
  return Paux;
}

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

int CheckKBrestart(WINDOW *w1){
    int b;
    int key;
    b=WCheckKBD(w1);

    if(b){
        key=WGetKey(w1);
        printf("\n%dKEY\n",key);
        if(key==29)return 1;
        else if(key==57)return 2;
    }
    return 0;
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

//passa as imagens da pasta intro com o backgroundno
void Intro (WINDOW* w1, PIC Pbkg)
{
  PIC Pf;
  int i, g;
  Sprite Intro = CriaIntro(w1);

  Pf = NewPic(w1, WW, WH);
  for (i=0; i<3; i++) {
    //passa as instruções uma por uma e coloca as coloca na posição (300,200)
    Pf = DrawOver (Pbkg, Intro->P[i], Intro->Aux[i], Intro->Msk[i], 300, 200);
    PutPic(w1, Pf,  0, 0, WW, WH, 0, 0);
    //checa se alguma tecla foi apertada, no caso apenas teclas do jogo (a, w, s, d e setas)
    g=0;
    do {
      g=CheckKB(w1);
    } while(!g);
    //para evitar que no mesmo apertar passe mais de uma instrução
    usleep(60000);
  }
  FreeSprite (Intro, 3);
}

//faz a animação da colisão e diz o vencedor
void GameOver (int win, WINDOW* w1, PIC Pf, PIC Pbkg,PIC Pplnt, Nave* Ns,double bounds[],Bullet* Bs,int nb,double t,double tb,int* dispclk,int pontos[])
{
  /*
  win = 3 naves colidiram, 1 player1 ganha, 2 player2 ganha
  */
  Sprite Explo = CriaExplo(w1);
  int i, j;
  PIC GO, GOA; //pic e mask do game over, e pic auxiliar
  MASK MGO;

  MGO = NewMask(w1, WW, WH);
  GOA = NewPic (w1, WW, WH);

/* Animação da explosão */
  if (win == 3) { //navexnave
    for (j=0; j<10; j++) {
      Pf=DrawShip(w1,Explo->Aux[j],Pbkg,Pplnt,Explo->P[j],Explo->Msk[j],(-bounds[3]+Ns[0]->c[0])*WW/(bounds[2]-bounds[3]),
      (-bounds[1]+Ns[0]->c[1])*WH/(bounds[0]-bounds[1]),bounds);

      Pf=DrawShip(w1,Explo->Aux[j],Pf,Pplnt,Explo->P[j],Explo->Msk[j],(-bounds[3]+Ns[1]->c[0])*WW/(bounds[2]-bounds[3]),
      (-bounds[1]+Ns[1]->c[1])*WH/(bounds[0]-bounds[1]),bounds);

      for(i=0;i<nb+2;i++){
        if((t<tb && i<nb) || (i>=nb && dispclk[i-nb])){
          Pf=DrawShip(w1,Bs[i]->spt->Aux[Bs[i]->o],Pf,Pplnt,Bs[i]->spt->P[Bs[i]->o],Bs[i]->spt->Msk[Bs[i]->o],(-bounds[3]+Bs[i]->c[0])*WW/(bounds[2]-bounds[3]),
          (-bounds[1]+Bs[i]->c[1])*WH/(bounds[0]-bounds[1]),bounds);
        }
      }
      PutPic(w1, Pf,  0, 0, WW, WH, 0, 0);
      usleep(20000);
    }
      GO = ReadPic(w1, "images/gameover/go3.xpm", MGO); //pic do planeta
      Pf = DrawOver (Pf, GO, GOA, MGO, 550, 100);
      puts ("Ambas as naves foram destruídas");
  }
  else if (win == 1) { //nave1 explodiu
    for (j=0; j<10; j++) {
      Pf=DrawShip(w1,Ns[0]->spt->Aux[Ns[0]->o],Pbkg,Pplnt,Ns[0]->spt->P[Ns[0]->o],Ns[0]->spt->Msk[Ns[0]->o],(-bounds[3]+Ns[0]->c[0])*WW/(bounds[2]-bounds[3]),
      (-bounds[1]+Ns[0]->c[1])*WH/(bounds[0]-bounds[1]),bounds);

      Pf=DrawShip(w1,Explo->Aux[j],Pbkg,Pplnt,Explo->P[j],Explo->Msk[j],(-bounds[3]+Ns[1]->c[0])*WW/(bounds[2]-bounds[3]),
      (-bounds[1]+Ns[1]->c[1])*WH/(bounds[0]-bounds[1]),bounds);

      for(i=0;i<nb+2;i++){
        if((t<tb && i<nb) || (i>=nb && dispclk[i-nb])){
          Pf=DrawShip(w1,Bs[i]->spt->Aux[Bs[i]->o],Pf,Pplnt,Bs[i]->spt->P[Bs[i]->o],Bs[i]->spt->Msk[Bs[i]->o],(-bounds[3]+Bs[i]->c[0])*WW/(bounds[2]-bounds[3]),
          (-bounds[1]+Bs[i]->c[1])*WH/(bounds[0]-bounds[1]),bounds);
        }
      }
      PutPic(w1, Pf,  0, 0, WW, WH, 0, 0);
      usleep(20000);
    }
    GO = ReadPic(w1, "images/gameover/go2.xpm", MGO); //pic do planeta
    Pf = DrawOver (Pf, GO, GOA, MGO, 450, 50);
    puts ("Player 2 venceu!");
    (pontos[1])++;

  }
  else if (win == 2) { //nave0 explodiu
    for (j=0; j<10; j++) {
      Pf=DrawShip(w1,Explo->Aux[j],Pbkg,Pplnt,Explo->P[j],Explo->Msk[j],(-bounds[3]+Ns[0]->c[0])*WW/(bounds[2]-bounds[3]),
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
      usleep(20000);
    }
    GO = ReadPic(w1, "images/gameover/go1.xpm", MGO); //pic do planeta
    Pf = DrawOver (Pf, GO, GOA, MGO, 450, 50);
    puts ("Player 1 venceu!");
    (pontos[0])++;
  }
  PutPic(w1, Pf,  0, 0, WW, WH, 0, 0);
  usleep(60000);

  FreeSprite (Explo,10);
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
                if(col)return 2;
                col=CheckCollision(Ns[1]->c,Bs[i]->c,0,0);
                if(col)return 1;
            }
            else if(dispclk[i-nb]){ //projéteis não atingem a nave que a disparou
                if(i!=nb)col=CheckCollision(Ns[0]->c,Bs[i]->c,1,dispclk[0]);
                if(col)return 2;
                if(i!=nb+1)col=CheckCollision(Ns[1]->c,Bs[i]->c,1,dispclk[1]);
                if(col)return 1;
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
    if(col)return 3; //naves colidiram entre si
    col=CheckCollisionPlnt(P,Ns[0]->c);
    if(col)return 2; //player2 ganhou
    col=CheckCollisionPlnt(P,Ns[1]->c);
    if(col)return 1; //player1 ganhou
    return col;
}

void Simulate(Bullet* Bs, Nave* Ns, Planet P, int nb, double dt, double tf, double tb,double bounds[],WINDOW* w1,PIC Pbkg,PIC Pplnt, PIC Pf,int pontos[]){
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
    int g;
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
            GameOver(col,w1,Pf,Pbkg,Pplnt,Ns,bounds,Bs,nb,t,tb,dispclk,pontos);
            break;
        }
    }

    //checa se alguma tecla foi apertada, no caso apenas teclas do jogo (a, w, s, d e setas)
    g=0;
    do {
      g=CheckKB(w1);
    } while(!g);

    printf("Simulação Finalizada!\n");
    #endif
}


/**************************************/
/*     Inicialização da rotina        */
/**************************************/

int main(int argc, char*argv[]){
    int i,j,nb,counter,pontos[2],stopcond=1,key;
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
    pontos[0]=pontos[1]=0;
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

    Intro (w1, Pbkg);

    while(stopcond){
        Ns = CriaNaves(DadosN, Nomes,w1);

        DefineBoundaries(Ns,bounds);

        Bs = CriaBullets(DadosB, nb,w1);

        for(i=0;i<nb;i++){
            Bs[i]->c=CheckLimits(Bs[i]->c,bounds);
        }

        PutPic(Pbkg, Pplnt,  0, 0,400,400,(-bounds[3])*WW/(bounds[2]-bounds[3])-200,
        (-bounds[1])*WH/(bounds[0]-bounds[1])-200);

        printf("Simulação prestes a ser iniciada...\n");
        printf("TimeStep: %.4lf\n",dt);
        printf("Número de projéteis: %d\n",nb);
        printf("Tempo de simulação: %.3lf\n",P->t);
        printf("Limite Superior do Mapa: %.3e\n",bounds[0]);
        printf("Limite Inferior do Mapa: %.3e\n",bounds[1]);
        printf("Limite Direito do Mapa: %.3e\n",bounds[2]);
        printf("Limite Esquerdo do Mapa: %.3e\n",bounds[3]);

        /*simulação do movimento*/
        Simulate(Bs, Ns, P, nb, dt, DadosP[2], tb,bounds,w1,Pbkg,Pplnt,Pf,pontos);
        while(1){
            key=CheckKBrestart(w1);
            if(key==1)
                break;
            else if (key==2){
                stopcond=0;
                break;
            }
        }
        printf("\nprontuacao: %d X %d\n",pontos[0],pontos[1]);
    }

    CloseGraph();
    WDestroy(w1);
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
