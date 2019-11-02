#include "spacewar.h"
#include "xwc.h"

#define WW 1500
#define WH 800

#define Fator 1 //quao efetivo é apertar os comandos

typedef struct sprite *Sprite;    /* struct do sprite */

struct sprite {
  PIC *P;
  PIC *Aux;
  MASK *Msk;
};

//norb 0 é bullet
Sprite CriaSprite(int norb, WINDOW* w1){
    Sprite S = malloc(sizeof(struct sprite));
    int i;
    char str[50],strint[10],strnave[2];
    S->P   = malloc(4*NumIma*(sizeof(PIC)));
    S->Aux = malloc(4*NumIma*(sizeof(PIC)));
    S->Msk = malloc(4*NumIma*(sizeof(MASK)));

    if(norb == 0){
        for(i = 0; i < 4*NumIma; i++){
            S->Msk[i] = NewMask(w1,WW,WH);
            strcpy (str,"images/bullet/tiro");
            sprintf(strint,"%d",i);
            strcat (str,strint);
            strcat (str,".xpm");
            S->P[i]   = ReadPic(w1, str, S->Msk[i]);
            S->Aux[i] = NewPic(w1, WW, WH);
            PrepareMask(w1, S->Aux[i], S->Msk[i]);
        }
        return S;
    }

    for(i = 0; i < 4*NumIma; i++){
        S->Msk[i] = NewMask(w1,WW,WH);
        strcpy (str,"images/nave/nave");
        sprintf(strnave,"%d",norb);
        sprintf(strint,"%d",i);
        strcat (str,strnave);
        strcat (str,"-");
        strcat (str,strint);
        strcat (str,".xpm");
        S->P[i]   = ReadPic(w1, str, S->Msk[i]);
        S->Aux[i] = NewPic(w1, WW, WH);
        PrepareMask(w1, S->Aux[i], S->Msk[i]);
    }
    return S;
}

//devolve nova aceleração após aplicar forca na orientação o, mod in
double* AplicaForca(double *a, double m, int o, double in)
{
  double *f;
  double alpha = o*(M_PI/(NumIma*2));
  f= malloc(sizeof(double)*2);
  f[1]=sin(alpha)*in;
  f[0]=cos(alpha)*in;
  return Add(a, Multk(f, (double)(1/m)));
}
//relação atrações
int Update(Bullet* Bs, Nave* Ns, Planet P, int nb, double dt, double t, double tb){
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
    double* aj = malloc(sizeof(double)*2);  /* aceleração da jogada */

    for(i = 0; i < 2; i++){
        Ns[i]->a[0] = 0;
        Ns[i]->a[1] = 0;
    }

    Ns[0]->a = Add(Ns[0]->a, Multk(gnn, (double)(-1/Ns[0]->m)));

    gnn = Multk(gnn, (double)(-Ns[0]->m));

    Ns[1]->a = Add(Ns[1]->a, Multk(gnn,(double)(1/Ns[1]->m)));
    Ns[0]->a = Add(Ns[0]->a, Multk(gpn1,(double)(1/Ns[0]->m)));
    Ns[1]->a = Add(Ns[1]->a, Multk(gpn2,(double)(1/Ns[1]->m)));

    /*adicionando a jogada
    Ns[0]->a = Add(Ns[0]->a, aj0);
    Ns[1]->a = Add(Ns[1]->a, aj1);

    /* liberação da memória alocada */
    free(gpn1);
    free(gpn2);
    free(gnn);
    free(aj);

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
            //Bs[i]->o = Orientacao(Bs[i]->a);

            /* Atualização da velocidade e posição do i+1-ésimo projétil */

            Bs[i]->v = Add(Bs[i]->v, Multk(Bs[i]->a,dt));
            Bs[i]->c = Add(Bs[i]->c, Multk(Bs[i]->v,dt));
            Bs[i]->v = Multk(Bs[i]->v,(double)(1/dt));

            free(gbn1);
            free(gbn2);
            free(gpb);
        }
    }

    /* atualiza velocidade e posição das naves */
    Ns[0]->v = Add(Ns[0]->v,Multk(Ns[0]->a,dt));
    Ns[1]->v = Add(Ns[1]->v,Multk(Ns[1]->a,dt));
    Ns[0]->c = Add(Ns[0]->c,Multk(Ns[0]->v,dt));
    Ns[1]->c = Add(Ns[1]->c,Multk(Ns[1]->v,dt));
    Ns[0]->v = Multk(Ns[0]->v,(double)1/dt);
    Ns[1]->v = Multk(Ns[1]->v,(double)1/dt);
}

void PrepareMask(WINDOW* w1, PIC P, MASK Msk){
   /*   recebe
    *       - window w1
    *       - pic auxiliar para aplicar mascara P
    *       - mascara a ser aplicada Msk
    *
    *   aplica Msk em P
    */

  PutPic(w1, P,  0, 0, WW, WH, 0, 0);
  SetMask(P,Msk);
}

PIC DrawShip(WINDOW* w1,PIC Paux,PIC Pbkg, PIC P,MASK Msk,int i,int j){
    /*  recebe
    *       - window w1
    *       - pic auxiliar Paux
    *       - pic do background Pbkg
    *       - pic do objeto P
    *       - mask do objeto Msk
    *
    *   desenha o objeto em Paux na posição (i,j) e retorna Paux
    */

    UnSetMask(Paux);
    PutPic(Paux, Pbkg,  0, 0, WW, WH, 0, 0);
    SetMask(Paux,Msk);
    PutPic(Paux, P,  0, 0, 100, 100, i, j);
    return Paux;
}

void DrawWindow(WINDOW* w1, PIC Pf, PIC Pbkg, Nave* Ns,double bounds[],Bullet* Bs,int nb,double t,double tb, Sprite Sn1, Sprite Sn2){
   /*  recebe
    *       - window w1
    *       - pic transitorio para sobrepor todos os objetos Pf
    *       - pic da nave1 P1
    *       - pic da nave2 P2
    *       - pic auxiliar da nave1 Paux1
    *       - pic do background Pbkg
    *       - pic auxiliar da nave2 Paux2
    *       - mask da nave1 Msk1
    *       - mask da nave2 Msk2
    *       - Vetor contendo as duas naves Ns
    *       - fronteiras do mapa bounds
    *       - vetor de pics auxiliares dos projeteis Paux
    *       - vetor de pics dos projeteis Pproj
    *       - vetor de masks dos projeteis MaskProj
    *       - Vetor contendo projeteis Bs
    *       - numero de projeteis nb
    *       - tempo atual t
    *       - tempo limite dos projeteis tb
    *
    *   desenha todos os objetos em w1
    *
    *   as posições calculadas são recalculadas para as dimensões da tela
    */

    int i;

    Pf = DrawShip(w1, Sn1->Aux[Ns[0]->o],
                  Pbkg,Sn1->P[Ns[0]->o],
                  Sn1->Msk[Ns[0]->o],
                  (-bounds[3]+Ns[0]->c[0])*WW/(bounds[2]-bounds[3]),
                  (-bounds[1]+Ns[0]->c[1])*WH/(bounds[0]-bounds[1]));

    Pf = DrawShip(w1,Sn2->Aux[Ns[1]->o],
                  Pf,Sn2->P[Ns[1]->o],
                  Sn2->Msk[Ns[1]->o],
                  (-bounds[3]+Ns[1]->c[0])*WW/(bounds[2]-bounds[3]),
                  (-bounds[1]+Ns[1]->c[1])*WH/(bounds[0]-bounds[1]));
  /*  if(t < tb){
        for(i = 0; i < nb; i++){
            Pf = DrawShip(w1,Bs[i]->spt->Aux[Bs[i]->o],
                          Pf,Bs[i]->spt->P[Bs[i]->o],
                          Bs[i]->spt->Msk[Bs[i]->o],
                          (-bounds[3]+Bs[i]->c[0])*WW/(bounds[2]-bounds[3]),
                          (-bounds[1]+Bs[i]->c[1])*WH/(bounds[0]-bounds[1]));
        }
    }*/
    PutPic(w1, Pf,  0, 0, WW, WH, 0, 0);
}

void Game(Bullet* Bs, Nave* Ns, Planet P, int nb, double dt, double tf, double tb,double bounds[],WINDOW* w1,PIC Pbkg, PIC Pf, Sprite Sn1, Sprite Sn2){
    /*  recebe
     *      - Vetor com projeteis Bs
     *      - vetor com naves Ns
     *      - Planeta P
     *      - numero de projeteis nb
     *      - timestep dt
     *      - tempo final tf
     *      - tempo maximo dos projeteis tb
     *      - vetor com fronteiras do mapa bounds
     *
     *  simula em t no intervalo [0,tf] em stepsizes dt,
     *  calculando as posições dos projéteis e desenhando em uma janela
    */
    #ifdef NOXPM
    puts("Este programa só funciona com a biblioteca Xpm!");
    #else
    double t  = 0;
    int i, bo = 1;
    int c = 0; // tecla
    printf("\nSimulação Iniciada!\n");

    for(t = 0; t < tf; t += dt){

        bo = ImprimeDados(bo,t,tb,nb,Ns,Bs);                      /* impressão dos dados*/

        Update(Bs, Ns, P, nb, dt, t, tb);                         /* atualiza valor das posições dos objetos*/

      //  for(i = 0; i < nb; i++){
        //    Bs[i]->c = CheckLimits(Bs[i]->c,bounds);              /* se a posição nova estiver fora da fronteira*/
        //}                                                         /* o objeto reaparece*/

        for(i = 0; i < 2; i++){                                   /*do outro lado do mapa*/
            Ns[i]->c=CheckLimits(Ns[i]->c,bounds);
        }


        DrawWindow(w1,Pf,Pbkg,Ns,bounds,Bs,nb,t,tb, Sn1, Sn2);     /*desenha objetos*/
        usleep(20000);
                                                        /*intervalo entre frames*/
    }

    getchar();
    CloseGraph();
    WDestroy(w1);

    printf("Simulação Finalizada!\n");
    #endif
}

//cria vetor aceleração da jogada
double* Jogada(double* f, int v)
{
  /*
  recebe o que foi a jogada (v)
  devolve vetor força sendo f[0]=horizontal
  */
  switch (v) {
    case 1: // cima
      f[0]=0;
      f[1]=1*Fator;
    break;

    case 2: // esquerda
      f[0]=-1*Fator;
      f[1]=0;
    break;

    case 3: // baixo
      f[0]=0;
      f[1]=-1*Fator;
    break;

    case 4: // direita
      f[0]=1*Fator;
      f[1]=0;
    break;

    case 5: // tiro
      f[0]=0;
      f[1]=0;
      puts("BANG!");
    break;
  }
  return f;
}

int main(int argc, char*argv[]){
    int i,j,nb;
    int c; //entrada
    int game_over=0;
    double dt = atof(argv[1]) /* timestep dado por linha de comando */,
           tb, bounds[4]      /* limites do mapa a serem definidos  */;
    double** DadosN = malloc(sizeof(double*)*2), *DadosP = malloc(sizeof(double*)*3), **DadosB;
    char**   Nomes  = malloc(sizeof(char*)*2);
    Nave*   Ns;
    Planet  P;
    Bullet* Bs;
  	WINDOW *w1;                           /* janela a ser utilizada*/
    PIC Pf,                               /* pic final utilizado para sobrepor os objetos 1 a 1*/
        Pplnt,                            /* pic do planeta*/
        Pbkg,                             /* pic do cenario*/
        Paux2;
      //  *Pproj=malloc(nb*sizeof(PIC)),    /*vetor de pics dos projeteis*/
      //  *Paux=malloc(nb*sizeof(PIC));     /*vetor de pics auxiliares dos projeteis*/
    MASK
        Mskplnt;                          /*mascara do planeta*/
      //  *MskProj=malloc(nb*sizeof(MASK)); /*mascara dos projeteis*/
    Sprite Sn1, Sn2; //sprite naves

    //PREPARAÇÕES INICIAIS
    //inicia janela
    w1 = InitGraph(WW,WH, "Arquivos");
    //mascara no planeta e coloca background
    Mskplnt = NewMask(w1,WW,WH);
    Pbkg = ReadPic(w1, "images/background.xpm", NULL); //pic do background
    Pplnt = ReadPic(w1, "images/planeta_rosa.xpm", Mskplnt); //pic do planeta
    //cria pic principal
    Pf = NewPic(w1, WW, WH);
    /* input para criação do Planeta: massa, tempo total e raio */
    for(i = 0; i < 3; i++) {
        if(!scanf("%lf", &DadosP[i])){
            printf("A leitura dos dados do planeta falhou!");
            return 1;
        }
    }
    P = CriaPlaneta(DadosP[0], DadosP[1], DadosP[2]);
    PrepareMask(w1,Pbkg,Mskplnt);
    PutPic(Pbkg, Pplnt, 0, 0, WW, WH, -WW/10, -WH/10);
    /* input para criação das Naves */
    for(i = 0; i < 2; i++) {
        DadosN[i] = malloc(sizeof(double)*5); /* Massa, pos_x, pos_y, vel_x, vel_y */
        Nomes[i]  = malloc(sizeof(char)*20);  /* Nome das naves */
        if(!scanf("%s", Nomes[i])) {
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
    //lista de naves
    Ns = CriaNaves(DadosN, Nomes);
    Sn1 = CriaSprite(1, w1);
    Sn2 = CriaSprite(2, w1);

    DefineBoundaries(Ns,bounds);
/*TIREI PROJETEIS POR ENQUANTO*/
    /*input para criação dos Projéteis*/
    /* número de projéteis */
    nb = 0;
    /*if(!scanf("%d", &nb)){
        printf("A leitura do numero de projeteis falhou");
        return 1;
    }

    /* tempo de simulação dos projéteis */
    /*if(!scanf("%lf", &tb)){
        printf("A leitura do tempo de simulação dos projéteis falhou");
        return 1;
    }
*/
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

    //lista de bullets de cada nave
    Bs = CriaBullets(DadosB, nb);

    for(i = 0;i < nb; i++){
        Bs[i]->c = CheckLimits(Bs[i]->c,bounds);
    }

    printf("Simulação prestes a ser iniciada...\n");
    printf("TimeStep: %.4lf\n",dt);
    printf("Número de projéteis: %d\n",nb);
    printf("Tempo de simulação: %.3lf\n",P->t);
    printf("Limite Superior do Mapa: %.3e\n",bounds[0]);
    printf("Limite Inferior do Mapa: %.3e\n",bounds[1]);
    printf("Limite Direito do Mapa: %.3e\n",bounds[2]);
    printf("Limite Esquerdo do Mapa: %.3e\n",bounds[3]);

    //set_conio_terminal_mode();
  /*
    while (!game_over){
      while (!kbhit()) {
        //roda enquanto não recebe nada*/
        Game(Bs, Ns, P, nb, dt, DadosP[2], tb,bounds,w1,Pbkg,Pf, Sn1, Sn2);

      c = getch();


    /* simulação do movimento */
    //Simulate(Bs, Ns, P, nb, dt, DadosP[2], tb,bounds,w1,Pbkg,Pf);


    /* desalocando os vetores */
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
