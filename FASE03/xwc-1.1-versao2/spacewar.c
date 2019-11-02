//BIBLIOTECA SPACEWAR (CODE)
#include "spacewar.h"

void reset_terminal_mode()
{
    tcsetattr(0, TCSANOW, &orig_termios);
}

void set_conio_terminal_mode()
{
    struct termios new_termios;

    /* take two copies - one for now, one for later */
    tcgetattr(0, &orig_termios);
    memcpy(&new_termios, &orig_termios, sizeof(new_termios));

    /* register cleanup handler, and set the new terminal mode */
    atexit(reset_terminal_mode);
    cfmakeraw(&new_termios);
    tcsetattr(0, TCSANOW, &new_termios);
}

int kbhit()
{
    struct timeval tv = { 0L, 0L };
    fd_set fds;
    FD_ZERO(&fds);
    FD_SET(0, &fds);
    return select(1, &fds, NULL, NULL, &tv);
}

int getch()
{
    int r;
    unsigned char c;
    if ((r = read(0, &c, sizeof(c))) < 0) {
        return r;
    } else {
        return c;
    }
}

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

Nave CriaNave(double m, double x, double y, double vx, double vy, char* nome){
    Nave N = malloc(sizeof(struct nave));
    double *c,*v,*a;

    c = malloc(sizeof(double)*2); /* vetor posição    */
    v = malloc(sizeof(double)*2); /* vetor velocidade */
    a = malloc(sizeof(double)*2); /* vetor aceleração */
    N->nome = nome;
    N->m = m;
    N->c = c;
    N->v = v;
    N->a = a;
    N->o = 0; //começa apontando para cima
    N->c[0] = x;
    N->c[1] = y;

    N->v[0] = vx;
    N->v[1] = vy;

    N->a[0] = 0;
    N->a[1] = 0;

    return N;
}

/*Criação de apenas 2 naves porém pode ser facilmente generalizada para n naves*/
Nave* CriaNaves(double**naves, char**nomes){
    int i = 0;
    Nave* NaveList = malloc(2*sizeof(Nave));
    for(i = 0; i < 2; i++){
        NaveList[i] = CriaNave(naves[i][0],naves[i][1],naves[i][2],naves[i][3],naves[i][4],nomes[i]);
    }
    return NaveList;
}

Planet CriaPlaneta(double r,double m,double t){
    Planet P = malloc(sizeof(struct planet));
    P->m = m; /* massa */
    P->t = t; /* tempo  de simulação*/
    P->r = r; /* raio  */
    return P;
}

Bullet CriaBullet(double m,double x, double y,double vx,double vy){
    Bullet B = malloc(sizeof(struct bullet));
    double *c,*v,*a;
    int i = 0; //contador
    struct dirent *de;  // Pointer for directory entry
    c = malloc(sizeof(double)*2); /* vetor posição    */
    v = malloc(sizeof(double)*2); /* vetor velocidade */
    a = malloc(sizeof(double)*2); /* vetor aceleração */
    B->m = m;
    B->c = c;
    B->v = v;
    B->a = a;
    B->o = 0;

    B->c[0] = x;
    B->c[1] = y;

    B->v[0] = vx;
    B->v[1] = vy;

    B->a[0] = 0;
    B->a[1] = 0;

    B->o = 0;

    return B;
}

/* Criação dos projéteis fornecidos na entrada padrão */
Bullet* CriaBullets(double**bullets, int nb){
    int i = 0;
    Bullet* BulletList = malloc(nb*sizeof(Bullet));
    for(i = 0; i < nb; i++){
        BulletList[i] = CriaBullet(bullets[i][0], bullets[i][1], bullets[i][2], bullets[i][3], bullets[i][4]);
    }
    return BulletList;
}

void DefineBoundaries(Nave* Ns,double bounds[]){
    /* Os limites são definidos de acordo com as posições das naves.
     * O limite superior, por exemplo, é definido como o máximo entre 0 e o dobro da maior das coordenadas y das naves
     * Analogamente, o mínimo é definido o mínimo entre 0 e o dobro da menor das coordenadas y das naves
     * O mesmo é feito para os limites direito e esquerdo
     */
    bounds[0] = max(0,2*max(Ns[0]->c[1],Ns[1]->c[1])); /*UpperBound*/
    bounds[1] = min(0,2*min(Ns[0]->c[1],Ns[1]->c[1])); /*LowerBound*/
    bounds[2] = max(0,2*max(Ns[0]->c[0],Ns[1]->c[0])); /*RightBound*/
    bounds[3] = min(0,2*min(Ns[0]->c[0],Ns[1]->c[0])); /*LeftBound */
}

int ImprimeDados(int bo,double t, double tb,int nb, Nave* Ns,Bullet* Bs){
    int i;
    if(bo && t>tb){
            bo=0;
            printf("Tempo Limite de simulação dos projéteis alcançado!\n\n");
        }

    /* resultados impressos na saída padrão */
    printf("Tempo %.3lf:\n",t);
    for(i = 0; i < 2; i++){
            printf("%5s:      x: %12.4e | y: %12.4e |",Ns[i]->nome, Ns[i]->c[0], Ns[i]->c[1]);
            printf(" vx: %12.4e | vy: %12.4e\n", Ns[i]->v[0], Ns[i]->v[1]);
    }
    if(bo)
            for(i = 0; i < nb; i++){
                printf("projétil%d:  x: %12.4e | y: %12.4e |", i + 1, Bs[i]->c[0], Bs[i]->c[1]);
                printf(" vx: %12.4e | vy: %12.4e \n", Bs[i]->v[0], Bs[i]->v[1]);
            }
    printf("\n");
    return bo;
}

/* Quadrado do módulo do vetor c */
double AbsSqrd(double* c){
    return (c[0])*(c[0])+(c[1])*(c[1]);
}

/* Multiplica o vetor c por uma constante k */
double* Multk(double* c, double k){
    double newx = c[0]*k, newy=c[1]*k;
    c[0] = newx;
    c[1] = newy;
    return c;
}

/* Adiciona o conteúdo de c2 em c */
double* Add(double* c,double* c2){
    c[0] = c[0] + c2[0];
    c[1] = c[1] + c2[1];
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
    double* c2 = malloc(sizeof(double)*2);
    c2[0] = c[0];
    c2[1] = c[1];
    return c2;
}

/* Força gravicional entre Planeta e Projétil */
double* gravityPB(Planet P, Bullet B){
    double* g =Copia(B->c);
    double dist = (double)(-G*P->m*B->m/(pow(AbsSqrd(B->c),3/2)));
    g=Multk(g,dist);
    return g;
}

/* Interação entre o Planeta e Nave */
double* gravityPN(Planet P, Nave N){
    double* g = Copia(N->c);
    double dist = (double)(-G* P->m * N->m/(pow(AbsSqrd(N->c),3/2)));

    g = Multk(g, dist);
    return g;
}

/* Interação entre as naves */
double* gravityNN(Nave N1,Nave N2){
    double* g    = Sub(Copia(N2->c), N1->c),
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

double AbsVal(double x){
    if(x >= 0)return x;
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
    if(x > y)return x;
    return y;
}

double min(double x, double y){
    if(x < y)return x;
    return y;
}

//0 aponta para cima, 4 para direita, 8 para baixo, 12 para esquerda
//só quem define a orientação é o controle
int Orientacao(int o, int n)
{
  /*
  Vetor de imagens com a seguinte organização (de 0 a 2n angulos a direita e 2n a 4n angulos a esquerda):
  0... |n  ...  |2n... |3n  .....|4n
  cima |direita |baixo |esquerda |cima
  Devolve uma int orientacao
  Só o que define a orientação é o controle
  Recebe n, se n>0, está girando para direita, n<0, para a esquerda
  */
  o = o + n;
  if (o>=16)
    return 0;
  return o;
}
//relação atrações
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

//recebe duas posições e calcula distancia entre elas
double Dist(double *d1, double *d2)
{
  double aux = sqrt((d1[0]-d2[0])*(d1[0]-d2[0])+(d1[1]-d2[1])*(d1[1]-d2[1]));
  return aux;
}

int Colisao(double *d1, double *d2, double r1, double r2)
{
  /*
  recebe posiçoes d1 e d2,
  raios r1 e r2;
  retorna 0 se não colide, 1 se colide
  */
  if (Dist(d1, d2)<=r1+r2)
    return 1;
  return 0;
}

//zera força
double* Zera(double* f)
{
  f[0]=0;
  f[1]=0;
  return f;
}

//checa se c pertence a s e devolve i+1; 0 se não pertence
int Valid(int c, char s[10]){
  /*
  Checa se c pertence a string s;
  devolve qual char da string é igual +1
  */
  int i;
  for (i=0; i < strlen(s); i++)
    if (c==s[i])
      return i+1;
  return 0;
}
