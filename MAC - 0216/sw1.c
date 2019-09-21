#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define G 1

typedef struct bullet *Bullet;

typedef struct nave *Nave;

typedef struct planet *Planet;

struct bullet{ /* Projéteis */
    double m;
    double *c,*v,*a;
};

struct nave{  
    char* nome;
    double m;
    double *c,*v,*a;
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


Nave CriaNave(double m, double x, double y, double vx, double vy, char* nome){
    Nave N = malloc(sizeof(struct nave));
    double *c,*v,*a;
    c = malloc(sizeof(double)*2); /* vetor posição    */
    v = malloc(sizeof(double)*2); /* vetor velocidade */
    a = malloc(sizeof(double)*2); /* vetor aceleração */
    N->nome=nome;
    N->m=m;
    N->c=c;
    N->v=v;
    N->a=a;
    N->c[0]=x;
    N->c[1]=y;

    N->v[0]=vx;
    N->v[1]=vy;
    
    N->a[0]=0;
    N->a[1]=0;
    
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
    P->m=m; /* massa */
    P->t=t; /* tempo  de simulação*/
    P->r=r; /* raio  */
    return P;
}

Bullet CriaBullet(double m,double x, double y,double vx,double vy){
    Bullet B = malloc(sizeof(struct nave));
    double *c,*v,*a;
    c = malloc(sizeof(double)*2); /* posições    */
    v = malloc(sizeof(double)*2); /* velocidade  */
    a = malloc(sizeof(double)*2); /* aceleração  */
    B->m=m;
    B->c=c;
    B->v=v;
    B->a=a;

    B->c[0]=x;
    B->c[1]=y;
    
    B->v[0]=vx;
    B->v[1]=vy;
    
    B->a[0]=0;
    B->a[1]=0;
    
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
    (c[0])+=(c2[0]);
    (c[1])+=(c2[1]);
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

void Update(Bullet* Bs, Nave* Ns, Planet P, int nb, double dt, double t, double tb){
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
            
            /* Atualização da velocidade e posição do i+1-ésimo projétil */

            Bs[i]->v = Add(Bs[i]->v, Multk(Bs[i]->a,dt));
            Bs[i]->c = Add(Bs[i]->c, Multk(Bs[i]->v,dt));
            
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
}

void Simulate(Bullet* Bs, Nave* Ns, Planet P, int nb, double dt, double tf, double tb){
    /*simulação em t no intervalo [0,tf] em stepsizes dt */
    double t = 0;
    for(t = 0; t < tf; t += dt){
        Update(Bs, Ns, P, nb, dt, t, tb);
        /* resultados jogados na saída padrão */
        printf("%s: %lf %lf / ", Ns[0]->nome, Ns[0]->c[0], Ns[0]->c[1]);
        printf("%s: %lf %lf / t=%lf\n", Ns[1]->nome, Ns[1]->c[0], Ns[1]->c[1], t);
    }
}


int main(){
    int i,j,nb;
    double dt = 1  /*stepsize pré-definido */, tb;
    double** DadosN = malloc(sizeof(double*)*2), *DadosP=malloc(sizeof(double*)*3), **DadosB;
    char** Nomes = malloc(sizeof(char*)*2);
    Nave* Ns;
    Planet P;
    Bullet* Bs;
    
    /*input para criação dos objetos: massa, tempo total e raio*/
    for(i = 0; i < 3; i++) {
        scanf("%lf", &DadosP[i]);
    }

    P = CriaPlaneta(DadosP[0], DadosP[1], DadosP[2]);
    
    for(i = 0; i < 2; i++) {
        DadosN[i] = malloc(sizeof(double)*5); /* Massa, pos_x, pos_y, vel_x, vel_y */
        Nomes[i]  = malloc(sizeof(char)*20);  /* Nome das naves */
        scanf("%s", Nomes[i]);
        
        for(j = 0; j < 5; j++){
            scanf("%lf", &DadosN[i][j]);
        }  
    }

    Ns = CriaNaves(DadosN, Nomes);

    scanf("%d", &nb);  /* número de projéteis */
    scanf("%lf", &tb); /* tempo de simulação dos projéteis */

    /* Dados dos projéteis {massa, pos_x, etc} */
    DadosB = malloc(sizeof(double*)*nb); 
    for(i = 0; i < nb ; i++){
        DadosB[i] = malloc(sizeof(double)*5);
        for(j = 0; j < 5; j++){
            scanf("%lf", &DadosB[i][j]);
        }
    }

    Bs = CriaBullets(DadosB, nb);

    /*simulação do movimento*/
    Simulate(Bs, Ns, P, nb, dt, DadosP[2], tb);  

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