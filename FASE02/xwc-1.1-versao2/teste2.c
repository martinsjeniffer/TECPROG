#include <stdio.h>
#include <unistd.h>
#include "xwc.h"
#include <math.h>

#define NumIma 5 //numero de imagens para cada lado rotacao
#define PI 3.14159
#define erro 0.001
int igual (double a, double b)
{
  if (a-b> -erro && a-b<erro)
    return 1;
  else
   return 0;
}

int ImaUpd (int io, int jo, int i, int j, int atual)
{
  /*
  Vetor de imagens com a seguinte organização:
  0... |n  ...  |2n... |3n
  cima |direita |baixo |esquerda
  Função definida para movimentos de 0 à 180 graus
  Recebe posição atual (io,jo), posição futura (i,j), int da pic atual
  Devolve uma int para ser colocada no vetor de imagens
  */
  double alpha; //angulo entre pontos
  int c;
  if ((i-io)*(i-io) + (j-jo)*(j-jo) == 0) // não se mexeu
  {
    return atual;
  }


  if (j>=jo && i>=io) {  //primeiro quadrante
    //calcula alpha
    //printf ("sin = %lf \n", (i-io)*1.0/sqrt((i-io)*(i-io) + (j-jo)*(j-jo)));
    alpha = asin((i-io)*1.0/sqrt((i-io)*(i-io) + (j-jo)*(j-jo)));
    //printf ("alpha = %lf \n", alpha);
    for (c=0; c <= NumIma; c++)
      if (alpha <= c*PI/(2*NumIma) || igual (alpha, c*PI/(2*NumIma)))
        return c;
    puts("Erro rotacao 1");
  }
  else if (j<jo && i>=io) {  //quarto quadrante
    alpha = asin((jo-j)*1.0/sqrt((i-io)*(i-io) + (j-jo)*(j-jo))*1.0);
    for (c=1; c <= NumIma; c++)
      if (alpha < c*PI/(2*NumIma) || igual (alpha, c*PI/(2*NumIma)))
        return NumIma+c;
   puts("Erro rotacao 4");
  }
  else if (j<=jo && i<io) {  //terceiro quadrante
    alpha = asin((io-i)*1.0/sqrt((i-io)*(i-io) + (j-jo)*(j-jo))*1.0);
    for (c=1; c <= NumIma; c++)
      if (alpha <= c*PI/(2*NumIma)  || igual (alpha, c*PI/(2*NumIma)))
        return 2*NumIma+c;
    puts("Erro rotacao 3");
  }
  else {  //segundo quadrante (j>jo e i<=io)
    //puts("entrou 2");
    //alpha = asin((j-jo)*1.0/sqrt((i-io)*(i-io) + (j-jo)*(j-jo))*1.0);
    printf ("alpha = %lf \n", alpha);
    for (c=1; c < NumIma; c++)
      if (alpha <= c*PI/(2*NumIma)  || igual (alpha, c*PI/(2*NumIma)))
        return 3*NumIma+c;
    return 0; //chegou a "cima" pela esquerda
  }
}

int main(int ac, char **av)
{
  PIC P1, P2, P3, Aux;
  PIC p[1+2*NumIma];
#ifdef NOXPM
  puts("Este programa só funciona com a biblioteca Xpm!");
#else
  WINDOW *w1;
  MASK msk;
  int i, io, j, jo, ant, ok;
  
  /*w1 = InitGraph(400,400, "Arquivos");
  P1 = ReadPic(w1, "images/nave.xpm", NULL);
  P2 = ReadPic(w1, "mascara.xpm", NULL);
  //P3 = ReadPic(w1, "Nose.xpm", NULL);
  //P2 = MountPic(w1, nose, NULL);
  
  for (int c = 0; c < NumIma; c++) {
     p[c] = P1;
     p[2*NumIma-c] = P2;
  }
  p[5] = P1;
 */
  //msk = NewMask(*p[1], 32, 32);
  //Aux = ReadPic(w1, "mascara.xpm", msk);
  ok=0;
  ant=0;
  while (!ok) {
   puts("digite io, jo, i, j");
   scanf ("%d %d %d %d", &io, &jo, &i, &j);
   ant = ImaUpd( io, jo, i, j, ant);
   printf("%d \n",  ant);
   puts("parar?");
   scanf ("%d", &ok);
  }
/*  puts("Tecle <enter>"); getchar();

  PutPic(w1, p[1], 0,0, 100, 100, 100, 0);
  puts("Tecle <enter>"); getchar();

  PutPic(w1, p[2], 0,0, 100, 100, 200, 0);
  puts("Tecle <enter>"); getchar();

  //PutPic(w1, P2, 0,0, 100, 100, 100, 0);

  /*puts("Agora  a figura  do arquivo mascara.xpm.");
  puts("Tecle <enter>"); getchar();
  PutPic(w1, Aux, 0,0, 100, 100, 200, 0);*/

  /*puts("Sobrepondo a última figura com a primeira\n"
	   "e usando sua  própria máscara.");
  SetMask(w1,msk);
  puts("Tecle <enter>"); getchar();
  PutPic(w1, Aux, 0,0, 100, 100, 0, 0);
  puts("Tecle <enter>"); getchar();

  puts("Gravando o narigudo em Nose.xpm.");
  WritePic(P2,"Nose.xpm", NULL);
  puts("Gravando samp.xpm mascarado em Tutti.xpm.");
  WritePic(P1,"Tutti.xpm", msk);
  puts("Tecle <enter>"); getchar();
  CloseGraph();*/
#endif
  return 0;
}
