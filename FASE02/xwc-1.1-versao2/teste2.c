#include <stdio.h>
#include <unistd.h>
#include "xwc.h"
#include <math.h>

#define NumIma 5 //numero de imagens para cada lado rotacao
#define PI 3.14159

PIC ImaUpd (int io, int jo, int i, int j, PIC imagens[1+2*NumIma], PIC atual)
{
  /*
  Função definida para movimentos de 0 à 180 graus
  Recebe posição atual (io,jo), posição futura (i,j) e vetor de imagens com a imagem na vertical a posiçao do meio e para cada lado apartir do meio, a nave deita um pouco (cada lado com NumIma de imagens), e pic atual
  Devolve uma PIC 
  */
  double alpha; //angulo entre pontos
  int c;
  if ((i-io)*(i-io) + (j-jo)*(j-jo) == 0) // não se mexeu
  {
    return atual;
  }
  //calcula alpha
  alpha = acos((i-io)*1.0/sqrt((i-io)*(i-io) + (j-jo)*(j-jo))*1.0);
  if (alpha < PI/2) //para direita  
  {
    for (c=1; c <= NumIma; c++)
      if (alpha < c*PI/2*NumIma)
        return imagens[2*NumIma - c];
   puts("Deu merda");
  }
  else  //para esquerda
  {
    for (c=1; c <= NumIma; c++)
      if (alpha < PI + c*PI/2*NumIma)
        return imagens[NumIma-c];
   puts("Deu merda2");
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
  
  w1 = InitGraph(400,400, "Arquivos");
  P1 = ReadPic(w1, "images/nave.xpm", NULL);
  P2 = ReadPic(w1, "mascara.xpm", NULL);
  //P3 = ReadPic(w1, "Nose.xpm", NULL);
  //P2 = MountPic(w1, nose, NULL);
  
  for (int c = 0; c < NumIma; c++) {
     p[c] = P1;
     p[2*NumIma-c] = P2;
  }
  p[5] = P1;
 
  //msk = NewMask(*p[1], 32, 32);
  //Aux = ReadPic(w1, "mascara.xpm", msk);

 // puts("Tecle <enter>"); getchar();
  PutPic(w1, ImaUpd( 0, 0, -1, 1, p, P1), 0,0, 100, 100, 0, 0);
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
  WritePic(P1,"Tutti.xpm", msk);*/
  puts("Tecle <enter>"); getchar();
  CloseGraph();
#endif
  return 0;
}
