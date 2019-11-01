#include <stdio.h>
#include <unistd.h>
#include "xwc.h"
#include <math.h>
#include <dirent.h>

#define NumIma 4 //numero de imagens para cada lado rotacao
#define erro 0.0001

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
   //ant = ImaUpd( io, jo, i, j, ant);
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
