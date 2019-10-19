#include <stdio.h>
#include <unistd.h>
#include "xwc.h"


void PrepareMask(WINDOW* w1, PIC P, MASK Msk){
  PutPic(w1, P,  0, 0, 1000, 1000, 0, 0);
  SetMask(P,Msk);
}

PIC DrawShip(WINDOW* w1,PIC P3,PIC P2, PIC P1,MASK Msk,int i,int j){
    UnSetMask(P3);
    PutPic(P3, P2,  0, 0, 1000, 1000, 0, 0);
    SetMask(P3,Msk);
    PutPic(P3, P1,  0, 0, 100, 100, i, j);
    return P3;
}


int main(int ac, char **av)
{
  PIC P,P1, P2, P3, P4, P5, P6, Aux;
#ifdef NOXPM
  puts("Este programa só funciona com a biblioteca Xpm!");
#else
  WINDOW *w1;
  MASK Msk1,Msk2,Msk3;
  int vezes = 2,ii,jj;

  w1 = InitGraph(1000,700, "Arquivos");

  Msk3 = NewMask(w1,400,400);
  Msk2 = NewMask(w1,100,100);
  Msk1 = NewMask(w1,100,100);
  P1 = ReadPic(w1, "images/nave.xpm", Msk1);
  P5 = ReadPic(w1, "images/nave2.xpm", Msk2);
  P2 = ReadPic(w1, "images/sky_pixelart.xpm", NULL);
  P3 = ReadPic(w1, "images/planeta_rosa.xpm", Msk3);
  P4 = NewPic(w1, 1000, 1000);
  P6 = NewPic(w1, 1000, 1000);
  P = NewPic(w1, 1000, 1000);
  
  /*SetMask(P2,Msk2);*/
  
  PrepareMask(w1,P2,Msk3);
  
  PutPic(P2, P3,  0, 0, 400, 400, -100, 500);
  PutPic(P6, P2,  0, 0, 1000, 1000, 0, 0);
  
  PrepareMask(w1,P6,Msk2);
  PrepareMask(w1,P4,Msk1);
  
	for (ii = 0; ii < 1000; ii++) {
    P=DrawShip(w1,P6,P2,P5,Msk2,500-ii,ii);
    P=DrawShip(w1,P4,P,P1,Msk1,ii,ii);
    PutPic(w1, P,  0, 0, 1000, 1000, 0, 0);
	  usleep(20000);
  }
  getchar();
  CloseGraph();
#endif
  return 0;
}



