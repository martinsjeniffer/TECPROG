#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "xwc.h"

int main(int ac, char **av)
{
  int ii, jj, vezes = 1;
  PIC P1, P2, P3;
  MASK Msk;
  WINDOW *w1;

  Color cara, olho, nariz, boca;

  if (ac > 1) vezes = atoi(av[1]);
  else puts(
		"Você pode dizer quantas vezes rodar uma animação\n"
		"pela linha de comando, por exemplo\n"
		"    teste2 5");

  w1 = InitGraph(400,400, "Teste animado");
  cara = WNamedColor("pink");
  olho = WNamedColor("blue");
  cara = WNamedColor("#f38080");
  boca = WNamedColor("#f80000");
  nariz = boca;


  /* Inicializa leitura do teclado */
  InitKBD(w1);

  P1 = NewPic(w1, 800, 100);
  /*  P1 = InitGraph( 800, 100, "P1");*/
  WFillRect(P1, 0, 0, 800, 100, 0);

  /* Rascunho */
  P2 = NewPic(w1, 100, 100);
  P3 = NewPic(w1, 400, 100);
  /*P3 = InitGraph(400, 100, "P3");*/

  Msk = NewMask(w1, 100, 100);
  WFillArc(Msk, 0, 0, 0, 360*64, 100, 100, 1);

  /* desenhos */

  WFillArc(P1,   0,  0, 0, 360*64, 100, 100, cara);
  WFillArc(P1,  25, 25, 0, 360*64,  10,  10, olho);
  WFillArc(P1,  75, 25, 0, 360*64,  10,  10, olho);
  WFillArc(P1,  25, 75, 0, 360*64,  50,   1, boca);
  WFillArc(P1,  45, 45, 0, 360*64,  10,  10, nariz);

  WFillArc(P1, 100,  0, 0, 360*64, 100, 100, cara);
  WFillArc(P1, 125, 25, 0, 360*64,  10,   8, olho);
  WFillArc(P1, 175, 25, 0, 360*64,  10,  10, olho);
  WFillArc(P1, 125, 75, 0, 360*64,  50,   5, boca);
  WFillArc(P1, 145, 45, 0, 360*64,  10,  10, nariz);

  WFillArc(P1, 200,  0, 0, 360*64, 100, 100, cara);
  WFillArc(P1, 225, 25, 0, 360*64,  10,   6, olho);
  WFillArc(P1, 275, 25, 0, 360*64,  10,  10, olho);
  WFillArc(P1, 225, 75, 0, 360*64,  50,   7, boca);
  WFillArc(P1, 245, 45, 0, 360*64,  10,  10, nariz);

  WFillArc(P1, 300,  0, 0, 360*64, 100, 100, cara);
  WFillArc(P1, 325, 25, 0, 360*64,  10,   4, olho);
  WFillArc(P1, 375, 25, 0, 360*64,  10,  10, olho);
  WFillArc(P1, 325, 75, 0, 360*64,  50,  10, boca);
  WFillArc(P1, 345, 45, 0, 360*64,  10,  10, nariz);

  WFillArc(P1, 400,  0, 0, 360*64, 100, 100, cara);
  WFillArc(P1, 425, 25, 0, 360*64,  10,   1, olho);
  WFillArc(P1, 475, 25, 0, 360*64,  10,  10, olho);
  WFillArc(P1, 425, 75, 0, 360*64,  50,   7, boca);
  WFillArc(P1, 445, 45, 0, 360*64,  10,  10, nariz);

  WFillArc(P1, 500,  0, 0, 360*64, 100, 100, cara);
  WFillArc(P1, 525, 25, 0, 360*64,  10,   4, olho);
  WFillArc(P1, 575, 25, 0, 360*64,  10,  10, olho);
  WFillArc(P1, 525, 75, 0, 360*64,  50,   4, boca);
  WFillArc(P1, 545, 45, 0, 360*64,  10,  10, nariz);

  WFillArc(P1, 600,  0, 0, 360*64, 100, 100, cara);
  WFillArc(P1, 625, 25, 0, 360*64,  10,   6, olho);
  WFillArc(P1, 675, 25, 0, 360*64,  10,  10, olho);
  WFillArc(P1, 625, 75, 0, 360*64,  50,   7, boca);
  WFillArc(P1, 645, 45, 0, 360*64,  10,  10, nariz);

  WFillArc(P1, 700,  0, 0, 360*64, 100, 100, cara);
  WFillArc(P1, 725, 25, 0, 360*64,  10,   8, olho);
  WFillArc(P1, 775, 25, 0, 360*64,  10,  10, olho);
  WFillArc(P1, 725, 75, 0, 360*64,  50,   1, boca);
  WFillArc(P1, 745, 45, 0, 360*64,  10,  10, nariz);

  /* Este é um exemplo de uma chamada direta ao X. Veja a página do */
  /* manual correspondente e olhe tempo xwc.c e xwc.h			 	*/
  XWriteBitmapFile(GetDisplay(), "caca", GetDraw(P1), 800,100, 0,0);

  for (ii = 0; ii < 400; ii+=4) 
	WLine(w1, 0, ii, 400, ii, (ii/4)*16);

  puts("Pronto, pressione uma tecla na janela");
  printf(">>> %d\n", WGetKey(w1));

  PutPic(P3, w1, 0, 100, 400, 100, 0, 0);

  SetMask(P3,Msk);
  for (jj = 0; jj < vezes; jj++) {
	for (ii = 0; ii < 100; ii++) {
	  PutPic(P2, P3,  50+ii, 0, 100, 100, 0, 0);
	  
	  PutPic(P3, P1, (ii%8)*100,   0, 100, 100, 50+ii, 0);
	  PutPic(w1, P3, 0, 0, 400, 100, 0, 100);
	  
	  PutPic(P3, P2,	0,   0, 100, 100, 50+ii, 0);
	  usleep(25000);
	}
  }

  puts("Pronto, pressione mais uma tecla na janela");
  printf(">>> %d\n", WGetKey(w1));

  PutPic(w1, P1,   0, 0, 400, 100,  0, 100);
  PutPic(w1, P1, 400, 0, 800, 100,  0, 300);

  puts("Fim, aperte <enter> aqui.");
  getchar();
  CloseGraph();
  return 0;
}
