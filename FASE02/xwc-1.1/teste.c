#include <stdio.h>
#include <math.h>
#include "xwc.h"

int main(int ac, char **av)
{
  int ii;
  Color c;
  double a;
  char cc[256];
  WINDOW *w;

  puts("Iniciando o sistema gráfico e abrindo uma janela");
  w = InitGraph(400,400, "Janelao");
  

  /* puts("Aperte <enter>"); */
  /* getchar(); */

  for (ii = 0; ii < 400; ii++) {
	Color cor = 256*256*255+ 256*(ii%256) + (ii%175);
	WPlot(w,ii,200, cor);
  }
  WLine(w,1,1,300,300,0xffd700);

  for (a = 0.; a < 1.75; a+=.1)
    WLine(w,0,0,100, 100*tan(a), 0x000495);

  WCor(w,c=WNamedColor("gold"));
  WPrint(w,10,10,"So isso");

  c = WNamedColor("Red");
  puts("Agora,  escreva vários nomes de cores,  ou uma cor no  formato\n"
	   "'#RRGGBB',  onde RR,  GG e BB correspondem às intensidades, em\n"
	   "hexadecimal,  das componentes vermelha,  verde e azul, respec.\n"
	   "Os nomes válidos de cores estão em /usr/X11R6/lib/X11/rgb.txt.\n"
	   "Termine a lista de cores com um 'q'.");

  for(;;) {
    int nada = scanf("%s", cc);
    if ((cc[0]=='q') && (cc[1]=='\x0')) break;
    c=WNamedColor(cc);
    WFillRect(w,20,20, 80, 230, c);
  }
  CloseGraph();
  return 0;
}


