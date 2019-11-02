//testes

#include "spacewar.h"

int game (int c){
  if (c == '.')
    return 0;
  return 1;
}

//devolve nova aceleração após aplicar forca na orientação o mod int
double* AplicaForca(double *a, double m, int o, double in)
{
  double *f;
  double alpha = o*(M_PI/(NumIma*2));
  f= malloc(sizeof(double)*2);
  f[1]=sin(alpha)*in;
  f[0]=cos(alpha)*in;
  return Add(a, Multk(f, (double)(1/m)));
}

int main ()
{
  int c=1;
  set_conio_terminal_mode();

  while (c){
    while (!kbhit()) {
        c=game(c);
    }
    c = getch(); /* consume the character */
    c = game (c);
  }
  puts("fim");
  return 0;
}
