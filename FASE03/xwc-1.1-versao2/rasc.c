#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define NumIma 4


int main (){

  int c;
  int v;
  char s[] = "ijklo";
  system ("/bin/stty raw");
  while((c=getchar())!= '.') {
    /* type a period to break out of the loop, since CTRL-D won't work raw*/
    if (v=Valid(c,s))
      printf("%d\n", v);
    putchar(c);
    puts("e?");
  }
  /* use system call to set terminal behaviour to more normal behaviour*/
  system ("/bin/stty cooked");
  return 0;
}
