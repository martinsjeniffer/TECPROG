#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define NumIma 4

int Valid(int c, char s[10]){
  /*
  Checa se c pertence a string s;
  */
  int i;
  for (i=0; i < strlen(s); i++)
    if (c==s[i])
      return i+1;
  return 0;
}


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
  }
  /* use system call to set terminal behaviour to more normal behaviour*/
  system ("/bin/stty cooked");
  return 0;
}
