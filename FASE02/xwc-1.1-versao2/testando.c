/***TESTANDO FUNÇÕES***/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <dirent.h>
#include <string.h>


int main(void)
{
  struct dirent *de;  // Pointer for directory entry
  char c[50];
  // opendir() returns a pointer of DIR type.
  DIR *dr = opendir("images/nave");

  if (dr == NULL)  // opendir returns NULL if couldn't open directory
  {
      printf("Could not open image directory" );
      return 0;
  }
  //vai ler os arquivos
  while ((de = readdir(dr)) != NULL) {
    strcpy(c, de->d_name);
    //if (!(strcmp(c,".")||(strcmp(c,"..")))) //não é diretório pai e tal
      printf("%s\n", c);
  }


  closedir(dr);
  return 0;
}
