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
  char e[] = "nave";
  char d[] = "images/";
  // opendir() returns a pointer of DIR type.
  DIR *dr = opendir(strcat(d,e));

  if (dr == NULL)  // opendir returns NULL if couldn't open directory
  {
      printf("Could not open image directory" );
      return 0;
  }
  //vai ler os arquivos
  while ((de = readdir(dr)) != NULL) {
    //strcpy(c, );
    //if (!(strcmp(c,".")||(strcmp(c,"..")))) //não é diretório pai e tal
      printf("%s\n", de->d_name);
  }


  closedir(dr);
  return 0;
}
