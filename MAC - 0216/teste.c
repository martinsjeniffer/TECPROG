#include <stdio.h>
#include <stdlib.h>

int* f(){
    int *a;
    a=malloc(sizeof(int)*2);
    a[0]=a[1]=2;
    return a;
}

int main(){
    int *b;
    b=f();
    printf("%d",b[1]);
    return 0;
}