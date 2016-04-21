#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ml6.h"
#include "display.h"
#include "draw.h"
#include "matrix.h"
#include "parser.h"

int main( int argc, char **argv ) {

  screen s;
  struct matrix *tr;
  struct matrix *pm = new_matrix(4,100);
  
  if ( argc == 2 )
    parse_file( argv[1],tr,pm, s );
  else
    parse_file( "stdin",tr,pm, s );
  
  return 0;
  
}  
