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
  struct matrix *a;
  struct matrix *b;
  
  if ( argc == 2 )
    parse_file( argv[1],a,b, s );
  else
    parse_file( "stdin",a,b, s );

  free_matrix(a);
  free_matrix(b);
  
  return 0;
  
}  
