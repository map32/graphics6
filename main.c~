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
  struct matrix edges;
  struct matrix transform;

  edges = *new_matrix(4, 4);
  transform = *new_matrix(4, 4);

  struct matrix *a = &edges;
  struct matrix *b = &transform;
  
  if ( argc == 2 )
    parse_file( argv[1], a,b, s );
  else
    parse_file( "stdin", a,b, s );
  print_matrix (b);
  free_matrix(a);
  free_matrix(b);
  
  return 0;
  
}  
