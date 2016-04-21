#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "ml6.h"
#include "display.h"
#include "draw.h"
#include "matrix.h"

/*======== void add_polygon() ==========
Inputs:   struct matrix *surfaces
         double x0
         double y0
         double z0
         double x1
         double y1
         double z1
         double x2
         double y2
         double z2  
Returns: 
Adds the vertices (x0, y0, z0), (x1, y1, z1)
and (x2, y2, z2) to the polygon matrix. They
define a single triangle surface.

04/16/13 13:05:59
jdyrlandweaver
====================*/
void add_polygon( struct matrix *polygons, 
		  double x0, double y0, double z0, 
		  double x1, double y1, double z1, 
		  double x2, double y2, double z2 ) {
	add_point(polygons,x0,y0,z0);
	add_point(polygons,x1,y1,z1);
	add_point(polygons,x2,y2,z2);
}

/*======== void draw_polygons() ==========
Inputs:   struct matrix *polygons
          screen s
          color c  
Returns: 
Goes through polygons 3 points at a time, drawing 
lines connecting each points to create bounding
triangles

04/16/13 13:13:27
jdyrlandweaver
====================*/
void draw_polygons( struct matrix *polygons, screen s, color c ) {
	int i=0;
	int nx,ny,nz;
	int ax,ay,az;
	int bx,by,bz;
	int dot;
	for(i;i+2<polygons->lastcol;i+=3){
	  ax = polygons->m[0][i+1]-polygons->m[0][i];
	  ay = polygons->m[1][i+1]-polygons->m[1][i];
	  az = polygons->m[2][i+1]-polygons->m[2][i];
	  bx = polygons->m[0][i+2]-polygons->m[0][i];
	  by = polygons->m[1][i+2]-polygons->m[1][i];
	  bz = polygons->m[2][i+2]-polygons->m[2][i];
	  nx = ay*bz-by*az;
	  ny = az*bx-ax*bz;
	  nz = ax*by-bx*ay;
	  dot = nx*0+ny*0+nz*-1;
	  if (dot<0){
	  draw_line(polygons->m[0][i],polygons->m[1][i],polygons->m[0][i+1],polygons->m[1][i+1],s,c);
	  draw_line(polygons->m[0][i+1],polygons->m[1][i+1],polygons->m[0][i+2],polygons->m[1][i+2],s,c);
	  draw_line(polygons->m[0][i+2],polygons->m[1][i+2],polygons->m[0][i],polygons->m[1][i],s,c);
	  }
	}
}

/*======== void add_sphere() ==========
  Inputs:   struct matrix * points
            double cx
	    double cy
	    double r
	    double step  
  Returns: 

  adds all the points for a sphere with center 
  (cx, cy) and radius r.

  should call generate_sphere to create the
  necessary points

  jdyrlandweaver
  ====================*/
void add_sphere( struct matrix * points, 
		 double cx, double cy, double r, 
		 double step ) {
  struct matrix * pts;
  pts = new_matrix(4,100);
  generate_sphere(pts,cx,cy,r,step);
  int i;
  int j;
  int n = (int)(1/step)+1;
  int i2 = n-1;
  int j2 = n-1;
  for(i=0;i<i2;i++){
    int k;
    for(j=0;j<j2;j++){
      k = i*n+j;
      
      add_polygon(points,pts->m[0][k],pts->m[1][k],pts->m[2][k],pts->m[0][k+n+1],pts->m[1][k+n+1],pts->m[2][k+n+1],pts->m[0][k+n],pts->m[1][k+n],pts->m[2][k+n]);
      add_polygon(points,pts->m[0][k],pts->m[1][k],pts->m[2][k],pts->m[0][k+1],pts->m[1][k+1],pts->m[2][k+1],pts->m[0][k+n+1],pts->m[1][k+n+1],pts->m[2][k+n+1]);
    }
  }
  free_matrix(pts);
}

/*======== void generate_sphere() ==========
  Inputs:   struct matrix * points
            double cx
	    double cy
	    double r
	    double step  
  Returns: 

  Generates all the points along the surface of a 
  sphere with center (cx, cy) and radius r

  Adds these points to the matrix parameter

  03/22/12 11:30:26
  jdyrlandweaver
  ====================*/
void generate_sphere( struct matrix * points, 
		      double cx, double cy, double r, 
		      double step ) {
  double a=0;
  double b=0;
  double x,y,z;
  for(a=0;a<1+step;a+=step){
    for(b=0;b<1+step;b+=step){
      x = r*cos(M_PI*a)+cx;
      y = r*sin(M_PI*a)*cos(M_PI*2*b)+cy;
      z = r*sin(M_PI*a)*sin(M_PI*2*b);
      add_point(points,x,y,z);
    }
  }
}    

/*======== void add_torus() ==========
  Inputs:   struct matrix * points
            double cx
	    double cy
	    double r1
	    double r2
	    double step  
  Returns: 

  adds all the points required to make a torus
  with center (cx, cy) and radii r1 and r2.

  should call generate_torus to create the
  necessary points

  03/22/12 13:34:03
  jdyrlandweaver
  ====================*/
void add_torus( struct matrix * points, 
		double cx, double cy, double r1, double r2, 
		double step ) {
  struct matrix * pts;
  pts = new_matrix(4,100);
  generate_torus(pts,cx,cy,r1,r2,step);
  /**int i=0;
  for(i=0;i<pts->lastcol-1;i+=2){
    add_edge(points,pts->m[0][i],pts->m[1][i],pts->m[2][i],
	     pts->m[0][i+1],pts->m[1][i+1],pts->m[2][i+1]);
  }
  free_matrix(pts);
  **/
  int i;
  int j;
  int n = (int)(1/step)+1;
  int i2 = n-1;
  int j2 = n-1;
  for(i=0;i<i2;i++){
    int k;
    for(j=0;j<j2;j++){
      k = i*n+j;
    add_polygon(points,pts->m[0][k],pts->m[1][k],pts->m[2][k],pts->m[0][k+n],pts->m[1][k+n],pts->m[2][k+n],pts->m[0][k+n+1],pts->m[1][k+n+1],pts->m[2][k+n+1]);
    add_polygon(points,pts->m[0][k+1],pts->m[1][k+1],pts->m[2][k+1],pts->m[0][k],pts->m[1][k],pts->m[2][k],pts->m[0][k+n+1],pts->m[1][k+n+1],pts->m[2][k+n+1]);
    }
  }
  free_matrix(pts);
}

/*======== void generate_torus() ==========
  Inputs:   struct matrix * points
            double cx
	    double cy
	    double r
	    double step  
  Returns: 

  Generates all the points along the surface of a 
  tarus with center (cx, cy) and radii r1 and r2

  Adds these points to the matrix parameter

  jdyrlandweaver
  ====================*/
void generate_torus( struct matrix * points, 
		     double cx, double cy, double r1, double r2, 
		     double step ) {
  double a=0;
  double b=0;
  double x,y,z;
  for(a=0;a<1+step;a+=step){
    for(b=0;b<1+step;b+=step){
      x = r1*cos(M_PI*2*a)+cx;
      y = cos(M_PI*2*b)*(r1*sin(M_PI*2*a)+r2)+cy;
      z = sin(M_PI*2*b)*(r1*sin(M_PI*2*a)+r2);
      add_point(points,x,y,z);
    }
  }
}

/*======== void add_box() ==========
  Inputs:   struct matrix * points
            double x
	    double y
	    double z
	    double width
	    double height
	    double depth
  Returns: 

  add the points for a rectagular prism whose 
  upper-left corner is (x, y, z) with width, 
  height and depth dimensions.

  jdyrlandweaver
  ====================*/
void add_box( struct matrix * polygons,
	      double x, double y, double z,
	      double width, double height, double depth ) {
  double x2, y2, z2;
  x2 = x + width;
  y2 = y - height;
  z2 = z - depth;
  //front
  add_polygon( polygons, 
	       x, y, z, 
	       x, y2, z,
	       x2, y2, z);
  add_polygon( polygons, 
	       x2, y2, z, 
	       x2, y, z,
	       x, y, z);
  //back
  add_polygon( polygons, 
	       x2, y, z2, 
	       x2, y2, z2,
	       x, y2, z2);
  add_polygon( polygons, 
	       x, y2, z2, 
	       x, y, z2,
	       x2, y, z2);
  //top
  add_polygon( polygons, 
	       x, y, z2, 
	       x, y, z,
	       x2, y, z);
  add_polygon( polygons, 
	       x2, y, z, 
	       x2, y, z2,
	       x, y, z2);
  //bottom
  add_polygon( polygons, 
	       x2, y2, z2, 
	       x2, y2, z,
	       x, y2, z);
  add_polygon( polygons, 
	       x, y2, z, 
	       x, y2, z2,
	       x2, y2, z2);
  //right side
  add_polygon( polygons, 
	       x2, y, z, 
	       x2, y2, z,
	       x2, y2, z2);
  add_polygon( polygons, 
	       x2, y2, z2, 
	       x2, y, z2,
	       x2, y, z);
  //left side
  add_polygon( polygons, 
	       x, y, z2, 
	       x, y2, z2,
	       x, y2, z);
  add_polygon( polygons, 
	       x, y2, z, 
	       x, y, z,
	       x, y, z2); 
  /**add_polygon(points,x,y,z,x,y-height,z,x+width,y-height,z);
  add_polygon(points,x+width,y-height,z,x+width,y,z,x,y,z);
  add_polygon(points,x+width,y,z-depth,x+width,y-height,z-depth,x,y-height,z-depth);
  add_polygon(points,x,y-height,z-depth,x,y,z-depth,x+width,y,z-depth);
  add_polygon(points,x,y,z-depth,x,y,z,x+width,y,z);
  add_polygon(points,x+width,y,z,x+width,y,z-depth,x,y,z-depth);
  add_polygon(points,x+width,y-height,z-depth,x+width,y-height,z,x,y-height,z);
  add_polygon(points,x,y-height,z,x,y-height,z-depth,x+width,y-height,z-depth);
  add_polygon(points,x,y,z-depth,x,y-height,z-depth,x,y-height,z);
  add_polygon(points,x,y-height,z,x,y,z,x,y,z-depth);
  add_polygon(points,x+width,y,z,x+width,y-height,z,x+width,y-height,z-depth);
  add_polygon(points,x+width,y-height,z-depth,x+width,y,z-depth,x+width,y,z);**/
  
}

/*======== void add_circle() ==========
  Inputs:   struct matrix * points
            double cx
	    double cy
	    double y
	    double step  
  Returns: 


  03/16/12 19:53:52
  jdyrlandweaver
  ====================*/
void add_circle( struct matrix * points, 
		 double cx, double cy, 
		 double r, double step ) {
  double x , y, t;
  double i;
  double x0 = cx;
  double y0 = r+cy;
  for(i=step;i<1;i+=step){
    t = 2*M_PI*i;
    x = r*sin(t)+cx;
    y = r*cos(t)+cy;
    add_edge(points,x0,y0,0,x,y,0);
    //printf("%lf %lf %lf %lf\n",x0,y0,x,y);
    x0 = x;
    y0 = y;
  }
  add_edge(points,x0,y0,0,cx,r+cy,0);
}

/*======== void add_curve() ==========
Inputs:   struct matrix *points
         double x0
         double y0
         double x1
         double y1
         double x2
         double y2
         double x3
         double y3
         double step
         int type  
Returns: 

Adds the curve bounded by the 4 points passsed as parameters
of type specified in type (see matrix.h for curve type constants)
to the matrix points

03/16/12 15:24:25
jdyrlandweaver
====================*/
void add_curve( struct matrix *points, 
		double x0, double y0, 
		double x1, double y1, 
		double x2, double y2, 
		double x3, double y3, 
		double step, int type ) {
  struct matrix *xcfs, *ycfs;
  double x,y,xs,ys;
  xs = x0;
  ys = y0;
  if(type==0){ //hermy
    double s1 = (y1-y0)/(x1-x0);
    double s2 = (y2-y3)/(x2-x3);
    if(x1==x0){
      if(y1>y0){
	s1 = y0-y1;
      } else {
	s1 = y1-y0;
      }
    }
    if(x3==x2){
      if(y3>y2){
	s2 = y3-y2;
      } else {
	s2 = y2-y3;
      }
    }
    xcfs = generate_curve_coefs(x0,x2,s1,s2,0);
    ycfs = generate_curve_coefs(y0,y2,s1,s2,0);
  } else { //bezzy
    xcfs = generate_curve_coefs(x0,x1,x2,x3,1);
    ycfs = generate_curve_coefs(y0,y1,y2,y3,1);
  }
  double i;
  for(i=step;i<1;i+=step){
      x = ((xcfs->m[0][0]*i+xcfs->m[1][0])*i+xcfs->m[2][0])*i+xcfs->m[3][0];
      y = ((ycfs->m[0][0]*i+ycfs->m[1][0])*i+ycfs->m[2][0])*i+ycfs->m[3][0];
      add_edge(points,xs,ys,0,x,y,0);
      xs = x;
      ys = y;
    }
}

/*======== void add_point() ==========
Inputs:   struct matrix * points
         int x
         int y
         int z 
Returns: 
adds point (x, y, z) to points and increment points.lastcol
if points is full, should call grow on points
====================*/
void add_point( struct matrix * points, double x, double y, double z) {
  if (points->lastcol == points->cols){
    grow_matrix(points, (points->cols)*2);
  }
  points->m[0][points->lastcol] = x;
  points->m[1][points->lastcol] = y;
  points->m[2][points->lastcol] = z;
  points->m[3][points->lastcol] = 1;
  points->lastcol++;
  //printf("%lf %lf %lf %lf %d\n",x,y,z,points->lastcol);
}

/*======== void add_edge() ==========
Inputs:   struct matrix * points
          int x0, int y0, int z0, int x1, int y1, int z1
Returns: 
add the line connecting (x0, y0, z0) to (x1, y1, z1) to points
should use add_point
====================*/
void add_edge( struct matrix * points, 
	       double x0, double y0, double z0, 
	       double x1, double y1, double z1) {
  add_point( points, x0, y0, z0 );
  add_point( points, x1, y1, z1 );
}

/*======== void draw_lines() ==========
Inputs:   struct matrix * points
         screen s
         color c 
Returns: 
Go through points 2 at a time and call draw_line to add that line
to the screen
====================*/

void draw_lines( struct matrix * points, screen s, color c) {
  int i=0;
  for(i;i<points->lastcol-1;i+=2){
    draw_line(points->m[0][i],points->m[1][i],points->m[0][i+1],points->m[1][i+1],s,c);
  }
}	       


void draw2(int always0, int conditional0, int always1, int flip, int dalways, int dconditional, int talways, int tconditional, screen s, color c){
  int D = talways+tconditional/2;
  //printf("%d %d %d",talways,tconditional,D);
  int always = always0;
  int conditional = conditional0;
  if(!flip){
    while(always!=always1){
      plot(s,c,always,conditional);
      if (D>0){
	conditional += dconditional;
	D += tconditional;
      }
      always += dalways;
      D += talways;
    }
    plot(s,c,always,conditional);
  } else {
    while(always!=always1){
      plot(s,c,conditional,always);
      if (D<0){
	conditional += dconditional;
	D += tconditional;
      }
      always += dalways;
      D += talways;
    }
    plot(s,c,conditional,always);
  }
}

//Insert your line algorithm here
void draw_line(int x0, int y0, int x1, int y1, screen s, color c) {
  int dx,dy;
  int A = 2*(y1-y0);
  int B = 2*(-(x1-x0));
  //printf("d\n");
  if (A<0){ //5-8 octants flipped into 1-4
    int t = x1;
    x1 = x0;
    x0 = t;
    t = y1;
    y1 = y0;
    y0 = t;
    //printf("%d %d %d %d\n",x0,x1,y0,y1);
    B = -B;
    A = -A;
  }
  if (A>B && A>-B){ // 2,3
    if (B<0){ // 2
      //printf("2");
      //draw2(x0,y0,x1,0,1,1,A,B,s,c);
      draw2(y0,x0,y1,1,1,1,B,A,s,c);
    } else { // 3
      //printf("3");
      draw2(y0,x0,y1,1,1,-1,-B,A,s,c);
    }
  } else { //1,4
    if (B<0){ // 1
      //printf("1\n");
      draw2(x0,y0,x1,0,1,1,A,B,s,c);
    } else { // 4
      //printf("4\n");
      draw2(x0,y0,x1,0,-1,1,A,-B,s,c);
    }
  }
}


