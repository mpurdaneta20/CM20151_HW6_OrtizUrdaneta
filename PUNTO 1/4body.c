#include <stdio.h>
#include <stdlib.h>
#include <math.h>
float * get_memory(int n_points);

/*Declaracion de las funciones que seran utilizadas en el metodo de RungeKutta4*/
float xDot_sol(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vx_s);
float vxDot_sol(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vx_s);
float yDot_sol(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vy_s);
float vyDot_sol(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vy_s);
float zDot_sol(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vz_s);
float vzDot_sol(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vz_s);

float xDot_tierra(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vx_t);
float vxDot_tierra(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vx_t);
float yDot_tierra(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vy_t);
float vyDot_tierra(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vy_t);
float zDot_tierra(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vz_t);
float vzDot_tierra(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vz_t);

float xDot_luna(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vx_l);
float vxDot_luna(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vx_l);
float yDot_luna(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vy_l);
float vyDot_luna(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vy_l);
float zDot_luna(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vz_l);
float vzDot_luna(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vz_l);

float xDot_ast(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vx_a);
float vxDot_ast(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vx_a);
float yDot_ast(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vy_a);
float vyDot_ast(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vy_a);
float zDot_ast(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vz_a);
float vzDot_ast(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vz_a);

// void initialize_pos(float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za);
// void initialize_vel(float *vx_s, float *vy_s, float *vz_s, float *vx_t, float *vy_t, float *vz_t, float *vx_l, float *vy_l, float *vz_l, float *vx_a, float *vy_a, float *vz_a);

float RungeKutta(int d, float h);
float *IC;
/*Posicion de todas las particulas*/
float *xs;
float *xl;
float *xt;
float *xa;
float *ys;
float *yl;
float *yt;
float *ya;
float *zs;
float *zl;
float *zt;
float *za;
/*Velocidad de todas las particulas*/
float *vx_s;
float *vy_s;
float *vz_s;
float *vx_l;
float *vy_l;
float *vz_l;
float *vx_t;
float *vy_t;
float *vz_t;
float *vx_a;
float *vy_a;
float *vz_a;
/*Aceleracion de todas las particulas*/
float *a_xs;
float *a_ys;
float *a_zs;
float *a_xt;
float *a_yt;
float *a_zt;
float *a_xl;
float *a_yl;
float *a_zl;
float *a_xa;
float *a_ya;
float *a_za;
/*Masa de todas las particulas*/
float ms;
float ml;
float ma;
float mt;
/*Constante gravitacional*/
float G;

int main(int argc, char **argv)
{
  if(argc != 4)
  {
    printf("El número de parámetros es diferente de 3");
    exit(1);
  }

  float h;
  float t;
  FILE *ifp;
  FILE *ofp;
  char *mode = "r";
  char *mode2= "w";
  
  /*Abrir el archivo ic.txt el cual contiene todas las condiciones iniciales del modelo*/
  ifp = fopen(argv[1], mode);
  h = atof(argv[2]);
  t = atoi(argv[3]); 

  if (ifp == NULL){
    fprintf(stderr, "El archivo no esta, el archivo se fue!\n");
    exit(1);
  }
 
  /*Crear .txt que guardara las posiciones finales en x,y de cada astro*/
  char na[20];
  sprintf(na, "orbitas-%dyrs.txt",atoi(argv[3]));
  ofp = fopen(na, mode2);

  if (ofp == NULL)
  {
    fprintf(stderr, "El archivo no esta, el archivo se fue!\n");
    exit(1);
  }
  /*Aumenta la memoria de los punteros*/
  int n_points=100000;
  float ictxt=29;
  IC = get_memory(ictxt);
  xs = get_memory(n_points);
  ys = get_memory(n_points);
  zs = get_memory(n_points);
  xl = get_memory(n_points);
  yl = get_memory(n_points);
  zl = get_memory(n_points);
  xt = get_memory(n_points);
  yt = get_memory(n_points);
  zt = get_memory(n_points);
  xa = get_memory(n_points);
  ya = get_memory(n_points);
  za = get_memory(n_points);

  a_xs = get_memory(n_points);
  a_ys = get_memory(n_points);
  a_zs = get_memory(n_points);
  a_xt = get_memory(n_points);
  a_yt = get_memory(n_points);
  a_zt = get_memory(n_points);
  a_xl = get_memory(n_points);
  a_yl = get_memory(n_points);
  a_zl = get_memory(n_points);
  a_xa = get_memory(n_points);
  a_ya = get_memory(n_points);
  a_za = get_memory(n_points);

  vx_s = get_memory(n_points);
  vy_s = get_memory(n_points);
  vz_s = get_memory(n_points);
  vx_l = get_memory(n_points);
  vy_l = get_memory(n_points);
  vz_l = get_memory(n_points);
  vx_t = get_memory(n_points);
  vy_t = get_memory(n_points);
  vz_t = get_memory(n_points);
  vx_a = get_memory(n_points);
  vy_a = get_memory(n_points);
  vz_a = get_memory(n_points);

  G= 39.42; 

  char line[20];

      int i;
      for(i=0; i<28; i++){
        fgets(line, sizeof(line), ifp);
        IC[i]=atof(line);
      } 

 /*Inicializacion de datos*/
   xs[0] = IC[0];
   xt[0] = IC[7];
   xl[0] = IC[14];
   xa[0] = IC[21];

   ys[0] = IC[1];
   yt[0] = IC[8];
   yl[0] = IC[15];
   ya[0] = IC[22];

   zs[0] = IC[2];
   zt[0] = IC[9];
   zl[0] = IC[16];
   za[0] = IC[23];

   ms = IC[6];
   mt = IC[13];
   ml = IC[20];
   ma = IC[27];

   vx_s[0]=IC[3];
   vy_s[0]=IC[4];
   vz_s[0]=IC[5];

   vx_t[0] = IC[10];
   vy_t[0] = IC[11]; 
   vz_t[0] = IC[12]; 

   vx_l[0] = IC[17];
   vy_l[0] = IC[18];
   vz_l[0] = IC[19];

   vx_a[0] = IC[24]; 
   vy_a[0] = IC[25]; 
   vz_a[0] = IC[26];

  int d;
  float N;
  N=(int)t/h;
  
  fprintf(ofp,"xs ys xt yt xl yl xa ya\n");
  fprintf(ofp,"%f %f %f %f %f %f %f %f\n", xs[0], ys[0], xt[0], yt[0], xl[0], yl[0], xa[0], ya[0]);

  for(d=1;d<N;d++)
  {
    RungeKutta(d,h);
    fprintf(ofp,"%f %f %f %f %f %f %f %f\n", xs[d], ys[d], xt[d], yt[d], xl[d], yl[d], xa[d], ya[d]);
  }

  fclose(ifp);
  fclose(ofp);
}

/*Funciones xDot y vDot del sol*/
float xDot_sol(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vx_s){
  return vx_s[i];
}

float yDot_sol(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vy_s){
  return vy_s[i];
}
float zDot_sol(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vz_s){
  return vz_s[i];
}

float vxDot_sol(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vx_s)
{
	float normt;
	float norml;
	float norma;
	float dirt;
	float dirl;
	float dira;
	normt=pow(sqrt(pow(((xs[i]) - (xt[i])),2.0) + pow((ys[i]) - (yt[i]),2.0) + pow((zs[i]) - (zt[i]),2.0)), 3.0);
	norml=pow(sqrt(pow(((xs[i]) - (xl[i])),2.0) + pow((ys[i]) - (yl[i]),2.0) + pow((zs[i]) - (zl[i]),2.0)), 3.0);
	norma=pow(sqrt(pow(((xs[i]) - (xa[i])),2.0) + pow((ys[i]) - (ya[i]),2.0) + pow((zs[i]) - (za[i]),2.0)), 3.0);
	dirt=(xs[i]-xt[i]);
	dirl=(xs[i]-xl[i]);
	dira=(xs[i]-xa[i]);
	a_xs[i]=-G*(((mt*dirt)/normt)+((ml*dirl)/norml)+((ma*dira)/norma));
	
	return a_xs[i];
}

float vyDot_sol(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vy_s)
{
	float normt_y;
	float norml_y;
	float norma_y;
	float dirt_y;
	float dirl_y;
	float dira_y;
	normt_y=pow(sqrt(pow(((xs[i]) - (xt[i])),2.0) + pow((ys[i]) - (yt[i]),2.0) + pow((zs[i]) - (zt[i]),2.0)), 3.0);
	norml_y=pow(sqrt(pow(((xs[i]) - (xl[i])),2.0) + pow((ys[i]) - (yl[i]),2.0) + pow((zs[i]) - (zl[i]),2.0)), 3.0);
	norma_y=pow(sqrt(pow(((xs[i]) - (xa[i])),2.0) + pow((ys[i]) - (ya[i]),2.0) + pow((zs[i]) - (za[i]),2.0)), 3.0);

	dirt_y=(ys[i]-yt[i]);
	dirl_y=(ys[i]-yl[i]);
	dira_y=(ys[i]-ya[i]);

	a_ys[i]=-G*(((mt*dirt_y)/normt_y)+((ml*dirl_y)/norml_y)+((ma*dira_y)/norma_y));
	return a_ys[i];
}

float vzDot_sol(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vz_s)
{
	float normt;
	float norml;
	float norma;
	float dirt;
	float dirl;
	float dira;

	normt=pow(sqrt(pow(((xs[i]) - (xt[i])),2.0) + pow((ys[i]) - (yt[i]),2.0) + pow((zs[i]) - (zt[i]),2.0)), 3.0);
	norml=pow(sqrt(pow(((xs[i]) - (xl[i])),2.0) + pow((ys[i]) - (yl[i]),2.0) + pow((zs[i]) - (zl[i]),2.0)), 3.0);
	norma=pow(sqrt(pow(((xs[i]) - (xa[i])),2.0) + pow((ys[i]) - (ya[i]),2.0) + pow((zs[i]) - (za[i]),2.0)), 3.0);
	dirt=(zs[i]-zt[i]);
	dirl=(zs[i]-zl[i]);
	dira=(zs[i]-za[i]);

	a_zs[i]=-G*(((mt*dirt)/normt)+((ml*dirl)/norml)+((ma*dira)/norma));
	return a_zs[i];
}

/*Funciones xDot y vDot de la tierra*/
float xDot_tierra(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vx_t){
  return vx_t[i];
}

float yDot_tierra(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vy_t){
  return vy_t[i];
}
float zDot_tierra(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vz_t){
  return vz_t[i];
}

float vxDot_tierra(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vx_t){

	float norms;
	float norml;
	float norma;
	float dirs;
	float dirl;
	float dira;

	norms=pow(sqrt(pow(((xt[i]) - (xs[i])),2.0) + pow((yt[i]) - (ys[i]),2.0) + pow((zt[i]) - (zs[i]),2.0)), 3.0);
	norml=pow(sqrt(pow(((xt[i]) - (xl[i])),2.0) + pow((yt[i]) - (yl[i]),2.0) + pow((zt[i]) - (zl[i]),2.0)), 3.0);
	norma=pow(sqrt(pow(((xt[i]) - (xa[i])),2.0) + pow((yt[i]) - (ya[i]),2.0) + pow((zt[i]) - (za[i]),2.0)), 3.0);

	dirs=(xt[i]-xs[i]);
	dirl=(xt[i]-xl[i]);
	dira=(xt[i]-xa[i]);

	a_xt[i]=-G*(((ms*dirs)/norms)+((ml*dirl)/norml)+((ma*dira)/norma));
	return a_xt[i];
	//printf("vxDot_tierra! %d %f \n", i, a_xt[i]);
}

float vyDot_tierra(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vx_t){

	float norms;
	float norml;
	float norma;
	float dirs;
	float dirl;
	float dira;

	norms=pow(sqrt(pow(((xt[i]) - (xs[i])),2.0) + pow((yt[i]) - (ys[i]),2.0) + pow((zt[i]) - (zs[i]),2.0)), 3.0);
	norml=pow(sqrt(pow(((xt[i]) - (xl[i])),2.0) + pow((yt[i]) - (yl[i]),2.0) + pow((zt[i]) - (zl[i]),2.0)), 3.0);
	norma=pow(sqrt(pow(((xt[i]) - (xa[i])),2.0) + pow((yt[i]) - (ya[i]),2.0) + pow((zt[i]) - (za[i]),2.0)), 3.0);

	dirs=(yt[i]-ys[i]);
	dirl=(yt[i]-yl[i]);
	dira=(yt[i]-ya[i]);

	a_yt[i]=-G*(((ms*dirs)/norms)+((ml*dirl)/norml)+((ma*dira)/norma));
	return a_yt[i];
}

float vzDot_tierra(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vx_t){

	float norms;
	float norml;
	float norma;
	float dirs;
	float dirl;
	float dira;

	norms=pow(sqrt(pow(((xt[i]) - (xs[i])),2.0) + pow((yt[i]) - (ys[i]),2.0) + pow((zt[i]) - (zs[i]),2.0)), 3.0);
	norml=pow(sqrt(pow(((xt[i]) - (xl[i])),2.0) + pow((yt[i]) - (yl[i]),2.0) + pow((zt[i]) - (zl[i]),2.0)), 3.0);
	norma=pow(sqrt(pow(((xt[i]) - (xa[i])),2.0) + pow((yt[i]) - (ya[i]),2.0) + pow((zt[i]) - (za[i]),2.0)), 3.0);

	dirs=(zt[i]-zs[i]);
	dirl=(zt[i]-zl[i]);
	dira=(zt[i]-za[i]);

	a_zt[i]=-G*(((ms*dirs)/norms)+((ml*dirl)/norml)+((ma*dira)/norma));
	return a_zt[i];
}

/*Funciones xDot y vDot de la luna*/
float xDot_luna(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vx_l){
  return vx_l[i];
}

float yDot_luna(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vy_l){
  return vy_l[i];
}
float zDot_luna(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vz_l){
  return vz_l[i];
}


float vxDot_luna(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vx_l){

	float norms;
	float normt;
	float norma;
	float dirs;
	float dirt;
	float dira;

	norms=pow(sqrt(pow(((xl[i]) - (xs[i])),2.0) + pow((yl[i]) - (ys[i]),2.0) + pow((zl[i]) - (zs[0]),2.0)), 3.0);
	normt=pow(sqrt(pow(((xl[i]) - (xt[i])),2.0) + pow((yl[i]) - (yt[i]),2.0) + pow((zl[i]) - (zt[0]),2.0)), 3.0);
	norma=pow(sqrt(pow(((xl[i]) - (xa[i])),2.0) + pow((yl[i]) - (ya[i]),2.0) + pow((zl[i]) - (za[0]),2.0)), 3.0);

	dirs=(xl[i]-xs[i]);
	dirt=(xl[i]-xt[i]);
	dira=(xl[i]-xa[i]);

	a_xl[i]=-G*(((ms*dirs)/norms)+((mt*dirt)/normt)+((ma*dira)/norma));
	return a_xl[i];
}

float vyDot_luna(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vy_l){

	float norms;
	float normt;
	float norma;
	float dirs;
	float dirt;
	float dira;

	norms=pow(sqrt(pow(((xl[i]) - (xs[i])),2.0) + pow((yl[i]) - (ys[i]),2.0) + pow((zl[i]) - (zs[i]),2.0)), 3.0);
	normt=pow(sqrt(pow(((xl[i]) - (xt[i])),2.0) + pow((yl[i]) - (yt[i]),2.0) + pow((zl[i]) - (zt[i]),2.0)), 3.0);
	norma=pow(sqrt(pow(((xl[i]) - (xa[i])),2.0) + pow((yl[i]) - (ya[i]),2.0) + pow((zl[i]) - (za[i]),2.0)), 3.0);

	dirs=(yl[i]-ys[i]);
	dirt=(yl[i]-yt[i]);
	dira=(yl[i]-ya[i]);

	a_yl[i]=-G*(((ms*dirs)/norms)+((mt*dirt)/normt)+((ma*dira)/norma));
	return a_yl[i];
}
float vzDot_luna(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vz_l){

	float norms;
	float normt;
	float norma;
	float dirs;
	float dirt;
	float dira;

	norms=pow(sqrt(pow(((xl[i]) - (xs[i])),2.0) + pow((yl[i]) - (ys[i]),2.0) + pow((zl[i]) - (zs[i]),2.0)), 3.0);
	normt=pow(sqrt(pow(((xl[i]) - (xt[i])),2.0) + pow((yl[i]) - (yt[i]),2.0) + pow((zl[i]) - (zt[i]),2.0)), 3.0);
	norma=pow(sqrt(pow(((xl[i]) - (xa[i])),2.0) + pow((yl[i]) - (ya[i]),2.0) + pow((zl[i]) - (za[i]),2.0)), 3.0);

	dirs=(zl[i]-zs[i]);
	dirt=(zl[i]-zt[i]);
	dira=(zl[i]-za[i]);

	a_zl[i]=-G*(((ms*dirs)/norms)+((mt*dirt)/normt)+((ma*dira)/norma));
	return a_zl[i];
}

float xDot_ast(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vx_a){
  return vx_a[i];
}

float yDot_ast(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vy_a){
  return vy_a[i];
}
float zDot_ast(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vz_a){
    return vz_a[i];
}

float vxDot_ast(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vx_a){

	float norms;
	float normt;
	float norml;
	float dirs;
	float dirt;
	float dirl;

	norms=pow(sqrt(pow(((xa[i]) - (xs[i])),2.0) + pow((ya[i]) - (ys[i]),2.0) + pow((za[i]) - (zs[i]),2.0)), 3.0);
	normt=pow(sqrt(pow(((xa[i]) - (xt[i])),2.0) + pow((ya[i]) - (yt[i]),2.0) + pow((za[i]) - (zt[i]),2.0)), 3.0);
	norml=pow(sqrt(pow(((xa[i]) - (xl[i])),2.0) + pow((ya[i]) - (yl[i]),2.0) + pow((za[i]) - (zl[i]),2.0)), 3.0);

	dirs=(xa[i]-xs[i]);
	dirt=(xa[i]-xt[i]);
	dirl=(xa[i]-xl[i]);

	a_xa[i]=-G*(((ms*dirs)/norms)+((mt*dirt)/normt)+((ma*dirl)/norml));
	return a_xa[i];
}

float vyDot_ast(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vy_a){

	float norms;
	float normt;
	float norml;
	float dirs;
	float dirt;
	float dirl;

	norms=pow(sqrt(pow(((xa[i]) - (xs[i])),2.0) + pow((ya[i]) - (ys[i]),2.0) + pow((za[i]) - (zs[i]),2.0)), 3.0);
	normt=pow(sqrt(pow(((xa[i]) - (xt[i])),2.0) + pow((ya[i]) - (yt[i]),2.0) + pow((za[i]) - (zt[i]),2.0)), 3.0);
	norml=pow(sqrt(pow(((xa[i]) - (xl[i])),2.0) + pow((ya[i]) - (yl[i]),2.0) + pow((za[i]) - (zl[i]),2.0)), 3.0);

	dirs=(ya[i]-ys[i]);
	dirt=(ya[i]-yt[i]);
	dirl=(ya[i]-yl[i]);

	a_ya[i]=-G*(((ms*dirs)/norms)+((mt*dirt)/normt)+((ma*dirl)/norml));
	return a_ya[i];
}

float vzDot_ast(int i, float *xs, float *ys, float *zs, float *xt, float *yt, float *zt, float *xl, float *yl, float *zl, float *xa, float *ya, float *za, float *vz_a){

	float norms;
	float normt;
	float norml;
	float dirs;
	float dirt;
	float dirl;

	norms=pow(sqrt(pow(((xa[i]) - (xs[i])),2.0) + pow((ya[i]) - (ys[i]),2.0) + pow((za[i]) - (zs[i]),2.0)), 3.0);
	normt=pow(sqrt(pow(((xa[i]) - (xt[i])),2.0) + pow((ya[i]) - (yt[i]),2.0) + pow((za[i]) - (zt[i]),2.0)), 3.0);
	norml=pow(sqrt(pow(((xa[i]) - (xl[i])),2.0) + pow((ya[i]) - (yl[i]),2.0) + pow((za[i]) - (zl[i]),2.0)), 3.0);

	dirs=(za[i]-zs[i]);
	dirt=(za[i]-zt[i]);
	dirl=(za[i]-zl[i]);

	a_za[i]=-G*(((ms*dirs)/norms)+((mt*dirt)/normt)+((ma*dirl)/norml));
	return a_za[i];
}

float RungeKutta(int d,float h)
{
  int i;
  i=d-1;

  float saltoxs, saltovxs, saltoys, saltovys, saltozs, saltovzs;
  float saltoxt, saltovxt, saltoyt, saltovyt, saltozt, saltovzt;
  float saltoxl, saltovxl, saltoyl, saltovyl, saltozl, saltovzl;
  float saltoxa, saltovxa, saltoya, saltovya, saltoza, saltovza;

  saltoxs=xs[i];
  saltovxs=vx_s[i];
  saltoys=ys[i];
  saltovys=vy_s[i];
  saltozs=zs[i];
  saltovzs=vz_s[i];

  saltoxt=xt[i];
  saltovxt=vx_t[i];
  saltoyt=yt[i];
  saltovyt=vy_t[i];
  saltozt=zt[i];
  saltovzt=vz_t[i];

  saltoxl=xl[i];
  saltovxl=vx_l[i];
  saltoyl=yl[i];
  saltovyl=vy_l[i];
  saltozl=zl[i];
  saltovzl=vz_l[i];

  saltoxa=xa[i];
  saltovxa=vx_a[i];
  saltoya=ya[i];
  saltovya=vy_a[i];
  saltoza=za[i];
  saltovza=vz_a[i];

  float K1xsol, K1ysol, K1zsol, K1xtierra, K1ytierra, K1ztierra, K1xluna, K1yluna, K1zluna, K1xast, K1yast, K1zast;
  float K1vxsol, K1vysol, K1vzsol, K1vxtierra, K1vytierra, K1vztierra, K1vxluna, K1vyluna, K1vzluna, K1vxast, K1vyast, K1vzast;

  float K2xsol, K2ysol, K2zsol, K2xtierra, K2ytierra, K2ztierra, K2xluna, K2yluna, K2zluna, K2xast, K2yast, K2zast;
  float K2vxsol, K2vysol, K2vzsol, K2vxtierra, K2vytierra, K2vztierra, K2vxluna, K2vyluna, K2vzluna, K2vxast, K2vyast, K2vzast;

  float K3xsol, K3ysol, K3zsol, K3xtierra, K3ytierra, K3ztierra, K3xluna, K3yluna, K3zluna, K3xast, K3yast, K3zast;
  float K3vxsol, K3vysol, K3vzsol, K3vxtierra, K3vytierra, K3vztierra, K3vxluna, K3vyluna, K3vzluna, K3vxast, K3vyast, K3vzast;

  float K4xsol, K4ysol, K4zsol, K4xtierra, K4ytierra, K4ztierra, K4xluna, K4yluna, K4zluna, K4xast, K4yast, K4zast;
  float K4vxsol, K4vysol, K4vzsol, K4vxtierra, K4vytierra, K4vztierra, K4vxluna, K4vyluna, K4vzluna, K4vxast, K4vyast, K4vzast;

	//K1
  K1xsol= xDot_sol(i,xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za,vx_s);
  K1ysol= yDot_sol(i,xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za,vy_s);
  K1zsol= zDot_sol(i,xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za,vz_s);

  K1xtierra= xDot_tierra(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_t);
  K1ytierra= yDot_tierra(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vy_t);
  K1ztierra= zDot_tierra(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vz_t);

  K1xluna= xDot_luna(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_l);
  K1yluna= yDot_luna(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vy_l);
  K1zluna= zDot_luna(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vz_l);

  K1xast= xDot_ast(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_a);
  K1yast= yDot_ast(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vy_a);
  K1zast= zDot_ast(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vz_a);
  
  K1vxsol= vxDot_sol(i,xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za,vx_s);
  K1vysol= vyDot_sol(i,xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za,vy_s);
  K1vzsol= vzDot_sol(i,xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za,vz_s);

  K1vxtierra= vxDot_tierra(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_t);
  K1vytierra= vyDot_tierra(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_t);
  K1vztierra= vzDot_tierra(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_t);

  K1vxluna= vxDot_luna(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_l);
  K1vyluna= vyDot_luna(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vy_l);
  K1vzluna= vzDot_luna(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vz_l);

  K1vxast= vxDot_ast(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_a);
  K1vyast= vyDot_ast(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vy_a);
  K1vzast= vzDot_ast(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vz_a);

  xs[i]=saltoxs+(K1xsol*h/2);
  ys[i]=saltoys+(K1ysol*h/2);
  zs[i]=saltozs+(K1zsol*h/2);
  xt[i]=saltoxt+(K1xtierra*h/2);
  yt[i]=saltoyt+(K1ytierra*h/2);
  zt[i]=saltozt+(K1ztierra*h/2);
  xl[i]=saltoxl+(K1xluna*h/2);
  yl[i]=saltoyl+(K1yluna*h/2);
  zl[i]=saltozl+(K1zluna*h/2);
  xa[i]=saltoxa+(K1xast*h/2);
  ya[i]=saltoya+(K1yast*h/2);
  za[i]=saltoza+(K1zast*h/2);

  vx_s[i]=saltovxs+(K1vxsol*h/2);
  vy_s[i]=saltovys+(K1vysol*h/2);
  vz_s[i]=saltovzs+(K1vzsol*h/2);
  vx_t[i]=saltovxt+(K1vxtierra*h/2);
  vy_t[i]=saltovyt+(K1vytierra*h/2);
  vz_t[i]=saltovzt+(K1vztierra*h/2);
  vx_l[i]=saltovxl+(K1vxluna*h/2);
  vy_l[i]=saltovyl+(K1vyluna*h/2);
  vz_l[i]=saltovzl+(K1vzluna*h/2);
  vx_a[i]=saltovxa+(K1vxast*h/2);
  vy_a[i]=saltovya+(K1vyast*h/2);
  vz_a[i]=saltovza+(K1vzast*h/2);

  //K2
  K2xsol= xDot_sol(i,xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za,vx_s);
  K2ysol= yDot_sol(i,xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za,vy_s);
  K2zsol= zDot_sol(i,xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za,vz_s);

  K2xtierra= xDot_tierra(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_t);
  K2ytierra= yDot_tierra(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vy_t);
  K2ztierra= zDot_tierra(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vz_t);

  K2xluna= xDot_luna(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_l);
  K2yluna= yDot_luna(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vy_l);
  K2zluna= zDot_luna(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vz_l);

  K2xast= xDot_ast(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_a);
  K2yast= yDot_ast(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vy_a);
  K2zast= zDot_ast(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vz_a);
  
  K2vxsol= vxDot_sol(i,xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za,vx_s);
  K2vysol= vyDot_sol(i,xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za,vy_s);
  K2vzsol= vzDot_sol(i,xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za,vz_s);

  K2vxtierra= vxDot_tierra(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_t);
  K2vytierra= vyDot_tierra(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_t);
  K2vztierra= vzDot_tierra(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_t);

  K2vxluna= vxDot_luna(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_l);
  K2vyluna= vyDot_luna(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vy_l);
  K2vzluna= vzDot_luna(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vz_l);

  K2vxast= vxDot_ast(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_a);
  K2vyast= vyDot_ast(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vy_a);
  K2vzast= vzDot_ast(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vy_a);

  xs[i]=saltoxs+(K2xsol*h/2);
  ys[i]=saltoys+(K2ysol*h/2);
  zs[i]=saltozs+(K2zsol*h/2);
  xt[i]=saltoxt+(K2xtierra*h/2);
  yt[i]=saltoyt+(K2ytierra*h/2);
  zt[i]=saltozt+(K2ztierra*h/2);
  xl[i]=saltoxl+(K2xluna*h/2);
  yl[i]=saltoyl+(K2yluna*h/2);
  zl[i]=saltozl+(K2zluna*h/2);
  xa[i]=saltoxa+(K2xast*h/2);
  ya[i]=saltoya+(K2yast*h/2);
  za[i]=saltoza+(K2zast*h/2);

  vx_s[i]=saltovxs+(K2vxsol*h/2);
  vy_s[i]=saltovys+(K2vysol*h/2);
  vz_s[i]=saltovzs+(K2vzsol*h/2);
  vx_t[i]=saltovxt+(K2vxtierra*h/2);
  vy_t[i]=saltovyt+(K2vytierra*h/2);
  vz_t[i]=saltovzt+(K2vztierra*h/2);
  vx_l[i]=saltovxl+(K2vxluna*h/2);
  vy_l[i]=saltovyl+(K2vyluna*h/2);
  vz_l[i]=saltovzl+(K2vzluna*h/2);
  vx_a[i]=saltovxa+(K2vxast*h/2);
  vy_a[i]=saltovya+(K2vyast*h/2);
  vz_a[i]=saltovza+(K2vzast*h/2);

  //K3
  K3xsol= xDot_sol(i,xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za,vx_s);
  K3ysol= yDot_sol(i,xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za,vy_s);
  K3zsol= zDot_sol(i,xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za,vz_s);

  K3xtierra= xDot_tierra(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_t);
  K3ytierra= yDot_tierra(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vy_t);
  K3ztierra= zDot_tierra(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vz_t);

  K3xluna= xDot_luna(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_l);
  K3yluna= yDot_luna(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vy_l);
  K3zluna= zDot_luna(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vz_l);

  K3xast= xDot_ast(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_a);
  K3yast= yDot_ast(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vy_a);
  K3zast= zDot_ast(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vz_a);
  
  K3vxsol= vxDot_sol(i,xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za,vx_s);
  K3vysol= vyDot_sol(i,xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za,vy_s);
  K3vzsol= vzDot_sol(i,xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za,vz_s);

  K3vxtierra= vxDot_tierra(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_t);
  K3vytierra= vyDot_tierra(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_t);
  K3vztierra= vzDot_tierra(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_t);

  K3vxluna= vxDot_luna(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_l);
  K3vyluna= vyDot_luna(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vy_l);
  K3vzluna= vzDot_luna(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vz_l);

  K3vxast= vxDot_ast(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_a);
  K3vyast= vyDot_ast(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vy_a);
  K3vzast= vzDot_ast(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vz_a);

  xs[i]=saltoxs+(K3xsol*h/2);
  ys[i]=saltoys+(K3ysol*h/2);
  zs[i]=saltozs+(K3zsol*h/2);
  xt[i]=saltoxt+(K3xtierra*h/2);
  yt[i]=saltoyt+(K3ytierra*h/2);
  zt[i]=saltozt+(K3ztierra*h/2);
  xl[i]=saltoxl+(K3xluna*h/2);
  yl[i]=saltoyl+(K3yluna*h/2);
  zl[i]=saltozl+(K3zluna*h/2);
  xa[i]=saltoxa+(K3xast*h/2);
  ya[i]=saltoya+(K3yast*h/2);
  za[i]=saltoza+(K3zast*h/2);
  
  vx_s[i]=saltovxs+(K3vxsol*h/2);
  vy_s[i]=saltovys+(K3vysol*h/2);
  vz_s[i]=saltovzs+(K3vzsol*h/2);
  vx_t[i]=saltovxt+(K3vxtierra*h/2);
  vy_t[i]=saltovyt+(K3vytierra*h/2);
  vz_t[i]=saltovzt+(K3vztierra*h/2);
  vx_l[i]=saltovxl+(K3vxluna*h/2);
  vy_l[i]=saltovyl+(K3vyluna*h/2);
  vz_l[i]=saltovzl+(K3vzluna*h/2);
  vx_a[i]=saltovxa+(K3vxast*h/2);
  vy_a[i]=saltovya+(K3vyast*h/2);
  vz_a[i]=saltovza+(K3vzast*h/2);

  //K4
  K4xsol= xDot_sol(i,xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za,vx_s);
  K4ysol= yDot_sol(i,xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za,vy_s);
  K4zsol= zDot_sol(i,xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za,vz_s);

  K4xtierra= xDot_tierra(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_t);
  K4ytierra= yDot_tierra(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vy_t);
  K4ztierra= zDot_tierra(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vz_t);

  K4xluna= xDot_luna(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_l);
  K4yluna= yDot_luna(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vy_l);
  K4zluna= zDot_luna(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vz_l);

  K4xast= xDot_ast(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_a);
  K4yast= yDot_ast(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vy_a);
  K4zast= zDot_ast(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vz_a);
  
  K4vxsol= vxDot_sol(i,xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za,vx_s);
  K4vysol= vyDot_sol(i,xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za,vy_s);
  K4vzsol= vzDot_sol(i,xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za,vz_s);

  K4vxtierra= vxDot_tierra(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_t);
  K4vytierra= vyDot_tierra(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_t);
  K4vztierra= vzDot_tierra(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_t);

  K4vxluna= vxDot_luna(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_l);
  K4vyluna= vyDot_luna(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vy_l);
  K4vzluna= vzDot_luna(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vz_l);

  K4vxast= vxDot_ast(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vx_a);
  K4vyast= vyDot_ast(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vy_a);
  K4vzast= vzDot_ast(i, xs, ys, zs, xt, yt, zt, xl, yl, zl, xa, ya, za, vz_a);
  
  /*Nuevas velocidades y posiciones en x*/
  xs[d]=saltoxs+(h*(K1xsol+2*K2xsol+2*K3xsol+K4xsol)/6);
  xt[d]=saltoxt+(h*(K1xtierra+2*K2xtierra+2*K3xtierra+K4xtierra)/6);
  xl[d]=saltoxl+(h*(K1xluna+2*K2xluna+2*K3xluna+K4xluna)/6);
  xa[d]=saltoxa+(h*(K1xast+2*K2xast+2*K3xast+K4xast)/6);
  vx_s[d]=saltovxs+(h*(K1vxsol+2*K2vxsol+2*K3vxsol+K4vxsol)/6);    
  vx_t[d]=saltovxt+(h*(K1vxtierra+2*K2vxtierra+2*K3vxtierra+K4vxtierra)/6);
  vx_l[d]=saltovxl+(h*(K1vxluna+2*K2vxluna+2*K3vxluna+K4vxluna)/6);
  vx_a[d]=saltovxa+(h*(K1vxast+2*K2vxast+2*K3vxast+K4vxast)/6); 

  /*Nuevas velocidades y posiciones en y*/
  ys[d]=saltoys+(h*(K1ysol+2*K2ysol+2*K3ysol+K4ysol)/6);
  yt[d]=saltoyt+(h*(K1ytierra+2*K2ytierra+2*K3ytierra+K4ytierra)/6);
  yl[d]=saltoyl+(h*(K1yluna+2*K2yluna+2*K3yluna+K4yluna)/6);
  ya[d]=saltoya+(h*(K1yast+2*K2yast+2*K3yast+K4yast)/6);
  vy_s[d]=saltovys+(h*(K1vysol+2*K2vysol+2*K3vysol+K4vysol)/6);    
  vy_t[d]=saltovyt+(h*(K1vytierra+2*K2vytierra+2*K3vytierra+K4vytierra)/6);
  vy_l[d]=saltovyl+(h*(K1vyluna+2*K2vyluna+2*K3vyluna+K4vyluna)/6);
  vy_a[d]=saltovya+(h*(K1vyast+2*K2vyast+2*K3vyast+K4vyast)/6);

  /*Nuevas velocidades y posiciones en z*/
  zs[d]=saltozs+(h*(K1zsol+2*K2zsol+2*K3zsol+K4zsol)/6);
  zt[d]=saltozt+(h*(K1ztierra+2*K2ztierra+2*K3ztierra+K4ztierra)/6);
  zl[d]=saltozl+(h*(K1zluna+2*K2zluna+2*K3zluna+K4zluna)/6);
  za[d]=saltoza+(h*(K1zast+2*K2zast+2*K3zast+K4zast)/6);
  vz_s[d]=saltovzs+(h*(K1vzsol+2*K2vzsol+2*K3vzsol+K4vzsol)/6);    
  vz_t[d]=saltovzt+(h*(K1vztierra+2*K2vztierra+2*K3vztierra+K4vztierra)/6);
  vz_l[d]=saltovzl+(h*(K1vzluna+2*K2vzluna+2*K3vzluna+K4vzluna)/6);
  vz_a[d]=saltovza+(h*(K1vzast+2*K2vzast+2*K3vzast+K4vzast)/6);

  return  xs[d], ys[d], xt[d], yt[d], xl[d], yl[d], xa[d], ya[d];
}

/*Memory allocation*/
float * get_memory(int n_points){
  float * x; 
  if(!(x = malloc(sizeof(float) * n_points))){
    printf("problem with memory allocation");
    exit(1);
  }
  return x;
}


