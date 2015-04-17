#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define G_GRAV 39.486
float * get_memory(int n_points);

void vel_i(float *vx, float *vy, float *vz, int n_points);
void pos_i(float *x, float *y, float *z, int n_points);
float get_norma2(int i, int j);
float get_ace(float m, float r_i, float r_j, int i , int j);
void get_ace_x();
void get_ace_y();
void get_ace_z();

//void initialize_acele(float *ax, float *ay, float *az, int n_points);
int n_points = 4;
float *r;
 /*positions of all particles*/
  float *x;
  float *y;
  float *z;
  
  /*velocities of all particles*/
  float *vx;
  float *vy;
  float *vz;

  /*accelerations of all particles*/
  float *a_x;
  float *a_y;
  float *a_z;

  float *mass;

int main(void){

  /*memory allocation*/
  r = get_memory(n_points);
  x = get_memory(n_points);
  y = get_memory(n_points);
  z = get_memory(n_points);

  a_x = get_memory(n_points);
  a_y = get_memory(n_points);
  a_z = get_memory(n_points);

  vx = get_memory(n_points);
  vy = get_memory(n_points);
  vz = get_memory(n_points);

  mass = get_memory(n_points);


  vel_i(vx,vy,vz, n_points);
  pos_i(x,y,z,n_points);
  // get_norma(x, y, z, n_points,r);

  get_ace_x();
  
//  int i;
//  for(i=0;i<10;i++){
 //   printf("%f ",r[i]);

}

void vel_i(float *vx, float *vy, float *vz, int n_points)
{
  int i; 
  //FLOAT delta_theta;
  //delta_theta = 2.0*PI/n_points;

  float vx_sol=0;
  float vy_sol=0;
  float vz_sol=0;
  float vx_luna=0;
  float vy_luna=1.023;
  float vz_luna=0;
  float vx_tierra=0;
  float vy_tierra=30;
  float vz_tierra=0;
  float vx_asteroide=1570;
  float vy_asteroide=0;
  float vz_asteroide=29959;

 vx[0]=vx_sol;
 vy[0]=vy_sol;
 vz[0]=vz_sol;

 vx[1] = vx_tierra;
 vy[1] = vy_tierra; 
 vz[1] = vz_tierra; 

 vx[2] = vx_luna;
 vy[2] = vy_luna;
 vz[2] = vz_luna;

 vx[3] = vx_asteroide; 
 vy[3] = vy_asteroide; 
 vz[3] = vz_asteroide;

 }
void pos_i(float *x, float *y, float *z, int n_points)
{
  int i; 
  //FLOAT delta_theta;
  //delta_theta = 2.0*PI/n_points;

  float x_sol=0;
  float y_sol=0;
  float z_sol=0;
  float x_luna=1.0026;
  float y_luna=0;
  float z_luna=0;
  float x_tierra=1;
  float y_tierra=0;
  float z_tierra=0;
  float x_asteroide=1;
  float y_asteroide=0.0963;
  float z_asteroide=0;

 x[0]=x_sol;
 y[0]=y_sol;
 z[0]=z_sol;

 x[1] = x_tierra;
 y[1] = y_tierra; 
 z[1] = z_tierra; 

 x[2] = x_luna;
 y[2] = y_luna;
 z[2] = z_luna;

 x[3] = x_asteroide; 
 y[3] = y_asteroide; 
 z[3] = z_asteroide;

 mass[0] = 1;
 mass[1] = 0.000003003;
 mass[2] = 0.00000003695;
 mass[3] = 0.00000000005;
 
 }

 float get_norma2(int i, int j){

 	return pow(sqrt(pow((x[i] - x[j]),2.0) + pow((y[i] - y[j]),2.0) + pow((z[i] - z[j]),2.0)), 3.0);
 	//return x[i];
 }

 float get_ace(float m, float r_i, float r_j, int i , int j){

 	float norm = get_norma2(i,j);

 	if(norm!=0){

 		return (G_GRAV * m * abs(r_i-r_j) ) / norm;

 	}else{

 		printf("Colision!");
 		return 0;
 	}
 }

void get_ace_x(){

	int i;
	int j;

	for(i=0;i<n_points;i++)
	{
		for(j=0;j<n_points;j++){

			if(i!=j){

				a_x[i] += get_ace(mass[j],x[i],x[j],i,j);
				//a_y[i] += get_ace(mass[j],y[i],y[j],i,j);
				//a_z[i] += get_ace(mass[j],z[i],z[j],i,j);
			}

		}
		printf("a_x! %i %f \n", i, a_x[i]);
		//printf("a_y! %i %f \n", i, a_y[i]);
		//printf("a_z! %i %f \n", i, a_z[i]);
	}

}

void get_ace_y(){

	int i;
	int j;

	for(i=0;i<n_points;i++)
	{
		for(j=0;j<n_points;j++){


			if(i!=j){

				a_y[i] += get_ace(mass[j],y[i],y[j],i,j);
			}
		}	
	}
}

void get_ace_z(){

	int i;
	int j;

	for(i=0;i<n_points;i++)
	{
		for(j=0;j<n_points;j++){


			if(i!=j){

				a_z[i] += get_ace(mass[j],z[i],z[j],i,j);
			
			}

		}	
	}

}

float * get_memory(int n_points){
  float * x; 
  if(!(x = malloc(sizeof(float) * n_points))){
    printf("problem with memory allocation");
    exit(1);
  }
  return x;
}


