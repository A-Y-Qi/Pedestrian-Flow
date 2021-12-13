#include "stdio.h"
#include "stdlib.h"
#include "math.h"

/* suppose all actual velocity are saved in an array */
/*suppose All desired direction e are saved in an array */

/* function declaration */
struct coordRep *Rep_F(double pedestrian_j_x, double pedestrian_j_y, double velo_j_x, double velo_j_y,
                    double pedestrian_i_x, double pedestrian_i_y, double des_j_x, double des_j_y, double dt, double Sigma);

//structure definition
struct coordRep
{
  double RepForce_x;
  double RepForce_y;
};

/*
int main()
{

    //setup constants

    //set up variables
    double X[pedestrian], Y[pedestrian];                   //position pedestrian
    double desirespeed0[pedestrian];                       //desired speed of the pedestrian (no direction)
    double act_velo_x[pedestrian], act_velo_y[pedestrian]; //actual velocity of the pedestrian
    double des_ang_x[pedestrian], des_ang_y[pedestrian];   //desired angle of the pedestrian
    double F_x, F_y;                                       //repulsive force
    struct coord *RepForce_poi = Rep_F(1.0, 2.0, 1.0, 2.0, 1.1, 2.2, 1 / (sqrt(5)), 2 / (sqrt(5)));
    F_x = RepForce_poi->RepForce_x;
    F_y = RepForce_poi->RepForce_y;
    printf("The result is (%f, %f)", F_x, F_y);
}*/

struct coordRep *Rep_F(double pedestrian_j_x, double pedestrian_j_y, double velo_j_x, double velo_j_y,
                    double pedestrian_i_x, double pedestrian_i_y, double des_j_x, double des_j_y, double dt, double Sigma)
{
  /* local variable declaration */
  double A, B, K, b2, b; //b is the semiminor axis of the ellipse
  double speed_j;
  double r_x, r_y;       //distance between i and j
  double j_x,j_y;        // distance between where j is going us current i
  double V_rep;          // repulsive potential
  double grad_x, grad_y; //find the gradient of the potential
  double norm_r,norm_j;
  double V=2.1;

  /*body function*/

  //find actual speed of j
  speed_j = sqrt(pow(velo_j_x, 2) + pow(velo_j_y, 2));
  

  //find the vector between i and j
  r_x = pedestrian_i_x - pedestrian_j_x;
  r_y = pedestrian_i_y - pedestrian_j_y;
  norm_r = sqrt(pow(r_x, 2) + pow(r_y, 2));


  //find the vector between desired j and current i
    j_x=r_x-speed_j*dt*des_j_x;
    j_y=r_y-speed_j*dt*des_j_y;
    norm_j=sqrt(pow(j_x, 2) + pow(j_y, 2));

    grad_x=1/(2*Sigma*sqrt(pow(norm_r+norm_j,2)-pow(speed_j*dt*des_j_x,2)))*(V*(norm_r+norm_j)\
    *(r_x/norm_r+(-2*speed_j*dt*des_j_x+2*r_x)/(2*norm_j))*exp(-(sqrt(pow(norm_r+norm_j,2)-pow(speed_j*dt*des_j_x,2))))/(2*Sigma));

    grad_y=1/(2*Sigma*sqrt(pow(norm_r+norm_j,2)-pow(speed_j*dt*des_j_y,2)))*(V*(norm_r+norm_j)\
    *(r_y/norm_r+(-2*speed_j*dt*des_j_y+2*r_y)/(2*norm_j))*exp(-(sqrt(pow(norm_r+norm_j,2)-pow(speed_j*dt*des_j_y,2))))/(2*Sigma));

  /*

  //find constant A B K b
  A = speed_j * dt * des_j_x;
  B = speed_j * dt * des_j_y;
  K = speed_j * dt;
  //b=1/2*(sqrt(pow(norm_r+sqrt(pow(r_x-A,2)+pow(r_y-B,2)),2)-K*K))
  b2 = (double)sqrt(pow(norm_r + sqrt(pow(r_x - A, 2) + pow(r_y - B, 2)), 2) - K * K);
  b = (double)0.5 * b2;
  //printf("The b is (%f)\n", b);
  if (b == 0)
  {
    grad_x = -r_y;
    grad_y = -r_x;
  }
  else
  {
    V_rep = 2.1 * exp(-b / Sigma);

    grad_x = -V_rep * (-1 / (2 * Sigma)) * ((r_x / norm_r) + (2 * r_x - 2 * A) / (2 * sqrt(pow(r_x - A, 2) + pow(r_y - B, 2)))) / (4 * b);
    grad_y = -V_rep * (-1 / (2 * Sigma)) * ((r_y / norm_r) + (2 * r_y - 2 * A) / (2 * sqrt(pow(r_x - A, 2) + pow(r_y - B, 2)))) / (4 * b);
  }

  */

  //test
  /*
    double test1, test2, test3;
    test1 = pow(norm_r + sqrt(pow(r_x - A, 2) + pow(r_y - B, 2)), 2);
    test2 = sqrt(pow(norm_r + sqrt(pow(r_x - A, 2) + pow(r_y - B, 2)), 2) - K * K);
    test3 = pow(norm_r + sqrt(pow(r_x - A, 2) + pow(r_y - B, 2)), 2) - K * K;
    printf("The test1 is (%f)\n", test1);
    printf("The test2 is (%f)\n", test2);
    printf("The test3 is (%f)\n", test3);
    */

  struct coordRep *RepForce_poi = malloc(sizeof(struct coordRep));
  RepForce_poi->RepForce_x = grad_x;
  RepForce_poi->RepForce_y = grad_y;

  return RepForce_poi;
}
