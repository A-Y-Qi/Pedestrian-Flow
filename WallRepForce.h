#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#define Uiw 10
#define R 0.2

/* suppose all actual velocity are saved in an array */
/*suppose All desired direction e are saved in an array */

/* function declaration */
struct coordWallF *WallRep_F(double wall_x, double wall_y, double pedestrian_i_x, double pedestrian_i_y);

//structure definition
struct coordWallF
{
  double WallRepForce_x;
  double WallRepForce_y;
};

struct coordWallF *WallRep_F(double wall_x, double wall_y, double pedestrian_i_x, double pedestrian_i_y)
{
  /* local variable declaration */
  double r_x, r_y; //distance between i and j
  double norm_r;
  double Fx_wall, Fy_wall;

  /*body function*/

  //find the vector between i and j
  r_x = pedestrian_i_x - wall_x;
  r_y = pedestrian_i_y - wall_y;
  norm_r = sqrt(pow(r_x, 2) + pow(r_y, 2));

  //Find F wall
  Fx_wall=Uiw*r_x*exp(-norm_r/R)*1/(norm_r*R);
  Fy_wall=Uiw*r_y*exp(-norm_r/R)*1/(norm_r*R);

  struct coordWallF *WallRepForce = malloc(sizeof(struct coordWallF));
  WallRepForce->WallRepForce_x = Fx_wall;
  WallRepForce->WallRepForce_y = Fy_wall;

  return WallRepForce;
}
