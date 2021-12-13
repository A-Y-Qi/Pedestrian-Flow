
#define PI 3.14159265

/* function declaration */
struct coord *Orig_F(double desired_speed, double rela_tau_i, double actual_velo_x, double actual_velo_y,
                     double pedestrian_i_x, double pedestrian_i_y, double temp_des_x, double temp_des_y);

//structure definition
struct coord
{
  double F_orig_x;
  double F_orig_y;
  double desdir_x;
  double desdir_y;
};

struct coord *Orig_F(double desired_speed, double rela_tau_i, double actual_velo_x, double actual_velo_y,
                     double pedestrian_i_x, double pedestrian_i_y, double temp_des_x, double temp_des_y)
{
  /* local variable declaration */
  double des_dir_x, des_dir_y;
  double norm;
  double F_x, F_y;

  /*body function*/

  //find desired direction
  norm = sqrt(pow(temp_des_x - pedestrian_i_x, 2) + pow(temp_des_y - pedestrian_i_y, 2));
  //printf("norm (%f)\n", norm);
  des_dir_x = (temp_des_x - pedestrian_i_x) / norm;
  des_dir_y = (temp_des_y - pedestrian_i_y) / norm;

  //find original force
  F_x = (des_dir_x * desired_speed - actual_velo_x) / rela_tau_i;
  F_y = (des_dir_y * desired_speed - actual_velo_y) / rela_tau_i;

  struct coord *F_orig = malloc(sizeof(struct coord));
  F_orig->F_orig_x = F_x;
  F_orig->F_orig_y = F_y;
  F_orig->desdir_x = des_dir_x;
  F_orig->desdir_y = des_dir_y;

  return F_orig;
}