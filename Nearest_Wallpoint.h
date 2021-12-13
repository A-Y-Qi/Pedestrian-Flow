#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#define PI 3.14159265

/* function declaration */
struct coordWallpoint *Nearest_Point(double wall_start_x, double wall_start_y,
                            double wall_end_x, double wall_end_y, double pedestrian_i_x, double pedestrian_i_y);

//structure definition
struct coordWallpoint
{
    double x_wall;
    double y_wall;
};

struct coordWallpoint *Nearest_Point(double wall_start_x, double wall_start_y, double wall_end_x, double wall_end_y, double pedestrian_i_x, double pedestrian_i_y)
{
    /* local variable declaration */
    double rela_ped_x, rela_ped_y;
    double rela_wall_x, rela_wall_y;
    double wall_length, ped_length;
    double dot, dot_norm;
    double nearest_x, nearest_y;

    /*body function*/
    //relative wall position(vector)
    rela_wall_x = wall_end_x - wall_start_x;
    rela_wall_y = wall_end_y - wall_start_y;
    //relative pedestrian position to the wall (vector)
    rela_ped_x = pedestrian_i_x - wall_start_x;
    rela_ped_y = pedestrian_i_y - wall_start_y;
    //length of wall vector and pedestrian vector
    wall_length = sqrt(pow(rela_wall_x, 2.0) + pow(rela_wall_y, 2.0));
    //printf("the wall length is %f\n",wall_length);
    ped_length = sqrt(pow(rela_ped_x, 2.0) + pow(rela_ped_y, 2.0));
    //dot product between the the ped vector and wall vec to determin the closest point
    dot = rela_wall_x * rela_ped_x + rela_wall_y * rela_ped_y;
    //printf("the wall vector is (%f,%f), the pedestrian vector is (%f,%f)\n",rela_wall_x,rela_wall_y,rela_ped_x,rela_ped_y);
    //printf("the dot product is %f\n",dot);
    //find the length of the nearist point on the wall wrt start point
    dot_norm = dot / wall_length;
    if (dot_norm < 0.0)
    {
        nearest_x = wall_start_x;
        nearest_y = wall_start_y;
    }
    else if (dot_norm > wall_length)
    {
        nearest_x = wall_end_x;
        nearest_y = wall_end_y;
    }
    else
    {
        nearest_x = rela_wall_x / wall_length * dot_norm + wall_start_x;
        nearest_y = rela_wall_y / wall_length * dot_norm + wall_start_y;
    }
    struct coordWallpoint *nearest_wall_point = malloc(sizeof(struct coord));
    nearest_wall_point->x_wall = nearest_x;
    nearest_wall_point->y_wall = nearest_y;

    return nearest_wall_point;
}
