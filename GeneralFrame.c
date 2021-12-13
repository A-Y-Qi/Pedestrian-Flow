#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "Original_Force.h"
#include "RepForce.h"
#include "Nearest_Wallpoint.h"
#include "WallRepForce.h"

#define pedestrian 70
#define PI 3.14159265
#define wall_number 2
#define dt 0.1
#define c 0.5
#define Vij 2.1

//index 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34

int main()
{
    /*let me start with no walls in the middle */
    // setup constants
    double v_max = 1.3;
    double sigma = 0.26, avg = 1.34; // gaussian distribution condition

    // set up variables
    float ran1, ran2, ran3, ran4;                                 // pseudo Random numbers with normal distribution
    float gr1, gr2, gr3;                                          // normal random numbers
    double v_0[pedestrian], v_0_max[pedestrian], tmax, max_speed; // desired speed and max speed
    double r0_x[pedestrian], r0_y[pedestrian];                    // destination
    double X[pedestrian], Y[pedestrian];                          // position pedestrian
    double e_x[pedestrian], e_y[pedestrian];                      // desired direction of pedestrian
    double v_x[pedestrian], v_y[pedestrian];                      // actual velocity

    // Forces
    // neigbour
    double F_rep_x[pedestrian], F_rep_y[pedestrian], F_rep; // repulsion force between neighbours
    double r_x, r_y, norm_r, j_x, j_y, norm_j, speed_j;
    double k = 0.1;
    // walls
    double F_wall1_x[pedestrian], F_wall1_y[pedestrian];            // Wall repulsion force
    double F_wall2_x[pedestrian], F_wall2_y[pedestrian];            // Wall repulsion force
    double Fsum_wall_x[pedestrian], Fsum_wall_y[pedestrian];        // Sum Wall repulsion force
    double F_0_x[pedestrian], F_0_y[pedestrian], F_0, dot_p, angle; // initial force with angle = 100
    double F_total_x[pedestrian], F_total_y[pedestrian];            // total force

    FILE *outwall1, *outwall2;
    outwall1 = fopen("General_wall1.xyz", "w");
    outwall2 = fopen("General_wall2.xyz", "w");

    // position of walls [x,y]
    double wall1_start[2] = {0, 0};
    double wall1_end[2] = {50, 0};
    double wall2_start[2] = {0, 10};
    double wall2_end[2] = {50, 10};
    double wall_len = 50, wall_hei = 10;
    double tau = 0.5;

    int i, j, kpp, tooclose; // for loop use

    //generate wall atoms
    fprintf(outwall1, "%i\n", 51);
    fprintf(outwall2, "%i\n", 51);

    //wall 1
    int w1x=0,w1y=0;
    fprintf(outwall1, "a%i %i  %i  0\n", 0, w1x, w1y);

    for ( i = 1; i <= 51; i++)
    {
        w1x++;

        fprintf(outwall1, "a%i %i  %i  0\n", i, w1x, w1y);
        
    }
    //wall2
    int w2x=0,w2y=10;
    fprintf(outwall2, "a%i %i  %i  0\n", 0, w2x, w2y);

    for ( i = 1; i <= 51; i++)
    {
        w2x++;

        fprintf(outwall2, "a%i %i  %i  0\n", i, w2x, w2y);
        
    }
    


    
   



    double W_x, W_y;         // nearist wall point
    double ei_x, ei_y;       // desired direction of pedestrian i
    

    double test1, test2, testmin = 100;
    int count = 0;

    // open up file to print

    FILE *outpos1, *outpos2;
    outpos1 = fopen("PedestrianPosition.xyz", "w");
    outpos2 = fopen("PedestrianPosition.dat", "w");
    // outx = fopen("rsGaussianBDtest2t45.dat", "w");
    // outtrj = fopen("rsGaussianBDtest2t45.xyz", "w");

    fprintf(outpos1, "%i\n", pedestrian - 1);

    // setup desired speed for each individual by box muller
    for (i = 0; i < pedestrian; i++)
    {
        ran1 = (float)rand() / RAND_MAX;
        ran2 = (float)rand() / RAND_MAX;
        gr1 = avg + sigma * (double)sqrt(-2 * log(ran1)) * cos(2 * 3.14159265 * ran2);
        v_0[i] = gr1;
        // printf("the gaussian speed is %f\n",gr1);
    }

    // setup desired max speed for each individual
    for (i = 0; i < pedestrian; i++)
    {
        v_0_max[i] = v_max * v_0[i];
    }

    // Find tmax
    tmax = wall_len / avg;
    float t;
    for (i = 0; i < pedestrian; i++)
    {
        t = wall_len * 0.5 / v_0_max[i];
        if (t < tmax)
        {
            tmax = t;
            max_speed = v_0_max[i];
        }
    }

    // printf("The max speed is j=%f", tmax);

    // set up the desired direction of the pedsestrian
    for (i = 0; i < pedestrian; i++)
    {
        if (i <= (int)floor(pedestrian / 2))
        {
            e_x[i] = 1;
            e_y[i] = 0;
        }
        else
        {
            e_x[i] = -1;
            e_y[i] = 0;
        }
    }

    // set up initial position of particles

    float distance;

    for (i = 0; i < pedestrian; i++)
    {
        kpp = 0;
        while (!kpp)
        {

            if (i <= (int)floor(pedestrian / 2))
            {
                ran3 = (float)rand() / RAND_MAX;
                ran4 = (float)rand() / RAND_MAX;
                X[i] = ran3 * wall_len * 0.5;
                Y[i] = ran4 * wall_hei;
            }
            else
            {
                ran3 = (float)rand() / RAND_MAX;
                ran4 = (float)rand() / RAND_MAX;
                X[i] = ran3 * wall_len * 0.5 + wall_len * 0.5;
                Y[i] = ran4 * wall_hei;
            }

            tooclose = 0;

            for (j = 0; j < i; j++)
            {
                distance = sqrt(pow(X[i] - X[j], 2.0) + pow(Y[i] - Y[j], 2.0));
                if (distance < dt * fmax(v_0[i], v_0[j]))
                // if (distance < 2.1)
                {
                    tooclose = 1;
                }
            }
            if (!tooclose)
            {
                kpp = 1;
            }
        }
        // printf("the initial position of %i pedestian is x=%f, y=%f\n", i,X[i],Y[i]);
        fprintf(outpos1, "a%i %f  %f  0\n", i, X[i], Y[i]);
        fprintf(outpos2, "%i %f  %f  0\n", i, X[i], Y[i]);
    }

    // setup initial velocity of the each individuals all start at the same speed =1.34
    for (i = 0; i < pedestrian; i++)
    {
        v_x[i] = e_x[i] * avg;
        v_y[i] = e_y[i] * avg;
    }

    // timestep
    double temp_des_x, temp_des_y;                     // temporary desired position
    double Fsum_nei_x, Fsum_nei_y, Fj_nei_x, Fj_nei_y; // temporary neighbour repulsion force

    for (t = 0; t < 100; t += dt)
    {
        count++;
        fprintf(outpos1, "%i\n", pedestrian - 1);
        for (i = 0; i < pedestrian; i++)
        {
            // original force if no obstacles
            temp_des_x = X[i] + e_x[i] * v_0[i] * dt;
            temp_des_y = Y[i] + e_y[i] * v_0[i] * dt;
            struct coord *orig_F_poi = Orig_F(v_0[i], tau, v_x[i], v_y[i],
                                              X[i], Y[i], temp_des_x, temp_des_y);
            F_0_x[i] = orig_F_poi->F_orig_x;
            F_0_y[i] = orig_F_poi->F_orig_y;

            F_0_x[i] = (v_0[i] * e_x[i] - v_x[i]) / tau;
            F_0_y[i] = (v_0[i] * e_y[i] - v_y[i]) / tau;
            // printf("The original force [%i] is x=%f,y=%f \n",i,F_0_x[i],F_0_y[i]);

            // neightbour repulsion Force
            Fsum_nei_x = 0;
            Fsum_nei_y = 0;

            for (j = 0; j < pedestrian; j++)
            {
                if (j != i)
                {
                    speed_j = sqrt(pow(v_x[j], 2) + pow(v_y[j], 2));
                    // find the vector between i and j
                    r_x = X[i] - X[j];
                    r_y = Y[i] - Y[j];
                    norm_r = sqrt(pow(r_x, 2) + pow(r_y, 2));
                    // printf("norm_r=%f",norm_r);

                    j_x = r_x - speed_j * dt * e_x[j];
                    j_y = r_y - speed_j * dt * e_y[j];
                    norm_j = sqrt(pow(j_x, 2) + pow(j_y, 2));

                    Fj_nei_x = 1 / (2 * k * sqrt(pow(norm_r + norm_j, 2) - pow(speed_j * dt * e_x[j], 2))) * (Vij * (norm_r + norm_j) * (r_x / norm_r + (-2 * speed_j * dt * e_x[j] + 2 * r_x) / (2 * norm_j)) * exp(-(sqrt(pow(norm_r + norm_j, 2) - pow(speed_j * dt * e_x[j], 2)))) / (2 * k));
                    Fj_nei_y = 1 / (2 * k * sqrt(pow(norm_r + norm_j, 2) - pow(speed_j * dt * e_y[j], 2))) * (Vij * (norm_r + norm_j) * (r_y / norm_r + (-2 * speed_j * dt * e_y[j] + 2 * r_y) / (2 * norm_j)) * exp(-(sqrt(pow(norm_r + norm_j, 2) - pow(speed_j * dt * e_y[j], 2)))) / (2 * k));

                    dot_p = e_x[i] * (-Fj_nei_x) + e_y[i] * (-Fj_nei_y);
                    angle = sqrt(Fj_nei_x * Fj_nei_x + Fj_nei_y * Fj_nei_y) * cos(5 / 9 * PI);
                    if (dot_p < angle)
                    {
                        Fj_nei_x = c * Fj_nei_x;
                        Fj_nei_y = c * Fj_nei_y;
                    }

                    Fsum_nei_x = Fsum_nei_x + Fj_nei_x;
                    Fsum_nei_y = Fsum_nei_y + Fj_nei_y;
                }
            }

            F_rep_x[i] = Fsum_nei_x;
            F_rep_y[i] = Fsum_nei_y;

            // wall repulson force

            double wall1_x = 0, wall1_y = 0, wall2_x = 0, wall2_y = 0;
            struct coordWallpoint *Wall1_poi = Nearest_Point(wall1_start[0], wall1_start[1],
                                                             wall1_end[0], wall1_end[1], X[i], Y[i]);
            wall1_x = Wall1_poi->x_wall;
            wall1_y = Wall1_poi->y_wall;
            // printf("The nearest for wall 1 of [%i] is x=%f,y=%f \n",i,wall1_x,wall1_y);

            struct coordWallpoint *Wall2_poi = Nearest_Point(wall2_start[0], wall2_start[1],
                                                             wall2_end[0], wall2_end[1], X[i], Y[i]);
            wall2_x = Wall2_poi->x_wall;
            wall2_y = Wall2_poi->y_wall;
            // printf("The nearest for wall 2 of [%i] is x=%f,y=%f \n",i,wall2_x,wall2_y);

            struct coordWallF *Wall1F_poi = WallRep_F(wall1_x, wall1_y, X[i], Y[i]);
            F_wall1_x[i] = Wall1F_poi->WallRepForce_x;
            F_wall1_y[i] = Wall1F_poi->WallRepForce_y;
            dot_p = e_x[i] * F_wall1_x[i] + e_y[i] * F_wall1_y[i];

            struct coordWallF *Wall2F_poi = WallRep_F(wall2_x, wall2_y, X[i], Y[i]);
            F_wall2_x[i] = Wall2F_poi->WallRepForce_x;
            F_wall2_y[i] = Wall2F_poi->WallRepForce_y;
            dot_p = e_x[i] * F_wall2_x[i] + e_y[i] * F_wall2_y[i];

            Fsum_wall_x[i] = F_wall1_x[i] + F_wall2_x[i];
            Fsum_wall_y[i] = F_wall1_y[i] + F_wall2_y[i];

            F_total_x[i] = F_0_x[i] + F_rep_x[i] + Fsum_wall_x[i];
            F_total_y[i] = F_0_y[i] + F_rep_y[i] + Fsum_wall_y[i];

            v_x[i] += F_total_x[i] * dt;
            v_y[i] += F_total_y[i] * dt;

            if (i == 1 && v_x[i] < 0)
            {
                printf("The velocity of pedestrian %i is vx=%f, vy=%f at t=%f \n", i, v_x[i], v_y[i], t);
            }

            if (sqrt(v_x[i] * v_x[i] + v_y[i] * v_y[i]) >= v_0_max[i])
            {
                v_x[i] = v_x[i] / sqrt(v_x[i] * v_x[i] + v_y[i] * v_y[i]) * v_0_max[i];
                v_y[i] = v_y[i] / sqrt(v_x[i] * v_x[i] + v_y[i] * v_y[i]) * v_0_max[i];
            }

            X[i] = X[i] + v_x[i] * dt;
            Y[i] = Y[i] + v_y[i] * dt;

            if (X[i] > 50)
            {
                X[i] = X[i] - 50;
            }
            if (X[i] < 0)
            {
                X[i] = X[i] + 50;
            }

            // printf("The new position of [%i] is x=%f,y=%f \n", i, X[i], Y[i]);

            fprintf(outpos1, "a%i %f  %f  0\n", i, X[i], Y[i]);
            fprintf(outpos2, "%i %f  %f  0\n", i, X[i], Y[i]);
        }
    }
}
