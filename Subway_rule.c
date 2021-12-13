#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "Original_Force.h"
#include "RepForce.h"
#include "Nearest_Wallpoint.h"
#include "WallRepForce.h"

#define pedestrian 60
#define PI 3.14159265
#define wall_number 7
#define dt 0.1
#define c 0.5
#define Vij 2.1

// index 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29

// the total time take is 12.950000 when <9

int main()
{
    printf("start");
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
    double F_wall3_x[pedestrian], F_wall3_y[pedestrian];            // Wall repulsion force
    double F_wall4_x[pedestrian], F_wall4_y[pedestrian];            // Wall repulsion force
    double F_wall5_x[pedestrian], F_wall5_y[pedestrian];            // Wall repulsion force
    double F_wall6_x[pedestrian], F_wall6_y[pedestrian];            // Wall repulsion force
    double F_wall7_x[pedestrian], F_wall7_y[pedestrian];            // Wall repulsion force
    double F_wall8_x[pedestrian], F_wall8_y[pedestrian];            // Wall repulsion force
    double F_wall9_x[pedestrian], F_wall9_y[pedestrian];            // Wall repulsion force
    double F_wall10_x[pedestrian], F_wall10_y[pedestrian];          // Wall repulsion force
    double F_wall11_x[pedestrian], F_wall11_y[pedestrian];          // Wall repulsion force
    double F_wall12_x[pedestrian], F_wall12_y[pedestrian];          // Wall repulsion force
    double F_wall13_x[pedestrian], F_wall13_y[pedestrian];          // Wall repulsion force
    double Fsum_wall_x[pedestrian], Fsum_wall_y[pedestrian];        // Sum Wall repulsion force
    double F_0_x[pedestrian], F_0_y[pedestrian], F_0, dot_p, angle; // initial force with angle = 100
    double F_total_x[pedestrian], F_total_y[pedestrian];            // total force

    // position of walls [x,y]
    double wall1_start[2] = {0, 0};
    double wall1_end[2] = {14, 0};
    double wall2_start[2] = {10, 0};
    double wall2_end[2] = {10, 4};
    double wall3_start[2] = {10, 7};
    double wall3_end[2] = {10, 10};
    double wall4_start[2] = {10, 13};
    double wall4_end[2] = {10, 16};
    double wall5_start[2] = {10, 19};
    double wall5_end[2] = {10, 23};
    double wall6_start[2] = {0, 23};
    double wall6_end[2] = {14, 23};
    double wall7_start[2] = {14, 23};
    double wall7_end[2] = {14, 0};

    // invisible walls
    double wall8_start[2] = {10, 4};
    double wall8_end[2] = {4, 4};
    double wall9_start[2] = {10, 7};
    double wall9_end[2] = {4, 7};
    double wall10_start[2] = {10, 10};
    double wall10_end[2] = {4, 10};
    double wall11_start[2] = {10, 13};
    double wall11_end[2] = {4, 13};
    double wall12_start[2] = {10, 16};
    double wall12_end[2] = {4, 16};
    double wall13_start[2] = {10, 19};
    double wall13_end[2] = {4, 19};

    double tau = 0.5, t = 0;

    double W_x, W_y;         // nearist wall point
    double ei_x, ei_y;       // desired direction of pedestrian i
    int i, j, kpp, tooclose; // for loop use

    FILE *outwall1, *outwall2, *outwall3, *outwall4, *outwall5, *outwall6, *outwall7;
    outwall1 = fopen("Subway_wall1.xyz", "w");
    outwall2 = fopen("Subway_wall2.xyz", "w");
    outwall3 = fopen("Subway_wall3.xyz", "w");
    outwall4 = fopen("Subway_wall4.xyz", "w");
    outwall5 = fopen("Subway_wall5.xyz", "w");
    outwall6 = fopen("Subway_wall6.xyz", "w");
    outwall7 = fopen("Subway_wall7.xyz", "w");

    // generate wall atoms
    fprintf(outwall1, "%i\n", 5);
    fprintf(outwall2, "%i\n", 5);
    fprintf(outwall3, "%i\n", 4);
    fprintf(outwall4, "%i\n", 4);
    fprintf(outwall5, "%i\n", 5);
    fprintf(outwall6, "%i\n", 5);
    fprintf(outwall7, "%i\n", 24);

    // wall 1
    int w1x = 10, w1y = 0;
    fprintf(outwall1, "a%i %i  %i  0\n", 0, w1x, w1y);
    fprintf(outwall1, "a%i %i  %i  0\n", 1, w1x, w1y);

    for (i = 1; i <= 5; i++)
    {
        w1x++;

        fprintf(outwall1, "a%i %i  %i  0\n", i, w1x, w1y);
    }
    // wall2
    int w2x = 10, w2y = 0;
    fprintf(outwall2, "a%i %i  %i  0\n", 0, w2x, w2y);
    fprintf(outwall2, "a%i %i  %i  0\n", 1, w2x, w2y);

    for (i = 1; i <= 5; i++)
    {
        w2y++;

        fprintf(outwall2, "a%i %i  %i  0\n", i, w2x, w2y);
    }

    // wall 3
    int w3x = 10, w3y = 7;
    fprintf(outwall3, "a%i %i  %i  0\n", 0, w3x, w3y);
    fprintf(outwall3, "a%i %i  %i  0\n", 1, w3x, w3y);

    for (i = 1; i <= 4; i++)
    {
        w3y++;

        fprintf(outwall3, "a%i %i  %i  0\n", i, w3x, w3y);
    }
    // wall4
    int w4x = 10, w4y = 13;
    fprintf(outwall4, "a%i %i  %i  0\n", 0, w4x, w4y);
    fprintf(outwall4, "a%i %i  %i  0\n", 1, w4x, w4y);

    for (i = 1; i <= 4; i++)
    {
        w4y++;

        fprintf(outwall4, "a%i %i  %i  0\n", i, w4x, w4y);
    }

    // wall 5
    int w5x = 10, w5y = 19;
    fprintf(outwall5, "a%i %i  %i  0\n", 0, w5x, w5y);
    fprintf(outwall5, "a%i %i  %i  0\n", 1, w5x, w5y);

    for (i = 1; i <= 5; i++)
    {
        w5y++;

        fprintf(outwall5, "a%i %i  %i  0\n", i, w5x, w5y);
    }
    // wall6
    int w6x = 10, w6y = 23;
    fprintf(outwall6, "a%i %i  %i  0\n", 0, w6x, w6y);
    fprintf(outwall6, "a%i %i  %i  0\n", 1, w6x, w6y);

    for (i = 1; i <= 5; i++)
    {
        w6x++;

        fprintf(outwall6, "a%i %i  %i  0\n", i, w6x, w6y);
    }

    // wall 7
    int w7x = 14, w7y = 0;
    fprintf(outwall7, "a%i %i  %i  0\n", 0, w7x, w7y);
    fprintf(outwall7, "a%i %i  %i  0\n", 1, w7x, w7y);

    for (i = 1; i <= 24; i++)
    {
        w7y++;

        fprintf(outwall7, "a%i %i  %i  0\n", i, w7x, w7y);
    }

    double test1, test2, testmin = 100;
    int count = 0;
    int count2 = 0;
    printf("constants and structure complete\n");

    // open up file to print

    FILE *outpos1, *outpos2;
    outpos1 = fopen("Subway_rule.xyz", "w");
    outpos2 = fopen("Subway_rule.dat", "w");
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
    printf("speed complete\n");

    // setup desired max speed for each individual
    for (i = 0; i < pedestrian; i++)
    {
        v_0_max[i] = v_max * v_0[i];
    }

    printf("max speed complete\n");

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
    printf("desire direction complete\n");

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
                X[i] = 10 * ran3;
                Y[i] = 23 * ran4;
            }
            else
            {
                ran3 = (float)rand() / RAND_MAX;
                ran4 = (float)rand() / RAND_MAX;
                X[i] = 10 + ran3 * 4;
                Y[i] = 23 * ran4;
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
    printf("initial position complete\n");

    // change the disired direction to prevent trapping and allow alignment for new riders

    for (i = 0; i < pedestrian; i++)
    {
        if (i <= (int)floor(pedestrian / 2))
        {
            if (Y[i] >= 3.9 && Y[i] <= 5.5)
            {
                e_x[i] = 10 - X[i];
                e_y[i] = 0 - Y[i];
                distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                e_x[i] = e_x[i] / distance;
                e_y[i] = e_y[i] / distance;
            }
            else if (Y[i] > 5.5 && Y[i] <= 7.1)
            {
                e_x[i] = 10 - X[i];
                e_y[i] = 10 - Y[i];
                distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                e_x[i] = e_x[i] / distance;
                e_y[i] = e_y[i] / distance;
            }
            else if (Y[i] >= 9.9 && Y[i] <= 11.5)
            {
                e_x[i] = 10 - X[i];
                e_y[i] = 7 - Y[i];
                distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                e_x[i] = e_x[i] / distance;
                e_y[i] = e_y[i] / distance;
            }
            else if (Y[i] > 11.5 && Y[i] <= 13.1)
            {
                e_x[i] = 10 - X[i];
                e_y[i] = 16 - Y[i];
                distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                e_x[i] = e_x[i] / distance;
                e_y[i] = e_y[i] / distance;
            }
            else if (Y[i] >= 15.9 && Y[i] <= 17.5)
            {
                e_x[i] = 10 - X[i];
                e_y[i] = 13 - Y[i];
                distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                e_x[i] = e_x[i] / distance;
                e_y[i] = e_y[i] / distance;
            }
            else if (Y[i] > 17.5 && Y[i] <= 19.1)
            {
                e_x[i] = 10 - X[i];
                e_y[i] = 23 - Y[i];
                distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                e_x[i] = e_x[i] / distance;
                e_y[i] = e_y[i] / distance;
            }
        }
        else
        {
            if (Y[i] >= 0 & Y[i] <= 4)
            {
                e_x[i] = 10 - X[i];
                e_y[i] = 7 - Y[i];
                distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                e_x[i] = e_x[i] / distance;
                e_y[i] = e_y[i] / distance;
            }
            else if (Y[i] >= 7 & Y[i] <= 8.5)
            {
                e_x[i] = 10 - X[i];
                e_y[i] = 4 - Y[i];
                distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                e_x[i] = e_x[i] / distance;
                e_y[i] = e_y[i] / distance;
            }
            else if (Y[i] > 8.5 & Y[i] <= 10)
            {
                e_x[i] = 10 - X[i];
                e_y[i] = 13 - Y[i];
                distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                e_x[i] = e_x[i] / distance;
                e_y[i] = e_y[i] / distance;
            }
            else if (Y[i] >= 13 & Y[i] <= 14.5)
            {
                e_x[i] = 10 - X[i];
                e_y[i] = 10 - Y[i];
                distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                e_x[i] = e_x[i] / distance;
                e_y[i] = e_y[i] / distance;
            }
            else if (Y[i] > 14.5 & Y[i] <= 16)
            {
                e_x[i] = 10 - X[i];
                e_y[i] = 19 - Y[i];
                distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                e_x[i] = e_x[i] / distance;
                e_y[i] = e_y[i] / distance;
            }
            else if (Y[i] >= 19 & Y[i] <= 23)
            {
                e_x[i] = 10 - X[i];
                e_y[i] = 16 - Y[i];
                distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                e_x[i] = e_x[i] / distance;
                e_y[i] = e_y[i] / distance;
            }
        }
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

    for (t = 0; t < 300; t += dt)
    {
        // determine whether all the passengers get out or not
        count = 0;
        count2 = 0;
        for (i = 0; i < pedestrian; i++)
        {
            if (i <= (int)floor(pedestrian / 2))
            {
                if (X[i] > 10.1)
                {
                    count2++;
                }
            }
            else
            {
                if (X[i] < 9)
                {
                    count++;
                }
            }
        }
        // stop when all the new riders come in

        if (count2 == 31 && count == 29)
        {
            printf("the total time take is %f\n", dt * t);
            break;
        }

        fprintf(outpos1, "%i\n", pedestrian - 1);
        for (i = 0; i < pedestrian; i++)
        {
            // setup the desired direction at different circumstances
            // start with people wanna get out
            if (i > (int)floor(pedestrian / 2))
            {
                if (X[i] > 10)
                {
                    if (Y[i] >= 0 & Y[i] <= 4)
                    {
                        e_x[i] = 10 - X[i];
                        e_y[i] = 7 - Y[i];
                        distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                        e_x[i] = e_x[i] / distance;
                        e_y[i] = e_y[i] / distance;
                    }
                    else if (Y[i] >= 7 & Y[i] <= 8.5)
                    {
                        e_x[i] = 10 - X[i];
                        e_y[i] = 4 - Y[i];
                        distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                        e_x[i] = e_x[i] / distance;
                        e_y[i] = e_y[i] / distance;
                    }
                    else if (Y[i] > 8.5 & Y[i] <= 10)
                    {
                        e_x[i] = 10 - X[i];
                        e_y[i] = 13 - Y[i];
                        distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                        e_x[i] = e_x[i] / distance;
                        e_y[i] = e_y[i] / distance;
                    }
                    else if (Y[i] >= 13 & Y[i] <= 14.5)
                    {
                        e_x[i] = 10 - X[i];
                        e_y[i] = 10 - Y[i];
                        distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                        e_x[i] = e_x[i] / distance;
                        e_y[i] = e_y[i] / distance;
                    }
                    else if (Y[i] > 14.5 & Y[i] <= 16)
                    {
                        e_x[i] = 10 - X[i];
                        e_y[i] = 19 - Y[i];
                        distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                        e_x[i] = e_x[i] / distance;
                        e_y[i] = e_y[i] / distance;
                    }
                    else if (Y[i] >= 19 & Y[i] <= 23)
                    {
                        e_x[i] = 10 - X[i];
                        e_y[i] = 16 - Y[i];
                        distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                        e_x[i] = e_x[i] / distance;
                        e_y[i] = e_y[i] / distance;
                    }
                    else
                    {
                        e_x[i] = -1;
                        e_y[i] = 0;
                    }
                }
                else
                {
                    e_x[i] = -1;
                    e_y[i] = 0;
                }
            }
            // then choose the people who wants to get in
            else
            {
                // when there is still passenger inside of the train
                // stick to the walls
                if (count != 29)
                {
                    if (3.9 <= Y[i] && 5.5 >= Y[i])
                    {
                        e_x[i] = 10 - X[i];
                        e_y[i] = 0 - Y[i];
                        distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                        e_x[i] = e_x[i] / distance;
                        e_y[i] = e_y[i] / distance;
                    }
                    else if (5.5 < Y[i] && 7.1 >= Y[i])
                    {
                        e_x[i] = 10 - X[i];
                        e_y[i] = 10 - Y[i];
                        distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                        e_x[i] = e_x[i] / distance;
                        e_y[i] = e_y[i] / distance;
                    }
                    else if (9.9 <= Y[i] && 11.5 >= Y[i])
                    {
                        e_x[i] = 10 - X[i];
                        e_y[i] = 7 - Y[i];
                        distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                        e_x[i] = e_x[i] / distance;
                        e_y[i] = e_y[i] / distance;
                    }
                    else if (11.5 < Y[i] && 13.1 >= Y[i])
                    {
                        e_x[i] = 10 - X[i];
                        e_y[i] = 16 - Y[i];
                        distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                        e_x[i] = e_x[i] / distance;
                        e_y[i] = e_y[i] / distance;
                    }
                    else if (15.9 <= Y[i] && 17.5 >= Y[i])
                    {
                        e_x[i] = 10 - X[i];
                        e_y[i] = 13 - Y[i];
                        distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                        e_x[i] = e_x[i] / distance;
                        e_y[i] = e_y[i] / distance;
                    }
                    else if (17.5 < Y[i] && 19.1 >= Y[i])
                    {
                        e_x[i] = 10 - X[i];
                        e_y[i] = 23 - Y[i];
                        distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                        e_x[i] = e_x[i] / distance;
                        e_y[i] = e_y[i] / distance;
                    }
                    else
                    {
                        e_x[i] = 1;
                        e_y[i] = 0;
                    }
                }
                // when there is no passengers in the train anymore
                else
                {
                    // when still outside of the train
                    if (X[i] < 10)
                    {
                        if (Y[i] > 4.1 && Y[i] < 6.9)
                        {
                            e_x[i] = 1;
                            e_y[i] = 0;
                        }
                        else if (Y[i] > 10.1 && Y[i] < 12.9)
                        {
                            e_x[i] = 1;
                            e_y[i] = 0;
                        }
                        else if (Y[i] > 16.1 && Y[i] < 18.9)
                        {
                            e_x[i] = 1;
                            e_y[i] = 0;
                        }
                        else if (Y[i] >= 0 && Y[i] <= 4.1)
                        {
                            e_x[i] = 10 - X[i];
                            e_y[i] = 7 - Y[i];
                            distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                            e_x[i] = e_x[i] / distance;
                            e_y[i] = e_y[i] / distance;

                            if (X[i] >= 7)
                            {
                                e_x[i] = 0;
                                e_y[i] = 1;
                            }
                        }
                        else if (Y[i] >= 6.9 && Y[i] <= 8.5)
                        {
                            e_x[i] = 10 - X[i];
                            e_y[i] = 4 - Y[i];
                            distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                            e_x[i] = e_x[i] / distance;
                            e_y[i] = e_y[i] / distance;
                            if (X[i] >= 7)
                            {
                                e_x[i] = 0;
                                e_y[i] = -1;
                            }
                        }
                        else if (Y[i] > 8.5 && Y[i] <= 10.1)
                        {
                            e_x[i] = 10 - X[i];
                            e_y[i] = 13 - Y[i];
                            distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                            e_x[i] = e_x[i] / distance;
                            e_y[i] = e_y[i] / distance;
                            if (X[i] >= 7)
                            {
                                e_x[i] = 0;
                                e_y[i] = 1;
                            }
                        }
                        else if (Y[i] >= 12.9 && Y[i] <= 14.5)
                        {
                            e_x[i] = 10 - X[i];
                            e_y[i] = 10 - Y[i];
                            distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                            e_x[i] = e_x[i] / distance;
                            e_y[i] = e_y[i] / distance;
                            if (X[i] >= 7)
                            {
                                e_x[i] = 0;
                                e_y[i] = -1;
                            }
                        }
                        else if (Y[i] > 14.5 && Y[i] <= 16.1)
                        {
                            e_x[i] = 10 - X[i];
                            e_y[i] = 19 - Y[i];
                            distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                            e_x[i] = e_x[i] / distance;
                            e_y[i] = e_y[i] / distance;
                            if (X[i] >= 7)
                            {
                                e_x[i] = 0;
                                e_y[i] = 1;
                            }
                        }
                        else if (Y[i] >= 18.9 && Y[i] <= 23)
                        {
                            e_x[i] = 10 - X[i];
                            e_y[i] = 16 - Y[i];
                            distance = sqrt(e_x[i] * e_x[i] + e_y[i] * e_y[i]);
                            e_x[i] = e_x[i] / distance;
                            e_y[i] = e_y[i] / distance;
                            if (X[i] >= 7)
                            {
                                e_x[i] = 0;
                                e_y[i] = -1;
                            }
                        }
                    }
                    // when inside of the train
                    else
                    {
                        e_x[i] = 1;
                        e_y[i] = 0;
                    }
                }
            }

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

                    dot_p = e_x[i] * (Fj_nei_x)+ e_y[i] * (Fj_nei_y);
                    angle = sqrt(Fj_nei_x * Fj_nei_x + Fj_nei_y * Fj_nei_y) * cos(5 / 9 * PI);
                    if (dot_p < angle)
                    {
                        Fj_nei_x = c * Fj_nei_x;
                        Fj_nei_y = c * Fj_nei_y;
                    }
                    if (i <= (int)floor(pedestrian / 2) && j <= (int)floor(pedestrian / 2))
                    {

                        if (X[i] <= 10 && X[j] > 10)
                        {
                            Fj_nei_x = -1 * Fj_nei_x;
                            Fj_nei_y = -1 * Fj_nei_y;
                        }
                    }
                    else if (i > (int)floor(pedestrian / 2) && j > (int)floor(pedestrian / 2))
                    {
                        if (X[i] >= 10 && X[j] < 10)
                        {
                            Fj_nei_x = -1 * Fj_nei_x;
                            Fj_nei_y = -1 * Fj_nei_y;
                        }
                    }

                    Fsum_nei_x = Fsum_nei_x + Fj_nei_x;
                    Fsum_nei_y = Fsum_nei_y + Fj_nei_y;
                }
            }

            F_rep_x[i] = Fsum_nei_x;
            F_rep_y[i] = Fsum_nei_y;

            // wall repulson force

            double wall1_x = 0, wall1_y = 0, wall2_x = 0, wall2_y = 0, wall3_x = 0, wall3_y = 0, wall4_x = 0, wall4_y = 0;
            double wall5_x = 0, wall5_y = 0, wall6_x = 0, wall6_y = 0, wall7_x = 0, wall7_y = 0;
            double wall8_x = 0, wall8_y = 0, wall9_x = 0, wall9_y = 0, wall10_x = 0, wall10_y = 0;
            double wall11_x = 0, wall11_y = 0, wall12_x = 0, wall12_y = 0, wall13_x = 0, wall13_y = 0;

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

            struct coordWallpoint *Wall3_poi = Nearest_Point(wall3_start[0], wall3_start[1],
                                                             wall3_end[0], wall3_end[1], X[i], Y[i]);
            wall3_x = Wall3_poi->x_wall;
            wall3_y = Wall3_poi->y_wall;

            struct coordWallpoint *Wall4_poi = Nearest_Point(wall4_start[0], wall4_start[1],
                                                             wall4_end[0], wall4_end[1], X[i], Y[i]);
            wall4_x = Wall4_poi->x_wall;
            wall4_y = Wall4_poi->y_wall;

            struct coordWallpoint *Wall5_poi = Nearest_Point(wall5_start[0], wall5_start[1],
                                                             wall5_end[0], wall5_end[1], X[i], Y[i]);
            wall5_x = Wall5_poi->x_wall;
            wall5_y = Wall5_poi->y_wall;

            struct coordWallpoint *Wall6_poi = Nearest_Point(wall6_start[0], wall6_start[1],
                                                             wall6_end[0], wall6_end[1], X[i], Y[i]);
            wall6_x = Wall6_poi->x_wall;
            wall6_y = Wall6_poi->y_wall;

            struct coordWallpoint *Wall7_poi = Nearest_Point(wall7_start[0], wall7_start[1],
                                                             wall7_end[0], wall7_end[1], X[i], Y[i]);
            wall7_x = Wall7_poi->x_wall;
            wall7_y = Wall7_poi->y_wall;

            struct coordWallpoint *Wall8_poi = Nearest_Point(wall8_start[0], wall8_start[1],
                                                             wall8_end[0], wall8_end[1], X[i], Y[i]);
            wall8_x = Wall8_poi->x_wall;
            wall8_y = Wall8_poi->y_wall;

            struct coordWallpoint *Wall9_poi = Nearest_Point(wall9_start[0], wall9_start[1],
                                                             wall9_end[0], wall9_end[1], X[i], Y[i]);
            wall9_x = Wall9_poi->x_wall;
            wall9_y = Wall9_poi->y_wall;

            struct coordWallpoint *Wall10_poi = Nearest_Point(wall10_start[0], wall10_start[1],
                                                              wall10_end[0], wall10_end[1], X[i], Y[i]);
            wall10_x = Wall10_poi->x_wall;
            wall10_y = Wall10_poi->y_wall;

            struct coordWallpoint *Wall11_poi = Nearest_Point(wall11_start[0], wall11_start[1],
                                                              wall11_end[0], wall11_end[1], X[i], Y[i]);
            wall11_x = Wall11_poi->x_wall;
            wall11_y = Wall11_poi->y_wall;

            struct coordWallpoint *Wall12_poi = Nearest_Point(wall12_start[0], wall12_start[1],
                                                              wall12_end[0], wall12_end[1], X[i], Y[i]);
            wall12_x = Wall12_poi->x_wall;
            wall12_y = Wall12_poi->y_wall;

            struct coordWallpoint *Wall13_poi = Nearest_Point(wall13_start[0], wall13_start[1],
                                                              wall13_end[0], wall13_end[1], X[i], Y[i]);
            wall13_x = Wall13_poi->x_wall;
            wall13_y = Wall13_poi->y_wall;

            struct coordWallF *Wall1F_poi = WallRep_F(wall1_x, wall1_y, X[i], Y[i]);
            F_wall1_x[i] = Wall1F_poi->WallRepForce_x;
            F_wall1_y[i] = Wall1F_poi->WallRepForce_y;
            dot_p = e_x[i] * F_wall1_x[i] + e_y[i] * F_wall1_y[i];

            struct coordWallF *Wall2F_poi = WallRep_F(wall2_x, wall2_y, X[i], Y[i]);
            F_wall2_x[i] = Wall2F_poi->WallRepForce_x;
            F_wall2_y[i] = Wall2F_poi->WallRepForce_y;
            dot_p = e_x[i] * F_wall2_x[i] + e_y[i] * F_wall2_y[i];

            struct coordWallF *Wall3F_poi = WallRep_F(wall3_x, wall3_y, X[i], Y[i]);
            F_wall3_x[i] = Wall3F_poi->WallRepForce_x;
            F_wall3_y[i] = Wall3F_poi->WallRepForce_y;
            dot_p = e_x[i] * F_wall3_x[i] + e_y[i] * F_wall3_y[i];

            struct coordWallF *Wall4F_poi = WallRep_F(wall4_x, wall4_y, X[i], Y[i]);
            F_wall4_x[i] = Wall4F_poi->WallRepForce_x;
            F_wall4_y[i] = Wall4F_poi->WallRepForce_y;
            dot_p = e_x[i] * F_wall4_x[i] + e_y[i] * F_wall4_y[i];

            struct coordWallF *Wall5F_poi = WallRep_F(wall5_x, wall5_y, X[i], Y[i]);
            F_wall5_x[i] = Wall5F_poi->WallRepForce_x;
            F_wall5_y[i] = Wall5F_poi->WallRepForce_y;
            dot_p = e_x[i] * F_wall5_x[i] + e_y[i] * F_wall5_y[i];

            struct coordWallF *Wall6F_poi = WallRep_F(wall6_x, wall6_y, X[i], Y[i]);
            F_wall6_x[i] = Wall6F_poi->WallRepForce_x;
            F_wall6_y[i] = Wall6F_poi->WallRepForce_y;
            dot_p = e_x[i] * F_wall6_x[i] + e_y[i] * F_wall6_y[i];

            struct coordWallF *Wall7F_poi = WallRep_F(wall7_x, wall7_y, X[i], Y[i]);
            F_wall7_x[i] = Wall7F_poi->WallRepForce_x;
            F_wall7_y[i] = Wall7F_poi->WallRepForce_y;
            dot_p = e_x[i] * F_wall7_x[i] + e_y[i] * F_wall7_y[i];

            // invisible walls set up
            F_wall8_x[i] = 0;
            F_wall8_y[i] = 0;
            F_wall9_x[i] = 0;
            F_wall9_y[i] = 0;
            F_wall10_x[i] = 0;
            F_wall10_y[i] = 0;
            F_wall11_x[i] = 0;
            F_wall11_y[i] = 0;
            F_wall12_x[i] = 0;
            F_wall12_y[i] = 0;
            F_wall13_x[i] = 0;
            F_wall13_y[i] = 0;

            if (i <= (int)floor(pedestrian / 2) && count != 29)
            {
                if ((Y[i] >= 0 && Y[i] <= 4) || (Y[i] >= 7 && Y[i] <= 10) || (Y[i] >= 13 && Y[i] <= 16) || (Y[i] >= 19 && Y[i] <= 23))
                {
                    struct coordWallF *Wall8F_poi = WallRep_F(wall8_x, wall8_y, X[i], Y[i]);
                    F_wall8_x[i] = Wall8F_poi->WallRepForce_x;
                    F_wall8_y[i] = Wall8F_poi->WallRepForce_y;
                    dot_p = e_x[i] * F_wall8_x[i] + e_y[i] * F_wall8_y[i];

                    struct coordWallF *Wall9F_poi = WallRep_F(wall9_x, wall9_y, X[i], Y[i]);
                    F_wall9_x[i] = Wall9F_poi->WallRepForce_x;
                    F_wall9_y[i] = Wall9F_poi->WallRepForce_y;
                    dot_p = e_x[i] * F_wall9_x[i] + e_y[i] * F_wall9_y[i];

                    struct coordWallF *Wall10F_poi = WallRep_F(wall10_x, wall10_y, X[i], Y[i]);
                    F_wall10_x[i] = Wall10F_poi->WallRepForce_x;
                    F_wall10_y[i] = Wall10F_poi->WallRepForce_y;
                    dot_p = e_x[i] * F_wall10_x[i] + e_y[i] * F_wall10_y[i];

                    struct coordWallF *Wall11F_poi = WallRep_F(wall11_x, wall11_y, X[i], Y[i]);
                    F_wall11_x[i] = Wall11F_poi->WallRepForce_x;
                    F_wall11_y[i] = Wall11F_poi->WallRepForce_y;
                    dot_p = e_x[i] * F_wall11_x[i] + e_y[i] * F_wall11_y[i];

                    struct coordWallF *Wall12F_poi = WallRep_F(wall12_x, wall12_y, X[i], Y[i]);
                    F_wall12_x[i] = Wall12F_poi->WallRepForce_x;
                    F_wall12_y[i] = Wall12F_poi->WallRepForce_y;
                    dot_p = e_x[i] * F_wall12_x[i] + e_y[i] * F_wall12_y[i];

                    struct coordWallF *Wall13F_poi = WallRep_F(wall13_x, wall13_y, X[i], Y[i]);
                    F_wall13_x[i] = Wall13F_poi->WallRepForce_x;
                    F_wall13_y[i] = Wall13F_poi->WallRepForce_y;
                    dot_p = e_x[i] * F_wall13_x[i] + e_y[i] * F_wall13_y[i];
                }
            }
            
            Fsum_wall_x[i] = F_wall1_x[i] + F_wall2_x[i] + F_wall3_x[i] + F_wall4_x[i] + F_wall5_x[i] + F_wall6_x[i] + F_wall7_x[i] + F_wall8_x[i] + F_wall9_x[i] + F_wall10_x[i] + F_wall11_x[i] + F_wall12_x[i] + F_wall13_x[i];
            Fsum_wall_y[i] = F_wall1_y[i] + F_wall2_y[i] + F_wall3_y[i] + F_wall4_y[i] + F_wall5_y[i] + F_wall6_y[i] + F_wall7_y[i] + F_wall8_y[i] + F_wall9_y[i] + F_wall10_y[i] + F_wall11_y[i] + F_wall12_y[i] + F_wall13_y[i];

            F_total_x[i] = F_0_x[i] + F_rep_x[i] + Fsum_wall_x[i];
            F_total_y[i] = F_0_y[i] + F_rep_y[i] + Fsum_wall_y[i];

            v_x[i] += F_total_x[i] * dt;
            v_y[i] += F_total_y[i] * dt;

            if (sqrt(v_x[i] * v_x[i] + v_y[i] * v_y[i]) >= v_0_max[i])
            {
                v_x[i] = v_x[i] / sqrt(v_x[i] * v_x[i] + v_y[i] * v_y[i]) * v_0_max[i];
                v_y[i] = v_y[i] / sqrt(v_x[i] * v_x[i] + v_y[i] * v_y[i]) * v_0_max[i];
            }

            X[i] = X[i] + v_x[i] * dt;
            Y[i] = Y[i] + v_y[i] * dt;

            // printf("The new position of [%i] is x=%f,y=%f \n", i, X[i], Y[i]);

            fprintf(outpos1, "a%i %f  %f  0\n", i, X[i], Y[i]);
            fprintf(outpos2, "%i %f  %f  0\n", i, X[i], Y[i]);
        }
    }
}
