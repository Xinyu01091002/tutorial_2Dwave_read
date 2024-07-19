# -*- coding: utf-8 -*-
# /**
# # Breaking wave
# 
# We solve the two-phase Navier--Stokes equations with surface tension
# and using a momentum-conserving transport of each phase. Gravity is
# taken into account using the "reduced gravity approach" and the
# results are visualised using Basilisk view. */
# 
# #include "navier-stokes/centered.h"
# #include "two-phase.h"
# #include "navier-stokes/conserving.h"
# #include "tension.h"
# #include "reduced.h"
# #include "view.h"
# #include "tag.h"
# 
# /**
# We log some profiling information. */
# 
# #include "navier-stokes/perfs.h"
# #include "profiling.h"
# 
# /**
# The primary parameters are the wave steepness $ak$, the Bond and
# Reynolds numbers. */
# 
# // double ak = 0.55;
# // double BO = 1000.;
# // double RE = 40000.;
# 
# /**
# The default maximum level of refinement depends on the dimension. */
# 
# int LEVEL = dimension == 2 ? 9 : 6;
# 
# /**
# The error on the components of the velocity field used for adaptive
# refinement. */
# 
# double uemax = 0.005;
# 
# /**
# The density and viscosity ratios are those of air and water. */
# 
# // #define RATIO (1. / 850.)
# // #define MURATIO (17.4e-6 / 8.9e-4)
# 
# /**
# Define if we want to use a Dirac viscous layer initialization. */
# 
# int DIRAC = 0;
# 
# /**
# The wave number, fluid depth and acceleration of gravity are set to
# these values. */
# 
# // #define k_ (2. * pi)
# // #define h_ 0.5
# // #define g_ 1.
# 
# /**
# The program takes optional arguments which are the level of
# refinement, steepness, Bond and Reynolds numbers, and optional Dirac
# initialisation. */
# 
# int main(int argc, char *argv[])
# {
#   // if (argc > 1)
#   //   LEVEL = atoi(argv[1]);
#   // if (argc > 2)
#   //   ak = atof(argv[2]);
#   // if (argc > 3)
#   //   BO = atof(argv[3]);
#   // if (argc > 4)
#   //   RE = atof(argv[4]);
#   // if (argc > 5)
#   //   DIRAC = atof(argv[5]);
# 
#   /**
#   The domain is a cubic box centered on the origin and of length
#   , periodic in the x- and z-directions. */
#   L0 = 1;
#   size(L0);
#   origin(-L0 / 2, -L0 / 2, -L0 / 2);
#   // origin (-L0/2, -L0/2);
#   //  size(L3);
#   //  origin(-L3/2,-L3/2,-L3/2);
#   periodic(right);
# #if dimension > 2
#   periodic(front);
# #endif
# 
#   /**
#   Here we set the densities and viscosities corresponding to the
#   parameters above. */
#   rho1 = 1000;
#   rho2 = 1.293;
#   mu1 = 0.001002;
#   mu2 = 1.825 * 1e-5;
#   f.sigma = 0.0728;
#   G.y = -9.8;
#   // rho1 = 1.;
#   // rho2 = RATIO;
#   // mu1 = 1.0 / RE; // using wavelength as length scale
#   // mu2 = 1.0 / RE * MURATIO;
#   // f.sigma = 1. / (BO * sq(k_));
#   // // f.sigma = 1. / (BO * sq(1));
#   // G.y = -g_/(2*pi);
#   // // G.y = 0;
#   // // Modify the vertical metric scale
# 
#   // /**
#   // When we use adaptive refinement, we start with a coarse mesh which
#   // will be refined as required when initialising the wave. */
# 
#   // #if TREE
#   //   N = 32;
#   // #else
#   //   N = 1 << LEVEL;
#   // #endif
#   N = 256;
#   // init_grid (N);
#   run();
# }
# 
# double my_eta(double x, double y, double *eta_s)
# {
#   int index_x;
#   double Delta = 0;
#   Delta = L0 / N;
#   index_x = round((x + L0 / 2) / Delta);
#   return *(eta_s + index_x) - y;
# }
# 
# event init(i = 0)
# {
# 
#   if (!restore("restart"))
#   {
#     FILE *fp;
#     float eta_posi, phi_posi;
# 
#     int max_n = N + 1;
#     double *eta_s = (double *)malloc(max_n * sizeof(double));
#     double *phi_s = (double *)malloc(max_n * sizeof(double));
# 
#     char filename[100];
#     sprintf(filename, "stoke3/Stoke3_%d.txt", N);
#     fp = fopen(filename, "r");
# 
#     // 读取文件内容到数组中
#     int read_index = 0;
#     while (fscanf(fp, "%f %f", &eta_posi, &phi_posi) != EOF)
#     {
#       if (read_index >= max_n)
#       {
#         fprintf(stderr, "Exceeded maximum number of elements: %d\n", max_n);
#         break;
#       }
#       eta_s[read_index] = (double)eta_posi;
#       phi_s[read_index] = (double)phi_posi;
#       read_index++;
#     }
# 
#     // 关闭文件
#     fclose(fp);
# 
#     // 处理完数据后，释放内存
# 
#     do
#     {
# 
#       // Try to initial the phase field with personalized input
# 
#       scalar my_f[];
#       scalar my_phi[];
#       foreach ()
#       {
#         static int index_x = 1;
#         index_x = round((x + L0 / 2) / Delta_x); // Remember to change when L0 is not 1
#         my_f[] = eta_s[index_x] - y - L0 / 2;
#         my_phi[] = phi_s[index_x];
#       }
#       foreach ()
#       {
#         static int index_x = 1;
#         index_x = round((x + L0 / 2) / Delta_x); // Remember to change when L0 is not 1
#         static int index_y = 1;
#         index_x = round((y + L0 / 2) / Delta_y); // Remember to change when L0 is not 1
#         static double printf_x = 0;
#         static double printf_y = 0;
#         printf_x = my_f[];
#       }
# 
#       // fraction(f, my_f[]);
#       fraction(f, my_eta(x, y, eta_s));
#       foreach ()
#       {
#         if (f[] == 1 || f[] == 0)
#         {
#           my_phi[] = 0;
#         }
#       }
#       boundary({my_phi});
#       my_phi[bottom] = dirichlet(0.);
#       my_phi[left] = neumann(0.);
#       my_phi[right] = neumann(0.);
#       scalar residual_J[];
#       foreach ()
#       {
#         residual_J[] = 0;
#       }
#       do
#       {
#         foreach ()
#         {
#           if (f[] ==1)
#           {
#             // my_phi[] = (my_phi[1, 0]*f[1,0] + my_phi[-1, 0]*f[-1,0] + my_phi[0, 1]*f[0,1] + my_phi[0, -1]*f[0,-1]) / 4;
#             my_phi[] = (my_phi[1, 0] + my_phi[-1, 0] + my_phi[0, 1] + my_phi[0, -1]) / 4;
#           }
#         }
#         foreach ()
#         {
#           if (f[] ==1)
#           {
#             residual_J[] = (my_phi[1, 0] + my_phi[-1, 0] + my_phi[0, 1] + my_phi[0, -1] - 4 * my_phi[0, 0]);
#           }
#         }
#       } while (normf(residual_J).max / normf(my_phi).max > 1e-3);
#       // Smooth the boundary phi
# 
#       foreach ()
#       {
#         foreach_dimension()
#         {
#           u.x[] = (my_phi[1] - my_phi[-1]) / (2.0 * Delta) * f[];
#           // u.y[] = (my_phi[0,1] - my_phi[0,-1]) / (2.0 * Delta) * f[];
#         }
#       }
#       foreach ()
#       // foreach_dimension()
#       {
#         if (f[] > 0 && f[] < 1 && f[0, -1] == 1)
#         {
#           u.x[] = u.x[0, -1] * f[];
#           u.y[] = u.y[0, -1] * f[];
#         }
#         if (f[0, -1] > 0 && f[0, -1] < 1)
#         {
#           u.x[] = u.x[0, -1] * f[];
#           u.y[] = u.y[0, -1] * f[];
#         }
#       }
#       foreach ()
#         foreach_dimension()
#         {
#           if (f[] > 0)
#           {
#             if (u.x[] > u.x[1, 0] && u.x[] > u.x[0, 1] && u.x[] > u.x[0, -1] && u.x[] > u.x[-1, 0])
#             {
#               u.x[] = (u.x[1, 0] + u.x[-1, 0] + u.x[0, 1] + u.x[0, -1]) / 4;
#             }
#           }
#         }
#       boundary((scalar *){u});
#       dump();
# 
#     }
# 
#     /**
#     On trees, we repeat this initialisation until mesh adaptation does
#     not refine the mesh anymore. */
# 
#     // #if TREE
#     //     while (adapt_wavelet({f, u},
#     //                          (double[]){0.02, uemax, uemax, uemax}, LEVEL, 5)
#     //                .nf);
#     // #else
#     //     while (0);
#     // #endif
#     while (0);
#   }
# }
# 
# event movies(t += 0.002)
# {
#   {
#     static FILE *fp = popen("ppm2mp4 f.mp4", "w");
#     output_ppm(f, fp, min = 0, max = 1, n = 512);
#   }
# 
# #if TREE
#   {
#     scalar l[];
#     foreach ()
#       l[] = level;
#     static FILE *fp = popen("ppm2mp4 level.mp4", "w");
#     output_ppm(l, fp, min = 5, max = LEVEL, n = 512);
#   }
# #endif
# }
# event snapshot(i += 200)
# {
#   char name[80];
#   sprintf(name, "dump-%03d", (int)i);
#   dump(file = name);
# }
# 
# /**
# ## End
# 
# The wave period is `k_/sqrt(g_*k_)`. We want to run up to 2
# (alternatively 4) periods. */
# 
# event end(t = 2)
# {
#   fprintf(fout, "i = %d t = %g\n", i, t);
#   dump("end");
# }
# 
# /**
# ## Mesh adaptation
# 
# On trees, we adapt the mesh according to the error on volume fraction
# and velocity. */
# 
# #if TREE
# event adapt(i++)
# {
#   adapt_wavelet({f, u}, (double[]){0.01, uemax, uemax, uemax}, LEVEL, 5);
# }
# #endif
# 
# /**
# ## Running in parallel
# 
# This file will work in 2D or 3D, either with parallel multigrid
# (without adaptivity), using for example:
# 
# ~~~bash
# qcc -source -D_MPI=1 -grid=multigrid3D wave.c
# scp _wave.c occigen.cines.fr:
# ~~~
# 
# and then following a recipe similar to that of the
# [atomisation](/src/examples/atomisation.c#on-occigen) example.
# 
# To use adaptivity, just do something like:
# 
# ~~~bash
# qcc -source -D_MPI=1 -grid=octree wave.c
# scp _wave.c occigen.cines.fr:
# ~~~
# */
