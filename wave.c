/**
# Breaking wave

We solve the two-phase Navier--Stokes equations with surface tension
and using a momentum-conserving transport of each phase. Gravity is
taken into account using the "reduced gravity approach" and the
results are visualised using Basilisk view. */
// #include "metric.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "view.h"
#include "tag.h"
#include "vtk.h"

/**
We log some profiling information. */

#include "navier-stokes/perfs.h"
#include "profiling.h"

/**
The primary parameters are the wave steepness $ak$, the Bond and
Reynolds numbers. */

double ak = 0.55;
double BO = 1000.;
double RE = 40000.;

/**
The default maximum level of refinement depends on the dimension. */

int LEVEL = dimension == 2 ? 9 : 6;

/**
The error on the components of the velocity field used for adaptive
refinement. */

double uemax = 0.005;

/**
The density and viscosity ratios are those of air and water. */

#define RATIO (1. / 850.)
#define MURATIO (17.4e-6 / 8.9e-4)

/**
Define if we want to use a Dirac viscous layer initialization. */

int DIRAC = 0;

/**
The wave number, fluid depth and acceleration of gravity are set to
these values. */

#define k_ (2. * pi)
#define h_ 0.5
#define g_ 1.
#define y_metric 1
/**
The program takes optional arguments which are the level of
refinement, steepness, Bond and Reynolds numbers, and optional Dirac
initialisation. */

int main(int argc, char *argv[])
{
  // if (argc > 1)
  //   LEVEL = atoi(argv[1]);
  // if (argc > 2)
  //   ak = atof(argv[2]);
  // if (argc > 3)
  //   BO = atof(argv[3]);
  // if (argc > 4)
  //   RE = atof(argv[4]);
  // if (argc > 5)
  //   DIRAC = atof(argv[5]);

  /**
  The domain is a cubic box centered on the origin and of length
  , periodic in the x- and z-directions. */
// Dimensionless group velocity 
double cp=1.56;
// Modify the metric factor
face vector fmv=fm;
scalar cmv=cm;
const face vector unityf1[] = {1.[0],y_metric[0],1.[0]};
const scalar unity1[] = y_metric[0];
fm.x=unityf1.x;
fm.y=unityf1.y;
cm=unity1;
  L0 = 2;
  size(L0);
  origin(-L0 / 2, -L0 / 2, -L0 / 2);
  // origin (-L0/2, -L0/2);
  //  size(L3);
  //  origin(-L3/2,-L3/2,-L3/2);
  periodic(right);
    //   u.n[left]=neumann(0.);
    // u.n[right]=neumann(0.);
    //   u.t[left]=neumann(0.);
    // u.t[right]=neumann(0.);
#if dimension > 2
  periodic(front);
#endif

  /**
  Here we set the densities and viscosities corresponding to the
  parameters above. */
  // rho1 = 1000;
  // rho2 = 1.17647;
  // // rho2 = 0.5;
  // mu1 = 17.4e-6 ;
  // mu2 = 8.9e-4;
  // // mu2=1e-10;
  // f.sigma = 0.0728;
  // G.y = -9.8;
  rho1 = 1.;
  rho2 = RATIO;
  mu1 = 1.0 / RE; // using wavelength as length scale
  mu2 = 1/ RE * MURATIO;
  f.sigma = 1. / (BO * sq(k_));
  G.y = -g_;
  // // Modify the vertical metric scale

  // /**
  // When we use adaptive refinement, we start with a coarse mesh which
  // will be refined as required when initialising the wave. */

  // #if TREE
  //   N = 32;
  // #else
  //   N = 1 << LEVEL;
  // #endif
  N = 512;

  run();
}

double my_eta(double x, double y, double *eta_s)
{
  int index_x;
  double Delta = 0;
  Delta = L0 / N;
  index_x = round((x + L0 / 2) / Delta);
  return *(eta_s + index_x) - y;
}

event init(i = 0)
{

  if (!restore("restart"))
  {
    FILE *fp;
    float eta_posi, phi_posi;

    int max_n = N + 1;
    double *eta_s = (double *)malloc(max_n * sizeof(double));
    double *phi_s = (double *)malloc(max_n * sizeof(double));

    char filename[100];
    sprintf(filename, "stoke3/Stoke3_%d.txt", N);
    // sprintf(filename, "OceanWave3D.init.512");
    // sprintf(filename, "OceanWave3D.init");
    fp = fopen(filename, "r");

    // 读取文件内容到数组中
    int read_index = 0;
    while (fscanf(fp, "%f %f", &eta_posi, &phi_posi) != EOF)
    {
      if (read_index >= max_n)
      {
        fprintf(stderr, "Exceeded maximum number of elements: %d\n", max_n);
        break;
      }
      eta_s[read_index] = (double)eta_posi;
      phi_s[read_index] = (double)phi_posi;
      read_index++;
    }

    // 关闭文件
    fclose(fp);

    // 处理完数据后，释放内存

    do
    {

      // Try to initial the phase field with personalized input

      scalar my_f[];
      scalar my_phi[];
      foreach ()
      {
        static int index_x = 1;
        index_x = round((x + L0 / 2) / Delta_x); // Remember to change when L0 is not 1
        my_f[] = eta_s[index_x]- y*fm.y[];
        my_phi[] = phi_s[index_x];
      }
      // foreach ()
      // {
      //   static int index_x = 1;
      //   index_x = round((x + L0 / 2) / Delta_x); // Remember to change when L0 is not 1
      //   static int index_y = 1;
      //   index_x = round((y + L0 / 2) / Delta_y); // Remember to change when L0 is not 1
      //   static double printf_x = 0;
      //   static double printf_y = 0;
      //   printf_x = my_f[];
      // }

      fraction(f, my_f[]);
      // fraction(f, my_eta(x, y, eta_s));
      foreach ()
      {
        if (f[] == 1 || f[] == 0)
        {
          my_phi[] = 0;
        }
      }
      boundary({my_phi});
      my_phi[bottom] = dirichlet(0.);
      scalar residual_J[];
      foreach ()
      {
        residual_J[] = 0;
      }
      do
      {
        foreach ()
        {
          if (f[] == 1)
          {
            // my_phi[] = (my_phi[1, 0]*f[1,0] + my_phi[-1, 0]*f[-1,0] + my_phi[0, 1]*f[0,1] + my_phi[0, -1]*f[0,-1]) / 4;
            my_phi[] = (my_phi[1, 0] + my_phi[-1, 0] + my_phi[0, 1] + my_phi[0, -1]) / 4;
          }
        }
        foreach ()
        {
          if (f[] == 1)
          {
            residual_J[] = (my_phi[1, 0] + my_phi[-1, 0] + my_phi[0, 1] + my_phi[0, -1] - 4 * my_phi[0, 0]);
          }
        }
      } while (normf(residual_J).max / normf(my_phi).max > 1e-6);
      // Smooth the boundary phi

      foreach ()
      {
        foreach_dimension()
        {
          u.x[] = (my_phi[1] - my_phi[-1]) / (2.0 * Delta*fm.x[]) * f[];
          // u.x[] = (my_phi[1] - my_phi[-1]) / (2.0 * Delta )* f[];
          // u.y[] = (my_phi[0,1] - my_phi[0,-1]) / (2.0 * Delta) * f[];
        }
      }
      foreach ()
      // foreach_dimension()
      {
        if (f[] > 0 && f[] < 1 && f[0, -1] == 1)
        {
          u.x[] = u.x[0, -1] * f[];
          u.y[] = u.y[0, -1] * f[];
        }
        if (f[0, -1] > 0 && f[0, -1] < 1)
        {
          u.x[] = u.x[0, -1] * f[];
          u.y[] = u.y[0, -1] * f[];
        }
      }
      foreach ()
        foreach_dimension()
        {
          if (f[] > 0)
          {
            if (u.x[] > u.x[1, 0] && u.x[] > u.x[0, 1] && u.x[] > u.x[0, -1] && u.x[] > u.x[-1, 0])
            {
              u.x[] = (u.x[1, 0] + u.x[-1, 0] + u.x[0, 1] + u.x[0, -1]) / 4;
            }
          }
        }
      boundary((scalar *){u});
      dump();

    }

    /**
    On trees, we repeat this initialisation until mesh adaptation does
    not refine the mesh anymore. */

    #if TREE
        while (adapt_wavelet({f, u},
                             (double[]){0.01, uemax, uemax/y_metric, uemax}, LEVEL, 5)
                   .nf);
    #else
        while (0);
    #endif
    // while (0);
  }
}

// event set_air_velocity_0 (i++)
// {
//   foreach(){
//     if (f[]==0)
//     {
//       u.x[]=0;
//       u.y[]=0;
//     }
//   }
// }

event movies(t+=0.05)
{
  {
    static FILE *fp = popen("ppm2mp4 f.mp4", "w");
    output_ppm(f, fp, min = 0, max = 1, n = 512);
  }

// #if TREE
//   {
//     scalar l[];
//     foreach ()
//       l[] = level;
//     static FILE *fp = popen("ppm2mp4 level.mp4", "w");
//     output_ppm(l, fp, min = 5, max = LEVEL, n = 512);
//   }
// #endif
}


event my_movie(t+=0.05){
  clear();
  view (quat = {0.000, 0.000, 0.000, 1.000},
      fov = 30, near = 0.01, far = 1000,
      tx = 0.051, ty = 0.139, tz = -2.607,
      width = 1024, height = 744,bg={0,0,0});
draw_vof (c = "f");

// squares (color = "u.x",spread=5,cbar=true);
squares("u.x", 
                   NULL, // 不指定 z 轴
                   -0.3, 0.3, // 自动计算最小和最大值
                   3, // 扩展因子
                   true, // 线性插值
                   cool_warm, // 颜色映射使用 jet
                   (float[]){0, 0, 1}, // 面的颜色为蓝色
                   (float[]){0, 0, 0}, // 线的颜色为黑色
                   false, // 不使用表达式
                   (coord){0, 0, 1}, // 法向量，适用于 3D
                   0, // 平面的 alpha 值
                   1, // 线宽
                   true, // 启用颜色条
                   15, // 颜色条的大小
                   (float[]){-0.2, -.85}, // 颜色条的位置
                   "u.x", // 标签
                   1, // 标签的缩放因子
                   true, // 水平
                   true, // 添加边框
                   false, // 不添加中间值标签
                   50, // 字体大小
                   "%g", // 标签格式
                   50 // 等级数量
    );
box ();
cells ();
 save ("ux.mp4");
}

event output_vtkfiles(t+=0.05){
char name[80];
    int scaled_value = (int)round(100 * t);
    const char* subdirectory = "vtk_output";  // Specify the subdirectory name

    // Create the subdirectory if it does not exist
    #ifdef _WIN32
    mkdir(subdirectory);  // Windows
    #else
    mkdir(subdirectory, 0777);  // Linux/Unix
    #endif

    // Update the file path to include the subdirectory
    sprintf(name, "%s/out_%03d.vtk", subdirectory, scaled_value);
    FILE *fp = fopen(name, "w");
// scalar * l = list_copy (all);
scalar * l = {f,u.x,u.y,p};
  output_vtk((scalar *) {l}, 1<<LEVEL, (FILE *) fp, true);
}


event snapshot(i += 200)
{
  char name[80];
  sprintf(name, "dump-%03d", (int)i);
  dump(file = name);
}

/**
## End

The wave period is `k_/sqrt(g_*k_)`. We want to run up to 2
(alternatively 4) periods. */

event end(t = 2*k_/sqrt(g_*k_))
{
  fprintf(fout, "i = %d t = %g\n", i, t);
  dump("end");
}


// event end(t=3)
// {
//   fprintf(fout, "i = %d t = %g\n", i, t);
//   dump("end");
// }


/**
## Mesh adaptation

On trees, we adapt the mesh according to the error on volume fraction
and velocity. */

#if TREE
event adapt(i++)
{
  adapt_wavelet({f, u}, (double[]){0.01, uemax, uemax/y_metric, uemax}, LEVEL, 6);
}
#endif

/**
## Running in parallel

This file will work in 2D or 3D, either with parallel multigrid
(without adaptivity), using for example:

~~~bash
qcc -source -D_MPI=1 -grid=multigrid3D wave.c
scp _wave.c occigen.cines.fr:
~~~

and then following a recipe similar to that of the
[atomisation](/src/examples/atomisation.c#on-occigen) example.

To use adaptivity, just do something like:

~~~bash
qcc -source -D_MPI=1 -grid=octree wave.c
scp _wave.c occigen.cines.fr:
~~~
*/
