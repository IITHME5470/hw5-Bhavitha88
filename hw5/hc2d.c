#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void grid(int nx, double xst, double xen, double *x, double *dx)
{
  int i;
  
  *dx = (xen-xst)/(double)(nx-1);

  for(i=0; i<nx; i++)
    x[i] = (double)i * (*dx); // ensure x[0] == 0.0 and x[nx-1] == 1.0
}

void enforce_bcs(int nx, int ny, double *x, double *y, double **T)
{
  int i, j;

  // left and right ends
  for(j=0; j<ny; j++)
  {
    T[0][j] = 0.0;    T[nx-1][j] = 0.0;
  }

  // top and bottom ends
  for(i=0; i<nx; i++)
  {
    T[i][0] = 0.0;    T[i][ny-1] = 0.0;
  }
}

void set_initial_condition(int nx, int ny, double *x, double *y, double **T, double dx, double dy)
{
  int i, j;
  double del=1.0;

  for(i=0; i<nx; i++)
    for(j=0; j<ny; j++)
    {
        T[i][j] = 0.25 * (tanh((x[i]-0.4)/(del*dx)) - tanh((x[i]-0.6)/(del*dx))) 
                       * (tanh((y[j]-0.4)/(del*dy)) - tanh((y[j]-0.6)/(del*dy)));
    //printf("\n%d %lf %lf %lf %lf", i, x[i], tanh((x[i]-0.45)/(del*dx)), tanh((x[i]-0.65)/(del*dx)), T[i]);
    }

  enforce_bcs(nx,ny,x,y,T); //ensure BCs are satisfied at t = 0
}

void timestep_FwdEuler(int nx, int ny, double dt, double dx, double dy, double kdiff, double *x, double *y, double **T, double **rhs)
{

  int i,j;
  double dxsq = dx*dx, dysq = dy*dy;
  // (Forward) Euler scheme
  for(i=1; i<nx-1; i++)
   for(j=1; j<ny-1; j++)
     rhs[i][j] = kdiff*(T[i+1][j]+T[i-1][j]-2.0*T[i][j])/dxsq +
           kdiff*(T[i][j+1]+T[i][j-1]-2.0*T[i][j])/dysq ;  // compute rhs

  for(i=1; i<nx-1; i++)
   for(j=1; j<ny-1; j++)
     T[i][j] = T[i][j] + dt*rhs[i][j];                           // update T^(it+1)[i]

  // set Dirichlet BCs
  enforce_bcs(nx,ny,x,y,T);

}

double get_error_norm_2d(int nx, int ny, double **arr1, double **arr2)
{
  double norm_diff = 0.0, local_diff;
  int i, j;

  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
     local_diff = arr1[i][j] - arr2[i][j];
     norm_diff += local_diff * local_diff;
   }
   norm_diff = sqrt(norm_diff/(double) (nx*ny));
   return norm_diff;
}

//void linsolve_hc2d_gs_adi(int nx, int ny, double rx, double ry, double **rhs, double **T, double **Tnew)
//{
//
//  // write the Alternating Direction Implicit method here
//  // The argument list may need to change to allow for a work array
//  // It will be similar to the Gauss-Seidel function, except 
//  // that the code-snippet labelled `update the solution' will be
//  // replaced by the following logic: In an outer loop, preform y and x sweeps.
//  // In the sweep along y-lines, for each j, solve a tridiagonal system along x
//  // In the sweep along x-lines, for each i, solve a tridiagonal system along y
//}
void linsolve_hc2d_gs_adi(int nx, int ny, double rx, double ry, double **rhs, double **T, double **Tnew) {
  // * X-direction sweep *
  for (int i = 1; i < nx - 1; i++) {
      for (int j = 1; j < ny - 1; j++) {
          Tnew[i][j] = (T[i][j] + rx * (T[i+1][j] + T[i-1][j])) / (1 + 2 * rx);
      }
  }

  // * Y-direction sweep *
  for (int i = 1; i < nx - 1; i++) {
      for (int j = 1; j < ny - 1; j++) {
          T[i][j] = (Tnew[i][j] + ry * (Tnew[i][j+1] + Tnew[i][j-1])) / (1 + 2 * ry);
      }
  }
}
//void linsolve_hc2d_gs_rb(int nx, int ny, double rx, double ry, double **rhs, double **T, double **Tnew)
//{
//
//  // write the red-black Gauss-Seidel method here
//  // It will be very similar to the original Gauss-Seidel, except 
//  // that the code-snippet labelled `update the solution' will be
//  // split into two portions
//}

void linsolve_hc2d_gs(int nx, int ny, double rx, double ry, double **rhs, double **T, double **Tnew)
{
  int i, j, k, max_iter;
  double tol, denom, local_diff, norm_diff;

  max_iter = 1000; tol = 1.0e-6;
  denom = 1.0 + 2.0*rx + 2.0*ry;

  for(k=0; k<max_iter;k++)
  {
    // update the solution
    for(i=1; i<nx-1; i++)
     for(j=1; j<ny-1; j++)
       Tnew[i][j] = (rhs[i][j] + rx*Tnew[i-1][j] + rx*T[i+1][j] + ry*Tnew[i][j-1] + ry*T[i][j+1]) /denom;

    // check for convergence
    norm_diff = get_error_norm_2d(nx, ny, T, Tnew);
    if(norm_diff < tol) break;

    // prepare for next iteration
    for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
       T[i][j] = Tnew[i][j];

  }
  printf("In linsolve_hc2d_gs: %d %e\n", k, norm_diff);
}

void linsolve_hc2d_jacobi(int nx, int ny, double rx, double ry, double **rhs, double **T, double **Tnew)
{
  int i, j, k, max_iter;
  double tol, denom, local_diff, norm_diff;

  max_iter = 1000; tol = 1.0e-6;
  denom = 1.0 + 2.0*rx + 2.0*ry;

  for(k=0; k<max_iter;k++)
  {
    // update the solution
    for(i=1; i<nx-1; i++)
     for(j=1; j<ny-1; j++)
       Tnew[i][j] = (rhs[i][j] + rx*T[i-1][j] + rx*T[i+1][j] + ry*T[i][j-1] + ry*T[i][j+1]) /denom;

    // check for convergence
    norm_diff = get_error_norm_2d(nx, ny, T, Tnew);
    if(norm_diff < tol) break;

    // prepare for next iteration
    for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
       T[i][j] = Tnew[i][j];

  }
  //printf("In linsolve_hc2d_jacobi: %d %e\n", k, norm_diff);

}

void timestep_BwdEuler(int nx, int ny, double dt, double dx, double dy, double kdiff, double *x, double *y, double **T, double **rhs, double **Tnew)
{

  int i,j;
  double rx, ry;

  // Backward (implicit) Euler scheme
  rx = kdiff*dt/(dx*dx);
  ry = kdiff*dt/(dy*dy);

  // initialize rhs to T at current time level
  for(i=1; i<nx-1; i++)
   for(j=1; j<ny-1; j++)
     rhs[i][j] = T[i][j];

  // boundaries: top and bottom
  for(i=0; i<nx; i++)
  {
    rhs[i][0] = 0.0;
    rhs[i][ny-1] = 0.0;
  }

  // boundaries: left and right
  for(j=0; j<ny; j++)
  {
    rhs[0][j] = 0.0;
    rhs[nx-1][j] = 0.0;
  }

  //// -- comment out all except one of the function calls below
  //linsolve_hc2d_jacobi(nx, ny, rx, ry, rhs, T, Tnew);
  // linsolve_hc2d_gs(nx, ny, rx, ry, rhs, T, Tnew);
  linsolve_hc2d_gs_adi(nx, ny, rx, ry, rhs, T, Tnew);
  //linsolve_hc2d_gs_rb(nx, ny, rx, ry, rhs, T, Tnew);

  // set Dirichlet BCs
  enforce_bcs(nx,ny,x,y,T);

}

void output_soln(int nx, int ny, int it, double tcurr, double *x, double *y, double **T)
{
  int i,j;
  FILE* fp;
  char fname[100];

  sprintf(fname, "T_x_y_%06d.dat", it);
  //printf("\n%s\n", fname);

  fp = fopen(fname, "w");
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
      fprintf(fp, "%lf %lf %lf\n", x[i], y[j], T[i][j]);
  fclose(fp);

  printf("Done writing solution for time step = %d, time level = %e\n", it, tcurr);
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main()
{
    int nx, ny;
    double *x, *y, **T, **rhs, tst, ten, xst, xen, yst, yen, dx, dy, dt, tcurr, kdiff;
    double min_dx_dy, **Tnew;
    int i, it, num_time_steps, it_print, j;
    FILE* fp;
    clock_t start_time, end_time;
    double total_time, step_time;

    // Read inputs
    fp = fopen("input2d1.in", "r");
    fscanf(fp, "%d %d\n", &nx, &ny);
    fscanf(fp, "%lf %lf %lf %lf\n", &xst, &xen, &yst, &yen);
    fscanf(fp, "%lf %lf\n", &tst, &ten);
    fscanf(fp, "%lf\n", &kdiff);
    fclose(fp);

    printf("Inputs are: %d %lf %lf %lf %lf %lf\n", nx, xst, xen, tst, ten, kdiff);
    printf("Inputs are: %d %lf %lf\n", ny, yst, yen);

    x = (double *)malloc(nx * sizeof(double));
    y = (double *)malloc(ny * sizeof(double));
    T = (double **)malloc(nx * sizeof(double *));
    for(i = 0; i < nx; i++)
        T[i] = (double *)malloc(ny * sizeof(double));
    rhs = (double **)malloc(nx * sizeof(double *));
    for(i = 0; i < nx; i++)
        rhs[i] = (double *)malloc(ny * sizeof(double));
    Tnew = (double **)malloc(nx * sizeof(double *));
    for(i = 0; i < nx; i++)
        Tnew[i] = (double *)malloc(ny * sizeof(double));

    grid(nx, xst, xen, x, &dx);  // Initialize the grid in x
    grid(ny, yst, yen, y, &dy);  // Initialize the grid in y

    set_initial_condition(nx, ny, x, y, T, dx, dy);  // Initial condition

    // Prepare for time loop
    min_dx_dy = fmin(dx, dy);
    dt = 0.25 / kdiff * (min_dx_dy * min_dx_dy);  // Ensure stability
    num_time_steps = (int)((ten - tst) / dt) + 1;
    it_print = num_time_steps / 5;

    start_time = clock();  // Start total time measurement

    // Start time stepping loop
    for(it = 0; it < num_time_steps; it++)
    {
        tcurr = tst + (double)(it + 1) * dt;

        clock_t step_start = clock();  // Start step time measurement

        // Forward (explicit) Euler
        timestep_FwdEuler(nx, ny, dt, dx, dy, kdiff, x, y, T, rhs);

        clock_t step_end = clock();  // End step time measurement
        step_time = ((double)(step_end - step_start)) / CLOCKS_PER_SEC;

        if(it % it_print == 0)
            printf("Time step %d: Step time = %f seconds\n", it, step_time);

        // Output solution every it_print time steps
        if (it == 9) {
            char filename[100];
            sprintf(filename, "serial_solution_t10.txt");
            FILE *fp = fopen(filename, "w");
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    fprintf(fp, "%0.15lf ", T[i][j]);
                }
                fprintf(fp, "\n");
            }
            fclose(fp);
        }

        if(it % it_print == 0)
            output_soln(nx, ny, it, tcurr, x, y, T);
    }
    
    // Output solution at the last time step
    output_soln(nx, ny, it, tcurr, x, y, T);

    end_time = clock();  // End total time measurement
    total_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

    printf("Total simulation time: %f seconds\n", total_time);
    printf("Average step time: %f seconds\n", total_time / num_time_steps);

    // Free allocated memory
    for(i = 0; i < nx; i++)
        free(T[i]);
    free(T);
    for(i = 0; i < nx; i++)
        free(rhs[i]);
    free(rhs);
    for(i = 0; i < nx; i++)
        free(Tnew[i]);
    free(Tnew);
    free(y);
    free(x);

    return 0;
}

