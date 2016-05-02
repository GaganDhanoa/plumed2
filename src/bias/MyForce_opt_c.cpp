/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Bias.h"
#include "ActionRegister.h"
#include <stdio.h>
#include <vector>
#include <math.h>
#include <string>
#include <mpi.h>
#include <time.h>

/* #define PI 3.1415926535897932*/
#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))



typedef unsigned int UINT_32;
using namespace std;

namespace PLMD{
namespace bias{

//+PLUMEDOC BIAS ReplicaForces 
/*
Replicas apply forces on each other to diverge.


To explore phase-space conformations of a protein we simulate it’s crystal structure in a MD simulation and wait
for protein’s trajectory to cross energy barriers. Once this transition state is passed protein now oscillates in
another free energy minima. This transformation changes protein conformation entirely and now the structure
of protein is different from what we started (crystal structure). In this hyper-dimensional phase-space there
are many such transitions possible and it is an area of interest to explore all such possibilities of conformation
structures of proteins.
There is a problem with this simple approach. The protein simulations tend to oscillate around free-energy
minimas. The simulation steps it takes for proteins to cross energy barriers is very high and hence it requires
a lot of computational power and clock time.
We tackle this problem by running multiple(n) simulations of same protein in parallel, and we try to push
these n simulations away from each other so that each one goes in a different direction rather that oscillating
around same position. For each protein simulation we create an external potential based on RMSD distance
between itself and rest of the n-1 simulations. This potential is inversely proportional to RMSD distance and
hence the net force due to the potential is pushing the simulation away from rest n-1 simulations. After 100
simulation steps each protein broadcasts its positional coordinates. Each simulation now calculates its forces
by
1. First optimally aligning itself to other proteins(canceling out translations and rotations)
2. Then by computing RMSD distances and potential
From this potential we calculate net force (by differentiation) and give the numbers back to MD simulator.

Usage of code:
Introduction:
Adds potential according to other (n-1) protein conformations running in parallel. The potential can be given as:
    

    Description of components:
    By default this Action calculates the following quantities.  (Output of program)
        Quantity     Description 
	bias         the instantaneous value of the bias potential 
	force2     the instantaneous value of the squared force due to this bias potential 
	rmsd        the instantaneous sum of all the rmsd due to this bias potential (optional)


Compulsory keywords: (Input of program)
ARG         
the input for this action is the scalar output from one or more other actions. The particular scalars that you will use are referenced using the label of the action. If the label appears on its own then it is assumed that the Action calculates a single scalar value. The value of this scalar is thus used as the input to this new action. If * or *.* appears the scalars calculated by all the proceding actions in the input file are taken. Some actions have multi-component outputs and each component of the output has a specific label. For example a DISTANCE action labelled dist may have three componets x, y and z. To take just the x component you should use dist.x, if you wish to take all three components then use dist.*.More information on the referencing of Actions can be found in the section of the manual on the PLUMED Getting started. Scalar values can also be referenced using POSIX regular expressions as detailed in the section on Regular Expressions. To use this feature you you must compile PLUMED with the appropriate flag. 
    ID            Specifies simulation id
        MAX_STEPS      specifies total number of steps
	FORCE_CONST     specifies Force Constant
	RMSD_THR        specifies Force Constant
	NUM_STRIDE    specifies the number of times it needed to be run
	WRITE_RMSD    Writes to RMSD val

Working Example:
For a two atom system we can apply our bias as:
p1: POSITION ATOM=1
p2: POSITION ATOM=2
myRestraint:    MyForce_opt ARG=p1.x,p1.y,p1.z,p2.x,p2.y,p2.z, ID=0 MAX_STEPS=150000000 
FORCE_CONST=15.0 NUM_STRIDE=1 RMSD_THR=0.0 WRITE_RMSD=0



*/
//+ENDPLUMEDOC

   class MyForce_opt : public Bias{
      UINT_32 id;
      UINT_32 num_sim;
      UINT_32 max_steps;
      UINT_32 num_args;
      double force_const;
      UINT_32 plumed_pos;
      std::vector<double> forces;
      std::vector<double> atoms_pos;
      std::vector<long double> temp_forces;
      UINT_32 stride;
      Value* valueBias;
      Value* valueForce2;
      Value* valuermsdAll;
      Value* valuetimeTaken;
      double rmsdAll;
      double timeTaken;
      double *othersim, *cursim;
      int rank, size;
      std::vector<double> cursimVec, othersimVec;
      FILE *fpRMSD; //*fpForce;
      FILE *fpCPU;
      int write_rmsd;
      long double mean_init[3], mean[3];
      void CenterCoords_INIT();
      void CenterCoords();
      void LoadMatrix(int);
      void LoadMatrix_INIT();
      void Compute_M();
      long double PYTHAG(long double, long double);
      int Compute_SVD();
      void Compute_R();
      void Compute_dM(int, int);
      void Compute_O();
      void Compute_dR();
      void Differenciate_Rotation(int, int);
      void Rotation(int);
      double getForce();
      void getData();
      std::vector<double> weight;
      int num_atoms;
      long double sqnum_atoms;
      std::vector<long double> coords1_xVec, coords1_yVec, coords1_zVec, coords2_xVec, coords2_yVec, coords2_zVec;
      long double *coords1_x, *coords1_y, *coords1_z, *coords2_x, *coords2_y, *coords2_z;
      long double *coords_init[3], *coords[3];
      long double R[3][3], U[3][3], V[3][3], S[3], M[3][3], dR[3][3], dM[3][3], O[3][3];
      std::vector<double> rmsd;
      double rmsdThres;
      public:
      MyForce_opt(const ActionOptions&);
      void calculate();
      static void registerKeywords(Keywords& keys);
   };


//Written by Gagan

   PLUMED_REGISTER_ACTION(MyForce_opt,"MyForce_opt")

   //Registers keyword that are required in plumed.dat file to specify user customisation of force applied
   void MyForce_opt::registerKeywords(Keywords& keys){
      Bias::registerKeywords(keys);
      keys.use("ARG");
      keys.add("compulsory","ID","0","specifies simulation ID");
      keys.add("compulsory","MAX_STEPS","2","specifies total number of steps");
      keys.add("compulsory","FORCE_CONST","1.0","specifies Force Constant");
      keys.add("compulsory","RMSD_THR","0.01","specifies Force Constant");
      keys.add("compulsory","NUM_STRIDE","1","specifies the number of times it needed to be run");
      keys.add("compulsory","WRITE_RMSD","1","Writes to RMSD val");
      componentsAreNotOptional(keys);
	  
	  //keywords for output 
      keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
      keys.addOutputComponent("force2","default","the instantaneous value of the squared force due to this bias potential");
      keys.addOutputComponent("rmsdAll","default","the instantaneous sum of all the rmsd due to this bias potential");
      keys.addOutputComponent("timeTaken","default","the instantaneous sum of all the rmsd due to this bias potential");
   }

   //constructor for class MyForce
   MyForce_opt::MyForce_opt(const ActionOptions&ao):
      PLUMED_BIAS_INIT(ao),
      id(0),
      max_steps(2),
      force_const(1.0),
      stride(1),
      rmsdThres(0.01),
      write_rmsd(1)
   {
   
	//automatic parsing of keywords specified by user in plumed.dat file
      parse("ID",id);
      parse("MAX_STEPS",max_steps);
      parse("FORCE_CONST", force_const);
      parse("NUM_STRIDE", stride);
      parse("RMSD_THR",rmsdThres);
      parse("WRITE_RMSD", write_rmsd);
      checkRead();

      addComponent("bias"); componentIsNotPeriodic("bias");
      addComponent("force2"); componentIsNotPeriodic("force2");
      addComponent("rmsdAll"); componentIsNotPeriodic("rmsdAll");
      addComponent("timeTaken"); componentIsNotPeriodic("timeTaken");
      valueBias=getPntrToComponent("bias");
      valueForce2=getPntrToComponent("force2");
      valuermsdAll=getPntrToComponent("rmsdAll");
      valuetimeTaken=getPntrToComponent("timeTaken");

      num_args = getNumberOfArguments();
      forces.resize(num_args);
      temp_forces.resize(num_args);
      atoms_pos.resize(num_args);
      
      rank = MPI::COMM_WORLD.Get_rank();
      size = MPI::COMM_WORLD.Get_size();
      std::fprintf(stderr, "rank - %d size - %d from id - %d\n", rank, size, id);

      num_sim = size;
      plumed_pos = 0;

      cursimVec.resize(num_args*num_sim);
      othersimVec.resize(num_args*num_sim);
      othersim = othersimVec.data();
      cursim = cursimVec.data();

      char fil_name[40]; 
      if (write_rmsd != 0 && force_const != 0){
         std::sprintf(fil_name,"%d.cpu",id); 
         fpRMSD = std::fopen(fil_name, "wb"); 
      }

      num_atoms = num_args/3;
      sqnum_atoms = sqrt((long double) num_atoms);
	  
	  //initializing vectors required for internal calculations
      coords1_xVec.resize(num_atoms); coords1_yVec.resize(num_atoms); coords1_zVec.resize(num_atoms);
      coords2_xVec.resize(num_atoms); coords2_yVec.resize(num_atoms); coords2_zVec.resize(num_atoms);
      coords1_x = coords1_xVec.data(); coords1_y = coords1_yVec.data(); coords1_z = coords1_zVec.data();
      coords2_x = coords2_xVec.data(); coords2_y = coords2_yVec.data(); coords2_z = coords2_zVec.data();
      coords_init[0] = coords1_x; coords_init[1] = coords1_y; coords_init[2] = coords1_z;
      coords[0] = coords2_x; coords[1] = coords2_y; coords[2] = coords2_z;

      //weight.resize(num_atoms,1.0);
      rmsd.resize(num_sim);
   }

//Written by Dilip


	//initialization of protein coordinate arrays
   inline void MyForce_opt::CenterCoords_INIT() 
   {
      int i;
      long double xsum, ysum, zsum;
      long double *x = coords1_x, *y = coords1_y, *z = coords1_z;

      xsum = ysum = zsum = 0.0;

      for (i = 0; i < num_atoms; ++i) 
      {
         xsum += x[i];
         ysum += y[i];
         zsum += z[i];
      }

      mean_init[0] = xsum / num_atoms;
      mean_init[1] = ysum / num_atoms;
      mean_init[2] = zsum / num_atoms;
      
      for (i = 0; i < num_atoms; ++i)
      {
         x[i] -= mean_init[0];
         y[i] -= mean_init[1];
         z[i] -= mean_init[2];
      }
   }

   inline void MyForce_opt::CenterCoords()
   {
      int i;
      long double xsum, ysum, zsum; 
      long double *x = coords2_x, *y = coords2_y, *z = coords2_z;

      xsum = ysum = zsum = 0.0;

      for (i = 0; i < num_atoms; ++i) 
      {
         xsum += x[i];
         ysum += y[i];
         zsum += z[i];
      }

      mean[0] = xsum / num_atoms;
      mean[1] = ysum / num_atoms;
      mean[2] = zsum / num_atoms;
      
      for (i = 0; i < num_atoms; ++i)
      {
         x[i] -= mean[0];
         y[i] -= mean[1];
         z[i] -= mean[2];
      }
   }
   //creates a matrix of atom positions of protein in coords1 vector
   void MyForce_opt::LoadMatrix_INIT() 
   {
      UINT_32 pos = 0;
      for (int i = 0; i < num_atoms; i++) 
      {
         coords1_x[i] = (long double) atoms_pos[pos++];
         coords1_y[i] = (long double) atoms_pos[pos++];
         coords1_z[i] = (long double) atoms_pos[pos++];
      }
   }
	//loads the atom-matrix previously created in coords2 vector
   void MyForce_opt::LoadMatrix(int sim_id) 
   {
      UINT_32 pos = sim_id*num_args;
      for (int i = 0; i < num_atoms; i++) 
      {
         coords2_x[i] = (long double) othersim[pos++];
         coords2_y[i] = (long double) othersim[pos++];
         coords2_z[i] = (long double) othersim[pos++];
      }
   }
	//computes rotation matrix
   void MyForce_opt::Compute_M()
   {
      int i,j,k;
      for (i = 0; i < 3; i++)
      {
         for (j = 0; j < 3; j++)
         {
            M[i][j] = 0.0;
            for (k = 0; k < num_atoms; k++)
               M[i][j] += coords[i][k]*coords_init[j][k];
            U[i][j] = M[i][j];
         }
      }
   }

   inline long double MyForce_opt::PYTHAG(long double a, long double b)
   {
       long double at = fabs(a), bt = fabs(b), ct, result;

       if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
       else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
       else result = 0.0;
       return(result);
   }


//Code taken from external source (for SVD calculation) 

//computation of SVD required for optimal aligned rotation
   int MyForce_opt::Compute_SVD()
   {
      int m = 3,n = 3;
      int flag, i, its, j, jj, k, l, nm;
      long double c, f, h, s, x, y, z;
      long double anorm = 0.0, g = 0.0, scale = 0.0;
      long double rv1[3];

      /* Householder reduction to bidiagonal form */
      for (i = 0; i < n; i++) 
      {
         /* left-hand reduction */
         l = i + 1;
         rv1[i] = scale * g;
         g = s = scale = 0.0;
         if (i < m) 
         {
            for (k = i; k < m; k++) 
               scale += fabs(U[k][i]);
            if (scale) 
            {
               for (k = i; k < m; k++) 
               {
                  U[k][i] = U[k][i]/scale;
                  s += U[k][i] * U[k][i];
               }
               f = U[i][i];
               g = -SIGN(sqrt(s), f);
               h = f * g - s;
               U[i][i] = f - g;
               if (i != n - 1) 
               {
                  for (j = l; j < n; j++) 
                  {
                     for (s = 0.0, k = i; k < m; k++) 
                        s += U[k][i] * U[k][j];
                     f = s / h;
                     for (k = i; k < m; k++) 
                        U[k][j] += f * U[k][i];
                  }
               }
               for (k = i; k < m; k++) 
                  U[k][i] = U[k][i]*scale;
            }
         }
         S[i] = scale * g;

         /* right-hand reduction */
         g = s = scale = 0.0;
         if (i < m && i != n - 1) 
         {
            for (k = l; k < n; k++) 
               scale += fabs(U[i][k]);
            if (scale) 
            {
               for (k = l; k < n; k++) 
               {
                  U[i][k] = U[i][k]/scale;
                  s += U[i][k] * U[i][k];
               }
               f = U[i][l];
               g = -SIGN(sqrt(s), f);
               h = f * g - s;
               U[i][l] = f - g;
               for (k = l; k < n; k++) 
                  rv1[k] = U[i][k] / h;
               if (i != m - 1) 
               {
                  for (j = l; j < m; j++) 
                  {
                     for (s = 0.0, k = l; k < n; k++) 
                        s += U[j][k] * U[i][k];
                     for (k = l; k < n; k++) 
                        U[j][k] += s * rv1[k];
                  }
               }
               for (k = l; k < n; k++) 
                  U[i][k] = U[i][k]*scale;
            }
         }
         anorm = MAX(anorm, (fabs(S[i]) + fabs(rv1[i])));
      }

      /* accumulate the right-hand transformation */
      for (i = n - 1; i >= 0; i--) 
      {
         if (i < n - 1) 
         {
            if (g) 
            {
               for (j = l; j < n; j++)
                  V[j][i] = (U[i][j] / U[i][l]) / g;
               /* long double division to avoid underflow */
               for (j = l; j < n; j++) 
               {
                  for (s = 0.0, k = l; k < n; k++) 
                     s += U[i][k] * V[k][j];
                  for (k = l; k < n; k++) 
                     V[k][j] += s * V[k][i];
               }
            }
            for (j = l; j < n; j++) 
               V[i][j] = V[j][i] = 0.0;
         }
         V[i][i] = 1.0;
         g = rv1[i];
         l = i;
      }

      /* accumulate the left-hand transformation */
      for (i = n - 1; i >= 0; i--) 
      {
         l = i + 1;
         g = S[i];
         if (i < n - 1) 
            for (j = l; j < n; j++) 
               U[i][j] = 0.0;
         if (g) 
         {
            g = 1.0 / g;
            if (i != n - 1) 
            {
               for (j = l; j < n; j++) 
               {
                  for (s = 0.0, k = l; k < m; k++) 
                     s += U[k][i] * U[k][j];
                  f = (s / U[i][i]) * g;
                  for (k = i; k < m; k++) 
                     U[k][j] += f * U[k][i];
               }
            }
            for (j = i; j < m; j++) 
               U[j][i] = U[j][i]*g;
         }
         else 
         {
            for (j = i; j < m; j++) 
               U[j][i] = 0.0;
         }
         ++U[i][i];
      }

      /* diagonalize the bidiagonal form */
      for (k = n - 1; k >= 0; k--) 
      {                             /* loop over singular values */
         for (its = 0; its < 30; its++) 
         {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--) 
            {                     /* test for splitting */
               nm = l - 1;
               if (fabs(rv1[l]) + anorm == anorm) 
               {
                  flag = 0;
                  break;
               }
               if (fabs(S[nm]) + anorm == anorm) 
                  break;
            }
            if (flag) 
            {
               c = 0.0;
               s = 1.0;
               for (i = l; i <= k; i++) 
               {
                  f = s * rv1[i];
                  if (fabs(f) + anorm != anorm) 
                  {
                     g = S[i];
                     h = PYTHAG(f, g);
                     S[i] = h; 
                     h = 1.0 / h;
                     c = g * h;
                     s = (- f * h);
                     for (j = 0; j < m; j++) 
                     {
                        y = U[j][nm];
                        z = U[j][i];
                        U[j][nm] = y * c + z * s;
                        U[j][i] = z * c - y * s;
                     }
                  }
               }
            }
            z = S[k];
            if (l == k) 
            {                  /* convergence */
               if (z < 0.0) 
               {              /* make singular value nonnegative */
                  S[k] = -z;
                  for (j = 0; j < n; j++) 
                     V[j][k] = (-V[j][k]);
               }
               break;
            }
            if (its >= 30) {
               //free((void*) rv1);
               fprintf(stderr, "No convergence after 30,000! iterations \n");
               return(0);
            }

            /* shift from bottom 2 x 2 minor */
            x = S[l];
            nm = k - 1;
            y = S[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++) 
            {
               i = j + 1;
               g = rv1[i];
               y = S[i];
               h = s * g;
               g = c * g;
               z = PYTHAG(f, h);
               rv1[j] = z;
               c = f / z;
               s = h / z;
               f = x * c + g * s;
               g = g * c - x * s;
               h = y * s;
               y = y * c;
               for (jj = 0; jj < n; jj++) 
               {
                  x = V[jj][j];
                  z = V[jj][i];
                  V[jj][j] = x * c + z * s;
                  V[jj][i] = z * c - x * s;
               }
               z = PYTHAG(f, h);
               S[j] = z;
               if (z) 
               {
                  z = 1.0 / z;
                  c = f * z;
                  s = h * z;
               }
               f = (c * g) + (s * y);
               x = (c * y) - (s * g);
               for (jj = 0; jj < m; jj++) 
               {
                  y = U[jj][j];
                  z = U[jj][i];
                  U[jj][j] = y * c + z * s;
                  U[jj][i] = z * c - y * s;
               }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            S[k] = x;
         }
      }
      return(1);
   } //Code from external source ends here


//computation of final rotation matrix
   inline void MyForce_opt::Compute_R()
   {
      //R = VxUT
      R[0][0] = V[0][0]*U[0][0]+V[0][1]*U[0][1]+V[0][2]*U[0][2];
      R[0][1] = V[0][0]*U[1][0]+V[0][1]*U[1][1]+V[0][2]*U[1][2];
      R[0][2] = V[0][0]*U[2][0]+V[0][1]*U[2][1]+V[0][2]*U[2][2];

      R[1][0] = V[1][0]*U[0][0]+V[1][1]*U[0][1]+V[1][2]*U[0][2];
      R[1][1] = V[1][0]*U[1][0]+V[1][1]*U[1][1]+V[1][2]*U[1][2];
      R[1][2] = V[1][0]*U[2][0]+V[1][1]*U[2][1]+V[1][2]*U[2][2];

      R[2][0] = V[2][0]*U[0][0]+V[2][1]*U[0][1]+V[2][2]*U[0][2];
      R[2][1] = V[2][0]*U[1][0]+V[2][1]*U[1][1]+V[2][2]*U[1][2];
      R[2][2] = V[2][0]*U[2][0]+V[2][1]*U[2][1]+V[2][2]*U[2][2];
   }

   inline void MyForce_opt::Compute_dM(int pos, int dir)
   {
      dM[0][dir] = coords_init[0][pos]; //-mean_init[0];
      dM[1][dir] = coords_init[1][pos]; //-mean_init[1];
      dM[2][dir] = coords_init[2][pos]; //-mean_init[2];
      int dir_1,dir_2;
      dir_1 = (dir+1)%3;
      dir_2 = (dir+2)%3;
      dM[0][dir_1] = dM[0][dir_2] = dM[1][dir_1] = dM[1][dir_2] = dM[2][dir_1] = dM[2][dir_2] = 0.0;
   }

   void MyForce_opt::Compute_O()
   {
      int p,q,m,n;
      for (p = 0; p < 3; p++)
      {
         for (q = 0; q < 3; q++)
         {
            O[p][q] = 0.0;
            for (m = 0; m < 3; m++)
            {
               for (n = 0; n < 3; n++)
               {
                  O[p][q] += (V[m][p]*dM[n][m]*U[n][q]) - (U[m][p]*dM[m][n]*V[n][q]);
               }
            } 
         }
      }
   }
	//computing differential of rotatoin matrix
   inline void MyForce_opt::Compute_dR()
   {
      long double OxUT[3][3];


      OxUT[0][0] = O[0][0]*U[0][0]+O[0][1]*U[0][1]+O[0][2]*U[0][2];
      OxUT[0][1] = O[0][0]*U[1][0]+O[0][1]*U[1][1]+O[0][2]*U[1][2];
      OxUT[0][2] = O[0][0]*U[2][0]+O[0][1]*U[2][1]+O[0][2]*U[2][2];

      OxUT[1][0] = O[1][0]*U[0][0]+O[1][1]*U[0][1]+O[1][2]*U[0][2];
      OxUT[1][1] = O[1][0]*U[1][0]+O[1][1]*U[1][1]+O[1][2]*U[1][2];
      OxUT[1][2] = O[1][0]*U[2][0]+O[1][1]*U[2][1]+O[1][2]*U[2][2];

      OxUT[2][0] = O[2][0]*U[0][0]+O[2][1]*U[0][1]+O[2][2]*U[0][2];
      OxUT[2][1] = O[2][0]*U[1][0]+O[2][1]*U[1][1]+O[2][2]*U[1][2];
      OxUT[2][2] = O[2][0]*U[2][0]+O[2][1]*U[2][1]+O[2][2]*U[2][2];

      dR[0][0] = V[0][0]*OxUT[0][0]+V[0][1]*OxUT[1][0]+V[0][2]*OxUT[2][0];
      dR[0][1] = V[0][0]*OxUT[0][1]+V[0][1]*OxUT[1][1]+V[0][2]*OxUT[2][1];
      dR[0][2] = V[0][0]*OxUT[0][2]+V[0][1]*OxUT[1][2]+V[0][2]*OxUT[2][2];

      dR[1][0] = V[1][0]*OxUT[0][0]+V[1][1]*OxUT[1][0]+V[1][2]*OxUT[2][0];
      dR[1][1] = V[1][0]*OxUT[0][1]+V[1][1]*OxUT[1][1]+V[1][2]*OxUT[2][1];
      dR[1][2] = V[1][0]*OxUT[0][2]+V[1][1]*OxUT[1][2]+V[1][2]*OxUT[2][2];

      dR[2][0] = V[2][0]*OxUT[0][0]+V[2][1]*OxUT[1][0]+V[2][2]*OxUT[2][0];
      dR[2][1] = V[2][0]*OxUT[0][1]+V[2][1]*OxUT[1][1]+V[2][2]*OxUT[2][1];
      dR[2][2] = V[2][0]*OxUT[0][2]+V[2][1]*OxUT[1][2]+V[2][2]*OxUT[2][2];
   }

   //applying optimal rotation to protein
   void MyForce_opt::Rotation(int sim_id)
   {
      LoadMatrix(sim_id);
      CenterCoords();
     
      Compute_M();
   
      Compute_SVD();
   
      Compute_R();

   }

   inline void MyForce_opt::Differenciate_Rotation(int pos, int dir)
   {
      Compute_dM(pos, dir);
      Compute_O();
      Compute_dR();
   }

	//computing final force 
   double MyForce_opt::getForce() 
   {
      double ene = 0;
      rmsdAll = 0.0;
      char rmsdcstr[100]; std::string rmsdstr;

      for (UINT_32 i = 0; i < num_args; i++) forces[i] = 0;

      long double dmu = 1/(1.0*num_atoms);

      for (UINT_32 sim_id = 0; sim_id < num_sim; sim_id++) 
      {
         if (sim_id != rank) 
         {
            Rotation(sim_id);

            long double rms = 0.0;
            long double dpq;
            for (int pos = 0; pos < num_atoms; pos++)
            {
               for (int dir = 0; dir < 3; dir++)
               {
                  int dir_1, dir_2;
                  dir_1 = (dir+1)%3;
                  dir_2 = (dir+2)%3;
                  int cur = pos*3+dir;

                  temp_forces[cur] = 0;
                  
                  // q = dir
                  dpq = coords_init[dir][pos]*R[dir][dir] + coords_init[dir_1][pos]*R[dir_1][dir] + coords_init[dir_2][pos]*R[dir_2][dir] - coords[dir][pos];
                  rms += dpq*dpq;
                  temp_forces[cur] += dpq*R[dir][dir];

                  //q = dir_1
                  dpq = coords_init[dir][pos]*R[dir][dir_1] + coords_init[dir_1][pos]*R[dir_1][dir_1] + coords_init[dir_2][pos]*R[dir_2][dir_1] - coords[dir_1][pos];
                  temp_forces[cur] += dpq*R[dir][dir_1];
                  
                  //q = dir_2
                  dpq = coords_init[dir][pos]*R[dir][dir_2] + coords_init[dir_1][pos]*R[dir_1][dir_2] + coords_init[dir_2][pos]*R[dir_2][dir_2] - coords[dir_2][pos];
                  temp_forces[cur] += dpq*R[dir][dir_2];

                  Differenciate_Rotation(pos,dir);

                  for (int p = 0; p < num_atoms; p++)
                  {
                     //q = dir
                     dpq = coords_init[dir][p]*R[dir][dir] + coords_init[dir_1][p]*R[dir_1][dir] + coords_init[dir_2][p]*R[dir_2][dir] - coords[dir][p];
                     temp_forces[cur] -= dpq*dmu*R[dir][dir];
                     temp_forces[cur] += dpq*(coords_init[dir][p]*dR[dir][dir] + coords_init[dir_1][p]*dR[dir_1][dir] + coords_init[dir_2][p]*dR[dir_2][dir]);

                     //q = dir_1
                     dpq = coords_init[dir][p]*R[dir][dir_1] + coords_init[dir_1][p]*R[dir_1][dir_1] + coords_init[dir_2][p]*R[dir_2][dir_1] - coords[dir_1][p];
                     temp_forces[cur] -= dpq*dmu*R[dir][dir_1];
                     temp_forces[cur] += dpq*(coords_init[dir][p]*dR[dir][dir_1] + coords_init[dir_1][p]*dR[dir_1][dir_1] + coords_init[dir_2][p]*dR[dir_2][dir_1]);

                     //q = dir_2
                     dpq = coords_init[dir][p]*R[dir][dir_2] + coords_init[dir_2][p]*R[dir_1][dir_2] + coords_init[dir_2][p]*R[dir_2][dir_2] - coords[dir_2][p];
                     temp_forces[cur] -= dpq*dmu*R[dir][dir_2];
                     temp_forces[cur] += dpq*(coords_init[dir][p]*dR[dir][dir_2] + coords_init[dir_1][p]*dR[dir_1][dir_2] + coords_init[dir_2][p]*dR[dir_2][dir_2]);
                 }   
               }
            }

            double rmsd_cur = (double) sqrt(rms);
            double rmsd_cub = rmsd_cur * (double) rms;
            
            for (int pos = 0; pos < num_args; pos++)
            {
               forces[pos] += (double)((-1.0*force_const*sqnum_atoms*temp_forces[pos]) / rmsd_cub);
            } 
            rmsd_cur /= sqnum_atoms;
            rmsdAll += rmsd_cur;
            rmsd[sim_id] = rmsd_cur;
	    if (write_rmsd != 0) {std::sprintf(rmsdcstr,"%lf ", rmsd_cur); rmsdstr.append(rmsdcstr);}

            ene += force_const/rmsd_cur;
         }
         else 
         {
            double rmsd_cur = rmsd[sim_id];
	    if (write_rmsd != 0) {std::sprintf(rmsdcstr,"%lf ", rmsd_cur); rmsdstr.append(rmsdcstr);}
         }
      }

      
      if (write_rmsd != 0)
      {
         std::sprintf(rmsdcstr,"%lf\n", rmsdAll); rmsdstr.append(rmsdcstr);
         std::fputs(rmsdstr.c_str(), fpRMSD);
      }
      return ene;
   }

   void MyForce_opt::getData()
   {
         for (UINT_32 i = 0, val = rank*num_args; i < num_args; i++) othersim[val++] = atoms_pos[i];

         for (UINT_32 i = 0; i < num_args; i++)
         {
            atoms_pos[i] = getArgument(i);
         } 
         
         LoadMatrix_INIT();
         CenterCoords_INIT();

         for(UINT_32 i = 0, val = 0; i < num_sim; i++)
         {
            for (UINT_32 j = 0; j < num_args; j++)
            {
               cursim[val+j] = atoms_pos[j];
            }
            val += num_args;
         }
	//broadcasting current simulations coordinates to all other simulations
         MPI::COMM_WORLD.Alltoall(cursim, num_args, MPI::DOUBLE, othersim, num_args, MPI::DOUBLE);
   }

   void MyForce_opt::calculate() 
   {
      double ene;
      double totf2;

      clock_t start, end;
      double cpu_time_used1;
      double cpu_time_used2;
      double cpu_time_used3;
      double cpu_time_used4;

      
      rmsdAll = 0;
      //calculate force according to stride 
      if (plumed_pos % stride == 0 && force_const != 0) 
      {
         start = 1000000*clock();
         getData();
         end = 1000000*clock();
      	 cpu_time_used1 = ((double) (end - start)) / CLOCKS_PER_SEC;
         start = 1000000*clock();
         ene=getForce();
         end = 1000000*clock();
      	 cpu_time_used2 = ((double) (end - start)) / CLOCKS_PER_SEC;
      }
      
      totf2 = 0;
      start = 1000000*clock();
      if (force_const != 0) {
         for(UINT_32 i = 0; i < num_args; ++i)
         {
            const double f_i = forces[i];
            setOutputForce(i,f_i);
            totf2+=f_i*f_i;
         }
      }
      end = 1000000*clock();
      cpu_time_used3 = ((double) (end - start)) / CLOCKS_PER_SEC;


      start = 1000000*clock();
	  
	  //giving final force values to gromacs
      valueBias->set(2.0*ene);
      valueForce2->set(totf2);
      valuermsdAll->set(rmsdAll);
      end = 1000000*clock();
      cpu_time_used4 = ((double) (end - start)) / CLOCKS_PER_SEC;

	

      if (plumed_pos == max_steps) 
      {
         if (write_rmsd != 0 && force_const != 0) std::fclose(fpRMSD);
      }
      plumed_pos++;
   }
}
}
