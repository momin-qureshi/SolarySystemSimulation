/* worker */
#include <cstdio>
#include <cmath>
#include <stddef.h>
#include <mpi.h>

#include "structs.h"


Force getForce(Planet* p, Node* root, double theta) {
   if (root == nullptr) {
      return Force();
   }
   Force force;
   force.x_force = 0;
   force.y_force = 0;

   double x_mid = (root->x1 + root->x2) / 2;
   double y_mid = (root->y1 + root->y2) / 2;

   if (root->p != nullptr && root->p->index != -1) {
      x_mid = root->p->x_pos;
      y_mid = root->p->y_pos;
   }

   double dx = x_mid - p->x_pos;
   double dy = y_mid - p->y_pos;
   double d = sqrt(dx * dx + dy * dy);

   double w = abs(root->x2 - root->x1);

   if (d < RLIMIT)
      d = RLIMIT;

   if (root->p != nullptr && root->p->index != -1 || (w / d) < theta) {
      if (root->p != nullptr && root->p->index == p->index) {
         return force;
      }
      force.x_force = G * p->mass * root->mass * dx / (d * d * d);
      force.y_force = G * p->mass * root->mass * dy / (d * d * d);
      return force;
   }

   for (int i = 0; i < 4; i++) {
      Force f = getForce(p, root->next[i], theta);
      force.x_force += f.x_force;
      force.y_force += f.y_force;
   }
   return force;
}


void updatePlanet(Planet *p, Node *root, double dt, double theta){

    Force force = getForce(p, root, theta);
    
    double ax = force.x_force / p->mass;
    double ay = force.y_force / p->mass;

    p->x_pos = p->x_pos + p->x_vel * dt + 0.5 * ax * (dt * dt);
    p->y_pos = p->y_pos + p->y_vel * dt + 0.5 * ay * (dt * dt);
    p->x_vel = p->x_vel + ax * dt;
    p->y_vel = p->y_vel + ay * dt;

    if (p->x_pos < 0 || p->x_pos > maxX) {
        p->mass = -1;
    }
    if (p->y_pos < 0 || p->y_pos > maxY) {
        p->mass = -1;
    }
}


int main(int argc, char* argv[]) {
   MPI_Comm parent;
   MPI_Init(&argc, &argv);

   int world_rank;
   int world_size;

   MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);


   MPI_Comm_get_parent(&parent);

   unsigned long long steps;
   int n;
   int k;
   double dt;
   double theta;
   MPI_Recv(&steps, 1, MPI_UNSIGNED_LONG_LONG, 0, 0, parent, MPI_STATUS_IGNORE);
   MPI_Recv(&n, 1, MPI_INT, 0, 0, parent, MPI_STATUS_IGNORE);
   MPI_Recv(&k, 1, MPI_INT, 0, 0, parent, MPI_STATUS_IGNORE);
   MPI_Recv(&dt, 1, MPI_DOUBLE, 0, 0, parent, MPI_STATUS_IGNORE);
   MPI_Recv(&theta, 1, MPI_DOUBLE, 0, 0, parent, MPI_STATUS_IGNORE);

   Planet* planets = new Planet[n];
   Planet** p = new Planet*[n];

   int *index = new int[n];
   double *mass = new double[n];
   double *x_pos = new double[n];
   double *y_pos = new double[n];
   double *x_vel = new double[n];
   double *y_vel = new double[n];

   for (int s = 0; s < steps; s++) {
      MPI_Recv(index, n, MPI_INT, 0, 0, parent, MPI_STATUS_IGNORE);
      MPI_Recv(mass, n, MPI_DOUBLE, 0, 0, parent, MPI_STATUS_IGNORE);
      MPI_Recv(x_pos, n, MPI_DOUBLE, 0, 0, parent, MPI_STATUS_IGNORE);
      MPI_Recv(y_pos, n, MPI_DOUBLE, 0, 0, parent, MPI_STATUS_IGNORE);
      MPI_Recv(x_vel, n, MPI_DOUBLE, 0, 0, parent, MPI_STATUS_IGNORE);
      MPI_Recv(y_vel, n, MPI_DOUBLE, 0, 0, parent, MPI_STATUS_IGNORE);

      for (int i = 0; i < n; i++) {
         planets[i].index = index[i];
         planets[i].mass = mass[i];
         planets[i].x_pos = x_pos[i];
         planets[i].y_pos = y_pos[i];
         planets[i].x_vel = x_vel[i];
         planets[i].y_vel = y_vel[i];
      }
      int z = 0;
      for (int i = 0; i < n; i++) {
         if(planets[i].mass != -1){
            p[z++] = &planets[i];
         }
      }
      Node *root = createTree(p, z, 0, maxX, 0, maxY);
      for (int i = 0; i < n; i++) {
         if (planets[i].mass != -1) {
            updatePlanet(&planets[i], root, dt, theta);
         }
      }

      MPI_Send(&planets[k].index, 1, MPI_INT, 0, 0, parent);
      MPI_Send(&planets[k].mass, 1, MPI_DOUBLE, 0, 0, parent);
      MPI_Send(&planets[k].x_pos, 1, MPI_DOUBLE, 0, 0, parent);
      MPI_Send(&planets[k].y_pos, 1, MPI_DOUBLE, 0, 0, parent);
      MPI_Send(&planets[k].x_vel, 1, MPI_DOUBLE, 0, 0, parent);
      MPI_Send(&planets[k].y_vel, 1, MPI_DOUBLE, 0, 0, parent);
   }
   MPI_Comm_disconnect(&parent);
   MPI_Finalize();

   delete [] planets;
   delete [] index;
   delete [] mass;
   delete [] x_pos;
   delete [] y_pos;
   return 0;
}
