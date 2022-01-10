#define GL_SILENCE_DEPRECATION

#include "structs.h"

#include <GLFW/glfw3.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <vector>
#include <mpi.h>


Planet* planetFromFile(char *filename, int& num_planets) {
    std::ifstream fin(filename);
    fin >> num_planets;

    Planet* planets = new Planet[num_planets];
    for (int i = 0; i < num_planets; i++) {
        fin >> planets[i].index >> planets[i].x_pos >> planets[i].y_pos >> planets[i].mass >> planets[i].x_vel >> planets[i].y_vel;
        planets[i].color[0] = (double)rand() / (double)RAND_MAX;
        planets[i].color[1] = (double)rand() / (double)RAND_MAX;
        planets[i].color[2] = (double)rand() / (double)RAND_MAX;
    }
    fin.close();
    return planets;
}


void planetsToFile(char *filename, Planet *planets, int num_planets) {
    std::ofstream fout(filename);
    fout << num_planets << std::endl;

    for (int i = 0; i < num_planets; i++) {
        fout << planets[i].index << "\t" << std::scientific << planets[i].x_pos << "\t" << planets[i].y_pos << "\t" << planets[i].mass << "\t" << planets[i].x_vel << "\t" << planets[i].y_vel << std::endl;
    }
    fout.close();
}


void updateParticles(int planet_count, Planet* planets, Node *tree, MPI_Comm *intercomm) {

    int n = planet_count;
    int *index = new int[n];
    double *mass = new double[n];
    double *x_pos = new double[n];
    double *y_pos = new double[n];
    double *x_vel = new double[n];
    double *y_vel = new double[n];   
    for (int i = 0; i < n; i++) {
        index[i] = planets[i].index;
        mass[i] = planets[i].mass;
        x_pos[i] = planets[i].x_pos;
        y_pos[i] = planets[i].y_pos;
        x_vel[i] = planets[i].x_vel;
        y_vel[i] = planets[i].y_vel;
    }
    for (int i = 0; i < n; i++) {
        MPI_Send(index, n, MPI_INT, 0, 0, intercomm[i]);
        MPI_Send(mass, n, MPI_DOUBLE, 0, 0, intercomm[i]);
        MPI_Send(x_pos, n, MPI_DOUBLE, 0, 0, intercomm[i]);
        MPI_Send(y_pos, n, MPI_DOUBLE, 0, 0, intercomm[i]);
        MPI_Send(x_vel, n, MPI_DOUBLE, 0, 0, intercomm[i]);
        MPI_Send(y_vel, n, MPI_DOUBLE, 0, 0, intercomm[i]);
    }
    for (int i = 0; i < n; i++) {
        MPI_Recv(&index[i], 1, MPI_INT, 0, MPI_ANY_TAG, intercomm[i], MPI_STATUS_IGNORE);
        MPI_Recv(&mass[i], 1, MPI_DOUBLE, 0, MPI_ANY_TAG, intercomm[i], MPI_STATUS_IGNORE);
        MPI_Recv(&x_pos[i], 1, MPI_DOUBLE, 0, MPI_ANY_TAG, intercomm[i], MPI_STATUS_IGNORE);
        MPI_Recv(&y_pos[i], 1, MPI_DOUBLE, 0, MPI_ANY_TAG, intercomm[i], MPI_STATUS_IGNORE);
        MPI_Recv(&x_vel[i], 1, MPI_DOUBLE, 0, MPI_ANY_TAG, intercomm[i], MPI_STATUS_IGNORE);
        MPI_Recv(&y_vel[i], 1, MPI_DOUBLE, 0, MPI_ANY_TAG, intercomm[i], MPI_STATUS_IGNORE);
    }
    for (int i = 0; i < n; i++) {
        planets[i].index = index[i];
        planets[i].mass = mass[i];
        planets[i].x_pos = x_pos[i];
        planets[i].y_pos = y_pos[i];
        planets[i].x_vel = x_vel[i];
        planets[i].y_vel = y_vel[i];
    }
}


double pointToPlane(const double &coord, const double &max) {
    return 2 * coord / max - 1;
}


void drawParticle2D(double x_coord, double y_coord,
    double radius,
    double* colors) {
    int k = 0;
    double angle = 0.0f;
    glBegin(GL_TRIANGLE_FAN);
    glColor3f(colors[0], colors[1], colors[2]);
    glVertex2f(x_coord, y_coord);
    for (k = 0;k < 20;k++) {
        angle = (double)(k) / 19 * 2 * 3.141592;
        glVertex2f(x_coord + radius * cos(angle), y_coord + radius * sin(angle));
    }
    glEnd();
}


void drawOctreeBounds2D(Node *node) {
    int i;
    if(node == nullptr)
        return;
    glBegin(GL_LINES);
    // set the color of lines to be white
    glColor3f(1.0f, 1.0f, 1.0f);
    double x_mid = pointToPlane((node->x1 + node->x2) / 2, maxX);
    double y_mid = pointToPlane((node->y1 + node->y2) / 2, maxY);
    // specify the start point's coordinates
    glVertex2f(pointToPlane(node->x1, maxX), y_mid);
    // specify the end point's coordinates
    glVertex2f(pointToPlane(node->x2, maxX), y_mid);
    // do the same for verticle line
    glVertex2f(x_mid, pointToPlane(node->y1, maxY));
    glVertex2f(x_mid, pointToPlane(node->y2, maxY));
    glEnd();
    for (int i = 0; i < 4; i++)
        drawOctreeBounds2D(node->next[i]);
}


int main(int argc, char* argv[]) {

    // Command line arguments
    /*
        -i inputfilename (char *): input file name
        -o outputfilename (char *): output file name
        -s steps (int): number of iterations
        -t \theta(double): threshold for MAC
        -d dt(double): timestep
        -V: (OPTIONAL, see below) flag to turn on visualization window
    */
    int ifile_index = -1;
    int ofile_index = -1;
    double dt = 0.005;
    double theta = 0;
    unsigned long long steps = 10;
    bool verbose = false;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-i") == 0) {
            ifile_index = i + 1;
            i++;
        }
        else if (strcmp(argv[i], "-o") == 0) {
            ofile_index = i + 1;
            i++;
        }
        else if (strcmp(argv[i], "-s") == 0) {
            steps = atoi(argv[i + 1]);
            i++;
        }
        else if (strcmp(argv[i], "-t") == 0) {
            theta = atof(argv[i + 1]);
            i++;
        }
        else if (strcmp(argv[i], "-d") == 0) {
            dt = atof(argv[i + 1]);
            i++;
        }
        else if (strcmp(argv[i], "-V") == 0) {
            verbose = true;
        }
    }


    // MPI SETUP
    int rank, size;
    int errcode;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &rank);

    srand(time(NULL));
    int n;
    Planet* planets = planetFromFile(argv[ifile_index], n);
    double maxMass = 0;
    for (int i = 0; i < n; i++) {
        if (planets[i].mass > maxMass)
            maxMass = planets[i].mass;
    }

    Planet **p = new Planet*[n];

    int *index = new int[n];
    int *mass = new int[n];
    double *x_pos = new double[n];
    double *y_pos = new double[n];
    double *x_vel = new double[n];
    double *y_vel = new double[n];

    MPI_Comm *intercomm = new MPI_Comm[n];

    for(int i = 0; i < n; i++)
        MPI_Comm_spawn("./worker", MPI_ARGV_NULL, 1, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &intercomm[i], MPI_ERRCODES_IGNORE);
    
    for (int i = 0; i < n; i++) {
        MPI_Send(&steps, 1, MPI_UNSIGNED_LONG_LONG, 0, 0, intercomm[i]);
    }
    for (int i = 0; i < n; i++) {
        MPI_Send(&n, 1, MPI_INT, 0, 0, intercomm[i]);
    }
    for (int i = 0; i < n; i++) {
        int k = i;
        MPI_Send(&k, 1, MPI_INT, 0, 0, intercomm[i]);
    }
    for (int i = 0; i < n; i++) {
        MPI_Send(&dt, 1, MPI_DOUBLE, 0, 0, intercomm[i]);
    }
    for (int i = 0; i < n; i++) {
        MPI_Send(&theta, 1, MPI_DOUBLE, 0, 0, intercomm[i]);
    }
    
    GLFWwindow* window;
    if (verbose) {

        /* OpenGL window dims */
        int width = 600, height = 600;
        
        if (!glfwInit()) {
            fprintf(stderr, "Failed to initialize GLFW\n");
            return -1;
        }
        // Open a window and create its OpenGL context
        window = glfwCreateWindow(width, height, "Simulation", NULL, NULL);
        if (window == NULL) {
            fprintf(stderr, "Failed to open GLFW window.\n");
            glfwTerminate();
            return -1;
        }
        glfwMakeContextCurrent(window); // Initialize GLEW

        // Ensure we can capture the escape key being pressed below
        glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
    }

    while (steps--) {

        int k = 0;
        for (int i = 0; i < n; i++) {
            if(planets[i].mass != -1){
                p[k++] = &planets[i];
            }
        }
        Node *tree = createTree(p, k, 0, maxX, 0, maxY);
        updateParticles(n, planets, tree, intercomm);
        
        if(verbose){
            glClear(GL_COLOR_BUFFER_BIT);
            drawOctreeBounds2D(tree);
            for (int p = 0; p < n; p++) {
                double x_window = 2 * planets[p].x_pos / maxX - 1;
                double y_window = 2 * planets[p].y_pos / maxY - 1;
                double radius = planets[p].mass / maxMass * PLANET_SIZE;
                drawParticle2D(x_window, y_window, radius, planets[p].color);
            }
            // Swap buffers
            glfwSwapBuffers(window);
            glfwPollEvents();
        }
    }

    for(int i = 0; i < n; i++)  
        MPI_Comm_disconnect(&intercomm[i]);
    MPI::Finalize();

    if(verbose){
        glfwDestroyWindow(window);

        glfwTerminate();

        planetsToFile(argv[ofile_index], planets, n);
    }
}
