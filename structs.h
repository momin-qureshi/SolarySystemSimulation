#define G 0.0001
#define RLIMIT 0.03
#define PLANET_SIZE 0.025
#define maxX 4
#define maxY 4


struct Planet {
    int index;
    double x_pos;
    double y_pos;
    double mass;
    double x_vel;
    double y_vel;
    double color[3];
};

struct Node {
    double x1, x2, y1, y2;
    double mass;
    Planet *p;

    Node *next[4];

    Node() {
        mass = 0;
        p = nullptr;
        for (int i = 0; i < 4; i++) {
            next[i] = nullptr;
        }
    }

    ~Node() {
        for (int i = 0; i < 4; i++) {
            if (next[i] != nullptr) {
                delete next[i];
            }
        }
    }

    double calcMass() {
        for (int i = 0; i < 4; i++) {
            if (next[i] != nullptr) {
                mass += next[i]->mass;
            }
        }
        return mass;
    }
};

struct Force {
    double x_force;
    double y_force;

    Force() {
        x_force = 0;
        y_force = 0;
    }
};


Node* createTree(Planet** planets, int n, double x1, double x2, double y1, double y2) {
   if (n == 0) {
      return nullptr;
   }

   Node* node = new Node();
   node->x1 = x1;
   node->x2 = x2;
   node->y1 = y1;
   node->y2 = y2;

   double x_mid = (x1 + x2) / 2;
   double y_mid = (y1 + y2) / 2;

   if (n == 1) {
      node->mass = planets[0]->mass;
      node->p = planets[0];
      return node;
   }

   Planet** quad1 = new Planet * [n];
   Planet** quad2 = new Planet * [n];
   Planet** quad3 = new Planet * [n];
   Planet** quad4 = new Planet * [n];

   int q1 = 0;
   int q2 = 0;
   int q3 = 0;
   int q4 = 0;

   for (int i = 0; i < n; i++) {
      if (planets[i]->mass != 0) {
         if (planets[i]->x_pos > x_mid && planets[i]->y_pos < y_mid) {
            quad1[q1] = planets[i];
            q1++;
         }
         else if (planets[i]->x_pos < x_mid && planets[i]->y_pos < y_mid) {
            quad2[q2] = planets[i];
            q2++;
         }
         else if (planets[i]->x_pos < x_mid && planets[i]->y_pos > y_mid) {
            quad3[q3] = planets[i];
            q3++;
         }
         else {
            quad4[q4] = planets[i];
            q4++;
         }
      }
   }

   node->next[0] = createTree(quad1, q1, x_mid, x2, y1, y_mid);
   node->next[1] = createTree(quad2, q2, x1, x_mid, y1, y_mid);
   node->next[2] = createTree(quad3, q3, x1, x_mid, y_mid, y2);
   node->next[3] = createTree(quad4, q4, x_mid, x2, y_mid, y2);

   delete[] quad1;
   delete[] quad2;
   delete[] quad3;
   delete[] quad4;

   node->mass = node->calcMass();
   return node;
}