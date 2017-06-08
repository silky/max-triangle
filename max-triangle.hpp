#include <gmpxx.h>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>

enum state_flag {status_ok, status_not_convex, status_no_interior, status_runtime_error, status_max_not_at_vertex, status_maxiter_exceeded};

/*struct state {
    mpq_class nx, ny;
    int ai,bi,ci;
    mpq_class bt, ct;
    state_flag status;
};*/

void print_status(state_flag status);
bool is_convex(std::vector<mpq_class> polygon);
//struct state anchored_triangle(std::vector<mpq_class> polygon, mpq_class nx, mpq_class ny);
void maximum_triangle(std::vector<mpq_class> polygon, unsigned int *ret, state_flag *status);
void brute_force_maximum_triangle(std::vector<mpq_class> polygon, unsigned int *ret, state_flag *status);
