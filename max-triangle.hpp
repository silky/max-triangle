#include <gmpxx.h>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>

enum state_flag {status_ok, status_not_convex, status_no_interior, status_runtime_error, status_max_not_at_vertex, status_maxiter_exceeded};

struct state {
    mpq_class nx, ny;
    int ai,bi,ci;
    mpq_class bt, ct;
    state_flag status;
};

void print_status(state_flag status);
bool is_convex(std::vector<mpq_class> polygon);
struct state anchored_triangle(std::vector<mpq_class> polygon, mpq_class nx, mpq_class ny);
struct state maximum_triangle(std::vector<mpq_class> polygon);
struct state naive_maximum_triangle(std::vector<mpq_class> polygon);
