/*
 * Genetic Algorithm.
 * version: 0.0
 * Author: Marco Arthur P. B. Silva
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>

#define NOT_IMPLEMENTED 0
#define MAX_INTERATION  0xffff
#define MAX_POP         0xffff
#define PC              0.65
#define PM              0.20
#define CRITERIA        0.0001
#define LE              2.71828
#define PI              3.14159
#define ALPHA           7.38904
#define GENES           fit_table[global_fit].gene_num

/* codes, check why on fit_table */
#define GP              0
#define HSK             1
#define EXP             2
#define ACK             3
#define CM              4

/* better random */
#define RAFFLE(p)       drand() < (p) ? 1 : 0
#define RS_SCALE        (1.0 / (1.0 + RAND_MAX))

struct individual_str {
    double *gene;
    double fitness;
    double per;
    int unsigned generation;
    int selected;
};

typedef double (*fitnessFunction)(struct individual_str * i);
typedef int (filterFunction)(const struct individual_str * i);

struct population_str {
    int unsigned size;
    int unsigned seed;
    struct individual_str list[MAX_POP];
    int filter[MAX_POP];
};

/* prototypes (our externables) */
void ga_crossover(struct population_str *, int unsigned);
void ga_mutation(struct population_str *, int unsigned);
double ga_fitness(const struct individual_str *);
double ga_selection(struct population_str *, int unsigned);

/* our hard workers */
struct population_str * ga_init(void);
void do_crossover(struct individual_str *, struct individual_str *);
void do_mutation(struct individual_str *);
double init_gene(const double, const double);
void summary(const struct population_str *);
void summary_ind( const struct individual_str *);
struct individual_str *find_best(struct population_str *);
double drand(void);
int get_domain(const int, double *, double *);

/* our fitness functions */
double gp_fitness(struct individual_str *);
double hsk_fitness(struct individual_str *);
double exp_fitness(struct individual_str *);
double ack_fitness(struct individual_str *);
double cm_fitness(struct individual_str *);

/* our crossover function */
void mean_cross(struct individual_str *, struct individual_str *);
void geom_cross(struct individual_str *, struct individual_str *);

/* table of implemented fitness functions */
static const
struct fit_table_str {
    int alg_code;
    char *description;
    fitnessFunction f;
    int gene_num;
    double best_fit;
    int min;
} fit_table[] = {
    {0, "the GP algorithm", &gp_fitness, 2, -3.0, 1},
    {1, "the HSK algorithm", &hsk_fitness, 2, -2.3458, 1},
    {2, "the EXP algorithm", &exp_fitness, 4, 1.0, 0},
    {3, "the ACK algorithm", &ack_fitness, 10, 0.0, 1},
    {4, "the CM algorithm", &cm_fitness, 4, 0.4, 0},
};

/*
 * global variable that chooses the fitness function,
 * so driving our entire execution
 */
static int global_fit = HSK;
