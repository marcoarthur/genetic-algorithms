/*
 * Genetic Algorithm.
 * version: 0.0
 * Author: Marco Arthur P. B. Silva
 *
 * TODO:
 *      - generalized mutation
 *      - decent UI ( nice command line selection for gp/ehc/ etc... )
 *      - save history of individuals (tag population, tag parent data)
 *      - a nice output with -verbose mode  also
 *      - a generalized match pair formation ( sex with two individuals)
 *      - optmize code (there is huge spaces for optmization)
 *
 */

#include "ag.h"

/* Implementation */
double drand (void) {
    double d;
    do {
        d = (((rand () * RS_SCALE) + rand ()) * RS_SCALE + rand ()) * RS_SCALE;
    } while (d >= 1); /* Round off */
    return d;
}

struct individual_str * find_best(struct population_str *p) {
    double max = p->list[0].fitness;
    struct individual_str ind;
    int i,best = 0, min = fit_table[global_fit].min;

    for(i = 0; i < p->size; i++) {
        ind = p->list[i];

        if ( min ) {
            if (ind.fitness < max) {
                max = ind.fitness;
                best = i;
            }
        }
        else {
            if (ind.fitness > max) {
                max = ind.fitness;
                best = i;
            }
        }
    }

    return &(p->list[best]);
}

double init_gene(const double inf, const double sup) {
    assert(inf < sup);
    double d = (double) rand() / RAND_MAX; // U[0,1]

    return d*(sup - inf) + inf; // U[inf, sup]
}

void do_crossover(struct individual_str *p1, struct individual_str *p2) {
    p1->generation++;
    p2->generation++;

    int cross = RAFFLE(PC);
    if ( ! cross )
        return;

    //mean_cross(p1, p2);
	geom_cross(p1, p2);
}

void mean_cross(struct individual_str *p1, struct individual_str *p2) {
    int i;
    for (i = 0; i < GENES; i++)
        p1->gene[i] = (p1->gene[i] + p2->gene[i])/ 2;
}

void geom_cross(struct individual_str *p1, struct individual_str *p2) {
    int i;
    for (i = 0; i < GENES; i++)
        p1->gene[i] = sqrt(p1->gene[i] * p2->gene[i]);
}

void blx_cross(struct individual_str *p1, struct individual_str *p2) {
    int i;
    double b = 2*drand()*ALPHA - ALPHA + 1;
    double inf, sup, c;

    for (i = 0; i < GENES; i++) {
        c = p1->gene[i] - b*(p2->gene[i] - p1->gene[i]);
        get_domain(i, &inf, &sup);

        /* gene outside range: fails */
        if (c < inf || c > sup)
            continue;

        p1->gene[i] = sqrt(p1->gene[i] * p2->gene[i]);
    }
}

void ga_crossover(struct population_str *p, int unsigned generation) {
    int i;

    if (generation == 0)
        return;

    struct individual_str *f = NULL;
    struct individual_str *m = NULL;

    for (i = 0; i < p->size; i++) {
        if (! p->list[i].selected)
            continue;

        if ( !m ) {
            m = &p->list[i];
            continue;
        }

        if ( !f ) {
            f = &p->list[i];
            do_crossover(f, m);
            f = m = NULL;
        }
    }
}

void do_mutation(struct individual_str *ind) {
    int mut = RAFFLE(PM);
    int i;

    if (!mut)
        return;

    double inf, sup;
    int direction;

    /* extreme mutation */
    for (i = 0; i < fit_table[global_fit].gene_num; i++) {
        get_domain(i, &inf, &sup);

        /* inf and sup has equal chances */
        direction = RAFFLE(0.5);
        ind->gene[i] = direction ? inf : sup;
    }
    /*linear mutation*/
    //TODO: assert(NOT_IMPLEMENTED);
}

void ga_mutation(struct population_str *p, int unsigned generation) {
    if (generation == 0)
       return;

    int i;
    int pm = RAFFLE(PM);

    if ( !pm )
        return;
    
    for (i = 0; i < p->size; i++) {
        if (! p->list[i].selected)
            continue;

        do_mutation(&p->list[i]);
    }
}

struct population_str * ga_init(void) {
    int i,j, n_dim = fit_table[global_fit].gene_num;
    struct population_str *p = malloc(sizeof(struct population_str));
    p->size = MAX_POP;
    p->seed = (time(0));
    srand(p->seed);
    double inf, sup;
    
    for ( i = 0; i < p->size; i++) {
        p->list[i].generation = 0;
        p->list[i].gene = calloc(n_dim, sizeof(double));

        for (j = 0; j < n_dim; j++) {
            get_domain(j, &inf, &sup);
            p->list[i].gene[j] = init_gene(inf, sup);
        }

        p->list[i].selected = 1;
        p->filter[i] = 0;
    }

    return p;
}

double gp_fit(const double x, const double y) {
    /* f(0,1) = - 3 best global*/

    double s1, s2;
    s1 = (double)(1.0 + pow(x + y + 1.0, 2.0) * (19.0 - 14.0*x + 3.0*pow(x, 2.0) - 14.0*y + 6.0 * x * y + 3.0 * pow(y,2.0)));
    s2 = (double)(30.0 + pow(2.0 * x - 3.0 * y,2.0) * (18.0 - 32.0*x + 12.0 * pow(x, 2.0) + 48.0 * y - 36.0 * x * y + 27.0 * pow(y,2.0)));

    return s1 * s2;
}

double hsk_fit(const double x, const double y) {
    return (1 - 8 * x + 7 * pow(x,2) - 7.0/3.0*pow(x,3) + 1.0/4.0*pow(x,4))*pow(y,2)*exp(-1*y);
}

double exp_fitness(struct individual_str *ind) {
    int i;
    double sum;

    for(i = 0; i < GENES ; i++)
        sum += ind->gene[i];

    ind->fitness = exp(-0.5*sum);
    return ind->fitness;
}

double ack_fitness(struct individual_str *ind) {
    int i;
    double sum_sr = 0, sum_cs = 0;

    for(i = 0; i < GENES; i++) {
        sum_sr += ind->gene[i]*ind->gene[i];
        sum_cs += cos(2*PI*ind->gene[i]);
    }

    ind->fitness = -20*exp(-0.02*sqrt(sum_sr/10)) - exp(sum_cs/10) + 20 + LE;
    return ind->fitness;
}

double hsk_fitness(struct individual_str *ind) {
    ind->fitness = hsk_fit(ind->gene[0], ind->gene[1]);
    return ind->fitness;
}

double gp_fitness(struct individual_str *ind) {
    ind->fitness = gp_fit(ind->gene[0], ind->gene[1]);
    return ind->fitness;
}

double cm_fitness(struct individual_str *ind) {
    int i;
    double sum_sr = 0, sum_cs = 0;

    for(i = 0; i < GENES; i++) {
        sum_sr += ind->gene[i]*ind->gene[i];
        sum_cs += cos(5*PI*ind->gene[i]);
    }

    ind->fitness = 0.1*sum_cs - sum_sr;
    return ind->fitness;
}

double ga_selection(struct population_str *p, int unsigned generation) {
    int i, cc = 0;
    struct individual_str *ind;
    double sum = 0, sum_sq, mean;
    fitnessFunction f = fit_table[global_fit].f;
    
    /* compute fitness of generation*/
    for (i = 0; i < p->size; i++) {
        ind = &p->list[i];

        /*a dead body*/
        if ( ! ind->selected )
            continue;

        sum += f( ind );
        cc++;
    }

    /* prepare the new generation */
    struct individual_str *best = find_best(p);

    mean = sum / cc;

    cc = 0;
    for (i = 0; i < p->size; i++) {
        ind = &p->list[i];
        
        ind->per = ind->fitness / sum;

        /* guarantee the ones near the best get in pool*/
        ind->selected = ind->fitness > best->fitness*0.7  ? 1 : 0;
        
        /* give a change again but based on pure luck */
        if ( ind->selected == 0 )
            ind->selected = RAFFLE(0.1);

        /* calculate standard deviation */
        if ( ind->selected == 1 ) {
            sum_sq = (ind->fitness - mean)*(ind->fitness - mean);
            cc++;
        }

    }

    /* guarantee the best */
    best->selected = 1;

    return sqrt(sum_sq/( cc - 1));
}

void summary_ind(const struct individual_str *p) {
    int i;
    printf("fitness:%3.4f\t", p->fitness);

    for (i = 0; i < GENES; i++)
        printf("x[%d]=%3.4f\t", i, p->gene[i]);

    printf("generation:%d\n", p->generation);
}

void summary(const struct population_str *p) {
    int i;

    for(i = 0; i < p->size; i++)
        summary_ind(&p->list[i]);
}

int filter_selection(struct population_str *p, filterFunction f) {
    int i;

    for (i = 0; i < p->size; i++)
        p->filter[i] = f(&p->list[i]);

    return 0;
}


int get_domain(const int gene, double *inf, double *sup) {
    *inf = *sup = 0;

    switch(global_fit) {
        case GP:
            if (gene == 0 || gene == 1) {
                *inf = -2.0;
                *sup = 2.0;
                return 0;
            }
            return -1;
        case HSK:
            if (gene == 0 || gene == 1) {
                *inf = 0;
                *sup = 5.0;
                return 0;
            }
        case EXP:
            assert(NOT_IMPLEMENTED);
        case ACK:
            /* for all genes are the same */
            *inf = -30.0;
            *sup = 30.0;
            return 0;
        case CM:
            /* all genes are the same */
            *inf = -1.0;
            *sup = 1.0;
            return 0;
    };

    return -1;
}

void summary_alg(void) {
    printf("Fitness: %s, known otimal solution: %3.4f\n",
            fit_table[global_fit].description,
            fit_table[global_fit].best_fit);
}

int main(int argc, const char *argv[]) {
    int i = 0;
    double std;
    struct individual_str *best;

    /* population */
    struct population_str *pop = ga_init();

    summary_alg();

    /* It's evolution babe :-) */
    do {
        ga_crossover(pop, i);
        ga_mutation(pop, i);
        std = ga_selection(pop, i++);
        best = find_best(pop);
        summary_ind(best);
    } while (i < MAX_INTERATION && std > CRITERIA);

    free(pop);
    return 0;
}
