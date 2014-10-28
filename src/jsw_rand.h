#ifndef JSW_RAND_H
#define JSW_RAND_H

/* Seed the RNG. Must be called first */
void          jsw_seed ( unsigned long s );

/* Return a 32-bit random number */
unsigned long jsw_rand ( void );

/* Seed with current system time */
unsigned      jsw_time_seed();


unsigned long jsw_random_integer(unsigned long a, unsigned long b);
double jsw_random_double(double a, double b);
float jsw_random_float(float a, float b);


#endif
