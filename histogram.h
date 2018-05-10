#include "stdlib.h"
#include "stdio.h"

typedef struct histogram histogram;

histogram *histo_build(size_t nbins, double min, double max);
void histo_update(histogram *h, double value);
void histo_update_multiple(histogram *h, unsigned *value, double *inc);
void histo_write(histogram *h, char *filename, unsigned normalize);
void histo_free(histogram *h);

  
