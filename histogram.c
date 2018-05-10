#include "histogram.h"

struct histogram{
  int    *data;
  size_t nbins;
  double min;
  double max;
  double range;
};


histogram *histo_build(size_t nbins, double min, double max) {
  histogram *h = malloc(sizeof(histogram));
  h->data = calloc(nbins + 1, sizeof(int));
  h->nbins = nbins;
  h->min   = min;
  h->max   = max;
  h->range = max - min;
  return h;
}

void histo_update(histogram *h, double value) {
  if (value <= h->max && value >= h->min) {
    h->data[(int) (((value - h->min) / h->range) * h->nbins)]++;
  }
}

void histo_update_multiple(histogram *h, unsigned *value, double *inc) {
  if (*value <= h->max && *value >= h->min) {
    h->data[(int) (((*value - h->min) / h->range) * h->nbins)] += *inc;
  }
}

void histo_write(histogram *h, char *filename, unsigned normalize) {
  FILE *outfile = fopen(filename, "w");
  double binsize = h->range / h->nbins;
  unsigned int bin;    

  double norm = 0.0;
  if (normalize) {
    for (bin=0; bin<h->nbins+1; bin++) {
      norm += h->data[bin];
    }
    norm *= binsize;
  }
  if (! (norm > 0) )
    norm = 1.0;

  /// Add maximum value to last bin
  h->data[h->nbins-1] += h->data[h->nbins];

  for (bin=0; bin<h->nbins; bin++) {
    fprintf(outfile, "%lf\t%lf\n", (bin + 0.5) * binsize, h->data[bin] / norm);
  }
  fclose(outfile);
}

void histo_free(histogram *h) {
  free(h->data);
  free(h);
}
  
