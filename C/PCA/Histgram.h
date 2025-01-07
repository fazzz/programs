

struct histgram {
  double *width;
  int *height;
};

int makeHistgram(double *data, int datanum, double width, double max, double min, struct histgram hist);
