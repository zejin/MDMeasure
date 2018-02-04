#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//
void Double_Center(int n, int p, double *X, double *XX);
double Inner_Prod(int n, double *XX, double *YY);
double Inner_Prod_Perm(int n, int *P, double *XX, double *YY);
void dist_sq(double *X, double *D, int nobs, int ndim, int ncomp, int *ICOMP);
double dCov_term0(double *D, int nobs, int ncomp, int start);
double dCov_term0_perm(double *D, int nobs, int ncomp, int start, int *IPERM);
void next_complete(int *index, int nobs, int ncomp);
void next_incomplete(int *index, int nobs, int ncomp);
double dCov_term1_complete(double *D, int nobs, int ncomp);
double dCov_term1_complete_perm(double *D, int nobs, int ncomp, int *IPERM);
double dCov_term1_complete_simple(double *D, int nobs, int ncomp);
double dCov_term1_complete_simple_perm(double *D, int nobs, int ncomp, int *IPERM);
double dCov_term1_asymmetric(double *D, int nobs, int ncomp, int start);
double dCov_term1_asymmetric_perm(double *D, int nobs, int ncomp, int start, int *IPERM);
double dCov_term1_asymmetric_simple(double *D, int nobs, int ncomp, int start);
double dCov_term1_asymmetric_simple_perm(double *D, int nobs, int ncomp, int start, int *IPERM);
double dCov_term1_symmetric(double *D, int nobs, int ncomp, int start);
double dCov_term1_symmetric_perm(double *D, int nobs, int ncomp, int start, int *IPERM);
double dCov_term1_symmetric_simple(double *D, int nobs, int ncomp, int start);
double dCov_term1_symmetric_simple_perm(double *D, int nobs, int ncomp, int start, int *IPERM);
double dCov_term2_complete(double *D, int nobs, int ncomp);
double dCov_term2_complete_perm(double *D, int nobs, int ncomp, int *IPERM);
double dCov_term2_complete_simple(double *D, int nobs, int ncomp);
double dCov_term2_complete_simple_perm(double *D, int nobs, int ncomp, int *IPERM);
double dCov_term2_asymmetric(double *D, int nobs, int ncomp, int start);
double dCov_term2_asymmetric_perm(double *D, int nobs, int ncomp, int start, int *IPERM);
double dCov_term2_asymmetric_simple(double *D, int nobs, int ncomp, int start);
double dCov_term2_asymmetric_simple_perm(double *D, int nobs, int ncomp, int start, int *IPERM);
double dCov_term2_symmetric(double *D, int nobs, int ncomp, int start);
double dCov_term2_symmetric_perm(double *D, int nobs, int ncomp, int start, int *IPERM);
double dCov_term2_symmetric_simple(double *D, int nobs, int ncomp, int start);
double dCov_term2_symmetric_simple_perm(double *D, int nobs, int ncomp, int start, int *IPERM);
double dCov_all_asymmetric(double *D, int nobs, int ncomp, int start);
double dCov_all_symmetric(double *D, int nobs, int ncomp, int start);
double dCov_all_asymmetric_perm(double *D, int nobs, int ncomp, int start, int *IPERM);
double dCov_all_symmetric_perm(double *D, int nobs, int ncomp, int start, int *IPERM);
void dCov(double *X, double *Y, double *XX, double *YY, double *Q, int *NOBS, int *NDIM);
void dCov_perm(double *XX, double *YY, double *Q, int *NOBS, int *IPERM);
void est_complete(double *X, double *D, double *Q, int *NOBS, int *NDIM, int *NCOMP, int *ICOMP);
void est_complete_perm(double *D, double *Q, int *NOBS, int *NCOMP, int *IPERM);
void est_complete_simple(double *X, double *D, double *Q, int *NOBS, int *NDIM, int *NCOMP, int *ICOMP);
void est_complete_simple_perm(double *D, double *Q, int *NOBS, int *NCOMP, int *IPERM);
void est_asymmetric(double *X, double *D, double *Q, int *NOBS, int *NDIM, int *NCOMP, int *ICOMP);
void est_asymmetric_perm(double *D, double *Q, int *NOBS, int *NCOMP, int *IPERM);
void dCov_asymmetric(double *X, double *D, double *Q, int *NOBS, int *NDIM, int *NCOMP, int *ICOMP);
void dCov_asymmetric_perm(double *D, double *Q, int *NOBS, int *NCOMP, int *IPERM);
void est_asymmetric_simple(double *X, double *D, double *Q, int *NOBS, int *NDIM, int *NCOMP, int *ICOMP);
void est_asymmetric_simple_perm(double *D, double *Q, int *NOBS, int *NCOMP, int *IPERM);
void est_symmetric(double *X, double *D, double *Q, int *NOBS, int *NDIM, int *NCOMP, int *ICOMP);
void est_symmetric_perm(double *D, double *Q, int *NOBS, int *NCOMP, int *IPERM);
void dCov_symmetric(double *X, double *D, double *Q, int *NOBS, int *NDIM, int *NCOMP, int *ICOMP);
void dCov_symmetric_perm(double *D, double *Q, int *NOBS, int *NCOMP, int *IPERM);
void est_symmetric_simple(double *X, double *D, double *Q, int *NOBS, int *NDIM, int *NCOMP, int *ICOMP);
void est_symmetric_simple_perm(double *D, double *Q, int *NOBS, int *NCOMP, int *IPERM);

// X is a pxn matrix, XX is a nxn matrix
void Double_Center(int n, int p, double *X, double *XX) {
  double* row_sum = (double*) calloc(n, sizeof(double));
  double* col_sum = (double*) calloc(n, sizeof(double));
  
  double total_sum = 0.0;
  double elem, part_sum;
  int i, j, k;

  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      if (i != j) {
        part_sum = 0.0;

        //XX[i, j] = |X[i, ] - X[j, ]|
        for (k = 0; k < p; ++k) {
          elem = X[i*p+k] - X[j*p+k];
          part_sum += elem * elem;
        }

        part_sum = sqrt(part_sum);

        XX[i+j*n] = part_sum;
        row_sum[i] += part_sum;
        col_sum[j] += part_sum;
        total_sum += part_sum;
      } else {
        XX[i+j*n] = 0.0;
      }
    }
  }

  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      XX[i+j*n] -= row_sum[i] / n + col_sum[j] / n - total_sum / n / n;
    }
  }

  free(row_sum);
  free(col_sum);
}

// XX is a nxn matrix, YY is a nxn matrix
double Inner_Prod(int n, double *XX, double *YY) {
  double sum = 0.0; 
  int i, j;

  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      // XX[i, j] * YY[i, j]
      sum += XX[i+j*n] * YY[i+j*n];
    }
  }

  return sum / n / n;
}

// XX is a nxn matrix, YY is a nxn matrix
double Inner_Prod_Perm(int n, int *P, double *XX, double *YY) {
  double sum = 0.0; 
  int i, j;

  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      // XX[i, j] * YY[P[i], P[j]]
      sum += XX[i+j*n] * YY[P[i]+P[j]*n];
    }
  }

  return sum / n / n;
}

void dist_sq(double *X, double *D, int nobs, int ndim, int ncomp, int *ICOMP) {
  int i, j, k, l;
  double diff, sumsq;

  for (j = 0; j < nobs; ++j) {
    for (i = 0; i < nobs; ++i) {
      for (k = 0; k < ncomp; ++k) {
        // calculate |X[i, k] - X[j, k]|^2
        sumsq = 0.0;

        if (i != j) {
          for (l = ICOMP[k]; l < ICOMP[k + 1]; ++l) { 
            diff = X[i * ndim + l] - X[j * ndim + l];
            sumsq += diff * diff;
          }
        }

        D[ncomp * (i + j * nobs) + k] = sumsq;
      }
    }
  } 

}

// double dCov_term0_mutual(double *D, int nobs, int ncomp) {
//   int i, j, k;
//   double total = 0.0; 
//   double sumsq;

//   for (j = 0; j < nobs; ++j) {
//     for (i = 0; i < nobs; ++i) {
//       sumsq = 0.0;

//       // calculate |X[i] - X[j]|^2
//       if (i != j) {
//         for (k = 0; k < ncomp; ++k) {
//           sumsq += D[ncomp * (i + j * nobs) + k];
//         }
//       }
        
//       total += sqrt(sumsq);
//     }
//   } 

//   return total / nobs / nobs;
// }

double dCov_term0(double *D, int nobs, int ncomp, int start) {
  int i, j, k;
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < nobs; ++j) {
    for (i = 0; i < nobs; ++i) {
      sumsq = 0.0;

      // calculate |X[i] - X[j]|^2
      if (i != j) {
        for (k = start; k < ncomp; ++k) {
          sumsq += D[ncomp * (i + j * nobs) + k];
        }
      }
        
      total += sqrt(sumsq);
    }
  } 

  return total / nobs / nobs;
}

double dCov_term0_perm(double *D, int nobs, int ncomp, int start, int *IPERM) {
  int i, j, k;
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < nobs; ++j) {
    for (i = 0; i < nobs; ++i) {
      sumsq = 0.0;

      // calculate |X[i] - X[j]|^2
      if (i != j) {
        for (k = start; k < ncomp; ++k) {
          sumsq += D[ncomp * (IPERM[ncomp * i + k] + IPERM[ncomp * j + k] * nobs) + k];
        }
      }
       
      total += sqrt(sumsq);
    }
  } 

  return total / nobs / nobs;
}

// int next_complete(int *index, int nobs, int ncomp) {
//   int i, j;
//   int flag = 0;

//   // find last index whose value is not nobs - 1
//   for (i = ncomp - 1; i >= 0; --i) {
//     if (index[i] != nobs - 1) {
//       index[i] += 1;

//       // update the index after to be 0
//       for (j = i + 1; j < ncomp; ++j) {
//         index[j] = 0;
//       }

//       flag = 1;
//       break;
//     }
//   }

//   return flag;
// }

void next_complete(int *index, int nobs, int ncomp) {
  int i, j;

  // find last index whose value is not nobs - 1
  for (i = ncomp - 1; i >= 0; --i) {
    if (index[i] != nobs - 1) {
      index[i] += 1;

      // update the index after to be 0
      for (j = i + 1; j < ncomp; ++j) {
        index[j] = 0;
      }

      break;
    }
  } 
}

void next_incomplete(int *index, int nobs, int ncomp) {
  int i, temp;

  for (i = 0; i < ncomp; ++i) {
    temp = index[i] + 1;
    index[i] = temp % nobs;
  }
}

double dCov_term1_complete(double *D, int nobs, int ncomp) {
  int i, j, k;
  int niter = pow(nobs, ncomp);
  // int index[ncomp];
  // memset(index, 0, ncomp * sizeof(int));
  int* index = (int*) calloc(ncomp, sizeof(int));
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < niter; ++j) {
    for (i = 0; i < nobs; ++i) {
      sumsq = 0.0;

      // calculate |X[i] - X[?]|^2
      for (k = 0; k < ncomp; ++k) {
        sumsq += D[ncomp * (i + index[k] * nobs) + k];
      }

      total += sqrt(sumsq);
    }

    next_complete(index, nobs, ncomp);
  } 

  free(index);
  return 2 * total / niter / nobs;
}

double dCov_term1_complete_perm(double *D, int nobs, int ncomp, int *IPERM) {
  int i, j, k;
  int niter = pow(nobs, ncomp);
  // int index[ncomp];
  // memset(index, 0, ncomp * sizeof(int));
  int* index = (int*) calloc(ncomp, sizeof(int));
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < niter; ++j) {
    for (i = 0; i < nobs; ++i) {
      sumsq = 0.0;

      // calculate |X[i] - X[?]|^2
      for (k = 0; k < ncomp; ++k) {
        sumsq += D[ncomp * (IPERM[ncomp * i + k] + IPERM[ncomp * index[k] + k] * nobs) + k];
      }

      total += sqrt(sumsq);
    }

    next_complete(index, nobs, ncomp);
  } 

  free(index);
  return 2 * total / niter / nobs;
}

double dCov_term1_complete_simple(double *D, int nobs, int ncomp) {
  int i, j, k, l;
  int* index = (int*) malloc(ncomp * sizeof(int));
  for (l = 0; l < ncomp; ++l) {
    index[l] = l;
  }
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < nobs; ++j) {
    for (i = 0; i < nobs; ++i) {
      sumsq = 0.0;

      // calculate |X[i] - X[?]|^2
      for (k = 0; k < ncomp; ++k) {
        sumsq += D[ncomp * (i + index[k] * nobs) + k];
      }

      total += sqrt(sumsq);
    }

    next_incomplete(index, nobs, ncomp);
  } 

  free(index);
  return 2 * total / nobs / nobs;
}

double dCov_term1_complete_simple_perm(double *D, int nobs, int ncomp, int *IPERM) {
  int i, j, k, l;
  int* index = (int*) malloc(ncomp * sizeof(int));
  for (l = 0; l < ncomp; ++l) {
    index[l] = l;
  }
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < nobs; ++j) {
    for (i = 0; i < nobs; ++i) {
      sumsq = 0.0;

      // calculate |X[i] - X[?]|^2
      for (k = 0; k < ncomp; ++k) {
        sumsq += D[ncomp * (IPERM[ncomp * i + k] + IPERM[ncomp * index[k] + k] * nobs) + k];
      }

      total += sqrt(sumsq);
    }

    next_incomplete(index, nobs, ncomp);
  } 

  free(index);
  return 2 * total / nobs / nobs;
}

double dCov_term1_asymmetric(double *D, int nobs, int ncomp, int start) {
  int i, j, k;
  int niter = pow(nobs, 2);
  int* index = (int*) calloc(2, sizeof(int));
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < niter; ++j) {
    for (i = 0; i < nobs; ++i) {
      sumsq = 0.0;

      // start
      sumsq += D[ncomp * (i + index[0] * nobs) + start];

      // start + 1, ..., end
      for (k = start + 1; k < ncomp; ++k) {
        sumsq += D[ncomp * (i + index[1] * nobs) + k];
      }

      total += sqrt(sumsq);
    }

    next_complete(index, nobs, 2);
  } 

  free(index);
  return 2 * total / niter / nobs;
}

double dCov_term1_asymmetric_perm(double *D, int nobs, int ncomp, int start, int *IPERM) {
  int i, j, k;
  int niter = pow(nobs, 2);
  int* index = (int*) calloc(2, sizeof(int));
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < niter; ++j) {
    for (i = 0; i < nobs; ++i) {
      sumsq = 0.0;

      // start
      sumsq += D[ncomp * (IPERM[ncomp * i + start] + IPERM[ncomp * index[0] + start] * nobs) + start];

      // start + 1, ..., end
      for (k = start + 1; k < ncomp; ++k) {
        sumsq += D[ncomp * (IPERM[ncomp * i + k] + IPERM[ncomp * index[1] + k] * nobs) + k];
      }

      total += sqrt(sumsq);
    }

    next_complete(index, nobs, 2);
  } 

  free(index);
  return 2 * total / niter / nobs;
}

double dCov_term1_asymmetric_simple(double *D, int nobs, int ncomp, int start) {
  int i, j, k, l;
  int* index = (int*) malloc(2 * sizeof(int));
  for (l = 0; l < 2; ++l) {
    index[l] = l;
  }
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < nobs; ++j) {
    for (i = 0; i < nobs; ++i) {
      sumsq = 0.0;

      // start
      sumsq += D[ncomp * (i + index[0] * nobs) + start];

      // start + 1, ..., end
      for (k = start + 1; k < ncomp; ++k) {
        sumsq += D[ncomp * (i + index[1] * nobs) + k];
      }

      total += sqrt(sumsq);
    }

    next_incomplete(index, nobs, 2);
  } 

  free(index);
  return 2 * total / nobs / nobs;
}

double dCov_term1_asymmetric_simple_perm(double *D, int nobs, int ncomp, int start, int *IPERM) {
  int i, j, k, l;
  int* index = (int*) malloc(2 * sizeof(int));
  for (l = 0; l < 2; ++l) {
    index[l] = l;
  }
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < nobs; ++j) {
    for (i = 0; i < nobs; ++i) {
      sumsq = 0.0;

      // start
      sumsq += D[ncomp * (IPERM[ncomp * i + start] + IPERM[ncomp * index[0] + start] * nobs) + start];

      // start + 1, ..., end
      for (k = start + 1; k < ncomp; ++k) {
        sumsq += D[ncomp * (IPERM[ncomp * i + k] + IPERM[ncomp * index[1] + k] * nobs) + k];
      }

      total += sqrt(sumsq);
    }

    next_incomplete(index, nobs, 2);
  } 

  free(index);
  return 2 * total / nobs / nobs;
}

double dCov_term1_symmetric(double *D, int nobs, int ncomp, int start) {
  int i, j, k;
  int niter = pow(nobs, 2);
  int* index = (int*) calloc(2, sizeof(int));
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < niter; ++j) {
    for (i = 0; i < nobs; ++i) {
      sumsq = 0.0;

      // start, ..., end
      for (k = 0; k < ncomp; ++k) {
        if (k == start) {
          sumsq += D[ncomp * (i + index[0] * nobs) + start];
        } else {
          sumsq += D[ncomp * (i + index[1] * nobs) + k];
        }
      }

      total += sqrt(sumsq);
    }

    next_complete(index, nobs, 2);
  } 

  free(index);
  return 2 * total / niter / nobs;
}

double dCov_term1_symmetric_perm(double *D, int nobs, int ncomp, int start, int *IPERM) {
  int i, j, k;
  int niter = pow(nobs, 2);
  int* index = (int*) calloc(2, sizeof(int));
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < niter; ++j) {
    for (i = 0; i < nobs; ++i) {
      sumsq = 0.0;

      // start, ..., end
      for (k = 0; k < ncomp; ++k) {
        if (k == start) {
          sumsq += D[ncomp * (IPERM[ncomp * i + start] + IPERM[ncomp * index[0] + start] * nobs) + start];
        } else {
          sumsq += D[ncomp * (IPERM[ncomp * i + k] + IPERM[ncomp * index[1] + k] * nobs) + k];
        }
      }

      total += sqrt(sumsq);
    }

    next_complete(index, nobs, 2);
  } 

  free(index);
  return 2 * total / niter / nobs;
}

double dCov_term1_symmetric_simple(double *D, int nobs, int ncomp, int start) {
  int i, j, k, l;
  int* index = (int*) malloc(2 * sizeof(int));
  for (l = 0; l < 2; ++l) {
    index[l] = l;
  }
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < nobs; ++j) {
    for (i = 0; i < nobs; ++i) {
      sumsq = 0.0;

      // start, ..., end
      for (k = 0; k < ncomp; ++k) {
        if (k == start) {
          sumsq += D[ncomp * (i + index[0] * nobs) + start];
        } else {
          sumsq += D[ncomp * (i + index[1] * nobs) + k];
        }
      }

      total += sqrt(sumsq);
    }

    next_incomplete(index, nobs, 2);
  } 

  free(index);
  return 2 * total / nobs / nobs;
}

double dCov_term1_symmetric_simple_perm(double *D, int nobs, int ncomp, int start, int *IPERM) {
  int i, j, k, l;
  int* index = (int*) malloc(2 * sizeof(int));
  for (l = 0; l < 2; ++l) {
    index[l] = l;
  }
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < nobs; ++j) {
    for (i = 0; i < nobs; ++i) {
      sumsq = 0.0;

      // start, ..., end
      for (k = 0; k < ncomp; ++k) {
        if (k == start) {
          sumsq += D[ncomp * (IPERM[ncomp * i + start] + IPERM[ncomp * index[0] + start] * nobs) + start];
        } else {
          sumsq += D[ncomp * (IPERM[ncomp * i + k] + IPERM[ncomp * index[1] + k] * nobs) + k];
        }
      }

      total += sqrt(sumsq);
    }

    next_incomplete(index, nobs, 2);
  } 

  free(index);
  return 2 * total / nobs / nobs;
}

double dCov_term2_complete(double *D, int nobs, int ncomp) {
  int i, j, k;
  int niter = pow(nobs, ncomp);
  // int left[ncomp];
  // memset(left, 0, ncomp * sizeof(int));
  // int right[ncomp];
  // memset(right, 0, ncomp * sizeof(int));
  int* left = (int*) calloc(ncomp, sizeof(int));
  int* right = (int*) calloc(ncomp, sizeof(int));
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < niter; ++j) {
    for (i = 0; i < niter; ++i) {
      sumsq = 0.0;

      // calculate |X[i] - X[?]|^2
      for (k = 0; k < ncomp; ++k) {
        sumsq += D[ncomp * (left[k] + right[k] * nobs) + k];
      }

      total += sqrt(sumsq);
      next_complete(left, nobs, ncomp);
    }

    next_complete(right, nobs, ncomp);
    memset(left, 0, ncomp * sizeof(int));
  } 

  free(left);
  free(right);
  return total / niter / niter;
}

double dCov_term2_complete_perm(double *D, int nobs, int ncomp, int *IPERM) {
  int i, j, k;
  int niter = pow(nobs, ncomp);
  // int left[ncomp];
  // memset(left, 0, ncomp * sizeof(int));
  // int right[ncomp];
  // memset(right, 0, ncomp * sizeof(int));
  int* left = (int*) calloc(ncomp, sizeof(int));
  int* right = (int*) calloc(ncomp, sizeof(int));
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < niter; ++j) {
    for (i = 0; i < niter; ++i) {
      sumsq = 0.0;

      // calculate |X[i] - X[?]|^2
      for (k = 0; k < ncomp; ++k) {
        sumsq += D[ncomp * (IPERM[ncomp * left[k] + k] + IPERM[ncomp * right[k] + k] * nobs) + k];
      }

      total += sqrt(sumsq);
      next_complete(left, nobs, ncomp);
    }

    next_complete(right, nobs, ncomp);
    memset(left, 0, ncomp * sizeof(int));
  } 

  free(left);
  free(right);
  return total / niter / niter;
}

double dCov_term2_complete_simple(double *D, int nobs, int ncomp) {
  int i, j, k, l;
  int* left = (int*) malloc(ncomp * sizeof(int));
  for (l = 0; l < ncomp; ++l) {
    left[l] = l;
  }
  int* right = (int*) malloc(ncomp * sizeof(int));
  for (l = 0; l < ncomp; ++l) {
    right[l] = l;
  }
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < nobs; ++j) {
    for (i = 0; i < nobs; ++i) {
      sumsq = 0.0;

      // calculate |X[i] - X[?]|^2
      for (k = 0; k < ncomp; ++k) {
        sumsq += D[ncomp * (left[k] + right[k] * nobs) + k];
      }

      total += sqrt(sumsq);
      next_incomplete(left, nobs, ncomp);
    }

    next_incomplete(right, nobs, ncomp);
    for (l = 0; l < ncomp; ++l) {
      left[l] = l;
    }
  } 

  free(left);
  free(right);
  return total / nobs / nobs;
}

double dCov_term2_complete_simple_perm(double *D, int nobs, int ncomp, int *IPERM) {
  int i, j, k, l;
  int* left = (int*) malloc(ncomp * sizeof(int));
  for (l = 0; l < ncomp; ++l) {
    left[l] = l;
  }
  int* right = (int*) malloc(ncomp * sizeof(int));
  for (l = 0; l < ncomp; ++l) {
    right[l] = l;
  }
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < nobs; ++j) {
    for (i = 0; i < nobs; ++i) {
      sumsq = 0.0;

      // calculate |X[i] - X[?]|^2
      for (k = 0; k < ncomp; ++k) {
        sumsq += D[ncomp * (IPERM[ncomp * left[k] + k] + IPERM[ncomp * right[k] + k] * nobs) + k];
      }

      total += sqrt(sumsq);
      next_incomplete(left, nobs, ncomp);
    }

    next_incomplete(right, nobs, ncomp);
    for (l = 0; l < ncomp; ++l) {
      left[l] = l;
    }
  } 

  free(left);
  free(right);
  return total / nobs / nobs;
}

double dCov_term2_asymmetric(double *D, int nobs, int ncomp, int start) {
  int i, j, k;
  int niter = pow(nobs, 2);
  int* left = (int*) calloc(2, sizeof(int));
  int* right = (int*) calloc(2, sizeof(int));
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < niter; ++j) {
    for (i = 0; i < niter; ++i) {
      sumsq = 0.0;

      // start
      sumsq += D[ncomp * (left[0] + right[0] * nobs) + start];

      // start + 1, ..., end
      for (k = start + 1; k < ncomp; ++k) {
        sumsq += D[ncomp * (left[1] + right[1] * nobs) + k];
      }

      total += sqrt(sumsq);
      next_complete(left, nobs, 2);
    }

    next_complete(right, nobs, 2);
    memset(left, 0, 2 * sizeof(int));
  } 

  free(left);
  free(right);
  return total / niter / niter;
}

double dCov_term2_asymmetric_perm(double *D, int nobs, int ncomp, int start, int *IPERM) {
  int i, j, k;
  int niter = pow(nobs, 2);
  int* left = (int*) calloc(2, sizeof(int));
  int* right = (int*) calloc(2, sizeof(int));
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < niter; ++j) {
    for (i = 0; i < niter; ++i) {
      sumsq = 0.0;

      // start
      sumsq += D[ncomp * (IPERM[ncomp * left[0] + start] + IPERM[ncomp * right[0] + start] * nobs) + start];

      // start + 1, ..., end
      for (k = start + 1; k < ncomp; ++k) {
        sumsq += D[ncomp * (IPERM[ncomp * left[1] + k] + IPERM[ncomp * right[1] + k] * nobs) + k];
      }

      total += sqrt(sumsq);
      next_complete(left, nobs, 2);
    }

    next_complete(right, nobs, 2);
    memset(left, 0, 2 * sizeof(int));
  } 

  free(left);
  free(right);
  return total / niter / niter;
}

double dCov_term2_asymmetric_simple(double *D, int nobs, int ncomp, int start) {
  int i, j, k, l;
  int* left = (int*) malloc(2 * sizeof(int));
  for (l = 0; l < 2; ++l) {
    left[l] = l;
  }
  int* right = (int*) malloc(2 * sizeof(int));
  for (l = 0; l < 2; ++l) {
    right[l] = l;
  }
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < nobs; ++j) {
    for (i = 0; i < nobs; ++i) {
      sumsq = 0.0;

      // start
      sumsq += D[ncomp * (left[0] + right[0] * nobs) + start];

      // start + 1, ..., end
      for (k = start + 1; k < ncomp; ++k) {
        sumsq += D[ncomp * (left[1] + right[1] * nobs) + k];
      }

      total += sqrt(sumsq);
      next_incomplete(left, nobs, 2);
    }

    next_incomplete(right, nobs, 2);
    for (l = 0; l < 2; ++l) {
      left[l] = l;
    }
  } 

  free(left);
  free(right);
  return total / nobs / nobs;
}

double dCov_term2_asymmetric_simple_perm(double *D, int nobs, int ncomp, int start, int *IPERM) {
  int i, j, k, l;
  int* left = (int*) malloc(2 * sizeof(int));
  for (l = 0; l < 2; ++l) {
    left[l] = l;
  }
  int* right = (int*) malloc(2 * sizeof(int));
  for (l = 0; l < 2; ++l) {
    right[l] = l;
  }
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < nobs; ++j) {
    for (i = 0; i < nobs; ++i) {
      sumsq = 0.0;

      // start
      sumsq += D[ncomp * (IPERM[ncomp * left[0] + start] + IPERM[ncomp * right[0] + start] * nobs) + start];

      // start + 1, ..., end
      for (k = start + 1; k < ncomp; ++k) {
        sumsq += D[ncomp * (IPERM[ncomp * left[1] + k] + IPERM[ncomp * right[1] + k] * nobs) + k];
      }

      total += sqrt(sumsq);
      next_incomplete(left, nobs, 2);
    }

    next_incomplete(right, nobs, 2);
    for (l = 0; l < 2; ++l) {
      left[l] = l;
    }
  } 

  free(left);
  free(right);
  return total / nobs / nobs;
}

double dCov_term2_symmetric(double *D, int nobs, int ncomp, int start) {
  int i, j, k;
  int niter = pow(nobs, 2);
  int* left = (int*) calloc(2, sizeof(int));
  int* right = (int*) calloc(2, sizeof(int));
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < niter; ++j) {
    for (i = 0; i < niter; ++i) {
      sumsq = 0.0;

      // start, ..., end
      for (k = 0; k < ncomp; ++k) {
        if (k == start) {
          sumsq += D[ncomp * (left[0] + right[0] * nobs) + start];
        } else {
          sumsq += D[ncomp * (left[1] + right[1] * nobs) + k];
        }
      }

      total += sqrt(sumsq);
      next_complete(left, nobs, 2);
    }

    next_complete(right, nobs, 2);
    memset(left, 0, 2 * sizeof(int));
  } 

  free(left);
  free(right);
  return total / niter / niter;
}

double dCov_term2_symmetric_perm(double *D, int nobs, int ncomp, int start, int *IPERM) {
  int i, j, k;
  int niter = pow(nobs, 2);
  int* left = (int*) calloc(2, sizeof(int));
  int* right = (int*) calloc(2, sizeof(int));
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < niter; ++j) {
    for (i = 0; i < niter; ++i) {
      sumsq = 0.0;

      // start, ..., end
      for (k = 0; k < ncomp; ++k) {
        if (k == start) {
          sumsq += D[ncomp * (IPERM[ncomp * left[0] + start] + IPERM[ncomp * right[0] + start] * nobs) + start];
        } else {
          sumsq += D[ncomp * (IPERM[ncomp * left[1] + k] + IPERM[ncomp * right[1] + k] * nobs) + k];
        }
      }

      total += sqrt(sumsq);
      next_complete(left, nobs, 2);
    }

    next_complete(right, nobs, 2);
    memset(left, 0, 2 * sizeof(int));
  } 

  free(left);
  free(right);
  return total / niter / niter;
}

double dCov_term2_symmetric_simple(double *D, int nobs, int ncomp, int start) {
  int i, j, k, l;
  int* left = (int*) malloc(2 * sizeof(int));
  for (l = 0; l < 2; ++l) {
    left[l] = l;
  }
  int* right = (int*) malloc(2 * sizeof(int));
  for (l = 0; l < 2; ++l) {
    right[l] = l;
  }
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < nobs; ++j) {
    for (i = 0; i < nobs; ++i) {
      sumsq = 0.0;

      // start, ..., end
      for (k = 0; k < ncomp; ++k) {
        if (k == start) {
          sumsq += D[ncomp * (left[0] + right[0] * nobs) + start];
        } else {
          sumsq += D[ncomp * (left[1] + right[1] * nobs) + k];
        }
      }

      total += sqrt(sumsq);
      next_incomplete(left, nobs, 2);
    }

    next_incomplete(right, nobs, 2);
    for (l = 0; l < 2; ++l) {
      left[l] = l;
    }
  } 

  free(left);
  free(right);
  return total / nobs / nobs;
}

double dCov_term2_symmetric_simple_perm(double *D, int nobs, int ncomp, int start, int *IPERM) {
  int i, j, k, l;
  int* left = (int*) malloc(2 * sizeof(int));
  for (l = 0; l < 2; ++l) {
    left[l] = l;
  }
  int* right = (int*) malloc(2 * sizeof(int));
  for (l = 0; l < 2; ++l) {
    right[l] = l;
  }
  double total = 0.0; 
  double sumsq;

  for (j = 0; j < nobs; ++j) {
    for (i = 0; i < nobs; ++i) {
      sumsq = 0.0;

      // start, ..., end
      for (k = 0; k < ncomp; ++k) {
        if (k == start) {
          sumsq += D[ncomp * (IPERM[ncomp * left[0] + start] + IPERM[ncomp * right[0] + start] * nobs) + start];
        } else {
          sumsq += D[ncomp * (IPERM[ncomp * left[1] + k] + IPERM[ncomp * right[1] + k] * nobs) + k];
        }
      }

      total += sqrt(sumsq);
      next_incomplete(left, nobs, 2);
    }

    next_incomplete(right, nobs, 2);
    for (l = 0; l < 2; ++l) {
      left[l] = l;
    }
  } 

  free(left);
  free(right);
  return total / nobs / nobs;
}

double dCov_all_asymmetric(double *D, int nobs, int ncomp, int start) {
  double* X_row_sum = (double*) calloc(nobs, sizeof(double));
  double* X_col_sum = (double*) calloc(nobs, sizeof(double));
  double* Y_row_sum = (double*) calloc(nobs, sizeof(double));
  double* Y_col_sum = (double*) calloc(nobs, sizeof(double));
  double* XX = (double*) calloc(nobs * nobs, sizeof(double));
  double* YY = (double*) calloc(nobs * nobs, sizeof(double));

  double X_total_sum = 0.0;
  double Y_total_sum = 0.0;
  double X_part_sum, Y_part_sum;
  int i, j, k;

  for (j = 0; j < nobs; ++j) {
    for (i = 0; i < nobs; ++i) {
      if (i != j) {
        X_part_sum = D[ncomp * (i + j * nobs) + start];
        Y_part_sum = 0.0;

        //XX[i, j] = |X[i, ] - X[j, ]|
        for (k = start + 1; k < ncomp; ++k) {
          Y_part_sum += D[ncomp * (i + j * nobs) + k];
        }

        X_part_sum = sqrt(X_part_sum);
        Y_part_sum = sqrt(Y_part_sum);

        XX[i+j*nobs] = X_part_sum;
        X_row_sum[i] += X_part_sum;
        X_col_sum[j] += X_part_sum;
        X_total_sum += X_part_sum;

        YY[i+j*nobs] = Y_part_sum;
        Y_row_sum[i] += Y_part_sum;
        Y_col_sum[j] += Y_part_sum;
        Y_total_sum += Y_part_sum;

      } else {
        XX[i+j*nobs] = 0.0;
        YY[i+j*nobs] = 0.0;
      }
    }
  }

  for (j = 0; j < nobs; ++j) {
    for (i = 0; i < nobs; ++i) {
      XX[i+j*nobs] -= X_row_sum[i] / nobs + X_col_sum[j] / nobs - X_total_sum / nobs / nobs;
      YY[i+j*nobs] -= Y_row_sum[i] / nobs + Y_col_sum[j] / nobs - Y_total_sum / nobs / nobs;
    }
  }

  free(X_row_sum);
  free(X_col_sum);
  free(Y_row_sum);
  free(Y_col_sum);

  // XX, YY
  double c0 = Inner_Prod(nobs, XX, YY);

  free(XX);
  free(YY);

  return c0;
}

double dCov_all_symmetric(double *D, int nobs, int ncomp, int start) {
  double* X_row_sum = (double*) calloc(nobs, sizeof(double));
  double* X_col_sum = (double*) calloc(nobs, sizeof(double));
  double* Y_row_sum = (double*) calloc(nobs, sizeof(double));
  double* Y_col_sum = (double*) calloc(nobs, sizeof(double));
  double* XX = (double*) calloc(nobs * nobs, sizeof(double));
  double* YY = (double*) calloc(nobs * nobs, sizeof(double));

  double X_total_sum = 0.0;
  double Y_total_sum = 0.0;
  double X_part_sum, Y_part_sum;
  int i, j, k;

  for (j = 0; j < nobs; ++j) {
    for (i = 0; i < nobs; ++i) {
      if (i != j) {
        X_part_sum = D[ncomp * (i + j * nobs) + start];
        Y_part_sum = 0.0;

        //XX[i, j] = |X[i, ] - X[j, ]|
        for (k = 0; k < ncomp; ++k) {
          if (k != start) {
            Y_part_sum += D[ncomp * (i + j * nobs) + k];
          }
        }

        X_part_sum = sqrt(X_part_sum);
        Y_part_sum = sqrt(Y_part_sum);

        XX[i+j*nobs] = X_part_sum;
        X_row_sum[i] += X_part_sum;
        X_col_sum[j] += X_part_sum;
        X_total_sum += X_part_sum;

        YY[i+j*nobs] = Y_part_sum;
        Y_row_sum[i] += Y_part_sum;
        Y_col_sum[j] += Y_part_sum;
        Y_total_sum += Y_part_sum;

      } else {
        XX[i+j*nobs] = 0.0;
        YY[i+j*nobs] = 0.0;
      }
    }
  }

  for (j = 0; j < nobs; ++j) {
    for (i = 0; i < nobs; ++i) {
      XX[i+j*nobs] -= X_row_sum[i] / nobs + X_col_sum[j] / nobs - X_total_sum / nobs / nobs;
      YY[i+j*nobs] -= Y_row_sum[i] / nobs + Y_col_sum[j] / nobs - Y_total_sum / nobs / nobs;
    }
  }

  free(X_row_sum);
  free(X_col_sum);
  free(Y_row_sum);
  free(Y_col_sum);

  // XX, YY
  double c0 = Inner_Prod(nobs, XX, YY);

  free(XX);
  free(YY);

  return c0;
}

double dCov_all_asymmetric_perm(double *D, int nobs, int ncomp, int start, int *IPERM) {
  double* X_row_sum = (double*) calloc(nobs, sizeof(double));
  double* X_col_sum = (double*) calloc(nobs, sizeof(double));
  double* Y_row_sum = (double*) calloc(nobs, sizeof(double));
  double* Y_col_sum = (double*) calloc(nobs, sizeof(double));
  double* XX = (double*) calloc(nobs * nobs, sizeof(double));
  double* YY = (double*) calloc(nobs * nobs, sizeof(double));

  double X_total_sum = 0.0;
  double Y_total_sum = 0.0;
  double X_part_sum, Y_part_sum;
  int i, j, k;

  for (j = 0; j < nobs; ++j) {
    for (i = 0; i < nobs; ++i) {
      if (i != j) {
        X_part_sum = D[ncomp * (IPERM[ncomp * i + start] + IPERM[ncomp * j + start] * nobs) + start];
        Y_part_sum = 0.0;

        //XX[i, j] = |X[i, ] - X[j, ]|
        for (k = start + 1; k < ncomp; ++k) {
          Y_part_sum += D[ncomp * (IPERM[ncomp * i + k] + IPERM[ncomp * j + k] * nobs) + k];
        }

        X_part_sum = sqrt(X_part_sum);
        Y_part_sum = sqrt(Y_part_sum);

        XX[i+j*nobs] = X_part_sum;
        X_row_sum[i] += X_part_sum;
        X_col_sum[j] += X_part_sum;
        X_total_sum += X_part_sum;

        YY[i+j*nobs] = Y_part_sum;
        Y_row_sum[i] += Y_part_sum;
        Y_col_sum[j] += Y_part_sum;
        Y_total_sum += Y_part_sum;

      } else {
        XX[i+j*nobs] = 0.0;
        YY[i+j*nobs] = 0.0;
      }
    }
  }

  for (j = 0; j < nobs; ++j) {
    for (i = 0; i < nobs; ++i) {
      XX[i+j*nobs] -= X_row_sum[i] / nobs + X_col_sum[j] / nobs - X_total_sum / nobs / nobs;
      YY[i+j*nobs] -= Y_row_sum[i] / nobs + Y_col_sum[j] / nobs - Y_total_sum / nobs / nobs;
    }
  }

  free(X_row_sum);
  free(X_col_sum);
  free(Y_row_sum);
  free(Y_col_sum);

  // XX, YY
  double c0 = Inner_Prod(nobs, XX, YY);

  free(XX);
  free(YY);

  return c0;
}

double dCov_all_symmetric_perm(double *D, int nobs, int ncomp, int start, int *IPERM) {
  double* X_row_sum = (double*) calloc(nobs, sizeof(double));
  double* X_col_sum = (double*) calloc(nobs, sizeof(double));
  double* Y_row_sum = (double*) calloc(nobs, sizeof(double));
  double* Y_col_sum = (double*) calloc(nobs, sizeof(double));
  double* XX = (double*) calloc(nobs * nobs, sizeof(double));
  double* YY = (double*) calloc(nobs * nobs, sizeof(double));

  double X_total_sum = 0.0;
  double Y_total_sum = 0.0;
  double X_part_sum, Y_part_sum;
  int i, j, k;

  for (j = 0; j < nobs; ++j) {
    for (i = 0; i < nobs; ++i) {
      if (i != j) {
        X_part_sum = D[ncomp * (IPERM[ncomp * i + start] + IPERM[ncomp * j + start] * nobs) + start];
        Y_part_sum = 0.0;

        //XX[i, j] = |X[i, ] - X[j, ]|
        for (k = 0; k < ncomp; ++k) {
          if (k != start) {
            Y_part_sum += D[ncomp * (IPERM[ncomp * i + k] + IPERM[ncomp * j + k] * nobs) + k];
          }
        }

        X_part_sum = sqrt(X_part_sum);
        Y_part_sum = sqrt(Y_part_sum);

        XX[i+j*nobs] = X_part_sum;
        X_row_sum[i] += X_part_sum;
        X_col_sum[j] += X_part_sum;
        X_total_sum += X_part_sum;

        YY[i+j*nobs] = Y_part_sum;
        Y_row_sum[i] += Y_part_sum;
        Y_col_sum[j] += Y_part_sum;
        Y_total_sum += Y_part_sum;

      } else {
        XX[i+j*nobs] = 0.0;
        YY[i+j*nobs] = 0.0;
      }
    }
  }

  for (j = 0; j < nobs; ++j) {
    for (i = 0; i < nobs; ++i) {
      XX[i+j*nobs] -= X_row_sum[i] / nobs + X_col_sum[j] / nobs - X_total_sum / nobs / nobs;
      YY[i+j*nobs] -= Y_row_sum[i] / nobs + Y_col_sum[j] / nobs - Y_total_sum / nobs / nobs;
    }
  }

  free(X_row_sum);
  free(X_col_sum);
  free(Y_row_sum);
  free(Y_col_sum);

  // XX, YY
  double c0 = Inner_Prod(nobs, XX, YY);

  free(XX);
  free(YY);

  return c0;
}

void dCov(double *X, double *Y, double *XX, double *YY, double *Q, int *NOBS, int *NDIM) {
  int nobs = NOBS[0];
  int ndim = NDIM[0];

  Double_Center(nobs, ndim, X, XX);
  Double_Center(nobs, ndim, Y, YY);

  Q[0] = Inner_Prod(nobs, XX, YY);
}

// 
void dCov_perm(double *XX, double *YY, double *Q, int *NOBS, int *IPERM) {
  int nobs = NOBS[0];

  Q[0] = Inner_Prod_Perm(nobs, IPERM, XX, YY);
}

void est_complete(double *X, double *D, double *Q, int *NOBS, int *NDIM, int *NCOMP, int *ICOMP) {
  int nobs = NOBS[0];
  int ndim = NDIM[0];
  int ncomp = NCOMP[0];

  // calculate D
  dist_sq(X, D, nobs, ndim, ncomp, ICOMP);

  // calculate Q
  double c0 = dCov_term0(D, nobs, ncomp, 0);
  double c1 = dCov_term1_complete(D, nobs, ncomp);
  double c2 = dCov_term2_complete(D, nobs, ncomp);

  // printf("========== complete mutual independence measure ==========\n");
  // printf("c0: %g\n", c0);  
  // printf("c1: %g\n", c1);
  // printf("c2: %g\n", c2);
 
  Q[0] = c1 - c0 - c2; 
  // printf("dCov: %g\n", Q[0]);
}

void est_complete_perm(double *D, double *Q, int *NOBS, int *NCOMP, int *IPERM) {
  int nobs = NOBS[0];
  int ncomp = NCOMP[0];

  // calculate Q
  double c0 = dCov_term0_perm(D, nobs, ncomp, 0, IPERM);
  double c1 = dCov_term1_complete_perm(D, nobs, ncomp, IPERM);
  double c2 = dCov_term2_complete_perm(D, nobs, ncomp, IPERM);
  // printf("c0: %g\n", c0);
  // printf("c1: %g\n", c1);
  // printf("c2: %g\n", c2);

  Q[0] = c1 - c0 - c2; 
  // printf("perm dCov: %g\n", Q[0]);
}

void est_complete_simple(double *X, double *D, double *Q, int *NOBS, int *NDIM, int *NCOMP, int *ICOMP) {
  int nobs = NOBS[0];
  int ndim = NDIM[0];
  int ncomp = NCOMP[0];

  // calculate D
  dist_sq(X, D, nobs, ndim, ncomp, ICOMP);

  // calculate Q
  double c0 = dCov_term0(D, nobs, ncomp, 0);
  double c1 = dCov_term1_complete_simple(D, nobs, ncomp);
  double c2 = dCov_term2_complete_simple(D, nobs, ncomp);

  // printf("========== incomplete mutual independence measure ==========\n");
  // printf("c0: %g\n", c0);  
  // printf("c1: %g\n", c1);
  // printf("c2: %g\n", c2);
 
  Q[0] = c1 - c0 - c2; 
  // printf("dCov: %g\n", Q[0]);
}

void est_complete_simple_perm(double *D, double *Q, int *NOBS, int *NCOMP, int *IPERM) {
  int nobs = NOBS[0];
  int ncomp = NCOMP[0];

  // calculate Q
  double c0 = dCov_term0_perm(D, nobs, ncomp, 0, IPERM);
  double c1 = dCov_term1_complete_simple_perm(D, nobs, ncomp, IPERM);
  double c2 = dCov_term2_complete_simple_perm(D, nobs, ncomp, IPERM);
  // printf("c0: %g\n", c0);  
  // printf("c1: %g\n", c1);
  // printf("c2: %g\n", c2);
 
  Q[0] = c1 - c0 - c2; 
  // printf("dCov: %g\n", Q[0]);
}

void est_asymmetric(double *X, double *D, double *Q, int *NOBS, int *NDIM, int *NCOMP, int *ICOMP) {
  int nobs = NOBS[0];
  int ndim = NDIM[0];
  int ncomp = NCOMP[0];
  int i;
  double sum = 0.0;
  double temp;

  // calculate D
  dist_sq(X, D, nobs, ndim, ncomp, ICOMP);

  // printf("========== asymmetric mutual independence measure ==========\n");

  // calculate Q
  for (i = 0; i < ncomp - 1; ++i) {
    double c0 = dCov_term0(D, nobs, ncomp, i);
    double c1 = dCov_term1_asymmetric(D, nobs, ncomp, i);
    double c2 = dCov_term2_asymmetric(D, nobs, ncomp, i);

    // printf("step %d\n", i);
    // printf("c0: %g\n", c0);  
    // printf("c1: %g\n", c1);
    // printf("c2: %g\n", c2);
   
    temp = c1 - c0 - c2; 
    // printf("dCov: %g\n", temp);
    sum += temp;
  }

  Q[0] = sum;
  // printf("total dCov: %g\n", sum);
}

void est_asymmetric_perm(double *D, double *Q, int *NOBS, int *NCOMP, int *IPERM) {
  int nobs = NOBS[0];
  int ncomp = NCOMP[0];
  int i;
  double sum = 0.0;
  double temp;

  // calculate Q
  for (i = 0; i < ncomp - 1; ++i) {
    double c0 = dCov_term0_perm(D, nobs, ncomp, i, IPERM);
    double c1 = dCov_term1_asymmetric_perm(D, nobs, ncomp, i, IPERM);
    double c2 = dCov_term2_asymmetric_perm(D, nobs, ncomp, i, IPERM);

    // printf("step %d\n", i);
    // printf("c0: %g\n", c0);  
    // printf("c1: %g\n", c1);
    // printf("c2: %g\n", c2);
   
    temp = c1 - c0 - c2; 
    // printf("dCov: %g\n", temp);
    sum += temp;
  }

  Q[0] = sum;
  // printf("total dCov: %g\n", sum);
}

void dCov_asymmetric(double *X, double *D, double *Q, int *NOBS, int *NDIM, int *NCOMP, int *ICOMP) {
  int nobs = NOBS[0];
  int ndim = NDIM[0];
  int ncomp = NCOMP[0];
  int i;
  double sum = 0.0;
  // double temp;

  // calculate D
  dist_sq(X, D, nobs, ndim, ncomp, ICOMP);

  // printf("========== asymmetric mutual independence measure ==========\n");

  // calculate Q
  for (i = 0; i < ncomp - 1; ++i) {
    double c0 = dCov_all_asymmetric(D, nobs, ncomp, i);
    // double c1 = dCov_term1_asymmetric(D, nobs, ncomp, i);
    // double c2 = dCov_term2_asymmetric(D, nobs, ncomp, i);

    // printf("step %d\n", i);
    // printf("c0: %g\n", c0);  
    // printf("c1: %g\n", c1);
    // printf("c2: %g\n", c2);
   
    // temp = c1 - c0 - c2; 
    // printf("dCov: %g\n", temp);
    sum += c0;
  }

  Q[0] = sum;
  // printf("total dCov: %g\n", sum);
}

void dCov_asymmetric_perm(double *D, double *Q, int *NOBS, int *NCOMP, int *IPERM) {
  int nobs = NOBS[0];
  int ncomp = NCOMP[0];
  int i;
  double sum = 0.0;
  // double temp;

  // calculate Q
  for (i = 0; i < ncomp - 1; ++i) {
    double c0 = dCov_all_asymmetric_perm(D, nobs, ncomp, i, IPERM);
    // double c1 = dCov_term1_asymmetric_perm(D, nobs, ncomp, i, IPERM);
    // double c2 = dCov_term2_asymmetric_perm(D, nobs, ncomp, i, IPERM);

    // printf("step %d\n", i);
    // printf("c0: %g\n", c0);  
    // printf("c1: %g\n", c1);
    // printf("c2: %g\n", c2);
   
    // temp = c1 - c0 - c2; 
    // printf("dCov: %g\n", temp);
    sum += c0;
  }

  Q[0] = sum;
  // printf("total dCov: %g\n", sum);
}

void est_asymmetric_simple(double *X, double *D, double *Q, int *NOBS, int *NDIM, int *NCOMP, int *ICOMP) {
  int nobs = NOBS[0];
  int ndim = NDIM[0];
  int ncomp = NCOMP[0];
  int i;
  double sum = 0.0;
  double temp;

  // calculate D
  dist_sq(X, D, nobs, ndim, ncomp, ICOMP);

  // printf("========== asymmetric mutual independence measure ==========\n");

  // calculate Q
  for (i = 0; i < ncomp - 1; ++i) {
    double c0 = dCov_term0(D, nobs, ncomp, i);
    double c1 = dCov_term1_asymmetric_simple(D, nobs, ncomp, i);
    double c2 = dCov_term2_asymmetric_simple(D, nobs, ncomp, i);

    // printf("step %d\n", i);
    // printf("c0: %g\n", c0);  
    // printf("c1: %g\n", c1);
    // printf("c2: %g\n", c2);
   
    temp = c1 - c0 - c2; 
    // printf("dCov: %g\n", temp);
    sum += temp;
  }

  Q[0] = sum;
  // printf("total dCov: %g\n", sum);
}

void est_asymmetric_simple_perm(double *D, double *Q, int *NOBS, int *NCOMP, int *IPERM) {
  int nobs = NOBS[0];
  int ncomp = NCOMP[0];
  int i;
  double sum = 0.0;
  double temp;

  // calculate Q
  for (i = 0; i < ncomp - 1; ++i) {
    double c0 = dCov_term0_perm(D, nobs, ncomp, i, IPERM);
    double c1 = dCov_term1_asymmetric_simple_perm(D, nobs, ncomp, i, IPERM);
    double c2 = dCov_term2_asymmetric_simple_perm(D, nobs, ncomp, i, IPERM);

    // printf("step %d\n", i);
    // printf("c0: %g\n", c0);  
    // printf("c1: %g\n", c1);
    // printf("c2: %g\n", c2);
   
    temp = c1 - c0 - c2; 
    // printf("dCov: %g\n", temp);
    sum += temp;
  }

  Q[0] = sum;
  // printf("total dCov: %g\n", sum);
}

void est_symmetric(double *X, double *D, double *Q, int *NOBS, int *NDIM, int *NCOMP, int *ICOMP) {
  int nobs = NOBS[0];
  int ndim = NDIM[0];
  int ncomp = NCOMP[0];
  int i;
  double sum = 0.0;
  double temp;

  // calculate D
  dist_sq(X, D, nobs, ndim, ncomp, ICOMP);

  // printf("========== symmetric mutual independence measure ==========\n");

  // calculate Q
  for (i = 0; i < ncomp; ++i) {
    double c0 = dCov_term0(D, nobs, ncomp, 0);
    double c1 = dCov_term1_symmetric(D, nobs, ncomp, i);
    double c2 = dCov_term2_symmetric(D, nobs, ncomp, i);

    // printf("step %d\n", i);
    // printf("c0: %g\n", c0);  
    // printf("c1: %g\n", c1);
    // printf("c2: %g\n", c2);
   
    temp = c1 - c0 - c2; 
    // printf("dCov: %g\n", temp);
    sum += temp;
  }

  Q[0] = sum;
  // printf("total dCov: %g\n", sum);
}

void est_symmetric_perm(double *D, double *Q, int *NOBS, int *NCOMP, int *IPERM) {
  int nobs = NOBS[0];
  int ncomp = NCOMP[0];
  int i;
  double sum = 0.0;
  double temp;

  // calculate Q
  for (i = 0; i < ncomp; ++i) {
    double c0 = dCov_term0_perm(D, nobs, ncomp, 0, IPERM);
    double c1 = dCov_term1_symmetric_perm(D, nobs, ncomp, i, IPERM);
    double c2 = dCov_term2_symmetric_perm(D, nobs, ncomp, i, IPERM);

    // printf("step %d\n", i);
    // printf("c0: %g\n", c0);  
    // printf("c1: %g\n", c1);
    // printf("c2: %g\n", c2);
   
    temp = c1 - c0 - c2; 
    // printf("dCov: %g\n", temp);
    sum += temp;
  }

  Q[0] = sum;
  // printf("total dCov: %g\n", sum);
}

void dCov_symmetric(double *X, double *D, double *Q, int *NOBS, int *NDIM, int *NCOMP, int *ICOMP) {
  int nobs = NOBS[0];
  int ndim = NDIM[0];
  int ncomp = NCOMP[0];
  int i;
  double sum = 0.0;
  // double temp;

  // calculate D
  dist_sq(X, D, nobs, ndim, ncomp, ICOMP);

  // printf("========== asymmetric mutual independence measure ==========\n");

  // calculate Q
  for (i = 0; i < ncomp - 1; ++i) {
    double c0 = dCov_all_symmetric(D, nobs, ncomp, i);
    // double c1 = dCov_term1_asymmetric(D, nobs, ncomp, i);
    // double c2 = dCov_term2_asymmetric(D, nobs, ncomp, i);

    // printf("step %d\n", i);
    // printf("c0: %g\n", c0);  
    // printf("c1: %g\n", c1);
    // printf("c2: %g\n", c2);
   
    // temp = c1 - c0 - c2; 
    // printf("dCov: %g\n", temp);
    sum += c0;
  }

  Q[0] = sum;
  // printf("total dCov: %g\n", sum);
}

void dCov_symmetric_perm(double *D, double *Q, int *NOBS, int *NCOMP, int *IPERM) {
  int nobs = NOBS[0];
  int ncomp = NCOMP[0];
  int i;
  double sum = 0.0;
  // double temp;

  // calculate Q
  for (i = 0; i < ncomp - 1; ++i) {
    double c0 = dCov_all_symmetric_perm(D, nobs, ncomp, i, IPERM);
    // double c1 = dCov_term1_asymmetric_perm(D, nobs, ncomp, i, IPERM);
    // double c2 = dCov_term2_asymmetric_perm(D, nobs, ncomp, i, IPERM);

    // printf("step %d\n", i);
    // printf("c0: %g\n", c0);  
    // printf("c1: %g\n", c1);
    // printf("c2: %g\n", c2);
   
    // temp = c1 - c0 - c2; 
    // printf("dCov: %g\n", temp);
    sum += c0;
  }

  Q[0] = sum;
  // printf("total dCov: %g\n", sum);
}

void est_symmetric_simple(double *X, double *D, double *Q, int *NOBS, int *NDIM, int *NCOMP, int *ICOMP) {
  int nobs = NOBS[0];
  int ndim = NDIM[0];
  int ncomp = NCOMP[0];
  int i;
  double sum = 0.0;
  double temp;

  // calculate D
  dist_sq(X, D, nobs, ndim, ncomp, ICOMP);

  // printf("========== symmetric mutual independence measure ==========\n");

  // calculate Q
  for (i = 0; i < ncomp; ++i) {
    double c0 = dCov_term0(D, nobs, ncomp, 0);
    double c1 = dCov_term1_symmetric_simple(D, nobs, ncomp, i);
    double c2 = dCov_term2_symmetric_simple(D, nobs, ncomp, i);

    // printf("step %d\n", i);
    // printf("c0: %g\n", c0);  
    // printf("c1: %g\n", c1);
    // printf("c2: %g\n", c2);
   
    temp = c1 - c0 - c2; 
    // printf("dCov: %g\n", temp);
    sum += temp;
  }

  Q[0] = sum;
  // printf("total dCov: %g\n", sum);
}

void est_symmetric_simple_perm(double *D, double *Q, int *NOBS, int *NCOMP, int *IPERM) {
  int nobs = NOBS[0];
  int ncomp = NCOMP[0];
  int i;
  double sum = 0.0;
  double temp;

  // calculate Q
  for (i = 0; i < ncomp; ++i) {
    double c0 = dCov_term0_perm(D, nobs, ncomp, 0, IPERM);
    double c1 = dCov_term1_symmetric_simple_perm(D, nobs, ncomp, i, IPERM);
    double c2 = dCov_term2_symmetric_simple_perm(D, nobs, ncomp, i, IPERM);

    // printf("step %d\n", i);
    // printf("c0: %g\n", c0);  
    // printf("c1: %g\n", c1);
    // printf("c2: %g\n", c2);
   
    temp = c1 - c0 - c2; 
    // printf("dCov: %g\n", temp);
    sum += temp;
  }

  Q[0] = sum;
  // printf("total dCov: %g\n", sum);
}