/* ----------------------------------------------------------------- */
/*             The Speech Signal Processing Toolkit (SPTK)           */
/*             developed by SPTK Working Group                       */
/*             http://sp-tk.sourceforge.net/                         */
/* ----------------------------------------------------------------- */
/*                                                                   */
/*  Copyright (c) 1984-2007  Tokyo Institute of Technology           */
/*                           Interdisciplinary Graduate School of    */
/*                           Science and Engineering                 */
/*                                                                   */
/*                1996-2016  Nagoya Institute of Technology          */
/*                           Department of Computer Science          */
/*                                                                   */
/* All rights reserved.                                              */
/*                                                                   */
/* Redistribution and use in source and binary forms, with or        */
/* without modification, are permitted provided that the following   */
/* conditions are met:                                               */
/*                                                                   */
/* - Redistributions of source code must retain the above copyright  */
/*   notice, this list of conditions and the following disclaimer.   */
/* - Redistributions in binary form must reproduce the above         */
/*   copyright notice, this list of conditions and the following     */
/*   disclaimer in the documentation and/or other materials provided */
/*   with the distribution.                                          */
/* - Neither the name of the SPTK working group nor the names of its */
/*   contributors may be used to endorse or promote products derived */
/*   from this software without specific prior written permission.   */
/*                                                                   */
/* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND            */
/* CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,       */
/* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF          */
/* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          */
/* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS */
/* BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,          */
/* EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED   */
/* TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,     */
/* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON */
/* ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,   */
/* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY    */
/* OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE           */
/* POSSIBILITY OF SUCH DAMAGE.                                       */
/* ----------------------------------------------------------------- */

/****************************************************************

    $Id: _gmm.c,v 1.27 2016/12/22 10:53:04 fjst15124 Exp $

    GMM output prob calculation functions

*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if defined(WIN32)
#include "SPTK.h"
#else
#include <SPTK.h>
#endif

int choleski(double **cov, double **S, const int L);

double cal_ldet(double **var, const int D)
{
   int i, j, l;
   double ldet = 0.0, **tri;

   tri = (double **) malloc(sizeof(double *) * D);
   for (l = 0; l < D; l++)
      tri[l] = dgetmem(D);

   for (i = 0; i < D; i++)
      for (j = 0; j < D; j++)
         tri[i][j] = 0.0;

   if (choleski(var, tri, D)) {
      for (l = 0; l < D; l++)
         ldet += log(tri[l][l]);

      for (l = 0; l < D; l++) {
         free(tri[l]);
      }
      free(tri);

      return (2.0 * ldet);
   } else {
      for (l = 0; l < D; l++) {
         free(tri[l]);
      }
      free(tri);

      return LZERO;
   }
}

double cal_gconst(double *var, const int D)
{
   int d;
   double gconst;

   gconst = D * log(M_2PI);
   for (d = 0; d < D; d++)
      gconst += log(var[d]);

   return (gconst);
}

double cal_gconstf(double **var, const int D)
{
   double gconst, tmp;

   tmp = cal_ldet(var, D);
   if (tmp == LZERO) {
      fprintf(stderr, "WARNING : det is 0!\n");
      return LZERO;
   }
   gconst = D * log(M_2PI);
   gconst += tmp;

   return (gconst);
}

void cal_tri_inv(double **S, double **S_inv, const int L)
{
   int i, j, k;

   for (i = 0; i < L; i++) {
      S_inv[i][i] = 1.0 / S[i][i];
   }
   for (i = 1; i < L; i++)
      for (j = i - 1; j >= 0; j--)
         for (k = j; k < i; k++) {
            S_inv[i][j] = S_inv[i][j] - S[i][k] * S_inv[k][j] / S[i][i];
         }
}

int choleski(double **cov, double **S, const int L)
{
   int i, j, k;
   double tmp;

   for (i = 0; i < L; i++) {
      for (j = 0; j < i; j++) {
         tmp = cov[i][j];
         for (k = 0; k < j; k++)
            tmp -= S[i][k] * S[j][k];
         S[i][j] = tmp / S[j][j];
      }
      tmp = cov[i][i];
      for (k = 0; k < i; k++)
         tmp -= S[i][k] * S[i][k];
      if (tmp <= 0) {
         return 0;
      }
      S[i][i] = sqrt(tmp);
   }
   return 1;
}

void cal_inv(double **cov, double **inv, const int L)
{
   int i, j, k;
   double **S, **S_inv;

   S = (double **) malloc(sizeof(double *) * L);
   S_inv = (double **) malloc(sizeof(double *) * L);

   for (i = 0; i < L; i++) {
      S[i] = dgetmem(L);
      S_inv[i] = dgetmem(L);
   }

   for (i = 0; i < L; i++) {
      for (j = 0; j < L; j++) {
         S[i][j] = 0.0;
         S_inv[i][j] = 0.0;
         inv[i][j] = 0.0;
      }
   }

   if (choleski(cov, S, L) == 0)
      return;
   cal_tri_inv(S, S_inv, L);

   for (i = 0; i < L; i++)
      for (j = 0; j < L; j++) {
         if (i > j)
            for (k = i; k < L; k++)
               inv[i][j] = inv[i][j] + S_inv[k][i] * S_inv[k][j];
         else
            for (k = j; k < L; k++)
               inv[i][j] = inv[i][j] + S_inv[k][i] * S_inv[k][j];
      }

   for (i = 0; i < L; i++) {
      free(S[i]);
      free(S_inv[i]);
   }
   free(S);
   free(S_inv);
}

void fillz_GMM(GMM * gmm)
{
   int m, l, ll;

   for (m = 0; m < gmm->nmix; m++) {
      gmm->weight[m] = 0.;
      if (gmm->full != TR) {
         for (l = 0; l < gmm->dim; l++) {
            gmm->gauss[m].mean[l] = 0.0;
            gmm->gauss[m].var[l] = 0.0;
         }
      } else {
         for (l = 0; l < gmm->dim; l++) {
            gmm->gauss[m].mean[l] = 0.0;
            for (ll = 0; ll < gmm->dim; ll++) {
               gmm->gauss[m].cov[l][ll] = 0.0;
               gmm->gauss[m].inv[l][ll] = 0.0;
            }
         }
      }
   }
}

void maskCov_GMM(GMM * gmm, const int *dim_list, const int cov_dim,
                 const Boolean block_full, const Boolean block_corr)
{

   int row, col, i, k, l, m, *offset;

   offset = (int *) malloc(sizeof(int) * cov_dim + 1);

   offset[0] = 0;
   for (i = 1; i < cov_dim + 1; i++) {
      offset[i] = offset[i - 1] + dim_list[i - 1];
   }

   for (m = 0; m < gmm->nmix; m++) {
      if (block_full == FA && block_corr == FA) {       /* without -c1 and -c2 */
         for (k = 0; k < gmm->dim; k++) {
            for (l = 0; l < gmm->dim; l++) {
               if (k != l) {
                  gmm->gauss[m].cov[k][l] = 0.0;
               }
            }
         }
      } else if (block_full == FA && block_corr == TR) {        /* with -c1 */
         for (row = 0; row < cov_dim; row++) {
            for (col = 0; col < cov_dim; col++) {
               for (k = offset[row]; k < offset[row] + dim_list[row]; k++) {
                  for (l = offset[col]; l < offset[col] + dim_list[col]; l++) {
                     if (dim_list[row] != dim_list[col]) {
                        gmm->gauss[m].cov[k][l] = 0.0;
                     } else {
                        if (offset[row + 1] - k != offset[col + 1] - l) {
                           gmm->gauss[m].cov[k][l] = 0.0;
                        }
                     }
                  }
               }
            }
         }
      } else if (block_full == TR && block_corr == FA) {        /* with -c2 */
         for (row = 0; row < cov_dim; row++) {
            for (col = 0; col < cov_dim; col++) {
               if (row != col) {
                  for (k = offset[row]; k < offset[row] + dim_list[row]; k++) {
                     for (l = offset[col]; l < offset[col] + dim_list[col]; l++) {
                        gmm->gauss[m].cov[k][l] = 0.0;
                     }
                  }
               }
            }
         }
      } else {                  /* with -c1 and -c2 */
         for (row = 0; row < cov_dim; row++) {
            for (col = 0; col < cov_dim; col++) {
               if (dim_list[row] != dim_list[col]) {
                  for (k = offset[row]; k < offset[row] + dim_list[row]; k++) {
                     for (l = offset[col]; l < offset[col] + dim_list[col]; l++) {
                        gmm->gauss[m].cov[k][l] = 0.0;
                     }
                  }
               }
            }
         }
      }
   }

   free(offset);

}


double log_wgd(const GMM * gmm, const int m, const int l1, const int l2,
               const double *dat)
{
   int l, ll;
   double sum, *diff = NULL, tmp, lwgd;

   sum = gmm->gauss[m].gconst;

   if (gmm->full != TR) {
      for (l = l1; l < l2; l++) {
         tmp = dat[l] - gmm->gauss[m].mean[l];
         sum += (tmp * tmp) / gmm->gauss[m].var[l];
      }
   } else {
      diff = dgetmem(l2);
      for (l = l1; l < l2; l++) {
         diff[l] = dat[l] - gmm->gauss[m].mean[l];
      }
      for (l = l1; l < l2; l++) {
         for (ll = l1, tmp = 0.0; ll < l2; ll++) {
            tmp += diff[ll] * gmm->gauss[m].inv[ll][l];
         }
         sum += tmp * diff[l];
      }
      free(diff);
   }

   lwgd = log(gmm->weight[m]) - 0.5 * sum;

   return (lwgd);
}

static double log_add_gmm(double logx, double logy)
{
   double swap, diff, minLogExp, z;

   if (logx < logy) {
      swap = logx;
      logx = logy;
      logy = swap;
   }

   diff = logy - logx;
   minLogExp = -log(-LZERO);

   if (diff < minLogExp)
      return ((logx < LSMALL) ? LZERO : logx);
   else {
      z = exp(diff);
      return (logx + log(1.0 + z));
   }
}

double log_outp(const GMM * gmm, const int l1, const int l2, const double *dat)
{
   int m;
   double logwgd, logb;

   for (m = 0, logb = LZERO; m < gmm->nmix; m++) {
      logwgd = log_wgd(gmm, m, l1, l2, dat);
      logb = log_add_gmm(logb, logwgd);
   }
   return (logb);
}

int alloc_GMM(GMM * gmm, const int M, const int L, const Boolean full)
{
   int m;
   gmm->nmix = M;
   gmm->dim = L;
   gmm->full = full;
   gmm->weight = dgetmem(M);
   gmm->gauss = (Gauss *) getmem(sizeof(Gauss), M);
   for (m = 0; m < M; m++) {
      gmm->gauss[m].mean = dgetmem(L);

      if (full != TR) {
         gmm->gauss[m].var = dgetmem(L);
      } else {
         gmm->gauss[m].cov = ddgetmem(L, L);
         gmm->gauss[m].inv = ddgetmem(L, L);
      }
   }

   return (0);
}

int load_GMM(GMM * gmm, FILE * fp)
{
   int m, l;

   freadf(gmm->weight, sizeof(*(gmm->weight)), gmm->nmix, fp);
   for (m = 0; m < gmm->nmix; m++) {
      freadf(gmm->gauss[m].mean, sizeof(*(gmm->gauss[m].mean)), gmm->dim, fp);

      if (gmm->full != TR) {
         freadf(gmm->gauss[m].var, sizeof(*(gmm->gauss[m].var)), gmm->dim, fp);
         gmm->gauss[m].gconst = cal_gconst(gmm->gauss[m].var, gmm->dim);
      } else {
         for (l = 0; l < gmm->dim; l++) {
            freadf(gmm->gauss[m].cov[l],
                   sizeof(*(gmm->gauss[m].cov[l])), gmm->dim, fp);
         }
      }
   }

   return (0);
}

int save_GMM(const GMM * gmm, FILE * fp)
{
   int m, i, j;

   fwritef(gmm->weight, sizeof(*(gmm->weight)), gmm->nmix, fp);
   for (m = 0; m < gmm->nmix; m++) {
      if (gmm->full != TR) {
         fwritef(gmm->gauss[m].mean, sizeof(*(gmm->gauss[m].mean)), gmm->dim,
                 fp);
         fwritef(gmm->gauss[m].var, sizeof(*(gmm->gauss[m].var)), gmm->dim, fp);
      } else {
         fwritef(gmm->gauss[m].mean, sizeof(*(gmm->gauss[m].mean)), gmm->dim,
                 fp);
         for (i = 0; i < gmm->dim; i++) {
            for (j = 0; j < i; j++) {
               gmm->gauss[m].cov[j][i] = gmm->gauss[m].cov[i][j];
            }
         }
         for (i = 0; i < gmm->dim; i++) {
            fwritef(gmm->gauss[m].cov[i],
                    sizeof(*(gmm->gauss[m].cov[i])), gmm->dim, fp);
         }
      }
   }

   return (0);
}

int prepareCovInv_GMM(GMM * gmm)
{
   int m;
   for (m = 0; m < gmm->nmix; m++) {
      cal_inv(gmm->gauss[m].cov, gmm->gauss[m].inv, gmm->dim);
   }

   return (0);
}

int prepareGconst_GMM(GMM * gmm)
{
   int m;

   for (m = 0; m < gmm->nmix; m++) {
      if (gmm->full == FA) {
         gmm->gauss[m].gconst = cal_gconst(gmm->gauss[m].var, gmm->dim);
      } else {
         gmm->gauss[m].gconst = cal_gconstf(gmm->gauss[m].cov, gmm->dim);
      }
      if (gmm->gauss[m].gconst == LZERO) {
         return -1;
      }
   }

   return (0);
}

int floorWeight_GMM(GMM * gmm, double floor)
{
   int m;
   double sum_w = 0.0, sum_floor = floor * gmm->nmix;

   for (m = 0; m < gmm->nmix; m++) {
      if (gmm->weight[m] < floor) {
         gmm->weight[m] = floor;
      }
      sum_w += gmm->weight[m];
   }
   if (sum_w != 1.0) {
      for (m = 0; m < gmm->nmix; m++) {
         gmm->weight[m] =
             (1.0 - sum_floor) / (sum_w - sum_floor) * (gmm->weight[m] -
                                                        floor) + floor;
      }
   }

   return (0);
}

int floorVar_GMM(GMM * gmm, double floor)
{
   int m, l;
   if (gmm->full == FA) {
      for (m = 0; m < gmm->nmix; m++) {
         for (l = 0; l < gmm->dim; l++) {
            if (gmm->gauss[m].var[l] < floor) {
               gmm->gauss[m].var[l] = floor;
            }
         }
      }
   } else {
      for (m = 0; m < gmm->nmix; m++) {
         for (l = 0; l < gmm->dim; l++) {
            if (gmm->gauss[m].cov[l][l] < floor) {
               gmm->gauss[m].cov[l][l] = floor;
            }
         }
      }
   }

   return (0);
}

int free_GMM(GMM * gmm)
{
   int m;

   for (m = 0; m < gmm->nmix; m++) {
      free(gmm->gauss[m].mean);

      if (gmm->full != TR) {
         free(gmm->gauss[m].var);
      } else {
         free(gmm->gauss[m].cov[0]);
         free(gmm->gauss[m].inv[0]);
         free(gmm->gauss[m].cov);
         free(gmm->gauss[m].inv);
      }
   }
   free(gmm->gauss);
   free(gmm->weight);
   gmm->nmix = 0;
   gmm->dim = 0;
   gmm->full = FA;
   gmm->weight = NULL;
   gmm->gauss = NULL;

   return (0);
}


int gmm_train(double *dat,int *dim_list,int L,int M,int T,int cov_dim,int S,int Imin,int Imax,double E,double V ,double W,double *out_weight,double *out_mean,double *out_cov)
{
	GMM gmm, tgmm;
	double *pd,*cb,*icb,*logwgd,logb,*sum,*xi=NULL,*eta=NULL,
		diff,ave_logp0=0.0,ave_logp1,change=1.0e10,tmp1,tmp2,mapt,test_sum;
	int l,m,N,t,full=1,i,j,k,*tindex,*cntcb,offset_row=0,offset_col=0,row=0,col=0;
	Boolean block_full=FA,block_corr=TR,multiple_dim=TR,full_cov=FA,init_success=0;
	
	logwgd = dgetmem(M);
	sum = dgetmem(M);
    
	/* Initialization of GMM parameters */
	alloc_GMM(&gmm, M, L, full);
	alloc_GMM(&tgmm, M, L, full);
	/* for VQ */
	N = 1;
	while (N < M)
		N *= 2;
	cb = dgetmem(N * L);
	icb = dgetmem(L);
	tindex = (int *) getmem(T, sizeof(int));
	cntcb = (int *) getmem(M, sizeof(int));


	fprintf(stderr, "T = %d  L = %d  M = %d\n", T, L, M);
	fprintf(stderr, "gmm_train:\n    start LBG\n");
	/* LBG */
	fprintf(stderr, "    start LBG vaverage()\n");
	vaverage(dat, L, T, icb);
	fprintf(stderr, "    start LBG lbg()\n");
	lbg(dat, L, T, icb, 1, cb, N, 1000, 1, S, 1, 0.0001, 0.0001);

	for (t = 0, pd = dat; t < T; t++, pd += L) {
		tindex[t] = vq(pd, cb, L, M);
		cntcb[tindex[t]]++;
	}

	for (m = 0; m < M; m++)
	{
		if (cntcb[m] == 0) {
			fprintf(stderr, "    Error: No data for mixture No.%d\n", m,i);
			free_GMM(&gmm);
			free_GMM(&tgmm);
			free(logwgd);
			free(sum);
			free(cb);
			free(icb);
			free(tindex);
			free(cntcb);
			return 1;
		}
	}



	/* weights */
	for (m = 0; m < M; m++) {
		gmm.weight[m] = (double) cntcb[m] / (double) T;
	}
	floorWeight_GMM(&gmm, W);

	/* mean */
	for (m = 0, pd = cb; m < M; m++, pd += L) {
		movem(pd, gmm.gauss[m].mean, sizeof(double), L);
	}

	/* variance */
	if (full != TR) {
		for (t = 0, pd = dat; t < T; t++, pd += L)
			for (l = 0; l < L; l++) {
				diff = gmm.gauss[tindex[t]].mean[l] - pd[l];
				gmm.gauss[tindex[t]].var[l] += diff * diff;
			}

		for (m = 0; m < M; m++) {
			for (l = 0; l < L; l++) {
				gmm.gauss[m].var[l] /= (double) cntcb[m];
			}
		}
	}
	/* full covariance */
	else {
		for (t = 0, pd = dat; t < T; t++, pd += L) {
			for (l = 0; l < L; l++) {
				for (i = 0; i <= l; i++) {
					diff =
						(gmm.gauss[tindex[t]].mean[l] -
						 pd[l]) * (gmm.gauss[tindex[t]].mean[i] - pd[i]);
					gmm.gauss[tindex[t]].cov[l][i] += diff;
				}
			}
		}

		for (m = 0; m < M; m++)
			for (l = 0; l < L; l++)
				for (i = 0; i <= l; i++) {
					gmm.gauss[m].cov[l][i] /= (double) cntcb[m];
				}

		/* masking */
		if (multiple_dim == TR) {
			maskCov_GMM(&gmm, dim_list, cov_dim, block_full, block_corr);
		}
	}
	floorVar_GMM(&gmm, V);

	/* end of initialization */

	fprintf(stderr, "    start EM:\n");
	/* EM training of GMM parameters */
	for (i = 0; (i <= Imax) && ((i <= Imin) || (fabs(change) > E)); i++) {
		fillz_GMM(&tgmm);
		fillz(sum, sizeof(double), M);

		if (full != TR) {
			for (m = 0; m < M; m++)
				gmm.gauss[m].gconst = cal_gconst(gmm.gauss[m].var, L);
		} else {
			for (m = 0; m < M; m++) {
				gmm.gauss[m].gconst = cal_gconstf(gmm.gauss[m].cov, L);
				if (gmm.gauss[m].gconst == LZERO) {
					fprintf(stderr, "    ERROR : Can't caluculate covdet\n");
					free_GMM(&gmm);
					free_GMM(&tgmm);
					free(logwgd);
					free(sum);
					free(cb);
					free(icb);
					free(tindex);
					free(cntcb);
					return 2;
				}
			}
		}
		if (full == TR) {
			prepareCovInv_GMM(&gmm);
		}

		for (t = 0, ave_logp1 = 0.0, pd = dat; t < T; t++, pd += L) {
			for (m = 0, logb = LZERO; m < M; m++) {
				logwgd[m] = log_wgd(&gmm, m, 0, L, pd);
				logb = log_add(logb, logwgd[m]);
			}
			ave_logp1 += logb;

			for (m = 0; m < M; m++) {
				tmp1 = exp(logwgd[m] - logb);
				sum[m] += tmp1;

				for (l = 0; l < L; l++) {
					tmp2 = tmp1 * pd[l];
					tgmm.gauss[m].mean[l] += tmp2;
					if (full != TR)
						tgmm.gauss[m].var[l] += tmp2 * pd[l];
					else {
						for (j = 0; j <= l; j++) {
							tgmm.gauss[m].cov[l][j] +=
								tmp1 * (pd[l] - gmm.gauss[m].mean[l]) * (pd[j] -
										gmm.
										gauss[m].mean
										[j]);
						}
					}
				}
			}
		}

		/* Output average log likelihood at each iteration */
		ave_logp1 /= (double) T;
		if (i == 1 && M == 1)
			ave_logp0 = ave_logp1;

		fprintf(stderr, "    iter %3d : ", i);
		fprintf(stderr, "ave_logprob = %g", ave_logp1);
		if (i) {
			change = ave_logp1 - ave_logp0;
			fprintf(stderr, "  change = %g", change);
		}
		fprintf(stderr, "\n");
		ave_logp0 = ave_logp1;

		/* Update perameters */
		/* weights */
		for (m = 0; m < M; m++) {
			gmm.weight[m] = sum[m] / (double) T;
		}

		floorWeight_GMM(&gmm, W);

		/* mean, variance */
		for (m = 0; m < M; m++) {
			for (l = 0; l < L; l++) {
				gmm.gauss[m].mean[l] = tgmm.gauss[m].mean[l] / sum[m];
			}

			if (multiple_dim == TR) {
				if (block_full == FA && block_corr == FA) {
					if (full_cov != TR) {    /* -f is not specified */
						for (l = 0; l < L; l++) {
							gmm.gauss[m].cov[l][l] = tgmm.gauss[m].cov[l][l] / sum[m];
						}
					} else {
						for (l = 0; l < L; l++) {
							for (j = 0; j <= l; j++) {
								gmm.gauss[m].cov[l][j] =
									tgmm.gauss[m].cov[l][j] / sum[m];
							}
						}
					}
				} else {
					/* for each block (lower triangle) */
					offset_row = 0;
					for (row = 0; row < cov_dim; row++) {    /* row block number */
						offset_col = 0;
						for (col = 0; col <= row; col++) {    /* column block number */
							if (dim_list[row] == dim_list[col]) {      /* block is square */
								if (block_full == FA && block_corr == TR) {
									/* blockwise diagonal */
									for (k = offset_row, l = offset_col;
											k < offset_row + dim_list[row]; k++, l++) {
										gmm.gauss[m].cov[k][l] =
											tgmm.gauss[m].cov[k][l] / sum[m];
									}
								} else {        /* block_full is TR */
									for (k = offset_row; k < offset_row + dim_list[row];
											k++) {
										for (l = offset_col;
												l < offset_col + dim_list[col]; l++) {
											if (row == col) {
												if (l <= k) {
													gmm.gauss[m].cov[k][l] =
														tgmm.gauss[m].cov[k][l] / sum[m];
												}
											} else {
												if (block_corr == TR) {
													gmm.gauss[m].cov[k][l] =
														tgmm.gauss[m].cov[k][l] / sum[m];
												}
											}
										}
									}
								}
							}
							offset_col += dim_list[col];
						}
						offset_row += dim_list[row];
					}
				}
			} else {
				if (full != TR) {
					for (l = 0; l < L; l++) {
						gmm.gauss[m].var[l] = tgmm.gauss[m].var[l] / sum[m]
							- gmm.gauss[m].mean[l] * gmm.gauss[m].mean[l];
					}
				} else {
					for (l = 0; l < L; l++) {
						for (j = 0; j <= l; j++) {
							gmm.gauss[m].cov[l][j] = tgmm.gauss[m].cov[l][j] / sum[m];
						}
					}
				}
			}
		}
		floorVar_GMM(&gmm, V);
	}
	/*
	test_sum=0;
	for (m=0;m<M;m++)
	{
		test_sum+=gmm.weight[m];
		fprintf(stderr,"the %dth weight is %lf\n",m,gmm.weight[m]);

	}
	fprintf(stderr,"the sum of weight is %lf\n",test_sum);
	*/
	memcpy(out_weight, gmm.weight, M*sizeof(*(gmm.weight)));
	for (m = 0; m < M; m++) {
			memcpy(out_mean + m * L, gmm.gauss[m].mean, L*sizeof(*(gmm.gauss[m].mean)));
			for (i = 0; i < L; i++) {
				for (j = 0; j < i; j++) {
					gmm.gauss[m].cov[j][i] = gmm.gauss[m].cov[i][j];
				}
			}
			for (i = 0; i < L; i++) {
				memcpy(out_cov + m*L*L + i*L,gmm.gauss[m].cov[i],L*sizeof(*(gmm.gauss[m].cov[i])));
			}
	}
	free_GMM(&gmm);
	free_GMM(&tgmm);
	free(logwgd);
	free(sum);
	free(cb);
	free(icb);
	free(tindex);
	free(cntcb);
	return 0;
}
