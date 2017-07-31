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

/*********************************************************************
*                                                                    *
*    $Id: _vc.c,v 1.10 2016/12/22 10:53:14 fjst15124 Exp $               *
*                                                                    *
*    module for vc command                                           *
*                                                                    *
**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if defined(WIN32)
#include "SPTK.h"
#else
#include <SPTK.h>
#endif

#include "hts_engine_API/HTS_engine.h"
#include "hts_engine_API/HTS_hidden.h"

/* perform conversion */
int vc(const GMM * gmm, const DELTAWINDOW * window, const size_t total_frame,
       const size_t source_vlen, const size_t target_vlen,
       const double *gv_mean, const double *gv_vari,
       const double *source, double *target)
{
   size_t t, i, j, k, max_num_mix = 0,
       src_vlen_dyn = source_vlen * window->win_size,
       tgt_vlen_dyn = target_vlen * window->win_size;
   int m, l, shift;
   double max_post_mix = 0.0, logoutp = LZERO, *input = NULL,
       *src_with_dyn = NULL, *logwgd = NULL,
       **cov_xx_inv = NULL, ***cov_yx_xx = NULL, *gv_weight = NULL,
       ***cond_mean = NULL, ***cond_vari = NULL, **cond_post_mix = NULL;
   GMM gmm_xx;
   HTS_SStreamSet sss;
   HTS_PStreamSet pss;

   /* append dynamic feature */
   src_with_dyn = dgetmem(total_frame * src_vlen_dyn);
   for (t = 0; t < total_frame; t++) {
      for (i = 0; i < window->win_size; i++) {
         j = window->win_size * source_vlen * t + source_vlen * i;
         for (shift = window->win_l_width[i];
              shift <= window->win_r_width[i]; shift++) {
            l = t + shift;
            if (l < 0) {
               l = 0;
            }
            if (!(l < (int) total_frame)) {
               l = total_frame - 1;
            }
            for (k = 0; k < source_vlen; k++) {
               src_with_dyn[j + k] += window->win_coefficient[i][shift]
                   * source[source_vlen * l + k];
            }
         }
      }
   }

   /* calculate mean and covariace of conditional distribution
      given source feature and mixture component */
   cond_post_mix = ddgetmem(total_frame, gmm->nmix);
   cond_mean = (double ***) getmem(gmm->nmix, sizeof(*(cond_mean)));
   for (m = 0; m < gmm->nmix; m++) {
      cond_mean[m] = ddgetmem(total_frame, tgt_vlen_dyn);
   }
   cond_vari = (double ***) getmem(gmm->nmix, sizeof(*(cond_vari)));
   for (m = 0; m < gmm->nmix; m++) {
      cond_vari[m] = ddgetmem(tgt_vlen_dyn, tgt_vlen_dyn);
   }
   cov_xx_inv = ddgetmem(src_vlen_dyn, src_vlen_dyn);
   cov_yx_xx = (double ***) getmem(gmm->nmix, sizeof(*(cov_yx_xx)));
   for (m = 0; m < gmm->nmix; m++) {
      cov_yx_xx[m] = ddgetmem(tgt_vlen_dyn, src_vlen_dyn);
   }
   for (m = 0; m < gmm->nmix; m++) {
      invert(gmm->gauss[m].cov, cov_xx_inv, src_vlen_dyn);
      for (i = 0; i < tgt_vlen_dyn; i++) {
         for (j = 0; j < src_vlen_dyn; j++) {
            for (k = 0; k < src_vlen_dyn; k++) {
               cov_yx_xx[m][i][j] += gmm->gauss[m].cov[src_vlen_dyn + i][k]
                   * cov_xx_inv[k][j];
            }
         }
      }
   }
   logwgd = dgetmem(gmm->nmix);
   input = dgetmem(src_vlen_dyn);
   alloc_GMM(&gmm_xx, gmm->nmix, src_vlen_dyn, gmm->full);
   for (i = 0; i < (size_t) gmm_xx.nmix; i++) {
      gmm_xx.weight[i] = gmm->weight[i];
      for (j = 0; j < (size_t) gmm_xx.dim; j++) {
         gmm_xx.gauss[i].mean[j] = gmm->gauss[i].mean[j];
         if (gmm_xx.full) {
            for (k = 0; k < (size_t) gmm_xx.dim; k++) {
               gmm_xx.gauss[i].cov[j][k] = gmm->gauss[i].cov[j][k];
            }
         } else {
            gmm_xx.gauss[i].var[j] = gmm->gauss[i].var[j];
         }
      }
   }
   for (i = 0; i < (size_t) gmm_xx.nmix; i++) {
      invert(gmm_xx.gauss[i].cov, gmm_xx.gauss[i].inv, src_vlen_dyn);
      gmm_xx.gauss[i].gconst = cal_gconstf(gmm_xx.gauss[i].cov, src_vlen_dyn);
   }
   for (t = 0; t < total_frame; t++) {
      for (i = 0; i < src_vlen_dyn; i++) {
         input[i] = src_with_dyn[t * src_vlen_dyn + i];
      }
      for (m = 0, logoutp = LZERO; m < gmm->nmix; m++) {
         logwgd[m] = log_wgd(&gmm_xx, m, 0, src_vlen_dyn, input);
         logoutp = log_add(logoutp, logwgd[m]);
      }
      for (m = 0; m < gmm->nmix; m++) {
         /* posterior probability of mixture component given source feature */
         cond_post_mix[t][m] = exp(logwgd[m] - logoutp);
         for (i = 0; i < tgt_vlen_dyn; i++) {
            for (j = 0; j < src_vlen_dyn; j++) {
               cond_mean[m][t][i] += cov_yx_xx[m][i][j]
                   * (input[j] - gmm->gauss[m].mean[j]);
            }
            cond_mean[m][t][i] += gmm->gauss[m].mean[src_vlen_dyn + i];
         }
      }
   }

   for (m = 0; m < gmm->nmix; m++) {
      for (i = 0; i < tgt_vlen_dyn; i++) {
         for (j = 0; j < tgt_vlen_dyn; j++) {
            for (k = 0; k < src_vlen_dyn; k++) {
               cond_vari[m][i][j] += cov_yx_xx[m][i][k]
                   * gmm->gauss[m].cov[k][src_vlen_dyn + j];
            }
            cond_vari[m][i][j] =
                gmm->gauss[m].cov[src_vlen_dyn + i][src_vlen_dyn + j]
                - cond_vari[m][i][j];
         }
      }
   }

   /* initialize parameter set of hts_engine */
   HTS_PStreamSet_initialize(&pss);
   sss.nstream = 1;
   sss.total_state = total_frame;
   sss.total_frame = total_frame;
   sss.duration = (size_t *) getmem(total_frame, sizeof(size_t));
   for (i = 0; i < total_frame; i++) {
      sss.duration[i] = 1;
   }
   sss.sstream = (HTS_SStream *) getmem(1, sizeof(HTS_SStream));
   sss.sstream->vector_length = target_vlen;
   sss.sstream->mean =
       (double **) getmem(sss.total_state, sizeof(*(sss.sstream->mean)));
   sss.sstream->vari =
       (double **) getmem(sss.total_state, sizeof(*(sss.sstream->vari)));
   for (i = 0; i < sss.total_state; i++) {
      sss.sstream->mean[i] = dgetmem(tgt_vlen_dyn);
      sss.sstream->vari[i] = dgetmem(tgt_vlen_dyn);
   }
   sss.sstream->msd = NULL;     /* no MSD */
   sss.sstream->win_size = window->win_size;
   sss.sstream->win_l_width =
       (int *) getmem(window->win_size, sizeof(*(sss.sstream->win_l_width)));
   sss.sstream->win_r_width =
       (int *) getmem(window->win_size, sizeof(*(sss.sstream->win_r_width)));
   sss.sstream->win_coefficient =
       (double **) getmem(window->win_size,
                          sizeof(*(sss.sstream->win_coefficient)));
   for (i = 0; i < window->win_size; i++) {
      sss.sstream->win_l_width[i] = window->win_l_width[i];
      sss.sstream->win_r_width[i] = window->win_r_width[i];
      if (sss.sstream->win_l_width[i] + sss.sstream->win_r_width[i] == 0) {
         sss.sstream->win_coefficient[i] =
             dgetmem(-2 * sss.sstream->win_l_width[i] + 1);
      } else {
         sss.sstream->win_coefficient[i] =
             dgetmem(-2 * sss.sstream->win_l_width[i]);
      }
      sss.sstream->win_coefficient[i] -= sss.sstream->win_l_width[i];
      for (shift = sss.sstream->win_l_width[i];
           shift <= sss.sstream->win_r_width[i]; shift++) {
         sss.sstream->win_coefficient[i][shift] =
             window->win_coefficient[i][shift];
      }
   }
   sss.sstream->win_max_width = window->win_max_width;
   if ((gv_mean != NULL) && (gv_vari != NULL)) {        /* set GV parameters */
      sss.sstream->gv_mean = dgetmem(sss.sstream->vector_length);
      sss.sstream->gv_vari = dgetmem(sss.sstream->vector_length);
      for (i = 0; i < sss.sstream->vector_length; i++) {
         sss.sstream->gv_mean[i] = gv_mean[i];
         sss.sstream->gv_vari[i] = gv_vari[i];
      }
   } else {
      sss.sstream->gv_mean = NULL;
      sss.sstream->gv_vari = NULL;
   }
   sss.sstream->gv_switch =
       (HTS_Boolean *) getmem(total_frame, sizeof(HTS_Boolean));
   for (i = 0; i < total_frame; i++) {
      sss.sstream->gv_switch[i] = TRUE;
   }
   gv_weight = dgetmem(tgt_vlen_dyn);
   for (i = 0; i < tgt_vlen_dyn; i++) {
      gv_weight[i] = 1.0;
   }

   /* initialize pdf sequence */
   for (t = 0; t < total_frame; t++) {
      max_post_mix = cond_post_mix[t][0];
      max_num_mix = 0;
      for (m = 1; m < gmm->nmix; m++) {
         if (max_post_mix < cond_post_mix[t][m]) {
            max_post_mix = cond_post_mix[t][m];
            max_num_mix = m;
         }
      }
      for (i = 0; i < tgt_vlen_dyn; i++) {
         sss.sstream->mean[t][i] = cond_mean[max_num_mix][t][i];
         sss.sstream->vari[t][i] = cond_vari[max_num_mix][i][i];
      }
   }

   /* parameter generation by hts_engine API */
   HTS_PStreamSet_create(&pss, &sss, NULL, gv_weight);
   for (t = 0; t < total_frame; t++) {
      k = t * target_vlen;
      for (i = 0; i < target_vlen; i++) {
         target[k + i] = pss.pstream->par[t][i];
      }
   }

   /* release memory */
   free(src_with_dyn);
   free(input);
   free(logwgd);
   free(cov_xx_inv[0]);
   free(cov_xx_inv);
   for (m = 0; m < gmm->nmix; m++) {
      free(cov_yx_xx[m][0]);
      free(cov_yx_xx[m]);
      free(cond_mean[m][0]);
      free(cond_mean[m]);
      free(cond_vari[m][0]);
      free(cond_vari[m]);
   }
   free(cov_yx_xx);
   free(cond_mean);
   free(cond_vari);
   free(cond_post_mix[0]);
   free(cond_post_mix);
   free(gv_weight);
   free_GMM(&gmm_xx);
   HTS_PStreamSet_clear(&pss);
   HTS_SStreamSet_clear(&sss);

   return (0);
}


int vc_run(double *source,double *gmm_weight,double *gmm_mean,double *gmm_cov,double *gv_mean_in,double *gv_vari_in,int source_vlen,int target_vlen,int total_frame,int num_mix,int delta_order,int gv_flag,double FLOOR,double *target)
{
	char *coef =NULL;
	int i,j,k,m,len_total,dw_num=1,win_max_width=0;
	double floor=FLOOR;
	double *gv_mean=NULL,*gv_vari=NULL;
	Boolean full =TR;
	GMM gmm;
	DELTAWINDOW window;
	
	dw_num+=delta_order;
	if(gv_flag!=0)
	{
		gv_mean=gv_mean_in;
		gv_vari=gv_vari_in;
	}
	len_total=(source_vlen+target_vlen)*dw_num;

	/* load GMM parameters */
	alloc_GMM(&gmm, num_mix, len_total, full);
	memcpy(gmm.weight, gmm_weight, num_mix*sizeof(*(gmm_weight)));
	for (m = 0; m < num_mix; m++) {
		memcpy(gmm.gauss[m].mean, gmm_mean + m * len_total, len_total*sizeof(*(gmm_mean)));
		for (i = 0; i < len_total; i++) {
			memcpy(gmm.gauss[m].cov[i],gmm_cov + m*len_total*len_total + i*len_total,len_total*sizeof(*(gmm.gauss[m].cov[i])));
		}
	}
	prepareCovInv_GMM(&gmm);
	prepareGconst_GMM(&gmm);

	/* flooring for diagonal component of covariance */
	if (floor != 0.0) {
		for (i = 0; i < num_mix; i++) {
			for (j = 0; j < (int) len_total; j++) {
				gmm.gauss[i].cov[j][j] += floor;
			}
		}
	}

	/* set window parameters*/
	window.win_size = dw_num;
	window.win_l_width =
	(int *) getmem(window.win_size, sizeof(*(window.win_l_width)));
	window.win_r_width =
	(int *) getmem(window.win_size, sizeof(*(window.win_r_width)));
	window.win_coefficient =
	(double **) getmem(window.win_size, sizeof(*(window.win_coefficient)));
	window.win_l_width[0] = 0;
	window.win_r_width[0] = 0;
	window.win_coefficient[0] = dgetmem(1);
	window.win_coefficient[0][0] = 1.0;
	if(delta_order>0)
	{
		int a0, a1, a2, dw_leng;
		for (i = 1; i < window.win_size; i++) {
			dw_leng = 1;
			window.win_l_width[i] = -dw_leng;
			window.win_r_width[i] = dw_leng;
			window.win_coefficient[i] = dgetmem(dw_leng * 2 + 1);
			window.win_coefficient[i] += dw_leng;
		}
		for (a1 = 0, j = -dw_leng; j <= dw_leng; a1 += j * j, j++);
		for (j = -dw_leng; j <= dw_leng; j++) {
			window.win_coefficient[1][j] = (double) j / (double) a1;
		}
	
		if (window.win_size > 2) {
			dw_leng = 1;
			for (a0 = a1 = a2 = 0, j = -dw_leng; j <= dw_leng;
				a0++, a1 += j * j, a2 += j * j * j * j, j++);
			for (j = -dw_leng; j <= dw_leng; j++) {
				window.win_coefficient[2][j]
				= 2 * ((double) (a0 * j * j - a1)) /
				((double) (a2 * a0 - a1 * a1));
			}
		}
	
	}
	win_max_width = window.win_r_width[0];       /* width of static window is 0 */
	for (i = 1; i < window.win_size; i++) {
		if (win_max_width < window.win_r_width[i]) {
			win_max_width = window.win_r_width[i];
		}
		if (win_max_width < -window.win_l_width[i]) {
			win_max_width = -window.win_l_width[i];
		}
	}
	window.win_max_width = win_max_width;
	
	for (m = 0; m < gmm.nmix; m++) {
		invert(gmm.gauss[m].cov, gmm.gauss[m].inv, source_vlen * window.win_size);
		if (full) {
			gmm.gauss[m].gconst =
				cal_gconstf(gmm.gauss[m].cov, source_vlen * window.win_size);
		} else {
			gmm.gauss[m].gconst =
				cal_gconst(gmm.gauss[m].var, source_vlen * window.win_size);
		}
	}

	/* perform conversion */
	vc(&gmm, &window, total_frame, source_vlen, target_vlen, gv_mean, gv_vari,source, target);

	free_GMM(&gmm);
	for (i = 0; i < window.win_size; i++) {
		free(window.win_coefficient[i] + window.win_l_width[i]);
	}
	free(window.win_l_width);
	free(window.win_r_width);

}





