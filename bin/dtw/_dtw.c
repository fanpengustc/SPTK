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

/*******************************************************************************
*                                                                              *
*    Dynamic Time Warping                                                      *
*                                                                              *
*                                      2011.12 Akira Tamamori                  *
*       usage:                                                                 *
*       		int dtw(double *x, double *y, int dim, int num_test, int num_ref, int path_type, int norm_type,int *outviterbi, double *z)
*       options:                                                               *
*               -m M          : order of vector                        [24]    *
*               -l L          : dimension of vector                    [m+1]   *
*               -t T          : number of test vectors                 [N/A]   *
*               -r R          : number of reference vectors            [N/A]   *
*               -n N          : type of norm used for calculation      [2]     *
*                               of local distance                              *
*                                 1 : L1-norm                                  *
*                                 2 : L2-norm                                  *
*               -p P          : local path constraint                  [5]     *
*               -s sfile      : output score of dynamic time warping   [FALSE] *
*               -v out_vfile  : output concatenated test/reference     [FALSE] *
*                               data sequence along the Vitebi path            *
*               -V in_vfile   : concatenate test and reference vectors [FALSE] *
*                               in accordance with the Vitebi path             *
*                               information written in in_vfile                *
*       infile:                                                                *
*              test vector sequence                                            *
*                  x_1(1), ..., x_1(L), x_2(1), ..., x_2(L), ...               *
*       reffile:                                                               *
*              reference vector sequence                                       *
*                  y_1(1), ..., y_1(L), y_2(1), ..., y_2(L), ...               *
*       stdout:                                                                *
*              concatenated test/reference vector sequence                     *
*              along the Vitebi path                                           *
*                  x_1(1), ..., x_1(L), y_1(1), ..., y_1(L), ...               *
*                                                                              *
*******************************************************************************/

static char *rcs_id = "$Id: dtw.c,v 1.14 2016/12/22 10:53:02 fjst15124 Exp $";

/*  Standard C Libraries  */
#include <stdio.h>

#ifdef HAVE_STRING_H
#include <string.h>
#else
#include <strings.h>
#ifndef HAVE_STRRCHR
#define strrchr rindex
#endif
#endif

#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#if defined(WIN32)
#include "SPTK.h"
#else
#include <SPTK.h>
#endif

/*  Default Values  */
#define  LENG  24

/*  Command Name  */

#define PATH_NG FA
#define PATH_OK TR

typedef struct _Dtw_Cell {
   double local;                /* local cost */
   double global;               /* global cost */
   int backptr[2];              /* back pointer for Viterbi path */
   Boolean is_region;
   int allow_path;
} Dtw_Cell;

typedef struct _Data {
   double *input;
   int total;                   /* total number of vectors */
   int dim;                     /* dimension of input vector */
   int *viterbi;                /* Viterbi path */
} Data;

typedef struct _Dtw_Table {
   Dtw_Cell **cell;
   Data data[2];                /* two comparative data */
   int vit_leng;                /* length of Viterbi path */
   int path;                    /* type of local constraint */
   int norm;                    /* type of norm for local cost */
   double *weight;              /* weight on the local path */
} Dtw_Table;

typedef struct _Float_List {
   float *f;
   struct _Float_List *next;
} Float_List;

static int round_up(double dat)
{
   return (int) (dat + 0.5);
}

Dtw_Cell **malloc_Dtw_Cell(int size1, int size2)
{
   Dtw_Cell **tmpcell, *tmpcell2;
   int i, j;

   tmpcell = (Dtw_Cell **) malloc(sizeof(Dtw_Cell *) * size1);
   if (tmpcell == NULL) {
      fprintf(stderr, "ERROR: Can't allocate memory !\n");
   }

   tmpcell2 = (Dtw_Cell *) malloc(sizeof(Dtw_Cell) * size1 * size2);
   if (tmpcell2 == NULL) {
      fprintf(stderr, "ERROR: Can't allocate memory !\n");
   }

   for (i = 0, j = 0; i < size1; i++, j += size2) {
      tmpcell[i] = tmpcell2 + j;
   }

   return (tmpcell);
}

int init_dtw(Dtw_Table * table, int leng, double *input1, double *input2,
              int total1, int total2, int path, int norm)
{

   int i, size[2] = { total1, total2 };

   if (path == 3 || path == 4) {
      if (total2 > total1) {
         fprintf(stderr, "Can't perform DTW !\n"
                 "The number of the reference vectors (= %d) must be less than "
                 "the number of the test vectors (= %d). \n", total2, total1);
      }
   } else if (path == 5 || path == 6 || path == 7) {
      if (total1 / 2 >= total2) {
         fprintf(stderr, "Can't perform DTW !\n"
                 "The number of the test vectors (= %d) must be less than "
                 "the twice of the reference vectors (= 2 * %d = %d). \n",
                 total1, total2, 2 * total2);
      } else if (total2 / 2 >= total1) {
         fprintf(stderr, "Can't perform DTW !\n"
                 "The number of the reference vectors (= %d) must be less than "
                 "the twice of the test vectors (= 2 * %d = %d). \n",
                 total2, total1, 2 * total1);
      }
   }

   table->cell = (Dtw_Cell **) malloc_Dtw_Cell(size[0], size[1]);

   table->data[0].input = input1;
   table->data[1].input = input2;
   for (i = 0; i < 2; i++) {
      table->data[i].total = size[i];
      table->data[i].dim = leng;
      table->data[i].viterbi =
          (int *) malloc(sizeof(int) * (size[0] + size[1]));
      if (table->data[i].viterbi == NULL) {
         fprintf(stderr, "ERROR: Can't allocate memory at init_dtw() !\n");
		 return 1;
      }
   }

   table->path = path;

   if (norm != 1 && norm != 2) {
      fprintf(stderr, "dtw : type of norm must be 1 or 2!\n");
   }
   table->norm = norm;

   if (path < 1 || path > 7) {
      fprintf(stderr, "dtw : local path constraint must be between 1 and 7!\n");
   }

   switch (path) {
   case 1:
      table->weight = dgetmem(2);
      table->weight[0] = 1.0;
      table->weight[1] = 1.0;
      break;
   case 2:
      table->weight = dgetmem(3);
      table->weight[0] = 1.0;
      table->weight[1] = 2.0;
      table->weight[2] = 1.0;
      break;
   case 3:
      table->weight = dgetmem(2);
      table->weight[0] = 1.0;
      table->weight[1] = 2.0;
      break;
   case 4:
      table->weight = dgetmem(3);
      table->weight[0] = 1.0;
      table->weight[1] = 2.0;
      table->weight[2] = 3.0;
      break;
   case 5:
      table->weight = dgetmem(5);
      table->weight[0] = 2.0;
      table->weight[1] = 1.0;
      table->weight[2] = 2.0;
      table->weight[3] = 2.0;
      table->weight[4] = 1.0;
      break;
   case 6:
      table->weight = dgetmem(3);
      table->weight[0] = 3.0;
      table->weight[1] = 2.0;
      table->weight[2] = 3.0;
      break;
   case 7:
      table->weight = dgetmem(6);
      table->weight[0] = 1.0;
      table->weight[1] = 1.0;
      table->weight[2] = 1.0;
      table->weight[3] = 1.0;
      table->weight[4] = 1.0;
      table->weight[5] = 1.0;
      break;
   default:
      break;
   }
   return 0;
}

void check_enabled_region_type_1(Dtw_Table * table)
{
   int i, j, tx = table->data[0].total, ty = table->data[1].total;
   for (j = 0; j < ty; j++) {
      for (i = 0; i < tx; i++) {
         table->cell[i][j].is_region = PATH_OK;
      }
   }
}

void check_enabled_region_type_2(Dtw_Table * table)
{
   int i, j, tx = table->data[0].total, ty = table->data[1].total;
   for (j = 0; j < ty; j++) {
      for (i = 0; i < tx; i++) {
         table->cell[i][j].is_region = PATH_OK;
      }
   }
}

void check_enabled_region_type_3(Dtw_Table * table)
{
   int i, j, tx = table->data[0].total, ty = table->data[1].total,
       range = tx - ty;
   Dtw_Cell **cell = table->cell;

   cell[0][0].is_region = PATH_OK;

   for (j = 1; j < ty; j++) {
      cell[0][j].is_region = PATH_NG;
   }
   for (i = 1; i < tx; i++) {
      cell[i][0].is_region = PATH_NG;
   }
   for (i = 1; i <= range; i++) {
      cell[i][0].is_region = PATH_OK;
   }

   for (j = 1; j < ty; j++) {
      for (i = 1; i < j; i++) {
         cell[i][j].is_region = PATH_NG;
      }
      for (i = j; i <= range + j; i++) {
         cell[i][j].is_region = PATH_OK;
      }
      for (; i < tx; i++) {
         cell[i][j].is_region = PATH_NG;
      }
   }

   for (j = 1; j < ty; j++) {
      for (i = 1; i < tx; i++) {
         if (cell[i][j].is_region == PATH_OK) {
            if (cell[i - 1][j].is_region == PATH_OK &&
                cell[i - 1][j - 1].is_region == PATH_OK) {
               cell[i][j].allow_path = 1;
            } else if (cell[i - 1][j].is_region == PATH_OK) {
               cell[i][j].allow_path = 2;
            } else if (cell[i - 1][j - 1].is_region == PATH_OK) {
               cell[i][j].allow_path = 3;
            }
         }
      }
   }
}

void check_enabled_region_type_4(Dtw_Table * table)
{
   int i, j, tx = table->data[0].total, ty = table->data[1].total;
   Dtw_Cell **cell = table->cell;

   cell[0][0].is_region = PATH_OK;

   for (j = 1; j < ty; j++) {
      cell[0][j].is_region = PATH_NG;
   }

   for (i = 1; i < tx - ty / 2; i++) {
      cell[i][0].is_region = PATH_OK;
   }
   for (; i < tx; i++) {
      cell[i][0].is_region = PATH_NG;
   }
   for (i = 1; i < tx - ty / 2 + 1; i++) {
      cell[i][1].is_region = PATH_OK;
   }
   for (; i < tx; i++) {
      cell[i][1].is_region = PATH_NG;
   }

   for (j = 2; j < ty; j++) {
      for (i = 0; i < round_up((double) j / 2); i++) {
         cell[i][j].is_region = PATH_NG;
      }
      for (; i < tx - ty / 2 + round_up((double) j / 2); i++) {
         cell[i][j].is_region = PATH_OK;
      }
      for (; i < tx; i++) {
         cell[i][j].is_region = PATH_NG;
      }
   }
   cell[tx - 1][ty - 1].is_region = PATH_OK;

   for (j = 2; j < ty; j++) {
      for (i = 2; i < tx; i++) {
         if (cell[i][j].is_region == PATH_OK) {
            if (cell[i - 1][j].is_region == PATH_OK &&
                cell[i - 1][j - 1].is_region == PATH_OK &&
                cell[i - 1][j - 2].is_region == PATH_OK) {
               cell[i][j].allow_path = 1;
            } else if (cell[i - 1][j].is_region == PATH_OK &&
                       cell[i - 1][j - 1].is_region == PATH_OK) {
               cell[i][j].allow_path = 2;
            } else if (cell[i - 1][j - 1].is_region == PATH_OK &&
                       cell[i - 1][j - 2].is_region == PATH_OK) {
               cell[i][j].allow_path = 3;
            } else if (cell[i - 1][j - 1].is_region == PATH_OK) {
               cell[i][j].allow_path = 4;
            } else if (cell[i - 1][j - 2].is_region == PATH_OK) {
               cell[i][j].allow_path = 5;
            }
         }
      }
   }
}

void check_enabled_region_type_5(Dtw_Table * table)
{
   int i, j, tx = table->data[0].total, ty = table->data[1].total;
   Dtw_Cell **cell = table->cell;

   cell[0][0].is_region = PATH_OK;

   for (j = 1; j < ty; j++) {
      for (i = 1; i < tx; i++) {
         cell[i][j].is_region = PATH_OK;
      }
   }

   for (i = 1; i < tx; i++) {
      cell[i][0].is_region = PATH_NG;
   }

   for (j = 1; j < ty; j++) {
      cell[0][j].is_region = PATH_NG;
   }

   for (j = 1; j < ty - 2; j++) {
      for (i = 2 * j + 1; i < tx; i++) {
         cell[i][j].is_region = PATH_NG;
      }
   }
   for (i = 1; i < tx - 2; i++) {
      for (j = 2 * i + 1; j < ty; j++) {
         cell[i][j].is_region = PATH_NG;
      }
   }

   /* Backward pruning */
   for (j = ty - 1; j > 0; j--) {
      for (i = 2 * (j - ty) + tx; i > 0; i--) {
         cell[i][j].is_region = PATH_NG;
      }
   }
   for (i = tx - 1; i > 0; i--) {
      for (j = 2 * (i - tx) + ty; j > 0; j--) {
         cell[i][j].is_region = PATH_NG;
      }
   }

   for (j = 2; j < ty - 1; j++) {
      for (i = 2; i < tx - 1; i++) {
         if (cell[i][j].is_region == PATH_OK) {
            if (cell[i - 2][j - 1].is_region == PATH_OK &&
                cell[i - 1][j - 1].is_region == PATH_OK &&
                cell[i - 1][j - 2].is_region == PATH_OK) {
               cell[i][j].allow_path = 1;
            } else if (cell[i - 2][j - 1].is_region == PATH_OK &&
                       cell[i - 1][j - 1].is_region == PATH_OK) {
               cell[i][j].allow_path = 2;
            } else if (cell[i - 1][j - 1].is_region == PATH_OK &&
                       cell[i - 1][j - 2].is_region == PATH_OK) {
               cell[i][j].allow_path = 3;
            } else if (cell[i - 2][j - 1].is_region == PATH_OK) {
               cell[i][j].allow_path = 4;
            } else if (cell[i - 1][j - 1].is_region == PATH_OK) {
               cell[i][j].allow_path = 5;
            } else if (cell[i - 1][j - 2].is_region == PATH_OK) {
               cell[i][j].allow_path = 6;
            }
         }
      }
   }
}

void check_enabled_region_type_6(Dtw_Table * table)
{
   int i, j, tx = table->data[0].total, ty = table->data[1].total;
   Dtw_Cell **cell = table->cell;

   cell[0][0].is_region = PATH_OK;

   for (j = 1; j < ty; j++) {
      for (i = 1; i < tx; i++) {
         cell[i][j].is_region = PATH_OK;
      }
   }

   for (i = 1; i < tx; i++) {
      cell[i][0].is_region = PATH_NG;
   }

   for (j = 1; j < ty; j++) {
      cell[0][j].is_region = PATH_NG;
   }

   for (j = 1; j < ty - 2; j++) {
      for (i = 2 * j + 1; i < tx; i++) {
         cell[i][j].is_region = PATH_NG;
      }
   }
   for (i = 1; i < tx - 2; i++) {
      for (j = 2 * i + 1; j < ty; j++) {
         cell[i][j].is_region = PATH_NG;
      }
   }

   /* Backward pruning */
   for (j = ty - 1; j > 0; j--) {
      for (i = tx - 1 + 2 * (j - ty) + 1; i > 0; i--) {
         cell[i][j].is_region = PATH_NG;
      }
   }
   for (i = tx - 1; i > 0; i--) {
      for (j = ty - 1 + 2 * (i - tx) + 1; j > 0; j--) {
         cell[i][j].is_region = PATH_NG;
      }
   }

   for (j = 2; j < ty; j++) {
      for (i = 2; i < tx; i++) {
         if (cell[i][j].is_region == PATH_OK) {
            if (cell[i - 2][j - 1].is_region == PATH_OK &&
                cell[i - 1][j - 1].is_region == PATH_OK &&
                cell[i - 1][j - 2].is_region == PATH_OK) {
               cell[i][j].allow_path = 1;
            } else if (cell[i - 2][j - 1].is_region == PATH_OK &&
                       cell[i - 1][j - 1].is_region == PATH_OK) {
               cell[i][j].allow_path = 2;
            } else if (cell[i - 1][j - 1].is_region == PATH_OK &&
                       cell[i - 1][j - 2].is_region == PATH_OK) {
               cell[i][j].allow_path = 3;
            } else if (cell[i - 2][j - 1].is_region == PATH_OK) {
               cell[i][j].allow_path = 4;
            } else if (cell[i - 1][j - 1].is_region == PATH_OK) {
               cell[i][j].allow_path = 5;
            } else if (cell[i - 1][j - 2].is_region == PATH_OK) {
               cell[i][j].allow_path = 6;
            }
         }
      }
   }
}

void check_enabled_region_type_7(Dtw_Table * table)
{
   int i, j, tx = table->data[0].total, ty = table->data[1].total;
   Dtw_Cell **cell = table->cell;

   cell[0][0].is_region = PATH_OK;

   for (j = 1; j < ty; j++) {
      for (i = 1; i < tx; i++) {
         cell[i][j].is_region = PATH_OK;
      }
   }

   for (i = 1; i < tx; i++) {
      cell[i][0].is_region = PATH_NG;
   }

   for (j = 1; j < ty; j++) {
      cell[0][j].is_region = PATH_NG;
   }

   for (j = 1; j < ty - 2; j++) {
      for (i = 2 * j + 1; i < tx; i++) {
         cell[i][j].is_region = PATH_NG;
      }
   }
   for (i = 1; i < tx - 2; i++) {
      for (j = 2 * i + 1; j < ty; j++) {
         cell[i][j].is_region = PATH_NG;
      }
   }

   /* Backward pruning */
   for (j = ty - 1; j > 0; j--) {
      for (i = tx - 1 + 2 * (j - ty) + 1; i > 0; i--) {
         cell[i][j].is_region = PATH_NG;
      }
   }
   for (i = tx - 1; i > 0; i--) {
      for (j = ty - 1 + 2 * (i - tx) + 1; j > 0; j--) {
         cell[i][j].is_region = PATH_NG;
      }
   }

   for (j = 2; j < ty; j++) {
      for (i = 2; i < tx; i++) {
         if (cell[i][j].is_region == PATH_OK) {
            if (cell[i - 2][j - 1].is_region == PATH_OK &&
                cell[i - 2][j - 2].is_region == PATH_OK &&
                cell[i - 1][j - 1].is_region == PATH_OK &&
                cell[i - 1][j - 2].is_region == PATH_OK) {
               cell[i][j].allow_path = 1;
            } else if (cell[i - 2][j - 2].is_region == PATH_OK &&
                       cell[i - 1][j - 1].is_region == PATH_OK &&
                       cell[i - 1][j - 2].is_region == PATH_OK) {
               cell[i][j].allow_path = 2;
            } else if (cell[i - 2][j - 1].is_region == PATH_OK &&
                       cell[i - 2][j - 2].is_region == PATH_OK &&
                       cell[i - 1][j - 1].is_region == PATH_OK) {
               cell[i][j].allow_path = 3;
            } else if (cell[i - 1][j - 1].is_region == PATH_OK &&
                       cell[i - 1][j - 2].is_region == PATH_OK) {
               cell[i][j].allow_path = 4;
            } else if (cell[i - 2][j - 2].is_region == PATH_OK &&
                       cell[i - 1][j - 1].is_region == PATH_OK) {
               cell[i][j].allow_path = 5;
            } else if (cell[i - 2][j - 1].is_region == PATH_OK &&
                       cell[i - 1][j - 1].is_region == PATH_OK) {
               cell[i][j].allow_path = 6;
            } else if (cell[i - 2][j - 1].is_region == PATH_OK) {
               cell[i][j].allow_path = 7;
            } else if (cell[i - 1][j - 1].is_region == PATH_OK) {
               cell[i][j].allow_path = 8;
            } else if (cell[i - 1][j - 2].is_region == PATH_OK) {
               cell[i][j].allow_path = 9;
            }
         }
      }
   }
}

void calc_global_cost_type_1(Dtw_Table * table)
{
   int i, j, tx = table->data[0].total, ty = table->data[1].total;
   double local, path1, path2;
   Dtw_Cell **cell = table->cell;
   double *weight = table->weight;

   for (i = 1; i < tx; i++) {
      cell[i][0].global = cell[i - 1][0].global +weight[0] * cell[i][0].local;
      cell[i][0].backptr[0] = i - 1;
      cell[i][0].backptr[1] = 0;
   }
   for (j = 1; j < ty; j++) {
      cell[0][j].global = cell[0][j - 1].global +weight[1] * cell[0][j].local;
      cell[0][j].backptr[0] = 0;
      cell[0][j].backptr[1] = j - 1;
   }
   for (i = 1; i < tx; i++) {
      for (j = 1; j < ty; j++) {
         local = cell[i][j].local;
         path1 = cell[i - 1][j].global +weight[0] * local;
         path2 = cell[i][j - 1].global +weight[1] * local;
         if (path1 < path2) {
            cell[i][j].global = path1;
            cell[i][j].backptr[0] = i - 1;
            cell[i][j].backptr[1] = j;
         } else {
            cell[i][j].global = path2;
            cell[i][j].backptr[0] = i;
            cell[i][j].backptr[1] = j - 1;
         }
      }
   }
   cell[tx - 1][ty - 1].global /=(tx + ty);
}

void calc_global_cost_type_2(Dtw_Table * table)
{
   int i, j, tx = table->data[0].total, ty = table->data[1].total;
   double local, min, path1, path2, path3;
   Dtw_Cell **cell = table->cell;
   double *weight = table->weight;

   for (i = 1; i < tx; i++) {
      if (cell[i][0].is_region == PATH_OK) {
         cell[i][0].global =
             cell[i - 1][0].global +weight[0] * cell[i][0].local;
         cell[i][0].backptr[0] = i - 1;
         cell[i][0].backptr[1] = 0;
      }
   }
   for (j = 1; j < ty; j++) {
      cell[0][j].global = cell[0][j - 1].global +weight[2] * cell[0][j].local;
      cell[0][j].backptr[0] = 0;
      cell[0][j].backptr[1] = j - 1;
   }
   for (j = 1; j < ty; j++) {
      for (i = 1; i < tx; i++) {
         local = cell[i][j].local;
         path1 = cell[i - 1][j].global +weight[0] * local;
         path2 = cell[i - 1][j - 1].global +weight[1] * local;
         path3 = cell[i][j - 1].global +weight[2] * local;
         cell[i][j].backptr[0] = i - 1;
         cell[i][j].backptr[1] = j;
         min = path1;
         if (min >= path2) {
            min = path2;
            cell[i][j].backptr[0] = i - 1;
            cell[i][j].backptr[1] = j - 1;
         }
         if (min >= path3) {
            min = path3;
            cell[i][j].backptr[0] = i;
            cell[i][j].backptr[1] = j - 1;
         }
         cell[i][j].global = min;
      }
   }
   cell[tx - 1][ty - 1].global /=(tx + ty);
}

void calc_global_cost_type_3(Dtw_Table * table)
{
   int i, j, tx = table->data[0].total, ty = table->data[1].total;
   double local, min = 0.0, path1, path2;
   Dtw_Cell **cell = table->cell;
   double *weight = table->weight;

   for (i = 1; i < tx; i++) {
      if (cell[i][0].is_region == PATH_OK) {
         cell[i][0].global =
             cell[i - 1][0].global +weight[0] * cell[i][0].local;
         cell[i][0].backptr[0] = i - 1;
         cell[i][0].backptr[1] = 0;
      }
   }
   for (j = 1; j < ty; j++) {
      for (i = 1; i < tx; i++) {
         local = cell[i][j].local;
         if (cell[i][j].is_region == PATH_OK) {
            path1 = cell[i - 1][j].global +weight[0] * local;
            path2 = cell[i - 1][j - 1].global +weight[1] * local;

            switch (cell[i][j].allow_path) {
            case 1:
               min = path1;
               cell[i][j].backptr[0] = i - 1;
               cell[i][j].backptr[1] = j;
               if (min >= path2) {
                  min = path2;
                  cell[i][j].backptr[0] = i - 1;
                  cell[i][j].backptr[1] = j - 1;
               }
               cell[i][j].global = min;
               break;
            case 2:
               cell[i][j].global = path1;
               cell[i][j].backptr[0] = i - 1;
               cell[i][j].backptr[1] = j;
               break;
            case 3:
               cell[i][j].global = path2;
               cell[i][j].backptr[0] = i - 1;
               cell[i][j].backptr[1] = j - 1;
            default:
               break;
            }
         }
      }
   }
   cell[tx - 1][ty - 1].global /=(tx + ty);
}

void calc_global_cost_type_4(Dtw_Table * table)
{
   int i, j, tx = table->data[0].total, ty = table->data[1].total;
   double local = 0.0, min = 0.0, path1, path2, path3;
   Dtw_Cell **cell = table->cell;
   double *weight = table->weight;

   for (i = 1; i < tx; i++) {
      if (cell[i][0].is_region == PATH_OK) {
         cell[i][0].global =
             cell[i - 1][0].global +weight[0] * cell[i][0].local;
         cell[i][0].backptr[0] = i - 1;
         cell[i][0].backptr[1] = 0;
      }
   }

   cell[1][1].global = cell[0][0].global +weight[1] * cell[1][1].local;
   cell[1][1].backptr[0] = 0;
   cell[1][1].backptr[1] = 0;

   for (i = 1; i < tx; i++) {
      if (cell[i][1].is_region == PATH_OK) {
         min = cell[i - 1][1].global +weight[0] * local;
         cell[i][1].backptr[0] = i - 1;
         cell[i][1].backptr[1] = 1;
         if (min >= cell[i - 1][0].global +weight[1] * local) {
            min = cell[i - 1][0].global +weight[1] * local;
            cell[i][1].backptr[0] = i - 1;
            cell[i][1].backptr[1] = 0;
         }
         cell[i][1].global = min;
      }
   }
   cell[1][2].global = cell[0][0].global +weight[2] * cell[1][2].local;
   cell[1][2].backptr[0] = 0;
   cell[1][2].backptr[1] = 0;

   for (j = 2; j < ty; j++) {
      for (i = 2; i < tx; i++) {
         local = cell[i][j].local;
         if (cell[i][j].is_region == PATH_OK) {
            path1 = cell[i - 1][j].global +weight[0] * local;
            path2 = cell[i - 1][j - 1].global +weight[1] * local;
            path3 = cell[i - 1][j - 2].global +weight[2] * local;

            switch (cell[i][j].allow_path) {
            case 1:
               min = path1;
               cell[i][j].backptr[0] = i - 1;
               cell[i][j].backptr[1] = j;
               if (min >= path2) {
                  min = path2;
                  cell[i][j].backptr[0] = i - 1;
                  cell[i][j].backptr[1] = j - 1;
               }
               if (min >= path3) {
                  min = path3;
                  cell[i][j].backptr[0] = i - 1;
                  cell[i][j].backptr[1] = j - 2;
               }
               cell[i][j].global = min;
               break;
            case 2:
               min = path1;
               cell[i][j].backptr[0] = i - 1;
               cell[i][j].backptr[1] = j;
               if (min >= path2) {
                  min = path2;
                  cell[i][j].backptr[0] = i - 1;
                  cell[i][j].backptr[1] = j - 1;
               }
               cell[i][j].global = min;
               break;
            case 3:
               min = path2;
               cell[i][j].backptr[0] = i - 1;
               cell[i][j].backptr[1] = j - 1;
               if (min >= path3) {
                  min = path3;
                  cell[i][j].backptr[0] = i - 1;
                  cell[i][j].backptr[1] = j - 2;
               }
               cell[i][j].global = min;
               break;
            case 4:
               cell[i][j].global = path2;
               cell[i][j].backptr[0] = i - 1;
               cell[i][j].backptr[1] = j - 1;
               break;
            case 5:
               cell[i][j].global = path3;
               cell[i][j].backptr[0] = i - 1;
               cell[i][j].backptr[1] = j - 2;
               break;
            default:
               break;
            }
         }
      }
   }
   cell[tx - 1][ty - 1].global /=(tx + ty);
}

void calc_global_cost_type_5(Dtw_Table * table)
{
   int i, j, tx = table->data[0].total, ty = table->data[1].total;
   double min = 0.0, path1, path2, path3;
   Dtw_Cell **cell = table->cell;
   double *weight = table->weight;

   cell[2][1].global = cell[0][0].global +weight[0] * cell[1][1].local +
       weight[1] * cell[2][1].local;
   cell[2][1].backptr[0] = 0;
   cell[2][1].backptr[1] = 0;

   cell[1][1].global = cell[0][0].global +weight[2] * cell[1][1].local;
   cell[1][1].backptr[0] = 0;
   cell[1][1].backptr[1] = 0;

   cell[1][2].global = cell[0][0].global +weight[3] * cell[1][1].local +
       weight[4] * cell[1][2].local;
   cell[1][2].backptr[0] = 0;
   cell[1][2].backptr[1] = 0;

   for (j = 2; j < ty - 1; j++) {
      for (i = 2; i < tx - 1; i++) {
         if (cell[i][j].is_region == PATH_OK) {
            path1 = cell[i - 2][j - 1].global +
                weight[0] * cell[i - 1][j].local + weight[1] * cell[i][j].local;
            path2 = cell[i - 1][j - 1].global +weight[2] * cell[i][j].local;
            path3 = cell[i - 1][j - 2].global +
                weight[3] * cell[i][j - 1].local + weight[4] * cell[i][j].local;

            switch (cell[i][j].allow_path) {
            case 1:
               min = path1;
               cell[i][j].backptr[0] = i - 2;
               cell[i][j].backptr[1] = j - 1;
               if (min > path2) {
                  min = path2;
                  cell[i][j].backptr[0] = i - 1;
                  cell[i][j].backptr[1] = j - 1;
               }
               if (min > path3) {
                  min = path3;
                  cell[i][j].backptr[0] = i - 1;
                  cell[i][j].backptr[1] = j - 2;
               }
               cell[i][j].global = min;
               break;
            case 2:
               min = path1;
               cell[i][j].backptr[0] = i - 2;
               cell[i][j].backptr[1] = j - 1;
               if (min > path2) {
                  min = path2;
                  cell[i][j].backptr[0] = i - 1;
                  cell[i][j].backptr[1] = j - 1;
               }
               cell[i][j].global = min;
               break;
            case 3:
               min = path2;
               cell[i][j].backptr[0] = i - 1;
               cell[i][j].backptr[1] = j - 1;
               if (min > path3) {
                  min = path3;
                  cell[i][j].backptr[0] = i - 1;
                  cell[i][j].backptr[1] = j - 2;
               }
               cell[i][j].global = min;
               break;
            case 4:
               cell[i][j].global = path1;
               cell[i][j].backptr[0] = i - 2;
               cell[i][j].backptr[1] = j - 1;
               break;
            case 5:
               cell[i][j].global = path2;
               cell[i][j].backptr[0] = i - 1;
               cell[i][j].backptr[1] = j - 1;
               break;
            case 6:
               cell[i][j].global = path3;
               cell[i][j].backptr[0] = i - 1;
               cell[i][j].backptr[1] = j - 2;
               break;
            default:
               break;
            }
         }
      }
   }
   cell[tx - 1][ty - 1].backptr[0] = tx - 2;
   cell[tx - 1][ty - 1].backptr[1] = ty - 2;
   cell[tx - 1][ty - 1].global = cell[tx - 2][ty - 2].global /(tx + ty);
}

void calc_global_cost_type_6(Dtw_Table * table)
{
   int i, j, tx = table->data[0].total, ty = table->data[1].total;
   double local = 0.0, min = 0.0, path1, path2, path3;
   Dtw_Cell **cell = table->cell;
   double *weight = table->weight;

   cell[2][1].global = cell[0][0].global +weight[0] * cell[2][1].local;
   cell[2][1].backptr[0] = 0;
   cell[2][1].backptr[1] = 0;

   cell[1][1].global = cell[0][0].global +weight[1] * cell[1][1].local;
   cell[1][1].backptr[0] = 0;
   cell[1][1].backptr[1] = 0;

   cell[1][2].global = cell[0][0].global +weight[2] * cell[1][2].local;
   cell[1][2].backptr[0] = 0;
   cell[1][2].backptr[1] = 0;

   for (j = 2; j < ty; j++) {
      for (i = 2; i < tx; i++) {
         local = cell[i][j].local;
         if (cell[i][j].is_region == PATH_OK) {
            path1 = cell[i - 2][j - 1].global +weight[0] * local;
            path2 = cell[i - 1][j - 1].global +weight[1] * local;
            path3 = cell[i - 1][j - 2].global +weight[2] * local;

            switch (cell[i][j].allow_path) {
            case 1:
               min = path1;
               cell[i][j].backptr[0] = i - 2;
               cell[i][j].backptr[1] = j - 1;
               if (min >= path2) {
                  min = path2;
                  cell[i][j].backptr[0] = i - 1;
                  cell[i][j].backptr[1] = j - 1;
               }
               if (min >= path3) {
                  min = path3;
                  cell[i][j].backptr[0] = i - 1;
                  cell[i][j].backptr[1] = j - 2;
               }
               cell[i][j].global = min;
               break;
            case 2:
               min = path1;
               cell[i][j].backptr[0] = i - 2;
               cell[i][j].backptr[1] = j - 1;
               if (min >= path2) {
                  min = path2;
                  cell[i][j].backptr[0] = i - 1;
                  cell[i][j].backptr[1] = j - 1;
               }
               cell[i][j].global = min;
               break;
            case 3:
               min = path2;
               cell[i][j].backptr[0] = i - 1;
               cell[i][j].backptr[1] = j - 1;
               if (min >= path3) {
                  min = path3;
                  cell[i][j].backptr[0] = i - 1;
                  cell[i][j].backptr[1] = j - 2;
               }
               cell[i][j].global = min;
               break;
            case 4:
               cell[i][j].global = path1;
               cell[i][j].backptr[0] = i - 2;
               cell[i][j].backptr[1] = j - 1;
               break;
            case 5:
               cell[i][j].global = path2;
               cell[i][j].backptr[0] = i - 1;
               cell[i][j].backptr[1] = j - 1;
               break;
            case 6:
               cell[i][j].global = path3;
               cell[i][j].backptr[0] = i - 1;
               cell[i][j].backptr[1] = j - 2;
               break;
            default:
               break;
            }
         }
      }
   }
   cell[tx - 1][ty - 1].global /=(tx + ty);
}

void calc_global_cost_type_7(Dtw_Table * table)
{
   int i, j, tx = table->data[0].total, ty = table->data[1].total;
   double local = 0.0, min = 0.0, path1, path2, path3, path4;
   Dtw_Cell **cell = table->cell;
   double *weight = table->weight;

   cell[1][1].global = cell[0][0].global +weight[4] * cell[1][1].local;
   cell[1][1].backptr[0] = 0;
   cell[1][1].backptr[1] = 0;

   cell[1][2].global = cell[0][0].global +weight[5] * cell[1][2].local;
   cell[1][2].backptr[0] = 0;
   cell[1][2].backptr[1] = 0;

   cell[2][1].global = cell[0][0].global +weight[0] * cell[1][1].local +
       weight[1] * cell[2][1].local;
   cell[2][1].backptr[0] = 0;
   cell[2][1].backptr[1] = 0;

   for (j = 2; j < ty - 1; j++) {
      for (i = 2; i < tx - 1; i++) {
         local = cell[i][j].local;
         if (cell[i][j].is_region == PATH_OK) {
            path1 = cell[i - 2][j - 1].global +
                weight[0] * cell[i - 1][j].local + weight[1] * local;
            path2 = cell[i - 2][j - 2].global +
                weight[2] * cell[i - 1][j].local + weight[3] * local;
            path3 = cell[i - 1][j - 1].global +weight[4] * local;
            path4 = cell[i - 1][j - 2].global +weight[5] * local;

            switch (cell[i][j].allow_path) {
            case 1:
               min = path1;
               cell[i][j].backptr[0] = i - 2;
               cell[i][j].backptr[1] = j - 1;
               if (min >= path2) {
                  min = path2;
                  cell[i][j].backptr[0] = i - 2;
                  cell[i][j].backptr[1] = j - 2;
               }
               if (min >= path3) {
                  min = path3;
                  cell[i][j].backptr[0] = i - 1;
                  cell[i][j].backptr[1] = j - 1;
               }
               if (min >= path4) {
                  min = path4;
                  cell[i][j].backptr[0] = i - 1;
                  cell[i][j].backptr[1] = j - 2;
               }
               cell[i][j].global = min;
               break;
            case 2:
               min = path2;
               cell[i][j].backptr[0] = i - 2;
               cell[i][j].backptr[1] = j - 2;
               if (min >= path3) {
                  min = path3;
                  cell[i][j].backptr[0] = i - 1;
                  cell[i][j].backptr[1] = j - 1;
               }
               if (min >= path4) {
                  min = path4;
                  cell[i][j].backptr[0] = i - 1;
                  cell[i][j].backptr[1] = j - 2;
               }
               cell[i][j].global = min;
               break;
            case 3:
               min = path1;
               cell[i][j].backptr[0] = i - 2;
               cell[i][j].backptr[1] = j - 1;
               if (min >= path2) {
                  min = path2;
                  cell[i][j].backptr[0] = i - 2;
                  cell[i][j].backptr[1] = j - 2;
               }
               if (min >= path3) {
                  min = path3;
                  cell[i][j].backptr[0] = i - 1;
                  cell[i][j].backptr[1] = j - 1;
               }
               cell[i][j].global = min;
               break;
            case 4:
               min = path3;
               cell[i][j].backptr[0] = i - 1;
               cell[i][j].backptr[1] = j - 1;
               if (min >= path4) {
                  min = path4;
                  cell[i][j].backptr[0] = i - 1;
                  cell[i][j].backptr[1] = j - 2;
               }
               cell[i][j].global = min;
               break;
            case 5:
               min = path2;
               cell[i][j].backptr[0] = i - 2;
               cell[i][j].backptr[1] = j - 2;
               if (min >= path4) {
                  min = path4;
                  cell[i][j].backptr[0] = i - 1;
                  cell[i][j].backptr[1] = j - 2;
               }
               cell[i][j].global = min;
               break;
            case 6:
               min = path1;
               cell[i][j].backptr[0] = i - 2;
               cell[i][j].backptr[1] = j - 1;
               if (min >= path3) {
                  min = path3;
                  cell[i][j].backptr[0] = i - 1;
                  cell[i][j].backptr[1] = j - 1;
               }
               cell[i][j].global = min;
               break;
            case 7:
               cell[i][j].global = path1;
               cell[i][j].backptr[0] = i - 2;
               cell[i][j].backptr[1] = j - 1;
               break;
            case 8:
               cell[i][j].global = path3;
               cell[i][j].backptr[0] = i - 1;
               cell[i][j].backptr[1] = j - 1;
               break;
            case 9:
               cell[i][j].global = path4;
               cell[i][j].backptr[0] = i - 1;
               cell[i][j].backptr[1] = j - 2;
               break;
            default:
               break;
            }
         }
      }
   }
   cell[tx - 1][ty - 1].backptr[0] = tx - 2;
   cell[tx - 1][ty - 1].backptr[1] = ty - 2;
   cell[tx - 1][ty - 1].global = cell[tx - 2][ty - 2].global /tx;
}

/* Check and mark region where global cost can be calculated */
void check_enabled_region(Dtw_Table * table)
{
   switch (table->path) {
   case 1:                     /* horizontal and vertical */
      check_enabled_region_type_1(table);
      break;
   case 2:                     /* horizontal, oblique and vertical */
      check_enabled_region_type_2(table);
      break;
   case 3:                     /* horizontal and oblique */
      check_enabled_region_type_3(table);
      break;
   case 4:                     /* horizontal, oblique1, oblique2 */
      check_enabled_region_type_4(table);
      break;
   case 5:                     /* default */
      check_enabled_region_type_5(table);
      break;
   case 6:
      check_enabled_region_type_6(table);
      break;
   case 7:
      check_enabled_region_type_7(table);
      break;
   default:
      break;
   }
}

/* Calculate local cost */
void calc_local_cost(Dtw_Table * table)
{
   int i, j, d, tdd = table->data[0].dim;
   double sum;
   Dtw_Cell **cell = table->cell;
   Data *data = table->data;
   int norm = table->norm;

   switch (norm) {
   case 1:
      for (i = 0; i < data[0].total; i++) {
         for (j = 0; j < data[1].total; j++) {
            if (cell[i][j].is_region == PATH_OK) {
               for (d = 0, sum = 0.0; d < tdd; d++) {
                  sum += fabs(data[0].input[i * tdd + d] -
                              data[1].input[j * tdd + d]);
               }
               cell[i][j].local = sum;
            }
         }
      }
      break;
   case 2:
      for (i = 0; i < data[0].total; i++) {
         for (j = 0; j < data[1].total; j++) {
            if (cell[i][j].is_region == PATH_OK) {
               for (d = 0, sum = 0.0; d < tdd; d++) {
                  sum += pow((data[0].input[i * tdd + d] -
                              data[1].input[j * tdd + d]), 2);
               }
               cell[i][j].local = sqrt(sum);
            }
         }
      }
      break;
   default:
      break;
   }
}

/* Calculate global cost recursively */
void calc_global_cost(Dtw_Table * table)
{
   table->cell[0][0].global = table->cell[0][0].local;
   table->cell[0][0].backptr[0] = -1;
   table->cell[0][0].backptr[1] = -1;

   switch (table->path) {
   case 1:
      calc_global_cost_type_1(table);
      break;
   case 2:
      calc_global_cost_type_2(table);
      break;
   case 3:
      calc_global_cost_type_3(table);
      break;
   case 4:
      calc_global_cost_type_4(table);
      break;
   case 5:
      calc_global_cost_type_5(table);
      break;
   case 6:
      calc_global_cost_type_6(table);
      break;
   case 7:
      calc_global_cost_type_7(table);
      break;
   default:
      break;
   }
}

/* Obtain Viterbi path */
void get_viterbi_path(Dtw_Table * table)
{
   int k, l, tx = table->data[0].total, ty = table->data[1].total,
       *back_x, *back_y, *phi_x, *phi_y;
   Dtw_Cell **cell = table->cell;
   Data *data = table->data;
   int path = table->path;

   back_x = (int *) malloc(sizeof(int) * (tx + ty));
   back_y = (int *) malloc(sizeof(int) * (tx + ty));

   phi_x = (int *) malloc(sizeof(int) * (tx + ty));
   phi_y = (int *) malloc(sizeof(int) * (tx + ty));

   back_x[0] = phi_x[0] = tx - 1;
   back_y[0] = phi_y[0] = ty - 1;
   k = l = 1;

   while (back_x[l - 1] != 0 && back_y[l - 1] != 0) {
      back_x[l] = cell[back_x[l - 1]][back_y[l - 1]].backptr[0];
      back_y[l] = cell[back_x[l - 1]][back_y[l - 1]].backptr[1];
      switch (path) {
      case 5:
         if (back_x[l - 1] - back_x[l] == 2 && back_y[l - 1] - back_y[l] == 1) {
            phi_x[k] = back_x[l - 1] - 1;
            phi_y[k] = back_y[l - 1];
            phi_x[k + 1] = back_x[l];
            phi_y[k + 1] = back_y[l];
            k += 2;
         } else if (back_x[l - 1] - back_x[l] == 1
                    && back_y[l - 1] - back_y[l] == 1) {
            phi_x[k] = back_x[l];
            phi_y[k] = back_y[l];
            k++;
         } else if (back_x[l - 1] - back_x[l] == 1
                    && back_y[l - 1] - back_y[l] == 2) {
            phi_x[k] = back_x[l - 1];
            phi_y[k] = back_y[l - 1] - 1;
            phi_x[k + 1] = back_x[l];
            phi_y[k + 1] = back_y[l];
            k += 2;
         }
         break;
      case 7:
         if (back_x[l - 1] - back_x[l] == 2 && back_y[l - 1] - back_y[l] == 1) {
            phi_x[k] = back_x[l - 1] - 1;
            phi_y[k] = back_y[l - 1];
            phi_x[k + 1] = back_x[l];
            phi_y[k + 1] = back_y[l];
            k += 2;
         } else if (back_x[l - 1] - back_x[l] == 2
                    && back_y[l - 1] - back_y[l] == 2) {
            phi_x[k] = back_x[l - 1] - 1;
            phi_y[k] = back_y[l - 1];
            phi_x[k + 1] = back_x[l];
            phi_y[k + 1] = back_y[l];
            k += 2;
         } else if (back_x[l - 1] - back_x[l] == 1
                    && back_y[l - 1] - back_y[l] == 1) {
            phi_x[k] = back_x[l];
            phi_y[k] = back_y[l];
            k++;
         } else if (back_x[l - 1] - back_x[l] == 1
                    && back_y[l - 1] - back_y[l] == 2) {
            phi_x[k] = back_x[l];
            phi_y[k] = back_y[l];
            k++;
         }
         break;
      default:
         phi_x[k] = back_x[l];
         phi_y[k] = back_y[l];
         k++;
         break;
      }
      l++;
   }

   table->vit_leng = k;
   for (k = 0; k < table->vit_leng; k++) {
      data[0].viterbi[k] = phi_x[table->vit_leng - k - 1];
      data[1].viterbi[k] = phi_y[table->vit_leng - k - 1];
   }

   free(back_x);
   free(back_y);
   free(phi_x);
   free(phi_y);
}

/* Concatenate two input vectors along Viterbi path */
double *concatenate(Dtw_Table * table)
{
   double *concat;
   int i, j, size = table->vit_leng,
       dim = table->data[0].dim + table->data[1].dim;
   Data *data = table->data;

   concat = dgetmem(size * dim);

   for (i = 0; i < size; i++) {
      for (j = 0; j < data[0].dim; j++) {
         concat[dim * i + j]
             = data[0].input[data[0].viterbi[i] * data[0].dim + j];
      }
      for (j = 0; j < data[1].dim; j++) {
         concat[dim * i + data[0].dim + j]
             = data[1].input[data[1].viterbi[i] * data[1].dim + j];
      }
   }

   return (concat);
}

/* Perform dynamic time warping */
void dtw_s(Dtw_Table * table, double **output)
{
   /* Check and mark region where global cost can be calculated */
   check_enabled_region(table);

   /* Calculate local cost */
   calc_local_cost(table);

   /* Calculate global cost recursively */
   calc_global_cost(table);

   /* Obtain Viterbi path */
   get_viterbi_path(table);

   /* Concatenate two input vectors along Viterbi path */
   *output = concatenate(table);
}

void dtw_free(Dtw_Table *table)
{
	int i=0;
	if(NULL != table->cell[0])
	{
		free(table->cell[0]);
		table->cell[0]=NULL;
	}
	if(NULL != table->cell)
	{
		free(table->cell);
		table->cell =NULL;
	}
	for(i=0;i<2;i++)
	{
		if(NULL != table->data[i].viterbi)
		{
			free(table->data[i].viterbi);
			table->data[i].viterbi=NULL;
		}
	}
	if(NULL != table->weight)
	{
		free(table->weight);
		table->weight=NULL;
	}
}

int dtw(double *x, double *y, int num_test, int num_ref, int dim, int path_type, int norm_type, int *out_addr,int *outviterbi_addr,long *length)
{
      /* Initialize */
	Dtw_Table table;
	int i=0,ret,*outviterbi=NULL;
	double *out=NULL,*z=NULL;
    ret=init_dtw(&table, dim, x, y, num_test, num_ref, path_type, norm_type);
	if(1==ret)
	{
		return 1;
	}
      /* Perform dynamic time warping */
    dtw_s(&table, &out);
	*length=table.vit_leng;

	z=(double *)malloc(2*sizeof(double)*table.vit_leng*dim);
	memcpy(z,out,2*sizeof(double)*table.vit_leng*dim);
	*out_addr=z;

	outviterbi=(int *)malloc(2*sizeof(int)*table.vit_leng);
    for (i = 0; i < table.vit_leng; i++) {
		outviterbi[2*i]=*(table.data[0].viterbi + i);
		outviterbi[2*i+1]=*(table.data[1].viterbi + i);
    }
	*outviterbi_addr=outviterbi;

	dtw_free(&table);
	return 0;
}

int dtw_result(double *out,int *viterbi_path,int out_addr,int outviterbi_addr,long length,int dim)
{
	int i=0,*ptr_viterbi=NULL;
	double *ptr_data=NULL;

	ptr_viterbi=outviterbi_addr;
	ptr_data=out_addr;
	for(i=0;i< (dim*length*2);i++)
	{
		out[i]=ptr_data[i];
	}
	for(i=0;i<length*2;i++)
	{
		viterbi_path[i]=ptr_viterbi[i];
	}

	free(ptr_data);
	free(ptr_viterbi);
	return 0;
}










