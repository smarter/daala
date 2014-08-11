/*Daala video codec
Copyright (c) 2002-2014 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS”
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vidinput.h"
#include "qr.h"
#include "../src/state.h"
#if defined(_WIN32)
# include <io.h>
# include <fcntl.h>
#endif
#include "getopt.h"
#include <math.h>

#define ORIGINAL_FILTERS 0

#define MIN(a,b) ((a)>(b)?(b):(a))

static void usage(char **_argv) {
  fprintf(stderr, "Usage: %s <input> <output>\n"
   "    <input> must be a YUV4MPEG file.\n\n", _argv[0]);
}

static const char *CHROMA_TAGS[4] = { " C420jpeg", "", " C422jpeg", " C444" };

#define MAX_TAPS 16

enum filters {
  EDGE_FILTER_1_H = 0,
  EDGE_FILTER_1_V,
  EDGE_FILTER_2_H,
  EDGE_FILTER_2_V,
  EDGE_FILTER_3_H,
  EDGE_FILTER_3_V,
  EDGE_FILTER_4_H,
  EDGE_FILTER_4_V,
  EDGE_FILTER_5_H,
  EDGE_FILTER_5_V,
  FALLBACK_FILTER_H,
  FALLBACK_FILTER_V,
  FILTERS_MAX
};

struct filter_regression {
  int taps;
  int *pred_col[MAX_TAPS];
  unsigned char *resp;
  int allocated_length;
  int length;
};

static void regression_init(struct filter_regression *reg, int taps) {
  int i;
  int length = 1024*1024;
  reg->taps = taps;
  for (i = 0; i < taps; i++) {
    reg->pred_col[i] = malloc(length*sizeof(**reg->pred_col));
  }
  reg->resp = malloc(length*sizeof(*reg->resp));
  reg->allocated_length = length;
  reg->length = 0;
}

static void regression_add_resp(struct filter_regression *reg, unsigned char value) {
  int i;
  reg->length++;
  if (reg->length > reg->allocated_length) {
    reg->allocated_length *= 2;
    for (i = 0; i < reg->taps; i++) {
      reg->pred_col[i] = realloc(reg->pred_col[i], reg->allocated_length*sizeof(**reg->pred_col));
    }
    reg->resp = realloc(reg->resp, reg->allocated_length*sizeof(*reg->resp));
  }
  reg->resp[reg->length - 1] = value;
}

static void regression_set_tap(struct filter_regression *reg, int tap, int value) {
  reg->pred_col[tap][reg->length - 1] = value;
}
static void regression_solve(struct filter_regression *reg) {
  int m = reg->length;
  int n = reg->taps;
  double *sol;
  double *mat;
  double d[MAX_TAPS];
  int rank;
  int i, j;
  double tap_sum = 0;
  sol = malloc(m*sizeof(*sol));
  for (i = 0; i < m; i++) {
    sol[i] = reg->resp[i];
  }
  mat = malloc(m*n*sizeof(*mat));
  for (i = 0; i < n; i++) {
    for (j = 0; j < m; j++) {
      mat[i*m + j] = reg->pred_col[i][j];
    }
  }
  rank = qrdecomp_hh(mat, m, d, NULL, 0, NULL, n, m);
  if (rank != n) {
    fprintf(stderr, "rank (%d) != n (%d)\n", rank, n);
    exit(1);
  }
  qrsolve_hh(mat, m, d, sol, m, n, m, 1);

  for (i = 0; i < reg->taps; i++) {
    printf("Tap %d: %f\n", i, sol[i]);
    tap_sum += sol[i];
  }
  printf("Tap sum: %f\n", tap_sum);
}

static int sinc_filter(const unsigned char *s, int stride) {
  return 20*s[0] - 5*s[1*stride] + s[2*stride];
}

static int reconstruct_v(const unsigned char *s, int xstride, int ystride,
 int a, int b, int c, int d)
{
  int x;

  x = sinc_filter(s + 0 - 3 * ystride, -xstride) * a;
  x += sinc_filter(s + 0 - 2 * ystride, -xstride) * b;
  x += sinc_filter(s + 0 - 1 * ystride, -xstride) * c;
  x += sinc_filter(s + 0 - 0 * ystride, -xstride) * d;
  x += sinc_filter(s + 1*xstride + 0 * ystride, xstride) * d;
  x += sinc_filter(s + 1*xstride + 1 * ystride, xstride) * c;
  x += sinc_filter(s + 1*xstride + 2 * ystride, xstride) * b;
  x += sinc_filter(s + 1*xstride + 3 * ystride, xstride) * a;
  return (x + 16*16) >> (5+4);
}

static int reconstruct_h(const unsigned char *s, int xstride, int ystride,
 int a, int b, int c, int d)
{
  int x;

  x = sinc_filter(s - 3*xstride, -ystride) * a;
  x += sinc_filter(s - 2*xstride, -ystride) * b;
  x += sinc_filter(s - 1*xstride, -ystride) * c;
  x += sinc_filter(s - 0*xstride, -ystride) * d;
  x += sinc_filter(s + ystride + 0*xstride, ystride) * d;
  x += sinc_filter(s + ystride + 1*xstride, ystride) * c;
  x += sinc_filter(s + ystride + 2*xstride, ystride) * b;
  x += sinc_filter(s + ystride + 3*xstride, ystride) * a;
  return (x + 16*16) >> (5+4);
}

static void train_upsample(od_state *state, od_img *dimg, const od_img *simg,
 struct filter_regression *regs, unsigned char *up_ref, int up_ref_stride) {
  int pli;
  for (pli = 0; pli < 1; pli++) {
    const od_img_plane *siplane;
    od_img_plane *diplane;
    int src_stride, dst_stride;
    const unsigned char *src;
    unsigned char *dst;
    const unsigned char *s;
    unsigned char *d;
    int xpad;
    int ypad;
    int w;
    int h;
    int x;
    int y;
    const int MARGIN = 3;

    siplane = simg->planes + pli;
    diplane = dimg->planes + pli;
    src_stride = siplane->ystride;
    dst_stride = diplane->ystride;
    xpad = OD_UMV_PADDING >> siplane->xdec;
    ypad = OD_UMV_PADDING >> siplane->ydec;
    w = simg->width >> siplane->xdec;
    h = simg->height >> siplane->ydec;
    src = siplane->data;
    dst = diplane->data;
    s = src;
    d = dst - dst_stride*(2*ypad);

    /* Padding, margin and source pixels copy */
    for (y = -ypad; y < h + ypad; y++) {
      memset(d - 2*xpad, s[0], 2*xpad);
      if (y < MARGIN || y >= h - MARGIN - 1) {
        d[0] = s[0];
        d[1] = (20*(s[0] + s[1]) - 5*(s[0] + s[2]) + s[0] + s[3] + 16) >> 5;
        d[2] = s[1];
        d[3] = (20*(s[1] + s[2]) - 5*(s[0] + s[3]) + s[0] + s[4] + 16) >> 5;
        for (x = 2; x < w - 3; x++) {
          d[2*x] = s[x];
          d[2*x + 1] = (20*(s[x] + s[x + 1])
           - 5*(s[x - 1] + s[x + 2]) + s[x - 2] + s[x + 3] + 16) >> 5;
        }
        d[2*x] = s[x];
        d[2*x + 1] = (20*(s[x] + s[x + 1])
         - 5*(s[x - 1] + s[x + 2]) + s[x - 2] + s[x + 2] + 16) >> 5;
        x++;
        d[2*x] = s[x];
        d[2*x + 1] = (20*(s[x] + s[x + 1])
         - 5*(s[x - 1] + s[x + 1]) + s[x - 2] + s[x + 1] + 16) >> 5;
        x++;
        d[2*x] = s[x];
        d[2*x + 1] = (36*s[x] - 5*s[x - 1] + s[x - 2] + 16) >> 5;
        x++;
      } else {
        for (x = 0; x < w; x++) {
          d[2*x] = s[x];
        }
      }
      memset(d + 2*w, s[w - 1], 2*xpad);
      if (y >= 0 && y < h - 1)
        s += src_stride;
      else
        OD_COPY(d + dst_stride - 2*xpad, d - 2*xpad, 2*(w + 2*xpad));
      d += 2*dst_stride;
    }
    /* Horizontal filtering */
    d = dst + dst_stride*2*MARGIN;
    for (y = MARGIN; y < h - MARGIN - 1; y++) {
      for (x = 0; x < 2*w; x += 2) {
        int dx, dy, dx2;
        int v;

        dx = -d[-2*dst_stride + x]
             - d[-2*dst_stride + x + 2]
             + d[2*dst_stride + x]
             + d[2*dst_stride + x + 2];
        dx *= 2;

        dy = -d[-2*dst_stride + x]
             - 2 * d[x]
             - d[2*dst_stride + x]
             + d[-2*dst_stride + x + 2]
             + 2 * d[x + 2]
             + d[2*dst_stride + x + 2];

        dx2 = -d[-2*dst_stride + x]
              + 2 * d[x]
              - d[2*dst_stride + x]
              - d[-2*dst_stride + x + 2]
              + 2 * d[x + 2]
              - d[2*dst_stride + x + 2];

        if (dy < 0) {
          dy = -dy;
          dx = -dx;
        }

        if (abs(dx) <= 4 * abs(dx2)) {
          struct filter_regression *r = regs + FALLBACK_FILTER_H;
#if ORIGINAL_FILTERS
          v = (20*(d[x] + d[x + 2])
               - 5*(d[x - 2] + d[x + 4]) + d[x - 4] + d[x + 6] + 16) >> 5;
#else
          v = (592*(d[x] + d[x + 2]) - 70*(d[x - 2] + d[x + 4])
            - 10*(d[x - 4] + d[x + 6]) + 512) >> 10;
#endif
          regression_add_resp(r , up_ref[2*y*up_ref_stride + (x + 1)]);
          regression_set_tap(r, 0, d[x - 4] + d[x + 6]);
          regression_set_tap(r, 1, d[x - 2] + d[x + 4]);
          regression_set_tap(r, 2, d[x] + d[x + 2]);
        } else if (dx < 0) {
          if (dx < -2 * dy) {
            v = reconstruct_v(d + x, 2, 2*dst_stride, 0, 0, 0, 16);
          } else if (dx < -dy) {
            v = reconstruct_v(d + x, 2, 2*dst_stride, 0, 0, 8, 8);
          } else if (2 * dx < -dy) {
            v = reconstruct_v(d + x, 2, 2*dst_stride, 0, 4, 8, 4);
          } else if (3 * dx < -dy) {
            v = reconstruct_v(d + x, 2, 2*dst_stride, 1, 7, 7, 1);
          } else {
            v = reconstruct_v(d + x, 2, 2*dst_stride, 4, 8, 4, 0);
          }
        } else {
          if (dx > 2 * dy) {
            v = reconstruct_v(d + x, 2, -2*dst_stride, 0, 0, 0, 16);
          } else if (dx > dy) {
            v = reconstruct_v(d + x, 2, -2*dst_stride, 0, 0, 8, 8);
          } else if (2 * dx > dy) {
            v = reconstruct_v(d + x, 2, -2*dst_stride, 0, 4, 8, 4);
          } else if (3 * dx > dy) {
            v = reconstruct_v(d + x, 2, -2*dst_stride, 1, 7, 7, 1);
          } else {
            v = reconstruct_v(d + x, 2, -2*dst_stride, 4, 8, 4, 0);
          }
        }
        d[x + 1] = OD_CLAMP255(v);
      }
      if (y >= 0 && y < h - 1)
        s += src_stride;
      d += 2*dst_stride;
    }
    /* Vertical filtering */
    d = dst;
    for (y = 0; y < h; y++) {
      unsigned char *d0 = d;
      unsigned char *d1 = d + dst_stride;
      unsigned char *d2 = d + 2*dst_stride;
      unsigned char *d4 = d + 4*dst_stride;
      unsigned char *d6 = d + 6*dst_stride;
      unsigned char *dm2 = d - 2*dst_stride;
      unsigned char *dm4 = d - 4*dst_stride;

      for (x = -2*xpad; x < 2*w + 2*xpad; x++) {
        if (x >= MARGIN && x < 2*w - MARGIN - 1) {
          int dx, dy;
          int dx2;
          int v;

          dx = -d0[x - 1]
              - d2[x - 1]
              + d0[x + 1]
              + d2[x + 1];
          dx *= 2;

          dy = -d0[x - 1]
              - 2 * d0[x]
              - d0[x + 1]
              + d2[x - 1]
              + 2 * d2[x]
              + d2[x + 1];

          dx2 = -d0[x - 1]
              + 2 * d0[x]
              - d0[x + 1]
              - d2[x - 1]
              + 2 * d2[x]
              - d2[x + 1];

          if (dy < 0) {
            dy = -dy;
            dx = -dx;
          }

          if (abs(dx) <= 4*abs(dx2)) {
            struct filter_regression *r = regs + FALLBACK_FILTER_V;
#if ORIGINAL_FILTERS
            v = (20*(d0[x] + d2[x])
             - 5*(dm2[x] + d4[x]) + dm4[x] + d6[x] + 16) >> 5;
#else
            v = (403*(d0[x] + d2[x])
                 + 42*(dm2[x] + d4[x]) + 67*(dm4[x] + d6[x]) + 512) >> 10;
#endif
            regression_add_resp(r , up_ref[2*(y + 1)*up_ref_stride + x]);
            regression_set_tap(r, 0, d[x - 4] + d[x + 6]);
            regression_set_tap(r, 1, d[x - 2] + d[x + 4]);
            regression_set_tap(r, 2, d[x] + d[x + 2]);
          } else if (dx < 0) {
            if (dx < -2 * dy) {
              v = reconstruct_h(d0 + x, 1, 2*dst_stride, 0, 0, 0, 16);
            } else if (dx < -dy) {
              v = reconstruct_h(d0 + x, 1, 2*dst_stride, 0, 0, 8, 8);
            } else if (2 * dx < -dy) {
              v = reconstruct_h(d0 + x, 1, 2*dst_stride, 0, 4, 8, 4);
            } else if (3 * dx < -dy) {
              v = reconstruct_h(d0 + x, 1, 2*dst_stride, 1, 7, 7, 1);
            } else {
              v = reconstruct_h(d0 + x, 1, 2*dst_stride, 4, 8, 4, 0);
            }
          } else {
            if (dx > 2 * dy) {
              v = reconstruct_h(d2 + x, 1, -2*dst_stride, 0, 0, 0, 16);
            } else if (dx > dy) {
              v = reconstruct_h(d2 + x, 1, -2*dst_stride, 0, 0, 8, 8);
            } else if (2 * dx > dy) {
              v = reconstruct_h(d2 + x, 1, -2*dst_stride, 0, 4, 8, 4);
            } else if (3 * dx > dy) {
              v = reconstruct_h(d2 + x, 1, -2*dst_stride, 1, 7, 7, 1);
            } else {
              v = reconstruct_h(d2 + x, 1, -2*dst_stride, 4, 8, 4, 0);
            }
          }
          d1[x] = OD_CLAMP255(v);
        } else {
          d1[x] = (20*(d0[x] + d2[x])
           - 5*(dm2[x] + d4[x]) + dm4[x] + d6[x] + 16) >> 5;
        }
      }
      d += 2*dst_stride;
    }
  }
  for (pli = 1; pli < state->io_imgs[OD_FRAME_REC].nplanes; pli++) {
    const od_img_plane *siplane;
    od_img_plane *diplane;
    const unsigned char *src;
    unsigned char *dst;
    int xpad;
    int ypad;
    int w;
    int h;
    int x;
    int y;
    siplane = simg->planes + pli;
    diplane = dimg->planes + pli;
    xpad = OD_UMV_PADDING >> siplane->xdec;
    ypad = OD_UMV_PADDING >> siplane->ydec;
    w = simg->width >> siplane->xdec;
    h = simg->height >> siplane->ydec;
    src = siplane->data;
    dst = diplane->data - (diplane->ystride << 1)*ypad;
    for (y = -ypad; y < h + ypad + 3; y++) {
      /*Horizontal filtering:*/
      if (y < h + ypad) {
        unsigned char *buf;
        buf = state->ref_line_buf[y & 7];
        memset(buf - (xpad << 1), src[0], (xpad - 2) << 1);
        /*for (x = -xpad; x < -2; x++) {
          *(buf + (x << 1)) = src[0];
          *(buf + (x << 1 | 1)) = src[0];
        }*/
        *(buf - 4) = src[0];
        *(buf - 3) = OD_CLAMP255((31*src[0] + src[1] + 16) >> 5);
        *(buf - 2) = src[0];
        *(buf - 1) = OD_CLAMP255((36*src[0] - 5*src[1] + src[1] + 16) >> 5);
        buf[0] = src[0];
        buf[1] = OD_CLAMP255((20*(src[0] + src[1])
         - 5*(src[0] + src[2]) + src[0] + src[3] + 16) >> 5);
        buf[2] = src[1];
        buf[3] = OD_CLAMP255((20*(src[1] + src[2])
         - 5*(src[0] + src[3]) + src[0] + src[4] + 16) >> 5);
        for (x = 2; x < w - 3; x++) {
          buf[x << 1] = src[x];
          buf[x << 1 | 1] = OD_CLAMP255((20*(src[x] + src[x + 1])
           - 5*(src[x - 1] + src[x + 2]) + src[x - 2] + src[x + 3] + 16) >> 5);
        }
        buf[x << 1] = src[x];
        buf[x << 1 | 1] = OD_CLAMP255((20*(src[x] + src[x + 1])
         - 5*(src[x - 1] + src[x + 2]) + src[x - 2] + src[x + 2] + 16) >> 5);
        x++;
        buf[x << 1] = src[x];
        buf[x << 1 | 1] = OD_CLAMP255((20*(src[x] + src[x + 1])
         - 5*(src[x - 1] + src[x + 1]) + src[x - 2] + src[x + 1] + 16) >> 5);
        x++;
        buf[x << 1] = src[x];
        buf[x << 1 | 1] =
         OD_CLAMP255((36*src[x] - 5*src[x - 1] + src[x - 2] + 16) >> 5);
        x++;
        buf[x << 1] = src[w - 1];
        buf[x << 1 | 1] = OD_CLAMP255((31*src[w - 1] + src[w - 2] + 16) >> 5);
        memset(buf + (++x << 1), src[w - 1], (xpad - 1) << 1);
        /*for (x++; x < w + xpad; x++) {
          buf[x << 1] = src[w - 1];
          buf[x << 1 | 1]=src[w - 1];
        }*/
        if (y >= 0 && y + 1 < h) src += siplane->ystride;
      }
      /*Vertical filtering:*/
      if (y >= -ypad + 3) {
        if (y < 1 || y > h + 3) {
          OD_COPY(dst - (xpad << 1),
           state->ref_line_buf[(y - 3) & 7] - (xpad << 1),
           (w + (xpad << 1)) << 1);
          /*fprintf(stderr, "%3i: ", (y - 3) << 1);
          for (x = -xpad << 1; x < (w + xpad) << 1; x++) {
            fprintf(stderr, "%02X", *(dst + x));
          }
          fprintf(stderr, "\n");*/
          dst += diplane->ystride;
          OD_COPY(dst - (xpad << 1),
           state->ref_line_buf[(y - 3) & 7] - (xpad << 1),
           (w + (xpad << 1)) << 1);
          /*fprintf(stderr, "%3i: ", (y - 3) << 1 | 1);
          for (x = -xpad << 1; x < (w + xpad) << 1; x++) {
            fprintf(stderr, "%02X", *(dst + x));
          }
          fprintf(stderr, "\n");*/
          dst += diplane->ystride;
        }
        else {
          unsigned char *buf[6];
          buf[0] = state->ref_line_buf[(y - 5) & 7];
          buf[1] = state->ref_line_buf[(y - 4) & 7];
          buf[2] = state->ref_line_buf[(y - 3) & 7];
          buf[3] = state->ref_line_buf[(y - 2) & 7];
          buf[4] = state->ref_line_buf[(y - 1) & 7];
          buf[5] = state->ref_line_buf[(y - 0) & 7];
          OD_COPY(dst - (xpad << 1),
           state->ref_line_buf[(y - 3) & 7] - (xpad << 1),
           (w + (xpad << 1)) << 1);
          /*fprintf(stderr, "%3i: ", (y - 3) << 1);
          for (x = -xpad << 1; x < (w + xpad) << 1; x++) {
            fprintf(stderr, "%02X", *(dst + x));
          }
          fprintf(stderr, "\n");*/
          dst += diplane->ystride;
          for (x = -xpad << 1; x < (w + xpad) << 1; x++) {
            *(dst + x) = OD_CLAMP255((20*(*(buf[2] + x) + *(buf[3] + x))
             - 5*(*(buf[1] + x) + *(buf[4] + x))
             + *(buf[0] + x) + *(buf[5] + x) + 16) >> 5);
          }
          /*fprintf(stderr, "%3i: ", (y - 3) << 1 | 1);
          for (x = -xpad << 1; x < (w + xpad) << 1; x++) {
            fprintf(stderr, "%02X", *(dst + x));
          }
          fprintf(stderr, "\n");*/
          dst += diplane->ystride;
        }
      }
    }
  }
}

int main(int _argc, char **_argv) {
  const char *optstring = "hv?";
  const struct option long_options[] = {
    { NULL, 0, NULL, 0 }
  };
  FILE *fin;
  FILE *fout;
  video_input vid1;
  video_input_info info1;
  int frameno;
  int pli;
  int xdec[3];
  int ydec[3];
  int w[3];
  int h[3];
  int limit;
  int long_option_index;
  int c;
  od_state state;
  daala_info info;
  struct filter_regression regs[FILTERS_MAX];

  limit = 0;
  while ((c = getopt_long(_argc, _argv, optstring, long_options,
   &long_option_index)) != EOF) {
    switch (c) {
      case 'v':
      case '?':
      case 'h':
      default:
      {
        usage(_argv);
        exit(EXIT_FAILURE);
      }
      break;
    }
  }
  if (optind+2 != _argc) {
    usage(_argv);
    exit(EXIT_FAILURE);
  }

  /*  1  2  3  4  3
      2  4  6  8  6 || *2 (a a b b ...)
      4  8 12 16 12 || *3 (sinc ...)*/
  regression_init(regs + EDGE_FILTER_1_H, 4);
  regression_init(regs + EDGE_FILTER_1_V, 4);
  regression_init(regs + EDGE_FILTER_2_H, 8);
  regression_init(regs + EDGE_FILTER_2_V, 8);
  regression_init(regs + EDGE_FILTER_3_H, 12);
  regression_init(regs + EDGE_FILTER_3_V, 12);
  regression_init(regs + EDGE_FILTER_4_H, 16);
  regression_init(regs + EDGE_FILTER_4_V, 16);
  regression_init(regs + EDGE_FILTER_5_H, 12);
  regression_init(regs + EDGE_FILTER_5_V, 12);
  regression_init(regs + FALLBACK_FILTER_H, 3);
  regression_init(regs + FALLBACK_FILTER_V, 3);

  /* TODO: for argv < argc ... */
  fin = strcmp(_argv[optind], "-") == 0 ? stdin : fopen(_argv[optind], "rb");
  if (fin == NULL) {
    fprintf(stderr, "Unable to open '%s' for extraction.\n", _argv[optind]);
    exit(EXIT_FAILURE);
  }
  fprintf(stderr, "Opening %s as input...\n", _argv[optind]);
  if (video_input_open(&vid1, fin) < 0) exit(EXIT_FAILURE);
  video_input_get_info(&vid1, &info1);

  daala_info_init(&info);
  for (pli = 0; pli < 3; pli++) {
    xdec[pli] = pli && !(info1.pixel_fmt&1);
    ydec[pli] = pli && !(info1.pixel_fmt&2);
    h[pli] = (info1.pic_h/2) >> ydec[pli];
    w[pli] = (info1.pic_w/2) >> xdec[pli];

    info.plane_info[pli].xdec = xdec[pli];
    info.plane_info[pli].ydec = ydec[pli];
  }
  info.nplanes = 3;
  info.pic_height = h[0];
  info.pic_width = w[0];

  od_state_init(&state, &info);

  fout = strcmp(_argv[optind+1], "-") == 0 ? stdout : fopen(_argv[optind+1],
   "wb");
  if (fout == NULL) {
    fprintf(stderr, "Error opening output file \"%s\".\n", _argv[optind+1]);
    return 1;
  }


  fprintf(fout, "YUV4MPEG2 W%i H%i F%i:%i Ip A%i:%i%s\n",
   w[0]*2, h[0]*2, (unsigned)info1.fps_n,
   (unsigned)info1.fps_d, info1.par_n, info1.par_d,
   CHROMA_TAGS[ydec[1] ? xdec[1] ? 0 : 2 : 3]);
  for (frameno = 0;; frameno++) {
    video_input_ycbcr in;
    int             ret1 = 0;
    char            tag1[5];
    od_img *simg = &state.io_imgs[OD_FRAME_REC];
    od_img *dimg = &state.ref_imgs[0];
    int x, y;
    if (!limit || frameno < limit) {
      ret1 = video_input_fetch_frame(&vid1, in, tag1);
    }
    if (ret1 == 0) break;
    for (pli = 0; pli < 3; pli++) {
      od_img_plane *siplane = simg->planes + pli;
      unsigned char *src = siplane->data;
      int src_stride = siplane->ystride;
      int plane_width = simg->width >> xdec[pli];
      int plane_height = simg->height >> ydec[pli];
      for (y = 0; y < h[pli]; y++) {
        for (x = 0; x < w[pli]; x++) {
          int cy = 2*y + (info1.pic_y >> ydec[pli]);
          int cx = 2*x + (info1.pic_x >> xdec[pli]);
          src[y*src_stride + x] = in[pli].data[cy*in[pli].stride + cx];
        }
      }
      /*From od_img_plane_copy_pad8*/
      /*Right side.*/
      for (x = w[pli]; x < plane_width; x++) {
        src = siplane->data + x - 1;
        for (y = 0; y < h[pli]; y++) {
          src[1] = (2*src[0] + (src - (src_stride & -(y > 0)))[0]
           + (src + (src_stride & -(y + 1 < h[pli])))[0] + 2) >> 2;
          src += src_stride;
        }
      }
      /*Bottom.*/
      src = siplane->data + src_stride*h[pli];
      for (y = h[pli]; y < plane_height; y++) {
        for (x = 0; x < plane_width; x++) {
          src[x] = (2*(src - src_stride)[x] + (src - src_stride)[x - (x > 0)]
           + (src - src_stride)[x + (x + 1 < plane_width)] + 2) >> 2;
        }
        src += src_stride;
      }
    }
    train_upsample(&state, dimg, simg, regs, in[0].data + info1.pic_y*in[0].stride + info1.pic_x,
                   in[0].stride);
    fprintf(fout, "FRAME\n");
    for (pli = 0; pli < 3; pli++) {
      od_img_plane *diplane = dimg->planes + pli;
      unsigned char *dst = diplane->data;
      for (y = 0; y < 2*h[pli]; y++) {
        if (fwrite(dst + diplane->ystride*y, 2*w[pli], 1, fout) < 1) {
          fprintf(stderr, "Error writing to output.\n");
          return EXIT_FAILURE;
        }
      }
    }
  }
  regression_solve(regs + FALLBACK_FILTER_H);
  /*regression_solve(regs + FALLBACK_FILTER_V);*/
  video_input_close(&vid1);
  /*regression_solve(&reg);*/
  if (fout != stdout) fclose(fout);
  return EXIT_SUCCESS;
}
