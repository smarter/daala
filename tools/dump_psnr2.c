#include "vidinput.h"
#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#if defined(_WIN32)
#include <io.h>
#include <fcntl.h>
#endif
#include "getopt.h"

#define OD_CLAMP255(x) \
  ((unsigned char)((((x) < 0) - 1) & ((x) | -((x) > 255))))

const char *optstring = "sa";
struct option options [] = {
  { "summary", no_argument, NULL, 's' },
  { "avg-ds", no_argument, NULL, 'a' },
  { NULL, 0, NULL, 0 }
};

static int summary_only;
static int avg_ds;

enum {
  INTERP_H = 0,
  INTERP_V,
  INTERP_HV,
  INTERP_MAX
};

enum {
  SRC_PIXEL = 0,
  FILTER_1,
  FILTER_2,
  FILTER_3,
  FILTER_4,
  FILTER_5,
  BILINEAR,
  NEAREST_NEIGHBOR,
  FILTERS_MAX
};

static void usage(char *_argv[]) {
  fprintf(stderr, "Usage: %s [options] <video1>\n"
   "    <video1> may be either YUV4MPEG or Ogg Theora.\n\n"
   "    Options:\n\n"
   "      -s --summary   Only output the summary line.\n"
   "      -a --avg-ds    Downsample using 2x2 average instead of decimating.\n",
   _argv[0]);
}

static int reconstruct_v(const unsigned char *src, int xstride, int ystride,
 int a, int b, int c, int d)
{
  int x;

  x = src[0 - 3 * ystride] * a;
  x += src[0 - 2 * ystride] * b;
  x += src[0 - 1 * ystride] * c;
  x += src[0 - 0 * ystride] * d;
  x += src[1*xstride + 0*ystride] * d;
  x += src[1*xstride + 1*ystride] * c;
  x += src[1*xstride + 2*ystride] * b;
  x += src[1*xstride + 3*ystride] * a;
  return (x + 16) >> 5;
}

static int reconstruct_h(const unsigned char *d1, const unsigned char *d2,
 int a, int b, int c, int d)
{
  int x;

  x = d1[-3] * a;
  x += d1[-2] * b;
  x += d1[-1] * c;
  x += d1[-0] * d;
  x += d2[0] * d;
  x += d2[1] * c;
  x += d2[2] * b;
  x += d2[3] * a;
  return (x + 16) >> 5;
}

void edi_hv(uint8_t *dst, int dst_stride, uint8_t *filter_buf, int filter_buf_stride,
            const uint8_t *src, int src_stride, int w, int h) {
  const uint8_t *s;
  uint8_t *d;
  uint8_t *f;
  int xpad, ypad;
  int x, y;
  const int MARGIN = 3;

  xpad = 0;
  ypad = 0;
  s = src;
  d = dst;
  f = filter_buf;

  /* Horizontal filtering */
  for (y = -ypad; y < h + ypad; y++) {
    if (y >= MARGIN && y < h - MARGIN - 1) {
      for (x = 0; x < w - 1; x++) {
        int dx, dy, dx2;
        int v;

        dx = -s[-src_stride + x]
             - s[-src_stride + x + 1]
             + s[src_stride + x]
             + s[src_stride + x + 1];
        dx *= 2;

        dy = -s[-src_stride + x]
             - 2 * s[x]
             - s[src_stride + x]
             + s[-src_stride + x + 1]
             + 2 * s[x + 1]
             + s[src_stride + x + 1];

        dx2 = -s[-src_stride + x]
              + 2 * s[x]
              - s[src_stride + x]
              - s[-src_stride + x + 1]
              + 2 * s[x + 1]
              - s[src_stride + x + 1];

        if (dy < 0) {
          dy = -dy;
          dx = -dx;
        }

        if (abs(dx) <= 4 * abs(dx2)) {
          v = (s[x] + s[x + 1] + 1) >> 1;
          f[x * 2 + 1] = BILINEAR;
        } else if (dx < 0) {
          if (dx < -2 * dy) {
            v = reconstruct_v(s + x, 1, src_stride, 0, 0, 0, 16);
            f[x * 2 + 1] = FILTER_1;
          } else if (dx < -dy) {
            v = reconstruct_v(s + x, 1, src_stride, 0, 0, 8, 8);
            f[x * 2 + 1] = FILTER_2;
          } else if (2 * dx < -dy) {
            v = reconstruct_v(s + x, 1, src_stride, 0, 4, 8, 4);
            f[x * 2 + 1] = FILTER_3;
          } else if (3 * dx < -dy) {
            v = reconstruct_v(s + x, 1, src_stride, 1, 7, 7, 1);
            f[x * 2 + 1] = FILTER_4;
          } else {
            v = reconstruct_v(s + x, 1, src_stride, 4, 8, 4, 0);
            f[x * 2 + 1] = FILTER_5;
          }
        } else {
          if (dx > 2 * dy) {
            v = reconstruct_v(s + x, 1, -src_stride, 0, 0, 0, 16);
            f[x * 2 + 1] = FILTER_1;
          } else if (dx > dy) {
            v = reconstruct_v(s + x, 1, -src_stride, 0, 0, 8, 8);
            f[x * 2 + 1] = FILTER_2;
          } else if (2 * dx > dy) {
            v = reconstruct_v(s + x, 1, -src_stride, 0, 4, 8, 4);
            f[x * 2 + 1] = FILTER_3;
          } else if (3 * dx > dy) {
            v = reconstruct_v(s + x, 1, -src_stride, 1, 7, 7, 1);
            f[x * 2 + 1] = FILTER_4;
          } else {
            v = reconstruct_v(s + x, 1, -src_stride, 4, 8, 4, 0);
            f[x * 2 + 1] = FILTER_5;
          }
        }
        d[x * 2] = s[x];
        d[x * 2 + 1] = OD_CLAMP255(v);
      }
      d[x * 2] = s[x];
      d[x * 2 + 1] = s[x];
      f[x * 2 + 1] = NEAREST_NEIGHBOR;
    } else {
      for (x = 0; x < w - 1; x++) {
        d[x * 2] = s[x];
        d[x * 2 + 1] = (s[x] + s[x + 1] + 1) >> 1;
        f[x * 2 + 1] = BILINEAR;
      }
      d[x * 2] = s[x];
      d[x * 2 + 1] = s[x];
      f[x * 2 + 1] = NEAREST_NEIGHBOR;
    }
    for (x = -xpad; x < 0; x++) {
      d[x * 2] = s[0];
      d[x * 2 + 1] = s[0];
      f[x * 2] = NEAREST_NEIGHBOR;
      f[x * 2 + 1] = NEAREST_NEIGHBOR;
    }
    for (x = w; x < w + xpad; x++) {
      d[x * 2] = s[w - 1];
      d[x * 2 + 1] = s[w - 1];
      f[x * 2] = NEAREST_NEIGHBOR;
      f[x * 2 + 1] = NEAREST_NEIGHBOR;
    }

    if (y >= 0 && y < h - 1)
      s += src_stride;
    d += 2*dst_stride;
    f += 2*dst_stride;
  }
  /* Vertical filtering */
  d = dst;
  f = filter_buf;
  for (y = -ypad; y < h + ypad - 1; y++) {
    unsigned char *d1 = d;
    unsigned char *d2 = d + dst_stride;
    unsigned char *d3 = d + 2*dst_stride;

    uint8_t *f2 = f + filter_buf_stride;

    for (x = -2*xpad; x < w*2  + xpad*2; x++) {
      if (x >= MARGIN && x < w * 2 - MARGIN - 1) {
        int dx, dy;
        int dx2;
        int v;

        dx = -d1[x - 1]
             - d3[x - 1]
             + d1[x + 1]
             + d3[x + 1];
        dx *= 2;

        dy = -d1[x - 1]
             - 2 * d1[x]
             - d1[x + 1]
             + d3[x - 1]
             + 2 * d3[x]
             + d3[x + 1];

        dx2 = -d1[x - 1]
              + 2 * d1[x]
              - d1[x + 1]
              - d3[x - 1]
              + 2 * d3[x]
              - d3[x + 1];

        if (dy < 0) {
          dy = -dy;
          dx = -dx;
        }

        if (abs(dx) <= 4*abs(dx2)) {
          v = (d1[x] + d3[x] + 1) >> 1;
          f2[x] = BILINEAR;
        } else if (dx < 0) {
          if (dx < -2 * dy) {
            v = reconstruct_h(d1 + x, d3 + x, 0, 0, 0, 16);
            f2[x] = FILTER_1;
          } else if (dx < -dy) {
            v = reconstruct_h(d1 + x, d3 + x, 0, 0, 8, 8);
            f2[x] = FILTER_2;
          } else if (2 * dx < -dy) {
            v = reconstruct_h(d1 + x, d3 + x, 0, 4, 8, 4);
            f2[x] = FILTER_3;
          } else if (3 * dx < -dy) {
            v = reconstruct_h(d1 + x, d3 + x, 1, 7, 7, 1);
            f2[x] = FILTER_4;
          } else {
            v = reconstruct_h(d1 + x, d3 + x, 4, 8, 4, 0);
            f2[x] = FILTER_5;
          }
        } else {
          if (dx > 2 * dy) {
            v = reconstruct_h(d3 + x, d1 + x, 0, 0, 0, 16);
            f2[x] = FILTER_1;
          } else if (dx > dy) {
            v = reconstruct_h(d3 + x, d1 + x, 0, 0, 8, 8);
            f2[x] = FILTER_2;
          } else if (2 * dx > dy) {
            v = reconstruct_h(d3 + x, d1 + x, 0, 4, 8, 4);
            f2[x] = FILTER_3;
          } else if (3 * dx > dy) {
            v = reconstruct_h(d3 + x, d1 + x, 1, 7, 7, 1);
            f2[x] = FILTER_4;
          } else {
            v = reconstruct_h(d3 + x, d1 + x, 4, 8, 4, 0);
            f2[x] = FILTER_5;
          }
        }
        d2[x] = OD_CLAMP255(v);
      } else {
        d2[x] = (d1[x] + d3[x] + 1) >> 1;
        f2[x] = BILINEAR;
      }
    }
    d += 2*dst_stride;
    f += 2*filter_buf_stride;
  }
  {
    unsigned char *d1 = d;
    unsigned char *d2 = d + dst_stride;
    uint8_t *f2 = f + filter_buf_stride;

    for (x = -xpad; x < w + xpad; x++) {
      d2[2 * x] = d1[x * 2];
      d2[2 * x + 1] = d1[x * 2];
      f2[2 * x] = NEAREST_NEIGHBOR;
      f2[2 * x + 1] = NEAREST_NEIGHBOR;
    }
    d += 2*dst_stride;
  }
}

int main(int _argc, char *_argv[]) {
  video_input  vid;
  video_input_info info;
  int64_t  gsqerr = 0;
  int64_t  gplsqerr[INTERP_MAX][FILTERS_MAX] = {{ 0 }};
  int64_t  gplnpixels[INTERP_MAX][FILTERS_MAX] = {{ 0 }};
  int64_t  gnpixels[INTERP_MAX] = { 0 };
  int          frameno;
  FILE *fin;
  int          long_option_index;
  int          c;
  int i, j;
#ifdef _WIN32
  _setmode(_fileno(stdin), _O_BINARY);
#endif
  /*Process option arguments.*/
  while ((c =
   getopt_long(_argc, _argv, optstring, options, &long_option_index)) != EOF) {
    switch (c) {
      case 's': summary_only = 1;
        break;
      case 'a': avg_ds = 1;
        break;
      default: usage(_argv);
        break;
    }
  }
  if (optind+1 != _argc) {
    usage(_argv);
    exit(1);
  }
  fin = strcmp(_argv[optind], "-") == 0 ? stdin : fopen(_argv[optind], "rb");
  if (fin == NULL) {
    fprintf(stderr, "Unable to open '%s' for extraction.\n", _argv[optind]);
    exit(1);
  }
  fprintf(stderr, "Opening %s...\n", _argv[optind]);
  if (video_input_open(&vid, fin) < 0) exit(1);
  video_input_get_info(&vid, &info);

  for (frameno = 0;; frameno++) {
    video_input_ycbcr f;
    int64_t plsqerr[INTERP_MAX][FILTERS_MAX] = {{ 0 }};
    long plnpixels[INTERP_MAX][FILTERS_MAX] = {{ 0 }};
    long npixels[INTERP_MAX] = { 0 };
    int64_t sqerr;
    int ret;

    int h = info.pic_h;
    int w = info.pic_w;
    uint8_t *ds_buf = malloc(w*h/4);
    int ds_buf_stride = w/2;
    uint8_t *edi_buf = malloc(w*h);
    int edi_buf_stride = w;
    uint8_t *filter_buf = malloc(w*h);
    int filter_buf_stride = w;

    int x, y;

    memset(filter_buf, SRC_PIXEL, w*h);

    ret = video_input_fetch_frame(&vid, f, NULL);
    if (ret <= 0) break;
    sqerr = 0;

    for (y = 0; y < h; y += 2) {
      for (x = 0; x < w; x += 2) {
        if (avg_ds) {
          ds_buf[(y/2)*ds_buf_stride + (x/2)] =
              (f[0].data[(info.pic_y + y)*f[0].stride + (info.pic_x + x)] +
               f[0].data[(info.pic_y + y)*f[0].stride + (info.pic_x + x + 1)] +
               f[0].data[(info.pic_y + y + 1)*f[0].stride + (info.pic_x + x)] +
               f[0].data[(info.pic_y + y + 1)*f[0].stride + (info.pic_x + x + 1)] + 1) >> 2;
        } else {
          ds_buf[(y/2)*ds_buf_stride + (x/2)] =
              f[0].data[(info.pic_y + y)*f[0].stride + (info.pic_x + x)];
        }
      }
    }
    edi_hv(edi_buf, edi_buf_stride, filter_buf, filter_buf_stride,
           ds_buf, ds_buf_stride, w/2, h/2);

    for (y = 0; y < h; y ++) {
      for (x = 0 ; x < w; x++) {
        uint8_t fi = filter_buf[y*filter_buf_stride + x];
        int d = edi_buf[y*edi_buf_stride + x] -
                f[0].data[(info.pic_y + y)*f[0].stride + (info.pic_x + x)];
        if (!(fi != SRC_PIXEL || (x & 1) == 0 && (y & 1) == 0))
          printf("%d, %d\n", x, y);
        assert(fi != SRC_PIXEL || (x & 1) == 0 && (y & 1) == 0);
        if ((x & 1) == 0 && (y & 1) == 0) {
          /*assert(d == 0);*/
        } else if ((x & 1) == 1 && (y & 1) == 0) {
          plsqerr[INTERP_H][fi] += d*d;
          plnpixels[INTERP_H][fi]++;
        } else if ((x & 1) == 0 && (y & 1) == 1) {
          plsqerr[INTERP_V][fi] += d*d;
          plnpixels[INTERP_V][fi]++;
        } else if ((x & 1) == 1 && (y & 1) == 1) {
          plsqerr[INTERP_HV][fi] += d*d;
          plnpixels[INTERP_HV][fi]++;
        }
      }
    }
    for (i = 0; i < INTERP_MAX; i++) {
      for (j = 0; j < FILTERS_MAX; j++) {
        assert(plnpixels[INTERP_H][SRC_PIXEL] == 0);
        gplsqerr[i][j] += plsqerr[i][j];
        gplnpixels[i][j] += plnpixels[i][j];
        npixels[i] += plnpixels[i][j];
        gnpixels[i] += plnpixels[i][j];
      }
    }
    if (!summary_only) {
      printf("%08i: ", frameno);
      for (j = 0; j < FILTERS_MAX; j++) {
        printf("%-7G ", (double)plnpixels[INTERP_H][j]/npixels[INTERP_H]);
      }
      printf("\n");
      printf("%08i: ", frameno);
      for (j = 0; j < FILTERS_MAX; j++) {
        printf("%-7G ",
             10*(log10(255*255)+log10(plnpixels[INTERP_H][j])-log10(plsqerr[INTERP_H][j])));
      }
      printf("\n");
    }
  }
  printf("Total: ");
  for (j = 0; j < FILTERS_MAX; j++) {
    printf("%-7G ", (double)gplnpixels[INTERP_H][j]/gnpixels[INTERP_H]);
  }
  printf("\n");
  printf("Total: ");
  for (j = 0; j < FILTERS_MAX; j++) {
    printf("%-7G ",
           10*(log10(255*255)+log10(gplnpixels[INTERP_H][j])-log10(gplsqerr[INTERP_H][j])));
  }
  printf("\n");
  video_input_close(&vid);
  return 0;
}
