/*Daala video codec
Copyright (c) 2006-2013 Daala project contributors.  All rights reserved.

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

#include <daala/codec.h>
#include <daala/daaladec.h>

#include "SDL.h"

#define ODS_NONE 0
#define ODS_HEADER 1
#define ODS_DATA 2

typedef struct {
  SDL_Surface *screen;
  daala_info di;
  daala_comment dc;
  ogg_sync_state oy;
  FILE *input;
  char *input_path;
  ogg_stream_state os;
  daala_dec_ctx *dctx;
  SDL_Surface *surf;
  od_img img;
  int width;
  int height;
  int done;
  int od_state;
  int paused;
  int restart;
  int slow;
  int loop;
  int step;
  int valid;
} player_example;

static void img_to_rgb(SDL_Surface *surf, const od_img *img);

int player_example_init(player_example *player);
player_example *player_example_create();
int player_example_clear(player_example *player);
int player_example_free(player_example *player);
int player_example_reset(player_example *player);
int player_example_play(player_example *player);
int player_example_restart(player_example *player);
int player_example_open(player_example *player, char *path);

int player_example_daala_stream_init(player_example *player, int serial) {
  int ret;
  if (player == NULL) return -1;
  ret = ogg_stream_init(&player->os, serial);
  if (ret != 0) return -1;
  daala_info_init(&player->di);
  daala_comment_init(&player->dc);
  player->od_state = ODS_HEADER;
  return 0;
}

int player_example_daala_stream_clear(player_example *player) {
  if (player == NULL) return -1;
  daala_info_clear(&player->di);
  daala_comment_clear(&player->dc);
  if (player->dctx != NULL) {
    daala_decode_free(player->dctx);
    player->dctx = NULL;
  }
  ogg_stream_clear(&player->os);
  player->od_state = ODS_NONE;
  return 0;
}

int player_example_init(player_example *player) {
  if (player == NULL) return -1;
  player->screen = NULL;
  player->surf = NULL;
  player->width = 0;
  player->height = 0;
  player->paused = 0;
  player->restart = 0;
  player->slow = 0;
  player->loop = 0;
  player->done = 0;
  player->step = 0;
  player->valid = 0;
  player->od_state = ODS_NONE;
  return 0;
}

int player_example_clear(player_example *player) {
  if (player == NULL) return -1;
  memset(player, 0, sizeof(player_example));
  return 0;
}

player_example *player_example_create() {
  int ret;
  player_example *player;
  player = malloc(sizeof(player_example));
  if (player == NULL) return NULL;
  ret = player_example_init(player);
  if (ret != 0) {
    free(player);
    return NULL;
  }
  return player;
}

int player_example_free(player_example *player) {
  int ret;
  if (player == NULL) return -1;
  ret = player_example_clear(player);
  if (ret == 0) {
    free(player);
    return 0;
  }
  return -1;
}

int player_example_reset(player_example *player) {
  int ret;
  if (player == NULL) return -1;
  ret = player_example_clear(player);
  if (ret != 0) return -1;
  ret = player_example_init(player);
  if (ret != 0) return -1;
  return 0;
}

int player_example_open_input(player_example *player, char *path) {
  if ((player == NULL) || ((path == NULL) || (path[0] == '\0'))) return -1;
  if ((path[0] == '-') && (path[1] == '\0')) {
    player->input = stdin;
  }
  else {
    player->input = fopen(path, "rb");
  }
  if (player->input == NULL) {
    player->input_path = "";
    return -1;
  }
  player->input_path = path;
  ogg_sync_init(&player->oy);
  return 0;
}

int player_example_close_input(player_example *player) {
  int ret;
  if (player == NULL) return -1;
  if ((player->input == stdin) || (player->input == NULL)) return -1;
  ret = fclose(player->input);
  player->input = NULL;
  ogg_sync_clear(&player->oy);
  if (ret != 0) return -1;
  return 0;
}

int player_example_input_restart(player_example *player) {
  int ret;
  if (player == NULL) return -1;
  if (player->input == stdin) return -1;
  ret = player_example_close_input(player);
  if (ret != 0) return -1;
  ret = player_example_open_input(player, player->input_path);
  return ret;
}

void player_example_handle_event(player_example *player, SDL_Event *event) {
  switch (event->type) {
    case SDL_QUIT: {
      player->done = 1;
      break;
    }
    case SDL_KEYDOWN: {
      switch (event->key.keysym.sym) {
        case SDLK_q: {
          player->done = 1;
          break;
        }
        case SDLK_s: {
          player->slow = !player->slow;
          break;
        }
        case SDLK_l: {
          player->loop = !player->loop;
          break;
        }
        case SDLK_ESCAPE: {
          player->done = 1;
          break;
        }
        case SDLK_SPACE: {
          player->paused = !player->paused;
          player->step = 0;
          break;
        }
        case SDLK_PERIOD: {
          player->step = 1;
          player->paused = 1;
          break;
        }
        case SDLK_r: {
          player->restart = 1;
          if (player->paused) {
            player->step = 1;
          }
          break;
        }
        default: break;
      }
      break;
    }
  }
}

void player_example_wait_user_input(player_example *player) {
  SDL_Event event;
  if (SDL_WaitEvent(&event)) {
    player_example_handle_event(player, &event);
  }
}

void player_example_check_user_input(player_example *player) {
  SDL_Event event;
  while (SDL_PollEvent(&event)) {
    player_example_handle_event(player, &event);
  }
}

int player_example_play(player_example *player) {
  size_t bytes;
  char *buffer;
  int ret;
  ogg_page page;
  ogg_packet packet;
  daala_setup_info *dsi;
  dsi = NULL;
  while (!player->done) {
    while (ogg_sync_pageout(&player->oy, &page) != 1) {
      buffer = ogg_sync_buffer(&player->oy, 4096);
      if (buffer == NULL) return -1;
      bytes = fread(buffer, 1, 4096, player->input);
      if (bytes > 0) {
        ret = ogg_sync_wrote(&player->oy, bytes);
        if (ret != 0) return -1;
      }
      else {
        if (!player->valid) {
          fprintf(stderr, "Invalid Ogg\n");
          exit(1);
        }
        if (player->od_state != ODS_NONE) {
          ret = player_example_daala_stream_clear(player);
          if (ret != 0) return -1;
        }
        if (player->input == stdin) {
          return 0;
        }
        if (player->loop == 1) {
          player_example_input_restart(player);
          continue;
        }
        for (;;) {
          player_example_wait_user_input(player);
          if (player->restart) {
            ret = player_example_input_restart(player);
            player->restart = 0;
            if (ret != 0) return -1;
            break;
          }
          if (player->done) {
            return 0;
          }
        }
      }
    }
    if (ogg_page_bos(&page)) {
      ret = player_example_daala_stream_init(player,
       ogg_page_serialno(&page));
      if (ret != 0) return -1;
    }
    ret = ogg_stream_pagein(&player->os, &page);
    if (ret != 0) return -1;
    while (ogg_stream_packetout(&player->os, &packet) == 1) {
      switch (player->od_state) {
        case ODS_HEADER: {
          ret = daala_decode_header_in(&player->di, &player->dc, &dsi, &packet);
          if (ret < 0) return -1;
          if (ret != 0) break;
          player->dctx = daala_decode_alloc(&player->di, dsi);
          if (player->dctx == NULL) return -1;
          daala_setup_free(dsi);
          dsi = NULL;
          player->od_state = ODS_DATA;
          /* Falling through */
        }
        case ODS_DATA: {
          if ((player->screen == NULL)
           || (player->width != player->di.pic_width)
           || (player->height != player->di.pic_height)) {
            player->width = player->di.pic_width;
            player->height = player->di.pic_height;
            player->screen = SDL_SetVideoMode(player->width, player->height,
             24,
             SDL_HWSURFACE | SDL_DOUBLEBUF);
            if (player->screen == NULL) return -1;
            player->surf = SDL_GetVideoSurface();
            if (player->surf == NULL) return -1;
          }
          ret = daala_decode_packet_in(player->dctx, &player->img, &packet);
          if (ret != 0) return -1;
          player->valid = 1;
          if ((player->slow) && (!player->step)) {
            SDL_Delay(420);
          }
          player_example_check_user_input(player);
          while ((player->paused) && (!player->done)) {
            if (player->restart) {
              break;
            }
            if (player->step) {
              player->step = 0;
              break;
            }
            player_example_wait_user_input(player);
          }
          if ((!player->restart) && (!player->done)) {
            SDL_LockSurface(player->surf);
            img_to_rgb(player->surf, &player->img);
            SDL_UnlockSurface(player->surf);
            SDL_Flip(player->screen);
          }
          break;
        }
      }
    }
    if ((player->restart) || (ogg_page_eos(&page))) {
      ret = player_example_daala_stream_clear(player);
      if (ret != 0) return -1;
    }
    if (player->restart) {
      ret = player_example_input_restart(player);
      player->restart = 0;
      if (ret != 0) return -1;
    }
  }
  if (player->od_state != ODS_NONE) {
    ret = player_example_daala_stream_clear(player);
    if (ret != 0) return -1;
  }
  return 0;
}

int main(int argc, char *argv[]) {
  int ret;
  char *input;
  int start_paused;
  player_example *player;
  daala_log_init();
  if ((argc == 3) && (memcmp(argv[1], "-p", 2) == 0)) {
    start_paused = 1;
    input = argv[2];
  }
  else {
    if ((argc != 2) || ((argc == 2)
     && ((memcmp(argv[1], "-h", 2) == 0)
     || (memcmp(argv[1] + 1, "-h", 2) == 0)))) {
      fprintf(stderr, "usage: %s input.ogg\n%s\n", argv[0],
       "\nProgram Options:\n-p to start paused\n- to read from stdin\n\n"
       "Playback Control: \n"
       "r to restart\nl to loop\ns for slow\n. to step\nspace to pause\n"
       "q to quit");
      exit(1);
    }
    start_paused = 0;
    input = argv[1];
  }
  if (SDL_Init(SDL_INIT_VIDEO) < 0) {
    fprintf(stderr, "error: unable to init SDL\n");
    exit(1);
  }
  atexit(SDL_Quit);
  SDL_EnableKeyRepeat(222, 100);

  player = player_example_create();
  if (player == NULL) {
    fprintf(stderr, "player example error: create player\n");
    return -1;
  }
  ret = player_example_open_input(player, input);
  if (ret != 0) {
    fprintf(stderr, "player example error: could not open: %s\n", input);
    player_example_free(player);
    return -1;
  }
  if (start_paused == 1) {
    player->step = 1;
    player->paused = 1;
  }
  ret = player_example_play(player);
  if (ret != 0) {
    fprintf(stderr, "player example error: playback error\n");
    exit(1);
  }
  ret = player_example_free(player);

  return ret;
}

#define OD_MINI(_a, _b)      ((_a) < (_b) ? (_a) : (_b))
#define OD_MAXI(_a, _b)      ((_a) > (_b) ? (_a) : (_b))
#define OD_CLAMPI(_a, _b, _c) (OD_MAXI(_a, OD_MINI(_b, _c)))
#define OD_SIGNMASK(_a)     (-((_a) < 0))
#define OD_FLIPSIGNI(_a, _b) (((_a) + OD_SIGNMASK(_b)) ^ OD_SIGNMASK(_b))
#define OD_DIV_ROUND(_x, _y) (((_x) + OD_FLIPSIGNI((_y)  >>  1, _x))/(_y))
#define OD_CLAMP255(x) \
  ((unsigned char)((((x) < 0) - 1) & ((x) | -((x) > 255))))

void img_to_rgb(SDL_Surface *surf, const od_img *img) {
  unsigned char *y_row;
  unsigned char *cb_row;
  unsigned char *cr_row;
  unsigned char *y;
  unsigned char *cb;
  unsigned char *cr;
  int            y_stride;
  int            cb_stride;
  int            cr_stride;
  int            width;
  int            height;
  int            xdec;
  int            ydec;
  int            i;
  int            j;
  unsigned char *pixels;
  int pitch;
  pixels = (unsigned char *)surf->pixels;
  pitch = surf->pitch;
  width = img->width;
  height = img->height;
  /*Assume both C planes are decimated.*/
  xdec = img->planes[1].xdec;
  ydec = img->planes[1].ydec;
  y_stride = img->planes[0].ystride;
  cb_stride = img->planes[1].ystride;
  cr_stride = img->planes[2].ystride;
  y_row = img->planes[0].data;
  cb_row = img->planes[1].data;
  cr_row = img->planes[2].data;
  /*Chroma up-sampling is just done with a box filter.
    This is very likely what will actually be used in practice on a real
     display, and also removes one more layer to search in for the source of
     artifacts.
    As an added bonus, it's dead simple.*/
  for (j = 0; j < height; j++) {
    int dc;
    y = y_row;
    cb = cb_row;
    cr = cr_row;
    for (i = 0; i < 3 * width;) {
      int64_t  yval;
      int64_t  cbval;
      int64_t  crval;
      unsigned rval;
      unsigned gval;
      unsigned bval;
      yval = *y - 16;
      cbval = *cb - 128;
      crval = *cr - 128;
      /*This is intentionally slow and very accurate.*/
      rval = OD_CLAMPI(0, (int32_t)OD_DIV_ROUND(
       2916394880000LL*yval + 4490222169144LL*crval, 9745792000LL), 65535);
      gval = OD_CLAMPI(0, (int32_t)OD_DIV_ROUND(
       2916394880000LL*yval - 534117096223LL*cbval - 1334761232047LL*crval,
       9745792000LL), 65535);
      bval = OD_CLAMPI(0, (int32_t)OD_DIV_ROUND(
       2916394880000LL*yval + 5290866304968LL*cbval, 9745792000LL), 65535);
      *(pixels + pitch*j + i++) = (unsigned char)(bval >> 8);
      *(pixels + pitch*j + i++) = (unsigned char)(gval >> 8);
      *(pixels + pitch*j + i++) = (unsigned char)(rval >> 8);
      dc = ((y - y_row) & 1) | (1 - xdec);
      y++;
      cb += dc;
      cr += dc;
    }
    y_row += y_stride;
    dc = -((j & 1) | (1 - ydec));
    cb_row += dc & cb_stride;
    cr_row += dc & cr_stride;
  }
}
