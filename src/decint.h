/*Daala video codec
Copyright (c) 2006-2013 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

#if !defined(_decint_H)
# define _decint_H (1)
# include "vp10/common/state.h"

typedef struct daala_dec_ctx od_dec_ctx;

/*Constants for the packet state machine specific to the decoder.*/
/*Next packet to read: Data packet.*/
# define OD_PACKET_DATA (0)

struct daala_dec_ctx {
  od_state state;
  oggbyte_buffer obb;
  od_ec_dec ec;
  int packet_state;
  /*User provided buffer for storing per frame block size information. These
   are set via daala_decode_ctl with OD_DECCTL_SET_BSIZE_BUFFER.*/
  unsigned char *user_bsize;
  int user_bstride;
  /*User provided buffer for storing the band flags per block per frame.  These
   are set via daala_decode_ctl with OD_DECCTL_SET_FLAGS_BUFFER.*/
  unsigned int *user_flags;
  int user_fstride;
  od_mv_grid_pt *user_mv_grid;
  daala_image *user_mc_img;
  od_output_queue out;
#if OD_ACCOUNTING
  int acct_enabled;
  od_accounting_internal acct;
#endif
  /*User provided buffer for storing the deringing filter flags per superblock.
    This is set via daala_decode_ctl with OD_DECCTL_SET_DERING_BUFFER.*/
  unsigned char *user_dering;
};

# if OD_ACCOUNTING
#  define OD_ACCOUNTING_SET_LOCATION(dec, layer, level, x, y) \
  od_accounting_set_location(&(dec)->acct, layer, level, x, y)
# else
#  define OD_ACCOUNTING_SET_LOCATION(dec, layer, level, x, y)
# endif

/*Stub for the daala_setup_info.*/
struct daala_setup_info {
  unsigned char dummy;
};

#endif
