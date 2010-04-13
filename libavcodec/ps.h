/*
 * MPEG-4 Parametric Stereo definitions and declarations
 * Copyright (c) 2010 Alex Converse <alex.converse@gmail.com>
 *
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * FFmpeg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FFmpeg; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#ifndef AVCODEC_PS_H
#define AVCODEC_PS_H

#include <stdint.h>

#define PS_MAX_NUM_ENV 4
#define PS_MAX_NR_IIDICC 34
#define PS_MAX_NR_IPDOPD 17

typedef struct {
    int    enable_iid;
    int    iid_mode;
    int    iid_quant;
    int    nr_iid_par;
    int    nr_ipdopd_par;
    int    enable_icc;
    int    icc_mode;
    int    nr_icc_par;
    int    enable_ext;
    int    frame_class;
    int    num_env_old;
    int    num_env;
    int    enable_ipdopd;
    int    border_position[PS_MAX_NUM_ENV];
    int8_t iid_dt [PS_MAX_NUM_ENV];
    int8_t iid_par[PS_MAX_NUM_ENV][PS_MAX_NR_IIDICC]; //<Inter-channel Intensity Difference Parameters
    int8_t icc_dt [PS_MAX_NUM_ENV];
    int8_t icc_par[PS_MAX_NUM_ENV][PS_MAX_NR_IIDICC]; //<Inter-Channel Coherence Parameters
    int8_t ipd_dt [PS_MAX_NUM_ENV];
    int8_t ipd_par[PS_MAX_NUM_ENV][PS_MAX_NR_IPDOPD]; //<Inter-channel Phase Difference Parameters
    int8_t opd_dt [PS_MAX_NUM_ENV];
    int8_t opd_par[PS_MAX_NUM_ENV][PS_MAX_NR_IPDOPD]; //<Overall Phase Difference Parameters
    int    is34bands;
    int    is34bands_old;

    float  in_buf[64][44][2];
} PSContext;

void ff_ps_init(void);
int ff_ps_data(GetBitContext *gb, PSContext *ps);
int ff_ps_apply(AVCodecContext *avctx, PSContext *ps, float L[2][38][64], float R[2][38][64]);

#endif /* AVCODEC_PS_H */
