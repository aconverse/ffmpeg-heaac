/*
 * Header file for hardcoded Parametric Stereo tables
 *
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

#ifndef PS_TABLEGEN_H
#define PS_TABLEGEN_H

#include <stdint.h>
#include <math.h>

#if CONFIG_HARDCODED_TABLES
#define ps_tableinit()
#include "libavcodec/ps_tables.h"
#else
#ifndef M_SQRT1_2
#define M_SQRT1_2      0.70710678118654752440  /* 1/sqrt(2) */
#endif
static float pd_re_smooth[8*8*8];
static float pd_im_smooth[8*8*8];

static void ps_tableinit(void)
{
    static const float ipdopd_sin[] = { 0, M_SQRT1_2, 1,  M_SQRT1_2,  0, -M_SQRT1_2, -1, -M_SQRT1_2 };
    static const float ipdopd_cos[] = { 1, M_SQRT1_2, 0, -M_SQRT1_2, -1, -M_SQRT1_2,  0,  M_SQRT1_2 };
    int pd0, pd1, pd2;
    for (pd0 = 0; pd0 < 8; pd0++) {
        float pd0_re = ipdopd_cos[pd0];
        float pd0_im = ipdopd_sin[pd0];
        for (pd1 = 0; pd1 < 8; pd1++) {
            float pd1_re = ipdopd_cos[pd1];
            float pd1_im = ipdopd_sin[pd1];
            for (pd2 = 0; pd2 < 8; pd2++) {
                float pd2_re = ipdopd_cos[pd2];
                float pd2_im = ipdopd_sin[pd2];
                float re_smooth = 0.25f * pd0_re + 0.5f * pd1_re + pd2_re;
                float im_smooth = 0.25f * pd0_im + 0.5f * pd1_im + pd2_im;
                float pd_mag = 1 / sqrt(im_smooth * im_smooth + re_smooth * re_smooth);
                pd_re_smooth[pd0*64+pd1*8+pd2] = re_smooth * pd_mag;
                pd_im_smooth[pd0*64+pd1*8+pd2] = im_smooth * pd_mag;
            }
        }
    }
}
#endif /* CONFIG_HARDCODED_TABLES */

#endif /* PS_TABLEGEN_H */
