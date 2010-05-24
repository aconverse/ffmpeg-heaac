/*
 * MPEG-4 Parametric Stereo decoding functions
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

#include <stdint.h>
#include "libavutil/mathematics.h"
#include "avcodec.h"
#include "get_bits.h"
#include "ps.h"
#include "psdata.c"

#define PS_BASELINE 0

#define numQMFSlots 32 //numTimeSlots * RATE

static const int8_t num_env_tab[2][4] = {
    { 0, 1, 2, 4, },
    { 1, 2, 3, 4, },
};

static const int8_t nr_iidicc_par_tab[] = {
    10, 20, 34, 10, 20, 34,
};

static const int8_t nr_iidopd_par_tab[] = {
     5, 11, 17,  5, 11, 17,
};

enum {
    huff_iid_df1,
    huff_iid_dt1,
    huff_iid_df0,
    huff_iid_dt0,
    huff_icc_df,
    huff_icc_dt,
    huff_ipd_df,
    huff_ipd_dt,
    huff_opd_df,
    huff_opd_dt,
};

static const int huff_iid[] = {
    huff_iid_df0,
    huff_iid_df1,
    huff_iid_dt0,
    huff_iid_dt1,
};

static VLC vlc_ps[10];

static int iid_data(AVCodecContext *avctx, GetBitContext *gb, PSContext *ps, int e, int dt)
{
    int b;
    int table_idx = huff_iid[2*dt+ps->iid_quant];
    VLC_TYPE (*vlc_table)[2] = vlc_ps[table_idx].table;
    if (dt) {
        int e_prev = e ? e - 1 : ps->num_env_old - 1;
        e_prev = FFMAX(e_prev, 0);
        for (b = 0; b < ps->nr_iid_par; b++) {
            ps->iid_par[e][b] = ps->iid_par[e_prev][b] +
                                get_vlc2(gb, vlc_table, 9, 3) -
                                huff_offset[table_idx];
            if (FFABS(ps->iid_par[e][b]) > 7 + 8 * ps->iid_quant) {
                av_log(avctx, AV_LOG_ERROR, "illegal iid\n");
                return -1;
            }
        }
    } else {
        int prev = 0;
        for (b = 0; b < ps->nr_iid_par; b++) {
            prev += get_vlc2(gb, vlc_table, 9, 3) -
                    huff_offset[table_idx];
            ps->iid_par[e][b] = prev;
            if (FFABS(ps->iid_par[e][b]) > 7 + 8 * ps->iid_quant) {
                av_log(avctx, AV_LOG_ERROR, "illegal iid\n");
                return -1;
            }
        }
    }
    return 0;
}

static int icc_data(AVCodecContext *avctx, GetBitContext *gb, PSContext *ps, int e, int dt)
{
    int b;
    int table_idx = dt ? huff_icc_dt : huff_icc_df;
    VLC_TYPE (*vlc_table)[2] = vlc_ps[table_idx].table;
    if (dt) {
        int e_prev = e ? e - 1 : ps->num_env_old - 1;
        e_prev = FFMAX(e_prev, 0);
        for (b = 0; b < ps->nr_icc_par; b++) {
            ps->icc_par[e][b] = ps->icc_par[e_prev][b] + get_vlc2(gb, vlc_table, 9, 3) - huff_offset[table_idx];
            if (ps->icc_par[e][b] > 7U) {
                av_log(avctx, AV_LOG_ERROR, "illegal icc\n");
                return -1;
            }
        }
    } else {
        int prev = 0;
        for (b = 0; b < ps->nr_icc_par; b++) {
            prev += get_vlc2(gb, vlc_table, 9, 3) - huff_offset[table_idx];
            ps->icc_par[e][b] = prev;
            if (ps->icc_par[e][b] > 7U) {
                av_log(avctx, AV_LOG_ERROR, "illegal icc\n");
                return -1;
            }
        }
    }
    return 0;
}

static void ipd_data(GetBitContext *gb, PSContext *ps, int e, int dt)
{
    int b;
    int table_idx = dt ? huff_ipd_dt : huff_ipd_df;
    VLC_TYPE (*vlc_table)[2] = vlc_ps[table_idx].table;
    if (dt) {
        int e_prev = e ? e - 1 : ps->num_env_old - 1;
        e_prev = FFMAX(e_prev, 0);
        for (b = 0; b < ps->nr_ipdopd_par; b++) {
            ps->ipd_par[e][b] = (ps->ipd_par[e_prev][b] + get_vlc2(gb, vlc_table, 9, 1)) & 0x07;
        }
    } else {
        int prev = 0;
        for (b = 0; b < ps->nr_ipdopd_par; b++) {
            prev += get_vlc2(gb, vlc_table, 9, 3);
            prev &= 0x07;
            ps->ipd_par[e][b] = prev;
        }
    }
}

static void opd_data(GetBitContext *gb, PSContext *ps, int e, int dt)
{
    int b;
    int table_idx = dt ? huff_opd_dt : huff_opd_df;
    VLC_TYPE (*vlc_table)[2] = vlc_ps[table_idx].table;
    if (dt) {
        int e_prev = e ? e - 1 : ps->num_env_old - 1;
        e_prev = FFMAX(e_prev, 0);
        for (b = 0; b < ps->nr_ipdopd_par; b++) {
            ps->opd_par[e][b] = (ps->opd_par[e_prev][b] + get_vlc2(gb, vlc_table, 9, 1)) & 0x07;
        }
    } else {
        int prev = 0;
        for (b = 0; b < ps->nr_ipdopd_par; b++) {
            prev += get_vlc2(gb, vlc_table, 9, 3);
            prev &= 0x07;
            ps->opd_par[e][b] = prev;
        }
    }
}

static int ps_extension(GetBitContext *gb, PSContext *ps, int ps_extension_id)
{
    int e;
    int count = get_bits_count(gb);
    if (!ps_extension_id) {
        ps->enable_ipdopd = get_bits1(gb);
        if (ps->enable_ipdopd) {
            for (e = 0; e < ps->num_env; e++) {
                int dt = get_bits1(gb);
                ipd_data(gb, ps, e, dt);
                dt = get_bits1(gb);
                opd_data(gb, ps, e, dt);
            }
        }
        skip_bits1(gb);      //reserved_ps
    }
    return get_bits_count(gb) - count;
}

static void ipdopd_reset(float (*opd_smooth)[2][2], float (*ipd_smooth)[2][2])
{
    int i;
    for (i = 0; i < PS_MAX_NR_IPDOPD; i++) {
        opd_smooth[i][0][0] = opd_smooth[i][1][0] = 1;
        ipd_smooth[i][0][0] = ipd_smooth[i][1][0] = 1;
        opd_smooth[i][0][1] = opd_smooth[i][1][1] = 0;
        ipd_smooth[i][0][1] = ipd_smooth[i][1][1] = 0;
    }
}

int ff_ps_read_data(AVCodecContext *avctx, GetBitContext *gb_host, PSContext *ps, int bits_left)
{
    int e;
    int bit_count_start = get_bits_count(gb_host);
    int header;
    int bits_consumed;
    GetBitContext gbc = *gb_host, *gb = &gbc;

    header = get_bits1(gb);
    if (header) {     //enable_ps_header
        ps->enable_iid = get_bits1(gb);
        if (ps->enable_iid) {
            ps->iid_mode = get_bits(gb, 3);
            if (ps->iid_mode > 5) {
                av_log(avctx, AV_LOG_ERROR, "iid_mode %d is reserved.\n",
                       ps->iid_mode);
                goto err;
            }
            ps->nr_iid_par    = nr_iidicc_par_tab[ps->iid_mode];
            ps->iid_quant     = ps->iid_mode > 2;
            ps->nr_ipdopd_par = nr_iidopd_par_tab[ps->iid_mode];
        }
        ps->enable_icc = get_bits1(gb);
        if (ps->enable_icc) {
            ps->icc_mode = get_bits(gb, 3);
            if (ps->icc_mode > 5) {
                av_log(avctx, AV_LOG_ERROR, "icc_mode %d is reserved.\n",
                       ps->icc_mode);
                goto err;
            }
            ps->nr_icc_par = nr_iidicc_par_tab[ps->icc_mode];
        }
        ps->enable_ext = get_bits1(gb);
    }

    ps->frame_class = get_bits1(gb);
    ps->num_env_old = ps->num_env;
    ps->num_env     = num_env_tab[ps->frame_class][get_bits(gb, 2)];

    ps->border_position[0] = -1;
    if (ps->frame_class) {
        for (e = 1; e <= ps->num_env; e++)
            ps->border_position[e] = get_bits(gb, 5);
    } else
        for (e = 1; e <= ps->num_env; e++)
            ps->border_position[e] = e * numQMFSlots / ps->num_env - 1;

    if (ps->enable_iid)
        for (e = 0; e < ps->num_env; e++) {
            int dt = get_bits1(gb);
            if (iid_data(avctx, gb, ps, e, dt))
                goto err;
        }
    else
        memset(ps->iid_par, 0, sizeof(ps->iid_par));

    if (ps->enable_icc)
        for (e = 0; e < ps->num_env; e++) {
            int dt = get_bits1(gb);
            if (icc_data(avctx, gb, ps, e, dt))
                goto err;
        }
    else
        memset(ps->icc_par, 0, sizeof(ps->icc_par));

    if (ps->enable_ext) {
        int cnt = get_bits(gb, 4);
        if (cnt == 15) {
            cnt += get_bits(gb, 8);
        }
        cnt *= 8;
        while (cnt > 7) {
            int ps_extension_id = get_bits(gb, 2);
            cnt -= 2 + ps_extension(gb, ps, ps_extension_id);
        }
        if (cnt < 0) {
            av_log(avctx, AV_LOG_ERROR, "ps extension overflow %d", cnt);
            goto err;
        }
        skip_bits(gb, cnt);
    }

    ps->enable_ipdopd &= !PS_BASELINE;

    //Fix up envelopes
    if (!ps->num_env) {
        ps->num_env = 1;
        ps->border_position[1] = 31;
        if (ps->enable_iid && ps->num_env_old > 1) {
            memcpy(ps->iid_par, ps->iid_par+ps->num_env_old-1, sizeof(ps->iid_par[0]));
        }
        if (ps->enable_icc && ps->num_env_old > 1) {
            memcpy(ps->icc_par, ps->icc_par+ps->num_env_old-1, sizeof(ps->icc_par[0]));
        }
        if (ps->enable_ipdopd && ps->num_env_old > 1) {
            memcpy(ps->ipd_par, ps->ipd_par+ps->num_env_old-1, sizeof(ps->ipd_par[0]));
            memcpy(ps->opd_par, ps->opd_par+ps->num_env_old-1, sizeof(ps->opd_par[0]));
        }
    } else if (ps->border_position[ps->num_env] < numQMFSlots - 1) {
        //Create a fake envelope
        if (ps->enable_iid && ps->num_env_old > 1) {
            memcpy(ps->iid_par+ps->num_env, ps->iid_par+ps->num_env-1, sizeof(ps->iid_par[0]));
        }
        if (ps->enable_icc && ps->num_env_old > 1) {
            memcpy(ps->icc_par+ps->num_env, ps->icc_par+ps->num_env-1, sizeof(ps->icc_par[0]));
        }
        if (ps->enable_ipdopd) {
            memcpy(ps->ipd_par+ps->num_env, ps->ipd_par+ps->num_env-1, sizeof(ps->ipd_par[0]));
            memcpy(ps->opd_par+ps->num_env, ps->opd_par+ps->num_env-1, sizeof(ps->opd_par[0]));
        }
        ps->num_env++;
        ps->border_position[ps->num_env] = numQMFSlots - 1;
    }


    ps->is34bands_old = ps->is34bands;
    if (!PS_BASELINE && (ps->enable_iid || ps->enable_icc))
        ps->is34bands = (ps->enable_iid && ps->nr_iid_par == 34) ||
                        (ps->enable_icc && ps->nr_icc_par == 34);

    //Baseline
    if (!ps->enable_ipdopd) {
        memset(ps->ipd_par, 0, sizeof(ps->ipd_par));
        memset(ps->opd_par, 0, sizeof(ps->opd_par));
    }

    if (header)
        ps->start = 1;

    bits_consumed = get_bits_count(gb) - bit_count_start;
    if (bits_consumed <= bits_left) {
        skip_bits_long(gb_host, bits_consumed);
        return bits_consumed;
    }
    av_log(avctx, AV_LOG_ERROR, "Expected to read %d PS bits actually read %d.\n", bits_left, bits_consumed);
err:
    ps->start = 0;
    skip_bits_long(gb_host, bits_left);
    return bits_left;
}

#if !CONFIG_HARDCODED_TABLES
static void make_filters_from_proto(float (*filter)[7][2], const float *proto, int bands)
{
    int q, n;
    //av_log(NULL, AV_LOG_ERROR, "{\n");
    for (q = 0; q < bands; q++) {
        //av_log(NULL, AV_LOG_ERROR, "    {\n        ");
        for (n = 0; n < 7; n++) {
            double theta = 2 * M_PI * (q + 0.5) * (n - 6) / bands;
            filter[q][n][0] = proto[n] *  cos(theta);
            filter[q][n][1] = proto[n] * -sin(theta);
            //av_log(NULL, AV_LOG_ERROR, "{ %13.10f, %13.10f  }, ", filter[q][n][0], filter[q][n][1]);
            //if ((n & 3) == 3) av_log(NULL, AV_LOG_ERROR, "\n        ");
        }
        //av_log(NULL, AV_LOG_ERROR, "\n    },\n");
    }
    //av_log(NULL, AV_LOG_ERROR, "};\n");
}
#endif

/** Split one subband into 2 subsubbands with a real filter */
static av_noinline void hybrid2_re(float (*in)[2], float (*out)[32][2], const float filter[7], int len, int reverse)
{
    int i, j;
    for (i = 0; i < len; i++) {
        float re_in = filter[6] * in[6+i][0];        //real inphase
        float re_op = 0.0f;                          //real out of phase
        float im_in = filter[6] * in[6+i][1];        //imag inphase
        float im_op = 0.0f;                          //imag out of phase
        for (j = 0; j < 6; j += 2) {
            re_in += filter[j  ] * (in[i+j  ][0] + in[12-j  +i][0]);
            im_in += filter[j  ] * (in[i+j  ][1] + in[12-j  +i][1]);
            re_op += filter[j+1] * (in[i+j+1][0] + in[12-j-1+i][0]);
            im_op += filter[j+1] * (in[i+j+1][1] + in[12-j-1+i][1]);
        }
        out[  reverse][i][0] = re_in + re_op;
        out[  reverse][i][1] = im_in + im_op;
        out[1-reverse][i][0] = re_in - re_op;
        out[1-reverse][i][1] = im_in - im_op;
    }
}

/** Split one subband into 6 subsubbands with a complex filter */
static av_noinline void hybrid6_cx(float (*in)[2], float (*out)[32][2], const float (*filter)[7][2], int len)
{
    int i, j, ssb;
    int N = 8;
    float temp[8][2];

    for (i = 0; i < len; i++) {
        for (ssb = 0; ssb < N; ssb++) {
            float sum_re = filter[ssb][6][0] * in[i+6][0], sum_im = filter[ssb][6][0] * in[i+6][1];
            for (j = 0; j < 6; j++) {
                float in0_re = in[i+j][0];
                float in0_im = in[i+j][1];
                float in1_re = in[i+12-j][0];
                float in1_im = in[i+12-j][1];
                sum_re += filter[ssb][j][0] * (in0_re + in1_re) - filter[ssb][j][1] * (in0_im - in1_im);
                sum_im += filter[ssb][j][0] * (in0_im + in1_im) + filter[ssb][j][1] * (in0_re - in1_re);
            }
            temp[ssb][0] = sum_re;
            temp[ssb][1] = sum_im;
        }
        out[0][i][0] = temp[6][0];
        out[0][i][1] = temp[6][1];
        out[1][i][0] = temp[7][0];
        out[1][i][1] = temp[7][1];
        out[2][i][0] = temp[0][0];
        out[2][i][1] = temp[0][1];
        out[3][i][0] = temp[1][0];
        out[3][i][1] = temp[1][1];
        out[4][i][0] = temp[2][0] + temp[5][0];
        out[4][i][1] = temp[2][1] + temp[5][1];
        out[5][i][0] = temp[3][0] + temp[4][0];
        out[5][i][1] = temp[3][1] + temp[4][1];
    }
}

static av_noinline void hybrid4_8_12_cx(float (*in)[2], float (*out)[32][2], const float (*filter)[7][2], int N, int len)
{
    int i, j, ssb;

    for (i = 0; i < len; i++) {
        for (ssb = 0; ssb < N; ssb++) {
            float sum_re = filter[ssb][6][0] * in[i+6][0], sum_im = filter[ssb][6][0] * in[i+6][1];
            for (j = 0; j < 6; j++) {
                float in0_re = in[i+j][0];
                float in0_im = in[i+j][1];
                float in1_re = in[i+12-j][0];
                float in1_im = in[i+12-j][1];
                sum_re += filter[ssb][j][0] * (in0_re + in1_re) - filter[ssb][j][1] * (in0_im - in1_im);
                sum_im += filter[ssb][j][0] * (in0_im + in1_im) + filter[ssb][j][1] * (in0_re - in1_re);
            }
            out[ssb][i][0] = sum_re;
            out[ssb][i][1] = sum_im;
        }
    }
}

static av_noinline void hybrid_analysis(float out[91][32][2], float in[64][44][2], int is34, int len)
{
    int i;
    if(is34) {
        hybrid4_8_12_cx(in[0], out,    f34_0_12, 12, len);
        hybrid4_8_12_cx(in[1], out+12, f34_1_8,   8, len);
        hybrid4_8_12_cx(in[2], out+20, f34_2_4,   4, len);
        hybrid4_8_12_cx(in[3], out+24, f34_2_4,   4, len);
        hybrid4_8_12_cx(in[4], out+28, f34_2_4,   4, len);
        for (i = 0; i < 59; i++) {
            memcpy(out[32 + i], in[5 + i]+6, len * sizeof(in[0][0]));
        }
    } else {
        hybrid6_cx(in[0], out, f20_0_8, len);
        hybrid2_re(in[1], out+6, g1_Q2, len, 1);
        hybrid2_re(in[2], out+8, g1_Q2, len, 0);
        for (i = 0; i < 61; i++) {
            memcpy(out[10 + i], in[3 + i]+6, len * sizeof(in[0][0]));
        }
    }
    //update in_buf
    for (i = 0; i < 5; i++) {
        memcpy(in[i], in[i]+32, 6 * sizeof(in[i][0]));
    }
}

static av_noinline void hybrid_synthesis(float out[64][32][2], float in[91][32][2], int is34, int len)
{
    int i, n;
    if(is34) {
        memset(out, 0, 5*sizeof(out[0]));
        for (n = 0; n < len; n++) {
            for(i = 0; i < 12; i++) {
                out[0][n][0] += in[   i][n][0];
                out[0][n][1] += in[   i][n][1];
            }
            for(i = 0; i < 8; i++) {
                out[1][n][0] += in[12+i][n][0];
                out[1][n][1] += in[12+i][n][1];
            }
            for(i = 0; i < 4; i++) {
                out[2][n][0] += in[20+i][n][0];
                out[2][n][1] += in[20+i][n][1];
                out[3][n][0] += in[24+i][n][0];
                out[3][n][1] += in[24+i][n][1];
                out[4][n][0] += in[28+i][n][0];
                out[4][n][1] += in[28+i][n][1];
            }
        }
        for (i = 0; i < 59; i++) {
            memcpy(out[5 + i], in[32 + i], len * sizeof(in[0][0]));
        }
    } else {
        for (n = 0; n < len; n++) {
            out[0][n][0] = in[0][n][0] + in[1][n][0] + in[2][n][0] +
                           in[3][n][0] + in[4][n][0] + in[5][n][0];
            out[0][n][1] = in[0][n][1] + in[1][n][1] + in[2][n][1] +
                           in[3][n][1] + in[4][n][1] + in[5][n][1];
            out[1][n][0] = in[6][n][0] + in[7][n][0];
            out[1][n][1] = in[6][n][1] + in[7][n][1];
            out[2][n][0] = in[8][n][0] + in[9][n][0];
            out[2][n][1] = in[8][n][1] + in[9][n][1];
        }
        for (i = 0; i < 61; i++) {
            memcpy(out[3 + i], in[10 + i], len * sizeof(in[0][0]));
        }
    }
}

/// All-pass filter decay slope
#define DECAY_SLOPE      0.05f
/// Number of frequency bands that can be addressed by the parameter index, b(k)
static const int   NR_PAR_BANDS[]      = { 20, 34 };
/// Number of frequency bands that can be addressed by the sub subband index, k
static const int   NR_BANDS[]          = { 71, 91 };
/// Start frequency band for the all-pass filter decay slope
static const int   DECAY_CUTOFF[]      = { 10, 32 };
#define NR_ALLPASS_BANDS20 30
#define NR_ALLPASS_BANDS34 50
/// Number of all-pass filer bands
static const int   NR_ALLPASS_BANDS[]  = { 30, 50 };
/// First stereo band using the short one sample delay
static const int   SHORT_DELAY_BAND[]  = { 42, 63 };

/** Table 8.46 */
#define MAP_GENERIC_10_TO_20(out, in, full) \
    int b;                                        \
    if (full)                                     \
        b = 9;                                    \
    else {                                        \
        b = 4;                                    \
        out[10] = 0;                              \
    }                                             \
    for (; b >= 0; b--) {                         \
        out[2*b+1] = out[2*b] = in[b];            \
    }

static void map_idx_10_to_20(int8_t *par_mapped, const int8_t *par, int full)
{
    MAP_GENERIC_10_TO_20(par_mapped, par, full)
}

#define MAP_GENERIC_34_TO_20(out, in, full) \
    out[ 0] = (2*in[ 0] +   in[ 1]) / 3;                               \
    out[ 1] = (  in[ 1] + 2*in[ 2]) / 3;                               \
    out[ 2] = (2*in[ 3] +   in[ 4]) / 3;                               \
    out[ 3] = (  in[ 4] + 2*in[ 5]) / 3;                               \
    out[ 4] = (  in[ 6] +   in[ 7]) / 2;                               \
    out[ 5] = (  in[ 8] +   in[ 9]) / 2;                               \
    out[ 6] =    in[10];                                               \
    out[ 7] =    in[11];                                               \
    out[ 8] = (  in[12] +   in[13]) / 2;                               \
    out[ 9] = (  in[14] +   in[15]) / 2;                               \
    out[10] =    in[16];                                               \
    if (full) {                                                        \
        out[11] =    in[17];                                           \
        out[12] =    in[18];                                           \
        out[13] =    in[19];                                           \
        out[14] = (  in[20] +   in[21]) / 2;                           \
        out[15] = (  in[22] +   in[23]) / 2;                           \
        out[16] = (  in[24] +   in[25]) / 2;                           \
        out[17] = (  in[26] +   in[27]) / 2;                           \
        out[18] = (  in[28] +   in[29] +   in[30] +   in[31]) / 4;     \
        out[19] = (  in[32] +   in[33]) / 2;                           \
    }

static void map_idx_34_to_20(int8_t *par_mapped, const int8_t *par, int full)
{
    MAP_GENERIC_34_TO_20(par_mapped, par, full)
}

static void map_val_34_to_20(float  par[PS_MAX_NUM_ENV][PS_MAX_NR_IIDICC], int e)
{
    MAP_GENERIC_34_TO_20(par[e], par[e], 1)
}

#define MAP_GENERIC_20_TO_34(out, in, full) \
    if (full) {                             \
        out[33] =  in[19];                 \
        out[32] =  in[19];                 \
        out[31] =  in[18];                 \
        out[30] =  in[18];                 \
        out[29] =  in[18];                 \
        out[28] =  in[18];                 \
        out[27] =  in[17];                 \
        out[26] =  in[17];                 \
        out[25] =  in[16];                 \
        out[24] =  in[16];                 \
        out[23] =  in[15];                 \
        out[22] =  in[15];                 \
        out[21] =  in[14];                 \
        out[20] =  in[14];                 \
        out[19] =  in[13];                 \
        out[18] =  in[12];                 \
        out[17] =  in[11];                 \
    }                                      \
    out[16] =  in[10];                     \
    out[15] =  in[ 9];                     \
    out[14] =  in[ 9];                     \
    out[13] =  in[ 8];                     \
    out[12] =  in[ 8];                     \
    out[11] =  in[ 7];                     \
    out[10] =  in[ 6];                     \
    out[ 9] =  in[ 5];                     \
    out[ 8] =  in[ 5];                     \
    out[ 7] =  in[ 4];                     \
    out[ 6] =  in[ 4];                     \
    out[ 5] =  in[ 3];                     \
    out[ 4] = (in[ 2] + in[ 3]) / 2;       \
    out[ 3] =  in[ 2];                     \
    out[ 2] =  in[ 1];                     \
    out[ 1] = (in[ 0] + in[ 1]) / 2;       \
    out[ 0] =  in[ 0];

static void map_idx_20_to_34(int8_t *par_mapped, const int8_t *par, int full)
{
    MAP_GENERIC_20_TO_34(par_mapped, par, full)
}

static void map_val_20_to_34(float  par[PS_MAX_NUM_ENV][PS_MAX_NR_IIDICC], int e)
{
    MAP_GENERIC_20_TO_34(par[e], par[e], 1)
}

static av_noinline void decorrelation(PSContext *ps, float (*out)[32][2], const float (*s)[32][2], int is34)
{
    float power[34][PS_QMF_TIME_SLOTS];
    float transient_gain[34][PS_QMF_TIME_SLOTS];
    float *peak_decay_nrg = ps->peak_decay_nrg;
    float *power_smooth = ps->power_smooth;
    float *peak_decay_diff_smooth = ps->peak_decay_diff_smooth;
    float (*delay)[PS_QMF_TIME_SLOTS + PS_MAX_DELAY][2] = ps->delay;
    float (*ap_delay)[PS_AP_LINKS][PS_QMF_TIME_SLOTS + PS_MAX_AP_DELAY][2] = ps->ap_delay;
    const int8_t *k_to_i = is34 ? k_to_i_34 : k_to_i_20;
    const float peak_decay_factor = 0.76592833836465f;
    const float transient_impact  = 1.5f;
    const float a_smooth          = 0.25f; //< Smoothing coefficient
    int i, k, m, n;
    int n0 = 0, nL = 32;
    static const int link_delay[] = { 3, 4, 5 };
    static const float a[] = { 0.65143905753106f,
                               0.56471812200776f,
                               0.48954165955695f };

    if (is34 != ps->is34bands_old) {
        memset(ps->peak_decay_nrg,         0, sizeof(ps->peak_decay_nrg));
        memset(ps->power_smooth,           0, sizeof(ps->power_smooth));
        memset(ps->peak_decay_diff_smooth, 0, sizeof(ps->peak_decay_diff_smooth));
        memset(ps->delay,                  0, sizeof(ps->delay));
        memset(ps->ap_delay,               0, sizeof(ps->ap_delay));
    }

    memset(power, 0, sizeof(power));
    for (n = n0; n < nL; n++) {
        for (k = 0; k < NR_BANDS[is34]; k++) {
            int i = k_to_i[k];
            power[i][n] += s[k][n][0] * s[k][n][0] + s[k][n][1] * s[k][n][1];
        }
    }

    //Transient detection
    for (i = 0; i < NR_PAR_BANDS[is34]; i++) {
        for (n = n0; n < nL; n++) {
            float decayed_peak = peak_decay_factor * peak_decay_nrg[i];
            peak_decay_nrg[i] = (decayed_peak < power[i][n]) ? power[i][n] : decayed_peak;
            power_smooth[i] = a_smooth * power[i][n] + (1.0f - a_smooth) * power_smooth[i];
            peak_decay_diff_smooth[i] = a_smooth * (peak_decay_nrg[i] - power[i][n]) +
                                         (1.0f - a_smooth) * peak_decay_diff_smooth[i];
            transient_gain[i][n]   = (transient_impact * peak_decay_diff_smooth[i] > power_smooth[i]) ?
                                         power_smooth[i] / (transient_impact * peak_decay_diff_smooth[i]) : 1.0f;
//av_log(NULL, AV_LOG_ERROR, "transient_gain[%2d][%2d] %f %f %f\n", i, n, transient_gain[i][n], peak_decay_diff_smooth[i], power_smooth[i]);
        }
    }

    //Decorrelation and transient reduction
    //                         PS_AP_LINKS - 1
    //                               -----
    //                                | |  Q_fract_allpass[k][m]*z^-link_delay[m] - a[m]*g_decay_slope[k]
    //H[k][z] = z^-2 * phi_fract[k] * | | ----------------------------------------------------------------
    //                                | | 1 - a[m]*g_decay_slope[k]*Q_fract_allpass[k][m]*z^-link_delay[m]
    //                               m = 0
    //d[k][z] (out) = transient_gain_mapped[k][z] * H[k][z] * s[k][z]
    for (k = 0; k < NR_ALLPASS_BANDS[is34]; k++) {
        int b = k_to_i[k];
        float g_decay_slope = 1.f - DECAY_SLOPE * (k - DECAY_CUTOFF[is34]);
        g_decay_slope = FFMIN(g_decay_slope, 1.f);
        g_decay_slope = FFMAX(g_decay_slope, 0.f);
        memcpy(delay[k], delay[k]+nL, PS_MAX_DELAY*sizeof(delay[k][0]));
        memcpy(delay[k]+PS_MAX_DELAY, s[k], numQMFSlots*sizeof(delay[k][0]));
        for (m = 0; m < PS_AP_LINKS; m++) {
            memcpy(ap_delay[k][m],   ap_delay[k][m]+numQMFSlots,           5*sizeof(ap_delay[k][m][0]));
        }
        for (n = n0; n < nL; n++) {
            float in_re = delay[k][n+PS_MAX_DELAY-2][0] * phi_fract[is34][k][0] -
                          delay[k][n+PS_MAX_DELAY-2][1] * phi_fract[is34][k][1];
            float in_im = delay[k][n+PS_MAX_DELAY-2][0] * phi_fract[is34][k][1] +
                          delay[k][n+PS_MAX_DELAY-2][1] * phi_fract[is34][k][0];
            for (m = 0; m < PS_AP_LINKS; m++) {
                float ag                  = a[m] * g_decay_slope;
                float a_re                = ag * in_re;
                float a_im                = ag * in_im;
                float link_delay_re       = ap_delay[k][m][n+5-link_delay[m]][0];
                float link_delay_im       = ap_delay[k][m][n+5-link_delay[m]][1];
                float fractional_delay_re = Q_fract_allpass[is34][k][m][0];
                float fractional_delay_im = Q_fract_allpass[is34][k][m][1];
                ap_delay[k][m][n+5][0] = in_re;
                ap_delay[k][m][n+5][1] = in_im;
                in_re = link_delay_re * fractional_delay_re - link_delay_im * fractional_delay_im - a_re;
                in_im = link_delay_re * fractional_delay_im + link_delay_im * fractional_delay_re - a_im;
                ap_delay[k][m][n+5][0] += ag * in_re;
                ap_delay[k][m][n+5][1] += ag * in_im;
            }
            out[k][n][0] = transient_gain[b][n] * in_re;
            out[k][n][1] = transient_gain[b][n] * in_im;
        }
    }
    for (; k < SHORT_DELAY_BAND[is34]; k++) {
        memcpy(delay[k], delay[k]+nL, PS_MAX_DELAY*sizeof(delay[k][0]));
        memcpy(delay[k]+PS_MAX_DELAY, s[k], numQMFSlots*sizeof(delay[k][0]));
        for (n = n0; n < nL; n++) {
            //H = delay 14
            out[k][n][0] = transient_gain[k_to_i[k]][n] * delay[k][n+PS_MAX_DELAY-14][0];
            out[k][n][1] = transient_gain[k_to_i[k]][n] * delay[k][n+PS_MAX_DELAY-14][1];
        }
    }
    for (; k < NR_BANDS[is34]; k++) {
        memcpy(delay[k], delay[k]+nL, PS_MAX_DELAY*sizeof(delay[k][0]));
        memcpy(delay[k]+PS_MAX_DELAY, s[k], numQMFSlots*sizeof(delay[k][0]));
        for (n = n0; n < nL; n++) {
            //H = delay 1
            out[k][n][0] = transient_gain[k_to_i[k]][n] * delay[k][n+PS_MAX_DELAY-1][0];
            out[k][n][1] = transient_gain[k_to_i[k]][n] * delay[k][n+PS_MAX_DELAY-1][1];
        }
    }
}

static av_noinline void stereo_processing(PSContext *ps, float (*l)[32][2], float (*r)[32][2], int is34)
{
    int e, b, k, n;

    float (*H11)[PS_MAX_NUM_ENV+1][PS_MAX_NR_IIDICC] = ps->H11;
    float (*H12)[PS_MAX_NUM_ENV+1][PS_MAX_NR_IIDICC] = ps->H12;
    float (*H21)[PS_MAX_NUM_ENV+1][PS_MAX_NR_IIDICC] = ps->H21;
    float (*H22)[PS_MAX_NUM_ENV+1][PS_MAX_NR_IIDICC] = ps->H22;
    float (*opd_smooth)[2][2] = ps->opd_smooth;
    float (*ipd_smooth)[2][2] = ps->ipd_smooth;
    int8_t iid_mapped_buf[PS_MAX_NUM_ENV][PS_MAX_NR_IIDICC];
    int8_t icc_mapped_buf[PS_MAX_NUM_ENV][PS_MAX_NR_IIDICC];
    int8_t ipd_mapped_buf[PS_MAX_NUM_ENV][PS_MAX_NR_IPDOPD];
    int8_t opd_mapped_buf[PS_MAX_NUM_ENV][PS_MAX_NR_IPDOPD];
    int8_t (*iid_mapped)[PS_MAX_NR_IIDICC] = iid_mapped_buf;
    int8_t (*icc_mapped)[PS_MAX_NR_IIDICC] = icc_mapped_buf;
    int8_t (*ipd_mapped)[PS_MAX_NR_IPDOPD] = ipd_mapped_buf;
    int8_t (*opd_mapped)[PS_MAX_NR_IPDOPD] = opd_mapped_buf;
    const int8_t *k_to_i = is34 ? k_to_i_34 : k_to_i_20;
    const float (*H_LUT)[8][4] = (PS_BASELINE || ps->icc_mode < 3) ? HA : HB;
    static const float ipdopd_sin[] = { 0, M_SQRT1_2, 1,  M_SQRT1_2,  0, -M_SQRT1_2, -1, -M_SQRT1_2 };
    static const float ipdopd_cos[] = { 1, M_SQRT1_2, 0, -M_SQRT1_2, -1, -M_SQRT1_2,  0,  M_SQRT1_2 };

    //Remapping
    for (b = 0; b < PS_MAX_NR_IIDICC; b++) {
        H11[0][0][b] = H11[0][ps->num_env_old][b];
        H12[0][0][b] = H12[0][ps->num_env_old][b];
        H21[0][0][b] = H21[0][ps->num_env_old][b];
        H22[0][0][b] = H22[0][ps->num_env_old][b];
        H11[1][0][b] = H11[1][ps->num_env_old][b];
        H12[1][0][b] = H12[1][ps->num_env_old][b];
        H21[1][0][b] = H21[1][ps->num_env_old][b];
        H22[1][0][b] = H22[1][ps->num_env_old][b];
    }
    if (is34) {
        for (e = 0; e < ps->num_env; e++) {
            if (ps->nr_icc_par == 20)
                map_idx_20_to_34(icc_mapped[e], ps->icc_par[e], 1);
            else if (ps->nr_icc_par == 10) {
                map_idx_10_to_20(icc_mapped[e], ps->icc_par[e], 1);
                map_idx_20_to_34(icc_mapped[e], icc_mapped[e], 1);
            } else
                icc_mapped = ps->icc_par;
            if (ps->nr_iid_par == 20)
                map_idx_20_to_34(iid_mapped[e], ps->iid_par[e], 1);
            else if (ps->nr_iid_par == 10) {
                map_idx_10_to_20(iid_mapped[e], ps->iid_par[e], 1);
                map_idx_20_to_34(iid_mapped[e], iid_mapped[e], 1);
            } else
                iid_mapped = ps->iid_par;
            if (ps->enable_ipdopd) {
                if (ps->nr_ipdopd_par == 11) {
                    map_idx_20_to_34(ipd_mapped[e], ps->ipd_par[e], 0);
                    map_idx_20_to_34(opd_mapped[e], ps->opd_par[e], 0);
                } else if (ps->nr_ipdopd_par == 5) {
                    map_idx_10_to_20(ipd_mapped[e], ps->ipd_par[e], 0);
                    map_idx_20_to_34(ipd_mapped[e], ipd_mapped[e], 0);
                    map_idx_10_to_20(opd_mapped[e], ps->opd_par[e], 0);
                    map_idx_20_to_34(opd_mapped[e], opd_mapped[e], 0);
                } else {
                    ipd_mapped = ps->ipd_par;
                    opd_mapped = ps->opd_par;
                }
            }
        }
        if (!ps->is34bands_old) {
            map_val_20_to_34(H11[0], 0);
            map_val_20_to_34(H11[1], 0);
            map_val_20_to_34(H12[0], 0);
            map_val_20_to_34(H12[1], 0);
            map_val_20_to_34(H21[0], 0);
            map_val_20_to_34(H21[1], 0);
            map_val_20_to_34(H22[0], 0);
            map_val_20_to_34(H22[1], 0);
            ipdopd_reset(ps->ipd_smooth, ps->opd_smooth);
        }
    } else {
        for (e = 0; e < ps->num_env; e++) {
            if (ps->nr_icc_par == 34)
                map_idx_34_to_20(icc_mapped[e], ps->icc_par[e], 1);
            else if (ps->nr_icc_par == 10)
                map_idx_10_to_20(icc_mapped[e], ps->icc_par[e], 1);
            else
                icc_mapped = ps->icc_par;
            if (ps->nr_iid_par == 34)
                map_idx_34_to_20(iid_mapped[e], ps->iid_par[e], 1);
            else if (ps->nr_iid_par == 10)
                map_idx_10_to_20(iid_mapped[e], ps->iid_par[e], 1);
            else
                iid_mapped = ps->iid_par;
            if (ps->enable_ipdopd) {
                if (ps->nr_ipdopd_par == 17) {
                    map_idx_34_to_20(ipd_mapped[e], ps->ipd_par[e], 0);
                    map_idx_34_to_20(opd_mapped[e], ps->opd_par[e], 0);
                } else if (ps->nr_ipdopd_par == 5) {
                    map_idx_10_to_20(ipd_mapped[e], ps->ipd_par[e], 0);
                    map_idx_10_to_20(opd_mapped[e], ps->opd_par[e], 0);
                } else {
                    ipd_mapped = ps->ipd_par;
                    opd_mapped = ps->opd_par;
                }
            }
        }
        if (ps->is34bands_old) {
            map_val_34_to_20(H11[0], 0);
            map_val_34_to_20(H11[1], 0);
            map_val_34_to_20(H12[0], 0);
            map_val_34_to_20(H12[1], 0);
            map_val_34_to_20(H21[0], 0);
            map_val_34_to_20(H21[1], 0);
            map_val_34_to_20(H22[0], 0);
            map_val_34_to_20(H22[1], 0);
            ipdopd_reset(ps->ipd_smooth, ps->opd_smooth);
        }
    }

    //Mixing
    for (e = 0; e < ps->num_env; e++) {
        for (b = 0; b < NR_PAR_BANDS[is34]; b++) {
            float h11, h12, h21, h22;
            h11 = H_LUT[iid_mapped[e][b] + 7 + 23 * ps->iid_quant][icc_mapped[e][b]][0];
            h12 = H_LUT[iid_mapped[e][b] + 7 + 23 * ps->iid_quant][icc_mapped[e][b]][1];
            h21 = H_LUT[iid_mapped[e][b] + 7 + 23 * ps->iid_quant][icc_mapped[e][b]][2];
            h22 = H_LUT[iid_mapped[e][b] + 7 + 23 * ps->iid_quant][icc_mapped[e][b]][3];
            if (!PS_BASELINE && ps->enable_ipdopd && b < ps->nr_ipdopd_par) {
                //The spec say says to only run this smoother when enable_ipdopd
                //is set but the reference decoder appears to run it constantly
                float h11i, h12i, h21i, h22i;
                float opd_mag, ipd_mag, ipd_adj_re, ipd_adj_im;
                float opd_re = ipdopd_cos[ps->opd_par[e][b]];
                float opd_im = ipdopd_sin[ps->opd_par[e][b]];
                float ipd_re = ipdopd_cos[ps->ipd_par[e][b]];
                float ipd_im = ipdopd_sin[ps->ipd_par[e][b]];
                float opd_im_smooth = 0.25f * opd_smooth[b][0][1] + 0.5f * opd_smooth[b][1][1] + opd_im;
                float opd_re_smooth = 0.25f * opd_smooth[b][0][0] + 0.5f * opd_smooth[b][1][0] + opd_re;
                float ipd_im_smooth = 0.25f * ipd_smooth[b][0][1] + 0.5f * ipd_smooth[b][1][1] + ipd_im;
                float ipd_re_smooth = 0.25f * ipd_smooth[b][0][0] + 0.5f * ipd_smooth[b][1][0] + ipd_re;
                opd_smooth[b][0][0] = opd_smooth[b][1][0];
                opd_smooth[b][0][1] = opd_smooth[b][1][1];
                opd_smooth[b][1][0] = opd_re;
                opd_smooth[b][1][1] = opd_im;
                ipd_smooth[b][0][0] = ipd_smooth[b][1][0];
                ipd_smooth[b][0][1] = ipd_smooth[b][1][1];
                ipd_smooth[b][1][0] = ipd_re;
                ipd_smooth[b][1][1] = ipd_im;
                opd_mag = 1 / sqrt(opd_im_smooth * opd_im_smooth + opd_re_smooth * opd_re_smooth);
                ipd_mag = 1 / sqrt(ipd_im_smooth * ipd_im_smooth + ipd_re_smooth * ipd_re_smooth);
                opd_re = opd_re_smooth * opd_mag;
                opd_im = opd_im_smooth * opd_mag;
                ipd_re = ipd_re_smooth * ipd_mag;
                ipd_im = ipd_im_smooth * ipd_mag;
                ipd_adj_re = opd_re*ipd_re + opd_im*ipd_im;
                ipd_adj_im = opd_im*ipd_re - opd_re*ipd_im;
                h11i = h11 * opd_im;
                h11  = h11 * opd_re;
                h12i = h12 * ipd_adj_im;
                h12  = h12 * ipd_adj_re;
                h21i = h21 * opd_im;
                h21  = h21 * opd_re;
                h22i = h22 * ipd_adj_im;
                h22  = h22 * ipd_adj_re;
                H11[1][e+1][b] = h11i;
                H12[1][e+1][b] = h12i;
                H21[1][e+1][b] = h21i;
                H22[1][e+1][b] = h22i;
            }
            H11[0][e+1][b] = h11;
            H12[0][e+1][b] = h12;
            H21[0][e+1][b] = h21;
            H22[0][e+1][b] = h22;
        }
        for (k = 0; k < NR_BANDS[is34]; k++) {
            float h11r, h12r, h21r, h22r;
            float h11i, h12i, h21i, h22i;
            float h11r_step, h12r_step, h21r_step, h22r_step;
            float h11i_step, h12i_step, h21i_step, h22i_step;
            int start = ps->border_position[e];
            int stop  = ps->border_position[e+1];
            float width = 1.f / (stop - start);
            b = k_to_i[k];
            h11r = H11[0][e][b];
            h12r = H12[0][e][b];
            h21r = H21[0][e][b];
            h22r = H22[0][e][b];
            if (!PS_BASELINE && ps->enable_ipdopd) {
            //Is this necessary? ps_04_new seems unchanged
            if ((is34 && k <= 13 && k >= 9) || (!is34 && k <= 1)) {
                h11i = -H11[1][e][b];
                h12i = -H12[1][e][b];
                h21i = -H21[1][e][b];
                h22i = -H22[1][e][b];
            } else {
                h11i = H11[1][e][b];
                h12i = H12[1][e][b];
                h21i = H21[1][e][b];
                h22i = H22[1][e][b];
            }
            }
            //Interpolation
            h11r_step = (H11[0][e+1][b] - h11r) * width;
            h12r_step = (H12[0][e+1][b] - h12r) * width;
            h21r_step = (H21[0][e+1][b] - h21r) * width;
            h22r_step = (H22[0][e+1][b] - h22r) * width;
            if (!PS_BASELINE && ps->enable_ipdopd) {
                h11i_step = (H11[1][e+1][b] - h11i) * width;
                h12i_step = (H12[1][e+1][b] - h12i) * width;
                h21i_step = (H21[1][e+1][b] - h21i) * width;
                h22i_step = (H22[1][e+1][b] - h22i) * width;
            }
            for (n = start + 1; n <= stop; n++) {
                //l is s, r is d
                float l_re = l[k][n][0];
                float l_im = l[k][n][1];
                float r_re = r[k][n][0];
                float r_im = r[k][n][1];
                h11r += h11r_step;
                h12r += h12r_step;
                h21r += h21r_step;
                h22r += h22r_step;
                if (!PS_BASELINE && ps->enable_ipdopd) {
                h11i += h11i_step;
                h12i += h12i_step;
                h21i += h21i_step;
                h22i += h22i_step;
#if 0
            h11r = 0;
            h12r = 1;
            h21r = 1;
            h22r = 0;
            h11i = h12i = h21i = h22i = 0;
#endif
                l[k][n][0] = h11r*l_re + h21r*r_re - h11i*l_im - h21i*r_im;
                l[k][n][1] = h11r*l_im + h21r*r_im + h11i*l_re + h21i*r_re;
                r[k][n][0] = h12r*l_re + h22r*r_re - h12i*l_im - h22i*r_im;
                r[k][n][1] = h12r*l_im + h22r*r_im + h12i*l_re + h22i*r_re;
                } else {
                l[k][n][0] = h11r*l_re + h21r*r_re;
                l[k][n][1] = h11r*l_im + h21r*r_im;
                r[k][n][0] = h12r*l_re + h22r*r_re;
                r[k][n][1] = h12r*l_im + h22r*r_im;
                }
            }
        }
    }
}

static void transpose_in(float Ltrans[64][44][2], float L[2][38][64])
{
    int i, j;
    for (i = 0; i < 64; i++) {
        for (j = 0; j < 38; j++) {
            Ltrans[i][j+6][0] = L[0][j][i];
            Ltrans[i][j+6][1] = L[1][j][i];
        }
    }
}

static void transpose_out(float in[64][32][2], float out[2][38][64])
{
    int i, j;
    for (i = 0; i < 64; i++) {
        for (j = 0; j < 32; j++) {
            out[0][j][i] = in[i][j][0];
            out[1][j][i] = in[i][j][1];
        }
    }
}

int ff_ps_apply(AVCodecContext *avctx, PSContext *ps, float L[2][38][64], float R[2][38][64], int top)
{
    float Lout[64][32][2];
    float Rout[64][32][2];
    float Lbuf[91][32][2];
    float Rbuf[91][32][2];
    const int len = 32;
    int is34 = ps->is34bands;

    top += NR_BANDS[is34] - 64;
    memset(ps->delay+top, 0, (NR_BANDS[is34] - top)*sizeof(ps->delay[0]));
    if (top < NR_ALLPASS_BANDS[is34])
        memset(ps->ap_delay + top, 0, (NR_ALLPASS_BANDS[is34] - top)*sizeof(ps->ap_delay[0]));

    transpose_in(ps->in_buf, L);

    memset(Lbuf, -1, sizeof(Lbuf));
    memset(Lout, -1, sizeof(Lout));
    hybrid_analysis(Lbuf, ps->in_buf, is34, len);
#if 1
    decorrelation(ps, Rbuf, Lbuf, is34);
    stereo_processing(ps, Lbuf, Rbuf, is34);
#endif
    hybrid_synthesis(Lout, Lbuf, is34, len);
    hybrid_synthesis(Rout, Rbuf, is34, len);

    transpose_out(Lout, L);
    transpose_out(Rout, R);

    return 0;
}

static av_cold void ps_init_dec()
{
#if !CONFIG_HARDCODED_TABLES
    int k, m;
    //TODO store these as int8 and divide them during initialization
    static const int8_t f_center_20[] = {
        -3, -1, 1, 3, 5, 7, 10, 14, 18, 22,
    };
    static const int8_t f_center_34[] = {
         2,  6, 10, 14, 18, 22, 26, 30,
        34,-10, -6, -2, 51, 57, 15, 21,
        27, 33, 39, 45, 54, 66, 78, 42,
       102, 66, 78, 90,102,114,126, 90,
    };
    static const float fractional_delay_links[] = { 0.43f, 0.75f, 0.347f };
    const float fractional_delay_gain = 0.39f;
    static const float icc_invq[] = {
        1, 0.937,      0.84118,    0.60092,    0.36764,   0,      -0.589,    -1
    };
    static const float acos_icc_invq[] = {
        0, 0.35685527, 0.57133466, 0.92614472, 1.1943263, M_PI/2, 2.2006171, M_PI
    };
    int iid, icc;
    for (k = 0; k < NR_ALLPASS_BANDS20; k++) {
        double f_center, theta;
        if (k < FF_ARRAY_ELEMS(f_center_20))
            f_center = f_center_20[k] * 0.125;
        else
            f_center = k - 6.5f;
        for (m = 0; m < PS_AP_LINKS; m++) {
            theta = -M_PI * fractional_delay_links[m] * f_center;
            Q_fract_allpass[0][k][m][0] = cos(theta);
            Q_fract_allpass[0][k][m][1] = sin(theta);
        }
        theta = -M_PI*fractional_delay_gain*f_center;
        phi_fract[0][k][0] = cos(theta);
        phi_fract[0][k][1] = sin(theta);
    }
    for (k = 0; k < NR_ALLPASS_BANDS34; k++) {
        double f_center, theta;
        if (k < FF_ARRAY_ELEMS(f_center_34))
            f_center = f_center_34[k] / 24.;
        else
            f_center = k - 26.5f;
        for (m = 0; m < PS_AP_LINKS; m++) {
            theta = -M_PI * fractional_delay_links[m] * f_center;
            Q_fract_allpass[1][k][m][0] = cos(theta);
            Q_fract_allpass[1][k][m][1] = sin(theta);
        }
        theta = -M_PI*fractional_delay_gain*f_center;
        phi_fract[1][k][0] = cos(theta);
        phi_fract[1][k][1] = sin(theta);
    }

    make_filters_from_proto(f20_0_8,  g0_Q8,   8);
    make_filters_from_proto(f34_0_12, g0_Q12, 12);
    make_filters_from_proto(f34_1_8,  g1_Q8,   8);
    make_filters_from_proto(f34_2_4,  g2_Q4,   4);

    for (iid = 0; iid < 46; iid++) {
        float c = iid_par_dequant[iid]; //<Linear Inter-channel Intensity Difference
        float c1 = (float)M_SQRT2 / sqrtf(1.0f + c*c);
        float c2 = c * c1;
        //av_log(NULL, AV_LOG_ERROR, "    {\n");
        for (icc = 0; icc < 8; icc++) {
            /*if (PS_BASELINE || ps->icc_mode < 3)*/ {
                float alpha = 0.5f * acos_icc_invq[icc];
                float beta  = alpha * (c1 - c2) * (float)M_SQRT1_2;
                HA[iid][icc][0] = c2 * cosf(beta + alpha);
                HA[iid][icc][1] = c1 * cosf(beta - alpha);
                HA[iid][icc][2] = c2 * sinf(beta + alpha);
                HA[iid][icc][3] = c1 * sinf(beta - alpha);
                //av_log(NULL, AV_LOG_ERROR, "        { %13.10f, %13.10f, %13.10f, %13.10f  },\n", HA[iid][icc][0], HA[iid][icc][1], HA[iid][icc][2], HA[iid][icc][3]);
            } /* else */ {
                float alpha, gamma, mu, rho;
                float alpha_c, alpha_s, gamma_c, gamma_s;
                rho = FFMAX(icc_invq[icc], 0.05f);
                alpha = 0.5f * atan2f(2.0f * c * rho, c*c - 1.0f);
                mu = c + 1.0f / c;
                mu = sqrtf(1 + (4 * rho * rho - 4)/(mu * mu));
                gamma = atanf(sqrtf((1.0f - mu)/(1.0f + mu)));
                if (alpha < 0) alpha += M_PI/2;
                alpha_c = cosf(alpha);
                alpha_s = sinf(alpha);
                gamma_c = cosf(gamma);
                gamma_s = sinf(gamma);
                HB[iid][icc][0] =  M_SQRT2 * alpha_c * gamma_c;
                HB[iid][icc][1] =  M_SQRT2 * alpha_s * gamma_c;
                HB[iid][icc][2] = -M_SQRT2 * alpha_s * gamma_s;
                HB[iid][icc][3] =  M_SQRT2 * alpha_c * gamma_s;
                //av_log(NULL, AV_LOG_ERROR, "        { %13.10f, %13.10f, %13.10f, %13.10f  },\n", HB[iid][icc][0], HB[iid][icc][1], HB[iid][icc][2], HB[iid][icc][3]);
            }
        }
        //av_log(NULL, AV_LOG_ERROR, "    },\n");
    }
#endif
}

#define PS_INIT_VLC_STATIC(num, size) \
    INIT_VLC_STATIC(&vlc_ps[num], 9, ps_tmp[num].table_size / ps_tmp[num].elem_size,    \
                    ps_tmp[num].ps_bits, 1, 1,                                          \
                    ps_tmp[num].ps_codes, ps_tmp[num].elem_size, ps_tmp[num].elem_size, \
                    size);

#define PS_VLC_ROW(name) \
    { name ## _codes, name ## _bits, sizeof(name ## _codes), sizeof(name ## _codes[0]) }

av_cold void ff_ps_init(void) {
    // Syntax initialization
    static const struct {
        const void *ps_codes, *ps_bits;
        const unsigned int table_size, elem_size;
    } ps_tmp[] = {
        PS_VLC_ROW(huff_iid_df1),
        PS_VLC_ROW(huff_iid_dt1),
        PS_VLC_ROW(huff_iid_df0),
        PS_VLC_ROW(huff_iid_dt0),
        PS_VLC_ROW(huff_icc_df),
        PS_VLC_ROW(huff_icc_dt),
        PS_VLC_ROW(huff_ipd_df),
        PS_VLC_ROW(huff_ipd_dt),
        PS_VLC_ROW(huff_opd_df),
        PS_VLC_ROW(huff_opd_dt),
    };

    PS_INIT_VLC_STATIC(0, 1544);
    PS_INIT_VLC_STATIC(1,  832);
    PS_INIT_VLC_STATIC(2, 1024);
    PS_INIT_VLC_STATIC(3, 1036);
    PS_INIT_VLC_STATIC(4,  544);
    PS_INIT_VLC_STATIC(5,  544);
    PS_INIT_VLC_STATIC(6,  512);
    PS_INIT_VLC_STATIC(7,  512);
    PS_INIT_VLC_STATIC(8,  512);
    PS_INIT_VLC_STATIC(9,  512);

    ps_init_dec();
}

av_cold void ff_ps_ctx_init(PSContext *ps)
{
    ipdopd_reset(ps->ipd_smooth, ps->opd_smooth);
}
