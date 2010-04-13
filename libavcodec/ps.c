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

#define DEBUG
#if defined(__llvm__) || !defined(DEBUG)
#define NO_OPT
#else
#define NO_OPT __attribute__((optimize(0)))
#endif

#define PS_BASELINE 0

#define numQMFSlots 32 //numTimeSlots * RATE

static int8_t num_env_tab[2][4] = {
    { 0, 1, 2, 4, },
    { 1, 2, 3, 4, },
};

static int8_t nr_iidicc_par_tab[] = {
    10, 20, 34, 10, 20, 34,
};

static int8_t nr_iidopd_par_tab[] = {
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
#define PS_INIT_VLC_STATIC(num, size) \
    INIT_VLC_STATIC(&vlc_ps[num], 9, ps_tmp[num].table_size / ps_tmp[num].elem_size,    \
                    ps_tmp[num].ps_bits, 1, 1,                                          \
                    ps_tmp[num].ps_codes, ps_tmp[num].elem_size, ps_tmp[num].elem_size, \
                    size);

#define PS_VLC_ROW(name) \
    { name ## _codes, name ## _bits, sizeof(name ## _codes), sizeof(name ## _codes[0]) }

static int iid_data(GetBitContext *gb, PSContext *ps, int e)
{
    int b;
    int table_idx = huff_iid[2*ps->iid_dt[e]+ps->iid_quant];
    VLC_TYPE (*vlc_table)[2] = vlc_ps[table_idx].table;
    if (ps->iid_dt[e]) {
        int e_prev = e ? e - 1 : ps->num_env_old - 1;
        e_prev = FFMAX(e_prev, 0); //TODO FIXME does this make sense for ps->num_env_old = 0
        for (b = 0; b < ps->nr_iid_par; b++) {
            ps->iid_par[e][b] = ps->iid_par[e_prev][b] +
                                get_vlc2(gb, vlc_table, 9, 3) -
                                huff_offset[table_idx];
            if (FFABS(ps->iid_par[e][b]) > 7 + 8 * ps->iid_quant) {
                av_log(NULL, AV_LOG_ERROR, "illegal iid\n");
                abort();
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
                av_log(NULL, AV_LOG_ERROR, "illegal iid\n");
                abort();
                return -1;
            }
        }
    }
    return 0;
}

static int icc_data(GetBitContext *gb, PSContext *ps, int e)
{
    int b;
    int table_idx = ps->icc_dt[e] ? huff_icc_dt : huff_icc_df;
    VLC_TYPE (*vlc_table)[2] = vlc_ps[table_idx].table;
    if (ps->icc_dt[e]) {
        int e_prev = e ? e - 1 : ps->num_env_old - 1;
        e_prev = FFMAX(e_prev, 0); //TODO FIXME does this make sense for ps->num_env_old = 0
        for (b = 0; b < ps->nr_icc_par; b++) {
            ps->icc_par[e][b] = ps->icc_par[e_prev][b] + get_vlc2(gb, vlc_table, 9, 3) - huff_offset[table_idx];
            if (ps->icc_par[e][b] > 7U) {
                av_log(NULL, AV_LOG_ERROR, "illegal icc\n");
                abort();
                return -1;
            }
        }
    } else {
        int prev = 0;
        for (b = 0; b < ps->nr_icc_par; b++) {
            prev += get_vlc2(gb, vlc_table, 9, 3) - huff_offset[table_idx];
            ps->icc_par[e][b] = prev;
            if (ps->icc_par[e][b] > 7U) {
                av_log(NULL, AV_LOG_ERROR, "illegal icc\n");
                abort();
                return -1;
            }
        }
    }
    return 0;
}

static void ipd_data(GetBitContext *gb, PSContext *ps, int e)
{
    int b;
    int table_idx = ps->ipd_dt[e] ? huff_ipd_dt : huff_ipd_df;
    VLC_TYPE (*vlc_table)[2] = vlc_ps[table_idx].table;
    for (b = 0; b < ps->nr_ipdopd_par; b++)
        ps->ipd_par[e][b] = get_vlc2(gb, vlc_table, 9, 1);
}

static void opd_data(GetBitContext *gb, PSContext *ps, int e)
{
    int b;
    int table_idx = ps->opd_dt[e] ? huff_opd_dt : huff_opd_df;
    VLC_TYPE (*vlc_table)[2] = vlc_ps[table_idx].table;
    for (b = 0; b < ps->nr_ipdopd_par; b++)
        ps->opd_par[e][b] = get_vlc2(gb, vlc_table, 9, 1);
}

static int ps_extension(GetBitContext *gb, PSContext *ps, int ps_extension_id)
{
    int e;
    int count = get_bits_count(gb);
    if (!ps_extension_id) {
        ps->enable_ipdopd = get_bits1(gb);
        if (ps->enable_ipdopd) {
            for (e = 0; e < ps->num_env; e++) {
                ps->ipd_dt[e] = get_bits1(gb);
                ipd_data(gb, ps, e);
                ps->opd_dt[e] = get_bits1(gb);
                opd_data(gb, ps, e);
            }
        }
        skip_bits1(gb);      //reserved_ps
    }
    return get_bits_count(gb) - count;
}

int ff_ps_data(GetBitContext *gb, PSContext *ps)
{
    int e;
    int bit_count_start = get_bits_count(gb);
    int header;

av_log(NULL, AV_LOG_ERROR, "ps_data\n");
    header = get_bits1(gb);
    if (header) {     //enable_ps_header
        ps->enable_iid = get_bits1(gb);
        if (ps->enable_iid) {
            ps->iid_mode = get_bits(gb, 3);
            if (ps->iid_mode > 5) {
                av_log(NULL, AV_LOG_ERROR, "iid_mode %d is reserved.\n",
                       ps->iid_mode);
                return -1;
            }
            ps->nr_iid_par    = nr_iidicc_par_tab[ps->iid_mode];
            ps->iid_quant     = ps->iid_mode > 2;
            ps->nr_ipdopd_par = nr_iidopd_par_tab[ps->iid_mode];
        }
        ps->enable_icc = get_bits1(gb);
        if (ps->enable_icc) {
            ps->icc_mode = get_bits(gb, 3);
            if (ps->icc_mode > 5) {
                av_log(NULL, AV_LOG_ERROR, "icc_mode %d is reserved.\n",
                       ps->icc_mode);
                return -1;
            }
            ps->nr_icc_par = nr_iidicc_par_tab[ps->icc_mode];
        }
        ps->enable_ext = get_bits1(gb);
    }
av_log(NULL, AV_LOG_ERROR, "header %d iid %d %d icc %d %d\n", header, ps->enable_iid, ps->iid_mode, ps->enable_icc, ps->icc_mode);

    ps->frame_class = get_bits1(gb);
    ps->num_env_old = ps->num_env;
    ps->num_env     = num_env_tab[ps->frame_class][get_bits(gb, 2)];

    if (ps->frame_class)
        for (e = 0; e < ps->num_env; e++)
            ps->border_position[e] = get_bits(gb, 5);
    else
        for (e = 0; e < ps->num_env; e++)
            ps->border_position[e] = (e + 1) * numQMFSlots / ps->num_env - 1;

    if (ps->enable_iid)
        for (e = 0; e < ps->num_env; e++) {
            ps->iid_dt[e] = get_bits1(gb);
            iid_data(gb, ps, e);
        }
    else
        memset(ps->iid_par, 0, sizeof(ps->iid_par));

    if (ps->enable_icc)
        for (e = 0; e < ps->num_env; e++) {
            ps->icc_dt[e] = get_bits1(gb);
            icc_data(gb, ps, e);
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
            av_log(NULL, AV_LOG_ERROR, "ps extension overflow %d", cnt);
            abort();
        }
        skip_bits(gb, cnt);
    }

    //Baseline
    ps->enable_ipdopd &= !PS_BASELINE;

//av_log(NULL, AV_LOG_ERROR, "bits consumed %d\n", get_bits_count(gb) - bit_count_start);
    return get_bits_count(gb) - bit_count_start;
}

static const float g0_Q8[] = {
    0.00746082949812f, 0.02270420949825f, 0.04546865930473f, 0.07266113929591f,
    0.09885108575264f, 0.11793710567217f, 0.125f,            0.11793710567217f,
    0.09885108575264f, 0.07266113929591f, 0.04546865930473f, 0.02270420949825f,
    0.00746082949812f
};

static const float g1_Q2[] = {
    0.0f,  0.01899487526049f, 0.0f, -0.07293139167538f,
    0.0f,  0.30596630545168f, 0.5f,  0.30596630545168f,
    0.0f, -0.07293139167538f, 0.0f,  0.01899487526049f,
    0.0f
};

static const float g0_Q12[] = {
    0.04081179924692f, 0.03812810994926f, 0.05144908135699f, 0.06399831151592f,
    0.07428313801106f, 0.08100347892914f, 0.08333333333333f, 0.08100347892914f,
    0.07428313801106f, 0.06399831151592f, 0.05144908135699f, 0.03812810994926f,
    0.04081179924692f
};

static const float g1_Q8[] = {
    0.01565675600122f, 0.03752716391991f, 0.05417891378782f, 0.08417044116767f,
    0.10307344158036f, 0.12222452249753f, 0.125f,            0.12222452249753f,
    0.10307344158036f, 0.08417044116767f, 0.05417891378782f, 0.03752716391991f,
    0.01565675600122f
};

static const float g2_Q4[] = {
    -0.05908211155639f, -0.04871498374946f, 0.0f,   0.07778723915851f,
     0.16486303567403f,  0.23279856662996f, 0.25f,  0.23279856662996f,
     0.16486303567403f,  0.07778723915851f, 0.0f,  -0.04871498374946f,
    -0.05908211155639f
};

static float f20_0_8[8][13][2];

static void make_filters_from_proto(float (*filter)[13][2], const float *proto, int bands)
{
    int q, n;
    for (q = 0; q < bands; q++) {
        for (n = 0; n < 13; n++) {
            float theta = 2 * M_PI * (q + 0.5) * (n - 6) / bands;
            filter[q][n][0] = proto[n] *  cosf(theta);
            filter[q][n][1] = proto[n] * -sinf(theta); //FIXME specbug? convolution?
        }
    }
}

/** Split one subband into 2 subsubbands with a real filter */
static void hybrid2_re(float (*in)[2], float (*out)[32][2], const float filter[13], int len, int reverse)
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
static void NO_OPT hybrid6_cx(float (*in)[2], float (*out)[32][2], const float (*filter)[13][2], int len)
{
    int i, j, ssb;
    int N = 8;
    float temp[8][2];

    for (i = 0; i < len; i++) {
        for (ssb = 0; ssb < N; ssb++) {
            //FIXME filter is conjugate symmetric
            //filter[6] is real
            float sum_re = 0.0f, sum_im = 0.0f;
            for (j = 0; j < 13; j++) {
                float in_re = in[i+j][0];
                float in_im = in[i+j][1];
                sum_re += filter[ssb][j][0] * in_re - filter[ssb][j][1] * in_im;
                sum_im += filter[ssb][j][0] * in_im + filter[ssb][j][1] * in_re;
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

static void NO_OPT hybrid_analysis(float out[91][32][2], float in[64][44][2], int is34, int len)
{
    int i;
    if(is34) {
        //XXX TODO
        av_log(NULL, AV_LOG_ERROR, "hybrid34!\n");
        abort();
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

static void hybrid_synthesis(float out[64][32][2], float in[91][32][2], int is34, int len)
{
    int i, n;
    if(is34) {
        //XXX TODO
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
/// Number of filter links for the all-pass filter
#define NR_ALLPASS_LINKS 3
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

static float Q_fract_allpass[2][NR_ALLPASS_BANDS34][NR_ALLPASS_LINKS][2];
static float phi_fract[2][NR_ALLPASS_BANDS34][2];

/// Inverse map i = b(k)

static int av_const map_k_to_i(int k, int is34)
{
    if (is34) {
        //TODO FIXME Table 8.49
        return k;
    } else {
        //Table 8.48
        if (k <= 1) {
            return 1 - k;
        } else if (k <= 16) {
            return k - 2;
        } else if (k <= 17) {
            return 14;
        } else if (k <= 20) {
            return 15;
        } else if (k <= 24) {
            return 16;
        } else if (k <= 29) {
            return 17;
        } else if (k <= 41) {
            return 18;
        } else {
            return 19;
        }
    }
}

/** Table 8.46 */
static void map_10_to_20(int8_t par[PS_MAX_NUM_ENV][PS_MAX_NR_IIDICC], int e)
{
    int b;
    for (b = 9; b >= 0; b--) {
        par[e][2*b+1] = par[e][2*b] = par[e][b];
    }
}
/** Table 8.46 */
static void map_34_to_20(int8_t par[PS_MAX_NUM_ENV][PS_MAX_NR_IIDICC], int e)
{
    par[e][ 0] = (2*par[e][ 0] +   par[e][ 1]) / 3;
    par[e][ 1] = (  par[e][ 1] + 2*par[e][ 2]) / 3;
    par[e][ 2] = (2*par[e][ 3] +   par[e][ 4]) / 3;
    par[e][ 3] = (  par[e][ 4] + 2*par[e][ 5]) / 3;
    par[e][ 4] = (  par[e][ 6] +   par[e][ 7]) / 2;
    par[e][ 5] = (  par[e][ 8] +   par[e][ 9]) / 2;
    par[e][ 6] =    par[e][10];
    par[e][ 7] =    par[e][11];
    par[e][ 8] = (  par[e][12] +   par[e][13]) / 2;
    par[e][ 9] = (  par[e][14] +   par[e][15]) / 2;
    par[e][10] =    par[e][16];
    par[e][11] =    par[e][17];
    par[e][12] =    par[e][18];
    par[e][13] =    par[e][19];
    par[e][14] = (  par[e][20] +   par[e][21]) / 2;
    par[e][15] = (  par[e][22] +   par[e][23]) / 2;
    par[e][16] = (  par[e][24] +   par[e][25]) / 2;
    par[e][17] = (  par[e][26] +   par[e][27]) / 2;
    par[e][18] = (  par[e][28] +   par[e][29] +   par[e][30] +   par[e][31]) / 4;
    par[e][19] = (  par[e][32] +   par[e][33]) / 2;
}

static void NO_OPT decorrelation(float (*out)[32][2], const float (*s)[32][2], int is34)
{
    static float power[34][32]; //[f][t]
    static float peak_decay_nrg[34][32];
    float transient_gain[34][32];
    const float peak_decay_factor = 0.76592833836465f;
    const float transient_impact  = 1.5f;
    const float a_smooth          = 0.25f; //< Smoothing coefficient
    int i, k, m, n;
    int n0 = 0, nL = 32; //ps->border_position[ps->num_env - 1]; //FIXME
    memset(power, 0, sizeof(power));
    for (n = n0; n < nL; n++) {
        for (k = 0; k < NR_BANDS[is34]; k++) { //TODO be careful about uninitialized bands in the 10,20 case
            int i = map_k_to_i(k, is34);
            power[i][n] += s[k][n][0] * s[k][n][0] + s[k][n][1] * s[k][n][1];
            if (!isfinite(power[i][n])) {
                av_log(NULL, AV_LOG_ERROR, "%d %d %d\n", i, k, n);
                av_log(NULL, AV_LOG_ERROR, "%f\n", power[i][n]);
                av_log(NULL, AV_LOG_ERROR, "%f %f\n", s[k][n][0], s[k][n][1]);
                abort();
            }
        }
    }

    //Transient detection
    for (i = 0; i < NR_PAR_BANDS[is34]; i++) {
        float decayed_peak = peak_decay_factor * peak_decay_nrg[i][nL - 1];
        peak_decay_nrg[i][n0] = (decayed_peak < power[i][n0]) ? power[i][n0] : decayed_peak;
        for (n = n0 + 1; n < nL; n++) {
            decayed_peak = peak_decay_factor * peak_decay_nrg[i][n - 1];
            peak_decay_nrg[i][n] = (decayed_peak < power[i][n]) ? power[i][n] : decayed_peak;
        }
    }
    for (i = 0; i < NR_PAR_BANDS[is34]; i++) {
        float power_smooth = 0.0f;
        float peak_decay_diff_smooth = 0.0f;
        for (n = n0; n < nL; n++) {
            power_smooth = a_smooth * power[i][n] + (1.0f - a_smooth) * power_smooth;
            if (!isfinite(power_smooth)) {
                av_log(NULL, AV_LOG_ERROR, "%d %d\n", i, n);
                av_log(NULL, AV_LOG_ERROR, "%f\n", power[i][n]);
                av_log(NULL, AV_LOG_ERROR, "%f\n", power_smooth);
                abort();
            }
            peak_decay_diff_smooth = a_smooth * (peak_decay_nrg[i][n] - power[i][n]) +
                                         (1.0f - a_smooth) * peak_decay_diff_smooth;
            transient_gain[i][n]   = (transient_impact * peak_decay_diff_smooth > power_smooth) ?
                                         power_smooth / (transient_impact * peak_decay_diff_smooth) : 1.0f;
//av_log(NULL, AV_LOG_ERROR, "transient_gain[%2d][%2d] %f\n", i, n, transient_gain[i][n]);
        }
    }

    //Decorrelation and transient reduction
    //g_decay_slope[k] = clip(1 - DECAY_SLOPE * (k - DECAY_CUTOFF), 0, 1)
    //                         NR_ALLPASS_LINKS - 1
    //                               -----
    //                                | |  Q_fract_allpass[k][m]*z^-link_delay[m] - a[m]*g_decay_slope[k]
    //H[k][z] = z^-2 * phi_fract[k] * | | ----------------------------------------------------------------
    //                                | | 1 - a[m]*g_decay_slope[k]*Q_fract_allpass[k][m]*z^-link_delay[m]
    //                               m = 0
    static const int link_delay[] = { 3, 4, 5 };
    static const float a[] = { 0.65143905753106f,
                               0.56471812200776f,
                               0.48954165955695f };
#define MAX_DELAY 14
    static float delay              [91 /*NR_BANDS[is34]*/][numQMFSlots+MAX_DELAY][2];
    static float all_pass_delay_buff[50 /*NR_ALLPASS_BANDS[is34]*/][ NR_ALLPASS_LINKS + 1][numQMFSlots+        5][2];
    //d[k][z] (out) = transient_gain_mapped[k][z] * H[k][z] * s[k][z]
    for (k = 0; k < NR_ALLPASS_BANDS[is34]; k++) {
        int b = map_k_to_i(k, is34);
        float g_decay_slope = 1.f - DECAY_SLOPE * (k - DECAY_CUTOFF[is34]);
        g_decay_slope = FFMIN(g_decay_slope, 1.f);
        g_decay_slope = FFMAX(g_decay_slope, 0.f);
        memcpy(delay[k], delay[k]+nL, MAX_DELAY*sizeof(delay[k][0]));
        memcpy(delay[k]+MAX_DELAY, s[k], numQMFSlots*sizeof(delay[k][0]));
        for (m = 0; m <= NR_ALLPASS_LINKS; m++) {
            memcpy(all_pass_delay_buff[k][m],   all_pass_delay_buff[k][m]+numQMFSlots,           5*sizeof(all_pass_delay_buff[k][m][0]));
            //memcpy(all_pass_delay_buff[k][m]+5, s[k],                                  numQMFSlots*sizeof(all_pass_delay_buff[k][m][0]));
        }
        for (n = n0; n < nL; n++) {
            float in_re = delay[k][n+MAX_DELAY-2][0] * phi_fract[is34][k][0] -
                          delay[k][n+MAX_DELAY-2][1] * phi_fract[is34][k][1];
            float in_im = delay[k][n+MAX_DELAY-2][0] * phi_fract[is34][k][1] +
                          delay[k][n+MAX_DELAY-2][1] * phi_fract[is34][k][0];
            for (m = 0; m < NR_ALLPASS_LINKS; m++) {
                float a_re                = a[m] * g_decay_slope * in_re;
                float a_im                = a[m] * g_decay_slope * in_im;
                float in_link_delay_re    = all_pass_delay_buff[k][m][n+5-link_delay[m]][0];
                float in_link_delay_im    = all_pass_delay_buff[k][m][n+5-link_delay[m]][1];
                float out_link_delay_re   = a[m] * g_decay_slope * all_pass_delay_buff[k][m+1][n+5-link_delay[m]][0];
                float out_link_delay_im   = a[m] * g_decay_slope * all_pass_delay_buff[k][m+1][n+5-link_delay[m]][1];
                float fractional_delay_re = Q_fract_allpass[is34][k][m][0];
                float fractional_delay_im = Q_fract_allpass[is34][k][m][1];
//av_log(NULL, AV_LOG_ERROR, "allpass stage %d in= %e %e\n", m, in_re, in_im);
                all_pass_delay_buff[k][m][n+5][0] = in_re;
                all_pass_delay_buff[k][m][n+5][1] = in_im;
                in_re  =  in_link_delay_re * fractional_delay_re -  in_link_delay_im * fractional_delay_im - a_re;
                in_im  =  in_link_delay_re * fractional_delay_im +  in_link_delay_im * fractional_delay_re - a_im;
                in_re += out_link_delay_re * fractional_delay_re - out_link_delay_im * fractional_delay_im;
                in_im += out_link_delay_re * fractional_delay_im + out_link_delay_im * fractional_delay_re;
            }
            all_pass_delay_buff[k][m][n+5][0] = in_re;
            all_pass_delay_buff[k][m][n+5][1] = in_im;
//av_log(NULL, AV_LOG_ERROR, "allpass[k=%2d][n=%2d] = %e %e ", k, n, in_re, in_im);
            out[k][n][0] = transient_gain[b][n] * in_re;
            out[k][n][1] = transient_gain[b][n] * in_im;
//av_log(NULL, AV_LOG_ERROR, "b %2d tr %f ", b, transient_gain[b][n]);
//av_log(NULL, AV_LOG_ERROR, "out = %e %e\n", out[k][n][0], out[k][n][1]);
        }
    }
#if 1
    for (; k < SHORT_DELAY_BAND[is34]; k++) {
        memcpy(delay[k], delay[k]+nL, MAX_DELAY*sizeof(delay[k][0]));
        memcpy(delay[k]+MAX_DELAY, s[k], numQMFSlots*sizeof(delay[k][0]));
        for (n = n0; n < nL; n++) {
            //H = delay 14
            out[k][n][0] = transient_gain[map_k_to_i(k, is34)][n] * delay[k][n+MAX_DELAY-14][0];
            out[k][n][1] = transient_gain[map_k_to_i(k, is34)][n] * delay[k][n+MAX_DELAY-14][1];
        }
    }
    for (; k < NR_BANDS[is34]; k++) {
        memcpy(delay[k], delay[k]+nL, MAX_DELAY*sizeof(delay[k][0]));
        memcpy(delay[k]+MAX_DELAY, s[k], numQMFSlots*sizeof(delay[k][0]));
        for (n = n0; n < nL; n++) {
            //H = delay 1
            out[k][n][0] = transient_gain[map_k_to_i(k, is34)][n] * delay[k][n+MAX_DELAY-1][0];
            out[k][n][1] = transient_gain[map_k_to_i(k, is34)][n] * delay[k][n+MAX_DELAY-1][1];
        }
    }
#else
memset(out+k, 0, (NR_BANDS[is34]-NR_ALLPASS_BANDS[is34])*sizeof(out[0]));
#endif
}

static void stereo_processing(PSContext *ps, float (*l)[32][2], float (*r)[32][2], int is34)
{
    int e, b, k, n;

    static float H11[PS_MAX_NR_IIDICC][numQMFSlots]; //TODO
    static float H12[PS_MAX_NR_IIDICC][numQMFSlots]; //make me a context var, or atleast my first row
    static float H21[PS_MAX_NR_IIDICC][numQMFSlots]; //deal 20-34 changes
    static float H22[PS_MAX_NR_IIDICC][numQMFSlots];
      //Table 8.28, Quantization grid for ICC
    static const float icc_invq[] = {
        1, 0.937,      0.84118,    0.60092,    0.36764,   0,      -0.589,    -1
    };
    static const float acos_icc_invq[] = {
        0, 0.35685527, 0.57133466, 0.92614472, 1.1943263, M_PI/2, 2.2006171, M_PI
    };

    int ne_prev = 0; //TODO value from the last frame

    for (b = 0; b < PS_MAX_NR_IIDICC; b++) {
        H11[b][0] = H11[b][63];
        H12[b][0] = H12[b][63];
        H21[b][0] = H21[b][63];
        H22[b][0] = H22[b][63];
    }
    //mixing
    //av_log(NULL, AV_LOG_ERROR, "num_env %d\n", ps->num_env);
    //av_log(NULL, AV_LOG_ERROR, "nr_iid_par %d\n", ps->nr_iid_par);
    //av_log(NULL, AV_LOG_ERROR, "nr_icc_par %d\n", ps->nr_icc_par);
    for (e = 0; e < ps->num_env; e++) {
        int ne = ps->border_position[e]; //TODO Spec says n[e+1] but that seems very dubious
av_log(NULL, AV_LOG_ERROR, "e %d border %d\n", e, ne);
        if (ps->nr_icc_par == 34 && !is34)
            map_34_to_20(ps->icc_par, e);
        else if (ps->nr_icc_par == 10 && !is34)
            map_10_to_20(ps->icc_par, e);
        if (ps->nr_iid_par == 34 && !is34)
            map_34_to_20(ps->iid_par, e);
        else if (ps->nr_iid_par == 10 && !is34)
            map_10_to_20(ps->iid_par, e);

        for (b = 0; b < NR_PAR_BANDS[is34]; b++) {
            float c = iid_par_dequant[ps->iid_quant][ps->iid_par[e][b] + 7 + 8 * ps->iid_quant]; //<Linear Inter-channel Intensity Difference
            float h11, h12, h21, h22;
            if (PS_BASELINE || ps->icc_mode < 3) {
                float c1 = (float)M_SQRT2 / sqrtf(1.0f + c*c);
                float c2 = c * c1;
                float alpha = 0.5f * acos_icc_invq[ps->icc_par[e][b]];
                float beta  = alpha * (c1 - c2) * (float)M_SQRT1_2;
                //av_log(NULL, AV_LOG_ERROR, "alpha %f beta %f c1 %f c2 %f\n", alpha, beta, c1, c2);
                h11 = c2 * cosf(beta + alpha);
                h12 = c1 * cosf(beta - alpha);
                h21 = c2 * sinf(beta + alpha);
                h22 = c1 * sinf(beta - alpha);
            } else {
                float rho = FFMAX(icc_invq[ps->icc_par[e][b]], 0.05f);
                float alpha = 0.5f * atan2f(2.0f * c * rho, c*c - 1.0f);
                float mu = c + 1.0f / c;
                mu = sqrtf(1 + (4 * rho * rho - 4)/(mu * mu));
                float gamma = atanf(sqrtf((1.0f - mu)/(1.0f + mu)));
                if (alpha < 0) alpha += M_PI/2;
                //av_log(NULL, AV_LOG_ERROR, "alpha %f gamma %f\n", alpha, gamma);
                h11 =  M_SQRT2 * cosf(alpha) * cosf(gamma);
                h12 =  M_SQRT2 * sinf(alpha) * cosf(gamma);
                h21 = -M_SQRT2 * sinf(alpha) * sinf(gamma);
                h22 =  M_SQRT2 * cosf(alpha) * sinf(gamma);
            }
            if (ps->enable_ipdopd) {
                //TODO FIXME
                abort();
            } else {
                //Interpolation
                float h11_step = (h11 - H11[b][ne_prev]) / (ne - ne_prev);
                float h12_step = (h12 - H12[b][ne_prev]) / (ne - ne_prev);
                float h21_step = (h21 - H21[b][ne_prev]) / (ne - ne_prev);
                float h22_step = (h22 - H22[b][ne_prev]) / (ne - ne_prev);
                for (n = ne_prev + 1; n < ne; n++) { //TODO optimize out the multiply with an iterative add
                    H11[b][n] = H11[b][ne_prev] + (n-ne_prev) * h11_step;
                    H12[b][n] = H12[b][ne_prev] + (n-ne_prev) * h12_step;
                    H21[b][n] = H21[b][ne_prev] + (n-ne_prev) * h21_step;
                    H22[b][n] = H22[b][ne_prev] + (n-ne_prev) * h22_step;
                }
                H11[b][ne] = h11;
                H12[b][ne] = h12;
                H21[b][ne] = h21;
                H22[b][ne] = h22;
            }
        }
        ne_prev = ne;
    }
    //fill out the rest of the frame
    for (b = 0; b < NR_PAR_BANDS[is34]; b++) {
        float h11 = H11[b][ne_prev];
        float h12 = H12[b][ne_prev];
        float h21 = H21[b][ne_prev];
        float h22 = H22[b][ne_prev];
        for (n = ne_prev + 1; n < numQMFSlots; n++) {
            H11[b][n] = h11;
            H12[b][n] = h12;
            H21[b][n] = h21;
            H22[b][n] = h22;
        }
    }
    for (n = 0; n < numQMFSlots; n++) {
        for (k = 0; k < NR_BANDS[is34]; k++) {
            //l is s, r is d
            float l_re = l[k][n][0];
            float l_im = l[k][n][1];
            float r_re = r[k][n][0];
            float r_im = r[k][n][1];
            float h11, h12, h21, h22;
            b = map_k_to_i(k, is34);
#if 0
            h11 = 0;
            h12 = 1;
            h21 = 1;
            h22 = 0;
#else
            h11 = H11[b][n];
            h12 = H12[b][n];
            h21 = H21[b][n];
            h22 = H22[b][n];
#endif
            l[k][n][0] = h11*l_re + h21*r_re;
            l[k][n][1] = h11*l_im + h21*r_im;
            r[k][n][0] = h12*l_re + h22*r_re;
            r[k][n][1] = h12*l_im + h22*r_im;
            //if (n==31) av_log(NULL, AV_LOG_ERROR, "ssb %2d parameter %2d h %f %f %f %f\n",
            //       k, b, h11, h21, h12, h22);
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

int NO_OPT ff_ps_apply(AVCodecContext *avctx, PSContext *ps, float L[2][38][64], float R[2][38][64])
{
   float Lout[64][32][2];
   float Rout[64][32][2];
   float Lbuf[91][32][2];
   float Rbuf[91][32][2];
   int is34 = !PS_BASELINE && (ps->nr_icc_par == 34 || ps->nr_iid_par == 34);
   const int len = 32;

av_log(NULL, AV_LOG_ERROR, "is34 %d\n", is34);
   transpose_in(ps->in_buf, L);

   memset(Lbuf, -1, sizeof(Lbuf));
   memset(Lout, -1, sizeof(Lout));
   hybrid_analysis(Lbuf, ps->in_buf, is34, len);
#if 1
   decorrelation(Rbuf, Lbuf, is34);
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
    int k, m;
    //TODO store these as int8 and divide them during initialization
    static const float f_center_20[] = {
        -3./8, -1./8, 1./8, 3./8, 5./8, 7./8, 5./4, 7./4, 9./4, 11./4,
    };
    static const float f_center_34[] = {
         1./12,  3./12,  5./12,  7./12,  9./12, 11./12, 13./12, 15./12,
        17./12, -5./12, -3./12, -1./12, 17./ 8, 19./8,   5./8,   7./8,
         9./8,  11./8,  13./8,  15./8,   9./4,  11./4,  13./4,   7./4,
        17./4,  11./4,  13./4,  15./4,  17./4,  19./4,  21./4,  15./4,
    };
    static const float fractional_delay_links[] = { 0.43f, 0.75f, 0.347f };
    const float fractional_delay_gain = 0.39f;
    for (k = 0; k < NR_ALLPASS_BANDS20; k++) {
        float f_center, theta;
        if (k < FF_ARRAY_ELEMS(f_center_20))
            f_center = f_center_20[k];
        else
            f_center = k - 7 + 0.5f;
        for (m = 0; m < NR_ALLPASS_LINKS; m++) {
            theta = -M_PI * fractional_delay_links[m] * f_center;
            Q_fract_allpass[0][k][m][0] = cosf(theta);
            Q_fract_allpass[0][k][m][1] = sinf(theta);
        }
        theta = -M_PI*fractional_delay_gain*f_center;
        phi_fract[0][k][0] = cosf(theta);
        phi_fract[0][k][1] = sinf(theta);
    }
    for (k = 0; k < NR_ALLPASS_BANDS34; k++) {
        float f_center, theta;
        if (k < FF_ARRAY_ELEMS(f_center_34))
            f_center = f_center_34[k];
        else
            f_center = k - 27 + 0.5f;
        for (m = 0; m < NR_ALLPASS_LINKS; m++) {
            Q_fract_allpass[1][k][m][0] = cosf(-M_PI * fractional_delay_links[m] * f_center);
            Q_fract_allpass[1][k][m][1] = sinf(-M_PI * fractional_delay_links[m] * f_center);
        }
        theta = -M_PI*fractional_delay_gain*f_center;
        phi_fract[1][k][0] = cosf(theta);
        phi_fract[1][k][1] = sinf(theta);
    }
    make_filters_from_proto(f20_0_8, g0_Q8, 8);
}

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
