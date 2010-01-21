/*
 * AAC Spectral Band Replication function delcarations
 * Copyright (c) 2008-2009 Robert Swain ( rob opendot cl )
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

/**
 * @file libavcodec/aacsbr_internal.h
 * AAC Spectral Band Replication function declarations
 * @author Robert Swain ( rob opendot cl )
 */

#ifndef AVCODEC_AACSBR_INTERNAL_H
#define AVCODEC_AACSBR_INTERNAL_H

#include "get_bits.h"
#include "aac.h"
#include "aacsbr.h"

av_cold void ff_aac_sbr_init(void);
av_cold void ff_aac_sbr_ctx_init(SpectralBandReplication *sbr);
int ff_decode_sbr_extension(AACContext *ac, SpectralBandReplication *sbr,
                            GetBitContext *gb, int crc, int cnt, int id_aac);
void ff_sbr_dequant(AACContext *ac, SpectralBandReplication *sbr, int id_aac);
void ff_sbr_apply(AACContext *ac, SpectralBandReplication *sbr, int ch, float* in, float* out);

#endif /* AVCODEC_AACSBR_INTERNAL_H */
