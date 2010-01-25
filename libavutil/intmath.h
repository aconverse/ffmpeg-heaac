/*
 * Copyright (c) 2010 Mans Rullgard <mans@mansr.com>
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

#ifndef AVUTIL_INTMATH_H
#define AVUTIL_INTMATH_H

#include "config.h"
#include "common.h"

extern const uint32_t ff_inverse[257];

#if   ARCH_ARM
#   include "arm/intmath.h"
#elif ARCH_X86
#   include "x86/intmath.h"
#endif

#if HAVE_FAST_CLZ && AV_GCC_VERSION_AT_LEAST(3,4)

#ifndef av_log2

#define av_log2(x) (31 - __builtin_clz((x)|1))

#ifndef av_log2_16bit
#define av_log2_16bit av_log2
#endif

#endif /* av_log2 */

#endif /* AV_GCC_VERSION_AT_LEAST(3,4) */

#ifndef FASTDIV

#if CONFIG_FASTDIV
#    define FASTDIV(a,b)   ((uint32_t)((((uint64_t)a) * ff_inverse[b]) >> 32))
#else
#    define FASTDIV(a,b)   ((a) / (b))
#endif

#endif /* FASTDIV */

#endif /* AVUTIL_INTMATH_H */
