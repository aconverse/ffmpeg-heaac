/*
 * Generate a header file for hardcoded Parametric Stereo tables
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

#include <stdlib.h>
#define CONFIG_HARDCODED_TABLES 0
#include "ps_tablegen.h"
#include "tableprint.h"

int main(void)
{
    ps_tableinit();

    write_fileheader();

    printf("static const float pd_re_smooth[8*8*8] = {\n");
    write_float_array(pd_re_smooth, 8*8*8);
    printf("};\n");
    printf("static const float pd_im_smooth[8*8*8] = {\n");
    write_float_array(pd_im_smooth, 8*8*8);
    printf("};\n");

    return 0;
}
