/*
 * Copyright (c) 2025 Jeff Boody
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 */

#ifndef atmo_uiInfoPanel_H
#define atmo_uiInfoPanel_H

#include "libvkk/vkk_ui.h"

typedef struct atmo_uiScreen_s atmo_uiScreen_t;

typedef struct atmo_uiInfoPanel_s
{
	vkk_uiInfoPanel_t  base;
	vkk_uiText_t*      heading_ctrl;
	vkk_uiText_t*      text_ctrl_h;
	vkk_uiText_t*      text_ctrl_phi;
	vkk_uiText_t*      text_ctrl_delta;
	vkk_uiText_t*      text_ctrl_omega;
	vkk_uiText_t*      text_ctrl_k;
	float              last_ctrl_h;
	float              last_ctrl_phi;
	float              last_ctrl_delta;
	float              last_ctrl_omega;
	float              last_ctrl_k;
} atmo_uiInfoPanel_t;

atmo_uiInfoPanel_t* atmo_uiInfoPanel_new(atmo_uiScreen_t* screen);
void                atmo_uiInfoPanel_delete(atmo_uiInfoPanel_t** _self);

#endif
