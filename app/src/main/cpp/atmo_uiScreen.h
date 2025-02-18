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

#ifndef atmo_uiScreen_H
#define atmo_uiScreen_H

#include "libvkk/vkk_ui.h"

typedef struct atmo_uiWindowHud_s atmo_uiWindowHud_t;

typedef struct atmo_uiScreen_s
{
	vkk_uiScreen_t      base;
	atmo_uiWindowHud_t* window_hud;
} atmo_uiScreen_t;

atmo_uiScreen_t* atmo_uiScreen_new(vkk_engine_t* engine);
void             atmo_uiScreen_delete(atmo_uiScreen_t** _self);
void             atmo_uiScreen_draw(atmo_uiScreen_t* self);
int              atmo_uiScreen_event(atmo_uiScreen_t* self,
                                     vkk_platformEvent_t* event);

#endif
