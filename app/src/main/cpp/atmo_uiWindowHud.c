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

#include <stdlib.h>

#define LOG_TAG "atmo"
#include "libcc/cc_log.h"
#include "atmo_renderer.h"
#include "atmo_solver.h"
#include "atmo_uiInfoPanel.h"
#include "atmo_uiScreen.h"
#include "atmo_uiWindowHud.h"

/***********************************************************
* private                                                  *
***********************************************************/

static void
atmo_uiWindowHud_rendererDraw(vkk_uiWidget_t* widget)
{
	ASSERT(widget);

	// widget is graphics_box

	atmo_uiScreen_t* screen;
	screen = (atmo_uiScreen_t*) widget->screen;

	atmo_uiWindowHud_t* window_hud = screen->window_hud;

	cc_rect1f_t* rect_draw = vkk_uiWidget_rectDraw(widget);

	atmo_renderer_draw(window_hud->renderer,
	                   window_hud->solver,
	                   (float) rect_draw->w,
	                   (float) rect_draw->h);
}

/***********************************************************
* public                                                   *
***********************************************************/

atmo_uiWindowHud_t*
atmo_uiWindowHud_new(atmo_uiScreen_t* screen)
{
	ASSERT(overlay);

	vkk_engine_t* engine = screen->base.engine;

	vkk_uiWindowFn_t wfn =
	{
		.priv = NULL,
	};

	uint32_t flags = VKK_UI_WINDOW_FLAG_LAYER0 |
	                 VKK_UI_WINDOW_FLAG_LAYER1 |
	                 VKK_UI_WINDOW_FLAG_TRANSPARENT;

	atmo_uiWindowHud_t* self;
	self = (atmo_uiWindowHud_t*)
	       vkk_uiWindow_new(&screen->base,
	                        sizeof(atmo_uiWindowHud_t),
	                        &wfn, flags);
	if(self == NULL)
	{
		return NULL;
	}

	vkk_uiGraphicsBoxFn_t gbfn =
	{
		.draw_fn = atmo_uiWindowHud_rendererDraw,
	};

	vkk_uiWidgetLayout_t gb_layout =
	{
		.wrapx    = VKK_UI_WIDGET_WRAP_STRETCH_PARENT,
		.wrapy    = VKK_UI_WIDGET_WRAP_STRETCH_PARENT,
		.stretchx = 1.0f,
		.stretchy = 1.0f,
	};

	cc_vec4f_t clear =
	{
		.a = 0.0f
	};

	self->graphics_box = vkk_uiGraphicsBox_new(&screen->base, 0,
	                                           &gbfn,
	                                           &gb_layout,
	                                           1, &clear);
	if(self->graphics_box == NULL)
	{
		goto failure;
	}

	self->info_panel = atmo_uiInfoPanel_new(screen);
	if(self->info_panel == NULL)
	{
		goto failure;
	}

	self->renderer = atmo_renderer_new(engine);
	if(self->renderer == NULL)
	{
		goto failure;
	}

	self->solver = atmo_solver_new(engine);
	if(self->solver == NULL)
	{
		goto failure;
	}

	vkk_uiWindow_t* window = (vkk_uiWindow_t*) self;
	vkk_uiLayer_t*  layer0 = vkk_uiWindow_layer0(window);
	vkk_uiLayer_t*  layer1 = vkk_uiWindow_layer1(window);
	vkk_uiLayer_add(layer0, (vkk_uiWidget_t*) self->graphics_box);
	vkk_uiLayer_add(layer1, (vkk_uiWidget_t*) self->info_panel);

	// success
	return self;

	// failure
	failure:
		atmo_uiWindowHud_delete(&self);
	return NULL;
}

void atmo_uiWindowHud_delete(atmo_uiWindowHud_t** _self)
{
	ASSERT(_self);

	atmo_uiWindowHud_t* self = *_self;
	if(self)
	{
		vkk_uiWindow_t* window = (vkk_uiWindow_t*) self;
		vkk_uiLayer_t*  layer0 = vkk_uiWindow_layer0(window);
		vkk_uiLayer_t*  layer1 = vkk_uiWindow_layer1(window);
		vkk_uiLayer_clear(layer0);
		vkk_uiLayer_clear(layer1);
		atmo_solver_delete(&self->solver);
		atmo_renderer_delete(&self->renderer);
		atmo_uiInfoPanel_delete(&self->info_panel);
		vkk_uiGraphicsBox_delete(&self->graphics_box);
		vkk_uiWindow_delete((vkk_uiWindow_t**) &self);
	}
}
