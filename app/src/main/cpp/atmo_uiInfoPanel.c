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
#include "libcc/math/cc_float.h"
#include "libcc/cc_log.h"
#include "atmo_renderer.h"
#include "atmo_solver.h"
#include "atmo_uiInfoPanel.h"
#include "atmo_uiScreen.h"
#include "atmo_uiWindowHud.h"

/***********************************************************
* private                                                  *
***********************************************************/

static int
atmo_uiInfoPanel_refresh(vkk_uiWidget_t* widget)
{
	ASSERT(widget);

	atmo_uiInfoPanel_t* self;
	self = (atmo_uiInfoPanel_t*) widget;

	atmo_uiScreen_t* screen;
	screen = (atmo_uiScreen_t*) widget->screen;

	atmo_uiWindowHud_t* window_hud = screen->window_hud;
	atmo_renderer_t*    renderer   = window_hud->renderer;
	atmo_solver_t*      solver     = window_hud->solver;
	atmo_solverParam_t* param      = &solver->param;

	if(self->last_ctrl_h != renderer->ctrl_h)
	{
		vkk_uiText_label(self->text_ctrl_h, "i/o: h=%0.1f",
		                 atmo_renderer_getH(renderer, param));
		self->last_ctrl_h = renderer->ctrl_h;
	}

	if(self->last_ctrl_phi != renderer->ctrl_phi)
	{
		vkk_uiText_label(self->text_ctrl_phi, "j/k: phi=%0.1f",
		                 cc_rad2deg(atmo_renderer_getPhi(renderer)));
		self->last_ctrl_phi = renderer->ctrl_phi;
	}

	if(self->last_ctrl_delta != renderer->ctrl_delta)
	{
		vkk_uiText_label(self->text_ctrl_delta, "w/s: delta=%0.1f",
		                 cc_rad2deg(atmo_renderer_getDelta(renderer)));
		self->last_ctrl_delta = renderer->ctrl_delta;
	}

	if(self->last_ctrl_omega != renderer->ctrl_omega)
	{
		vkk_uiText_label(self->text_ctrl_omega, "a/d: omega=%0.1f",
		                 cc_rad2deg(atmo_renderer_getOmega(renderer)));
		self->last_ctrl_omega = renderer->ctrl_omega;
	}

	if(self->last_ctrl_k != renderer->ctrl_k)
	{
		vkk_uiText_label(self->text_ctrl_k, "123: k=%u",
		                 atmo_renderer_getK(renderer));
		self->last_ctrl_k = renderer->ctrl_k;
	}

	return 0;
}

/***********************************************************
* public                                                   *
***********************************************************/

atmo_uiInfoPanel_t*
atmo_uiInfoPanel_new(atmo_uiScreen_t* screen)
{
	ASSERT(overlay);

	vkk_uiInfoPanelFn_t ipfn =
	{
		.refresh_fn = atmo_uiInfoPanel_refresh
	};

	atmo_uiInfoPanel_t* self;
	self = (atmo_uiInfoPanel_t*)
	       vkk_uiInfoPanel_new(&screen->base,
	                           sizeof(atmo_uiInfoPanel_t),
	                           &ipfn);
	if(self == NULL)
	{
		return NULL;
	}

	self->heading_ctrl = vkk_uiText_newInfoHeading(&screen->base);
	if(self->heading_ctrl == NULL)
	{
		goto failure;
	}

	self->text_ctrl_h = vkk_uiText_newInfoItem(&screen->base);
	if(self->text_ctrl_h == NULL)
	{
		goto failure;
	}

	self->text_ctrl_phi = vkk_uiText_newInfoItem(&screen->base);
	if(self->text_ctrl_phi == NULL)
	{
		goto failure;
	}

	self->text_ctrl_delta = vkk_uiText_newInfoItem(&screen->base);
	if(self->text_ctrl_delta == NULL)
	{
		goto failure;
	}

	self->text_ctrl_omega = vkk_uiText_newInfoItem(&screen->base);
	if(self->text_ctrl_omega == NULL)
	{
		goto failure;
	}

	self->text_ctrl_k = vkk_uiText_newInfoItem(&screen->base);
	if(self->text_ctrl_k == NULL)
	{
		goto failure;
	}

	vkk_uiText_label(self->heading_ctrl, "%s", "Controls");
	vkk_uiText_label(self->text_ctrl_h,     "i/o: h=%0.1f", 0.0f);
	vkk_uiText_label(self->text_ctrl_phi,   "j/k: phi=%0.1f", 0.0f);
	vkk_uiText_label(self->text_ctrl_delta, "w/s: delta=%0.1f", 0.0f);
	vkk_uiText_label(self->text_ctrl_omega, "a/d: omega=%0.1f", 0.0f);
	vkk_uiText_label(self->text_ctrl_k,     "123: k=%u", 0);

	vkk_uiInfoPanel_t* ip = &self->base;
	vkk_uiInfoPanel_add(ip, (vkk_uiWidget_t*) self->heading_ctrl);
	vkk_uiInfoPanel_add(ip, (vkk_uiWidget_t*) self->text_ctrl_h);
	vkk_uiInfoPanel_add(ip, (vkk_uiWidget_t*) self->text_ctrl_phi);
	vkk_uiInfoPanel_add(ip, (vkk_uiWidget_t*) self->text_ctrl_delta);
	vkk_uiInfoPanel_add(ip, (vkk_uiWidget_t*) self->text_ctrl_omega);
	vkk_uiInfoPanel_add(ip, (vkk_uiWidget_t*) self->text_ctrl_k);

	// success
	return self;

	// failure
	failure:
		atmo_uiInfoPanel_delete(&self);
	return NULL;
}

void atmo_uiInfoPanel_delete(atmo_uiInfoPanel_t** _self)
{
	ASSERT(_self);

	atmo_uiInfoPanel_t* self = *_self;
	if(self)
	{
		vkk_uiInfoPanel_clear(&self->base);
		vkk_uiText_delete(&self->text_ctrl_k);
		vkk_uiText_delete(&self->text_ctrl_omega);
		vkk_uiText_delete(&self->text_ctrl_delta);
		vkk_uiText_delete(&self->text_ctrl_phi);
		vkk_uiText_delete(&self->text_ctrl_h);
		vkk_uiText_delete(&self->heading_ctrl);
		vkk_uiInfoPanel_delete((vkk_uiInfoPanel_t**) _self);
	}
}
