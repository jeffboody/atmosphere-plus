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
#include "atmo_uiScreen.h"
#include "atmo_uiWindowHud.h"

/***********************************************************
* public                                                   *
***********************************************************/

atmo_uiScreen_t* atmo_uiScreen_new(vkk_engine_t* engine)
{
	ASSERT(engine);

	vkk_renderer_t* rend;
	rend = vkk_engine_defaultRenderer(engine);

	char resource[256];
	snprintf(resource, 256, "%s/resource.bfs",
	         vkk_engine_internalPath(engine));

	vkk_uiWidgetStyle_t widget_style =
	{
		.color_primary =
		{
			.r=0.306f,
			.g=0.8f,
			.b=0.639f,
			.a=1.0f
		},
		.color_secondary =
		{
			.r=0.224f,
			.g=0.243f,
			.b=0.275f,
			.a=1.0f
		},
		.color_text =
		{
			.r=1.0f,
			.g=1.0f,
			.b=1.0f,
			.a=1.0f
		},
		.color_background =
		{
			.r=0.122f,
			.g=0.122f,
			.b=0.122f,
			.a=1.0f
		},
	};

	atmo_uiScreen_t* self;
	self = (atmo_uiScreen_t*)
	       vkk_uiScreen_new(sizeof(atmo_uiScreen_t), engine,
	                        rend, resource, &widget_style);
	if(self == NULL)
	{
		return NULL;
	}

	self->window_hud = atmo_uiWindowHud_new(self);
	if(self->window_hud == NULL)
	{
		goto failure;
	}

	vkk_uiScreen_windowReset(&self->base,
	                         &self->window_hud->base);

	// success
	return self;

	// failure
	failure:
		atmo_uiScreen_delete(&self);
	return NULL;
}

void atmo_uiScreen_delete(atmo_uiScreen_t** _self)
{
	ASSERT(_self);

	atmo_uiScreen_t* self = *_self;
	if(self)
	{
		atmo_uiWindowHud_delete(&self->window_hud);
		vkk_uiScreen_delete((vkk_uiScreen_t**) _self);
	}
}

void atmo_uiScreen_draw(atmo_uiScreen_t* self)
{
	ASSERT(self);

	vkk_renderer_t* rend;
	rend = vkk_engine_defaultRenderer(self->base.engine);

	float clear_color[4] =
	{
		0.0f, 0.0f, 0.0f, 1.0f
	};
	if(vkk_renderer_beginDefault(rend,
	                             VKK_RENDERER_MODE_DRAW,
	                             clear_color) == 0)
	{
		return;
	}
	vkk_uiScreen_draw(&self->base);
	vkk_renderer_end(rend);
}

int atmo_uiScreen_event(atmo_uiScreen_t* self,
                        vkk_platformEvent_t* event)
{
	ASSERT(self);
	ASSERT(event);

	if((event->type == VKK_PLATFORM_EVENTTYPE_ACTION_DOWN) ||
	   (event->type == VKK_PLATFORM_EVENTTYPE_ACTION_MOVE) ||
	   (event->type == VKK_PLATFORM_EVENTTYPE_ACTION_UP))
	{
		vkk_uiScreen_eventAction(&self->base, event);
	}
	else if(event->type == VKK_PLATFORM_EVENTTYPE_DENSITY)
	{
		vkk_uiScreen_eventDensity(&self->base, event->density);
	}
	else if((event->type == VKK_PLATFORM_EVENTTYPE_KEY_UP) ||
	        ((event->type == VKK_PLATFORM_EVENTTYPE_KEY_DOWN) &&
	         (event->key.repeat)))
	{
		atmo_renderer_t* renderer = self->window_hud->renderer;
		return atmo_renderer_event(renderer, event) ||
		       vkk_uiScreen_eventKey(&self->base, &event->key);
	}
	else if(event->type == VKK_PLATFORM_EVENTTYPE_CONTENT_RECT)
	{
		vkk_uiScreen_eventContentRect(&self->base,
		                              &event->content_rect);
	}

	return 1;
}
