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
#include "libcc/cc_memory.h"
#include "atmo_engine.h"

/***********************************************************
* public                                                   *
***********************************************************/

atmo_engine_t* atmo_engine_new(vkk_engine_t* engine)
{
	ASSERT(engine);

	atmo_engine_t* self;
	self = (atmo_engine_t*) CALLOC(1, sizeof(atmo_engine_t));
	if(self == NULL)
	{
		LOGE("CALLOC failed");
		return NULL;
	}

	self->engine = engine;

	return self;
}

void atmo_engine_delete(atmo_engine_t** _self)
{
	ASSERT(_self);

	atmo_engine_t* self = *_self;
	if(self)
	{
		FREE(self);
		*_self = NULL;
	}
}

void atmo_engine_draw(atmo_engine_t* self)
{
	ASSERT(self);

	vkk_renderer_t* rend;
	rend = vkk_engine_defaultRenderer(self->engine);

	float clear_color[4] =
	{
		1.0f, 0.0f, 1.0f, 1.0f,
	};

	if(vkk_renderer_beginDefault(rend, VKK_RENDERER_MODE_DRAW,
	                             clear_color) == 0)
	{
		return;
	}

	vkk_renderer_end(rend);
}

int atmo_engine_event(atmo_engine_t* self,
                      vkk_platformEvent_t* event)
{
	ASSERT(self);
	ASSERT(event);

	vkk_platformEventKey_t* e = &event->key;

	if(e->keycode == VKK_PLATFORM_KEYCODE_ESCAPE)
	{
		return 0;
	}

	return 1;
}
