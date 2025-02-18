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
#include "atmo_uiScreen.h"

/***********************************************************
* callbacks                                                *
***********************************************************/

void* atmo_onCreate(vkk_engine_t* engine)
{
	ASSERT(engine);

	return (void*) atmo_uiScreen_new(engine);
}

void atmo_onDestroy(void** _priv)
{
	ASSERT(_priv);

	atmo_uiScreen_delete((atmo_uiScreen_t**) _priv);
}

void atmo_onDraw(void* priv)
{
	ASSERT(priv);

	atmo_uiScreen_draw((atmo_uiScreen_t*) priv);
}

void atmo_onPause(void* priv)
{
	ASSERT(priv);

	// ignore
}

int atmo_onEvent(void* priv, vkk_platformEvent_t* event)
{
	ASSERT(priv);
	ASSERT(event);

	return atmo_uiScreen_event((atmo_uiScreen_t*) priv, event);
}

vkk_platformInfo_t VKK_PLATFORM_INFO =
{
	.app_name    = "atmosphere-elek",
	.app_version =
	{
		.major = 1,
		.minor = 0,
		.patch = 0,
	},
	.app_dir     = "atmosphere-elek",
	.onCreate    = atmo_onCreate,
	.onDestroy   = atmo_onDestroy,
	.onPause     = atmo_onPause,
	.onDraw      = atmo_onDraw,
	.onEvent     = atmo_onEvent,
};
