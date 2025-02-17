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

#ifndef atmo_engine_H
#define atmo_engine_H

#include "libvkk/vkk_platform.h"
#include "libvkk/vkk.h"

typedef struct atmo_engine_s
{
	vkk_engine_t* engine;

	float ctrl_h;
	float ctrl_phi;
	float ctrl_delta;
	float ctrl_omega;

	uint32_t        sphere_ic;
	vkk_indexType_e sphere_it;
	vkk_buffer_t*   sphere_ib;
	vkk_buffer_t*   sphere_vb;

	vkk_uniformSetFactory_t* scene_usf0;
	vkk_pipelineLayout_t*    scene_pl;
	vkk_buffer_t*            scene_ub000_mvp;
	vkk_buffer_t*            scene_ub001_L;
	vkk_uniformSet_t*        scene_us0;

	vkk_graphicsPipeline_t*  planet_gp;
	vkk_graphicsPipeline_t*  sky_gp;
} atmo_engine_t;

atmo_engine_t* atmo_engine_new(vkk_engine_t* engine);
void           atmo_engine_delete(atmo_engine_t** _self);
void           atmo_engine_draw(atmo_engine_t* self);
int            atmo_engine_event(atmo_engine_t* self,
                                 vkk_platformEvent_t* event);

#endif
