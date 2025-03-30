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

#ifndef atmo_renderer_H
#define atmo_renderer_H

#include "libvkk/vkk_platform.h"
#include "libvkk/vkk.h"

typedef struct atmo_solver_s      atmo_solver_t;
typedef struct atmo_solverParam_s atmo_solverParam_t;

typedef struct atmo_renderer_s
{
	vkk_engine_t* engine;

	float    ctrl_h;
	float    ctrl_phi;
	float    ctrl_delta;
	float    ctrl_omega;
	uint32_t ctrl_k;

	vkk_buffer_t* vb_vertex;
	vkk_buffer_t* vb_V;

	vkk_uniformSetFactory_t* scene_usf0;
	vkk_uniformSetFactory_t* scene_usf1;
	vkk_pipelineLayout_t*    scene_pl;
	vkk_buffer_t*            scene_ub000_mvp;
	vkk_buffer_t*            scene_ub001_RaRp; // Ra, Rp
	vkk_buffer_t*            scene_ub002_L4;
	vkk_buffer_t*            scene_ub003_P0H; // P0, H
	vkk_buffer_t*            scene_ub100_Unused;
	vkk_buffer_t*            scene_ub101_Zenith4;
	vkk_buffer_t*            scene_ub102_II4;
	vkk_buffer_t*            scene_ub103_phase_g_mie;
	vkk_uniformSet_t*        scene_us0;
	vkk_uniformSet_t*        scene_us1;

	vkk_graphicsPipeline_t*  sky_flat_gp;
	vkk_graphicsPipeline_t*  sky_atmo_gp;
} atmo_renderer_t;

atmo_renderer_t* atmo_renderer_new(vkk_engine_t* engine);
void             atmo_renderer_delete(atmo_renderer_t** _self);
void             atmo_renderer_draw(atmo_renderer_t* self,
                                    atmo_solver_t* solver,
                                    float width, float height);
int              atmo_renderer_event(atmo_renderer_t* self,
                                     vkk_platformEvent_t* event);
float            atmo_renderer_getH(atmo_renderer_t* self,
                                    atmo_solverParam_t* param);
float            atmo_renderer_getPhi(atmo_renderer_t* self);
float            atmo_renderer_getDelta(atmo_renderer_t* self);
float            atmo_renderer_getOmega(atmo_renderer_t* self);
uint32_t         atmo_renderer_getK(atmo_renderer_t* self);

#endif
