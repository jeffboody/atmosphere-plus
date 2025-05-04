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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LOG_TAG "atmo"
#include "libbfs/bfs_file.h"
#include "libcc/math/cc_float.h"
#include "libcc/math/cc_mat4f.h"
#include "libcc/math/cc_vec2f.h"
#include "libcc/math/cc_vec3f.h"
#include "libcc/cc_log.h"
#include "libcc/cc_memory.h"
#include "atmo_renderer.h"
#include "atmo_solver.h"

/***********************************************************
* private                                                  *
***********************************************************/

static void atmo_renderer_resetCtrl(atmo_renderer_t* self)
{
	ASSERT(self);

	self->ctrl_h     = 0.1f;
	self->ctrl_phi   = 0.5f;
	self->ctrl_delta = 0.0f;
	self->ctrl_omega = 0.0f;
	self->ctrl_k     = 1;
}

/***********************************************************
* public                                                   *
***********************************************************/

atmo_renderer_t* atmo_renderer_new(vkk_engine_t* engine)
{
	ASSERT(engine);

	vkk_renderer_t* rend;
	rend = vkk_engine_defaultRenderer(engine);

	vkk_updateMode_e um = vkk_renderer_updateMode(rend);

	atmo_renderer_t* self;
	self = (atmo_renderer_t*)
	       CALLOC(1, sizeof(atmo_renderer_t));
	if(self == NULL)
	{
		LOGE("CALLOC failed");
		return NULL;
	}

	self->engine = engine;

	atmo_renderer_resetCtrl(self);

	cc_vec3f_t vertices[] =
	{
		{ .x=-1.0f, .y= 1.0f, .z=-1.0f, },
		{ .x=-1.0f, .y=-1.0f, .z=-1.0f, },
		{ .x= 1.0f, .y= 1.0f, .z=-1.0f, },
		{ .x= 1.0f, .y=-1.0f, .z=-1.0f, },
	};
	self->vb_vertex = vkk_buffer_new(self->engine,
	                                 VKK_UPDATE_MODE_STATIC,
	                                 VKK_BUFFER_USAGE_VERTEX,
	                                 sizeof(vertices),
	                                 vertices);
	if(self->vb_vertex == NULL)
	{
		goto failure;
	}

	self->vb_V = vkk_buffer_new(self->engine, um,
	                            VKK_BUFFER_USAGE_VERTEX,
	                            sizeof(vertices),
	                            NULL);
	if(self->vb_V == NULL)
	{
		goto failure;
	}

	vkk_uniformBinding_t scene_ub0_array[] =
	{
		// ub000_mvp
		{
			.binding = 0,
			.type    = VKK_UNIFORM_TYPE_BUFFER,
			.stage   = VKK_STAGE_VS,
		},
		// ub001_RaRp
		{
			.binding = 1,
			.type    = VKK_UNIFORM_TYPE_BUFFER,
			.stage   = VKK_STAGE_VSFS,
		},
		// ub002_L4
		{
			.binding = 2,
			.type    = VKK_UNIFORM_TYPE_BUFFER,
			.stage   = VKK_STAGE_FS,
		},
		// ub003_P0H
		{
			.binding = 3,
			.type    = VKK_UNIFORM_TYPE_BUFFER,
			.stage   = VKK_STAGE_VSFS,
		},
	};

	self->scene_usf0 = vkk_uniformSetFactory_new(engine, um, 4,
	                                             scene_ub0_array);
	if(self->scene_usf0 == NULL)
	{
		goto failure;
	}

	vkk_uniformBinding_t scene_ub1_array[] =
	{
		// ub100_Unused
		{
			.binding = 0,
			.type    = VKK_UNIFORM_TYPE_BUFFER,
			.stage   = VKK_STAGE_FS,
		},
		// ub101_Zenith4
		{
			.binding = 1,
			.type    = VKK_UNIFORM_TYPE_BUFFER,
			.stage   = VKK_STAGE_FS,
		},
		// ub102_II4
		{
			.binding = 2,
			.type    = VKK_UNIFORM_TYPE_BUFFER,
			.stage   = VKK_STAGE_FS,
		},
		// ub103_phase_g_mie
		{
			.binding = 3,
			.type    = VKK_UNIFORM_TYPE_BUFFER,
			.stage   = VKK_STAGE_FS,
		},
		// sampler104_fIS
		{
			.binding = 4,
			.type    = VKK_UNIFORM_TYPE_IMAGE_REF,
			.stage   = VKK_STAGE_FS,
			.si      =
			{
				.min_filter     = VKK_SAMPLER_FILTER_LINEAR,
				.mag_filter     = VKK_SAMPLER_FILTER_LINEAR,
				.mipmap_mode    = VKK_SAMPLER_MIPMAP_MODE_NEAREST,
				.anisotropy     = 0,
				.max_anisotropy = 0.0f,
			},
		},
	};

	self->scene_usf1 = vkk_uniformSetFactory_new(engine, um, 5,
	                                             scene_ub1_array);
	if(self->scene_usf1 == NULL)
	{
		goto failure;
	}

	vkk_uniformSetFactory_t* usf_array[] =
	{
		self->scene_usf0,
		self->scene_usf1,
	};
	self->scene_pl = vkk_pipelineLayout_new(engine, 2,
	                                        usf_array);
	if(self->scene_pl == NULL)
	{
		goto failure;
	}

	self->scene_ub000_mvp = vkk_buffer_new(engine, um,
	                                       VKK_BUFFER_USAGE_UNIFORM,
	                                       sizeof(cc_mat4f_t),
	                                       NULL);
	if(self->scene_ub000_mvp == NULL)
	{
		goto failure;
	}

	self->scene_ub001_RaRp = vkk_buffer_new(engine, um,
	                                        VKK_BUFFER_USAGE_UNIFORM,
	                                        sizeof(cc_vec2f_t),
	                                        NULL);
	if(self->scene_ub001_RaRp == NULL)
	{
		goto failure;
	}

	self->scene_ub002_L4 = vkk_buffer_new(engine, um,
	                                      VKK_BUFFER_USAGE_UNIFORM,
	                                      sizeof(cc_vec4f_t),
	                                      NULL);
	if(self->scene_ub002_L4 == NULL)
	{
		goto failure;
	}

	self->scene_ub003_P0H = vkk_buffer_new(engine, um,
	                                       VKK_BUFFER_USAGE_UNIFORM,
	                                       sizeof(cc_vec4f_t),
	                                       NULL);
	if(self->scene_ub003_P0H == NULL)
	{
		goto failure;
	}

	self->scene_ub100_Unused = vkk_buffer_new(engine, um,
	                                          VKK_BUFFER_USAGE_UNIFORM,
	                                          sizeof(cc_vec4f_t),
	                                          NULL);
	if(self->scene_ub100_Unused == NULL)
	{
		goto failure;
	}

	self->scene_ub101_Zenith4 = vkk_buffer_new(engine, um,
	                                           VKK_BUFFER_USAGE_UNIFORM,
	                                           sizeof(cc_vec4f_t),
	                                           NULL);
	if(self->scene_ub101_Zenith4 == NULL)
	{
		goto failure;
	}

	self->scene_ub102_II4 = vkk_buffer_new(engine, um,
	                                       VKK_BUFFER_USAGE_UNIFORM,
	                                       sizeof(cc_vec4f_t),
	                                       NULL);
	if(self->scene_ub102_II4 == NULL)
	{
		goto failure;
	}

	self->scene_ub103_phase_g_mie = vkk_buffer_new(engine, um,
	                                               VKK_BUFFER_USAGE_UNIFORM,
	                                               sizeof(float),
	                                               NULL);
	if(self->scene_ub103_phase_g_mie == NULL)
	{
		goto failure;
	}

	vkk_uniformAttachment_t scene_ua0_array[] =
	{
		// ub000_mvp
		{
			.binding = 0,
			.type    = VKK_UNIFORM_TYPE_BUFFER,
			.buffer  = self->scene_ub000_mvp,
		},
		// ub001_RaRp
		{
			.binding = 1,
			.type    = VKK_UNIFORM_TYPE_BUFFER,
			.buffer  = self->scene_ub001_RaRp,
		},
		// ub002_L4
		{
			.binding = 2,
			.type    = VKK_UNIFORM_TYPE_BUFFER,
			.buffer  = self->scene_ub002_L4,
		},
		// ub003_P0H
		{
			.binding = 3,
			.type    = VKK_UNIFORM_TYPE_BUFFER,
			.buffer  = self->scene_ub003_P0H,
		},
	};
	self->scene_us0 = vkk_uniformSet_new(engine, 0, 4,
	                                     scene_ua0_array,
	                                     self->scene_usf0);
	if(self->scene_us0 == NULL)
	{
		goto failure;
	}

	vkk_uniformAttachment_t scene_ua1_array[] =
	{
		// ub100_Unused
		{
			.binding = 0,
			.type    = VKK_UNIFORM_TYPE_BUFFER,
			.buffer  = self->scene_ub100_Unused,
		},
		// ub101_Zenith4
		{
			.binding = 1,
			.type    = VKK_UNIFORM_TYPE_BUFFER,
			.buffer  = self->scene_ub101_Zenith4,
		},
		// ub102_II4
		{
			.binding = 2,
			.type    = VKK_UNIFORM_TYPE_BUFFER,
			.buffer  = self->scene_ub102_II4,
		},
		// ub103_phase_g_mie
		{
			.binding = 3,
			.type    = VKK_UNIFORM_TYPE_BUFFER,
			.buffer  = self->scene_ub103_phase_g_mie,
		},
	};
	self->scene_us1 = vkk_uniformSet_new(engine, 1, 4,
	                                     scene_ua1_array,
	                                     self->scene_usf1);
	if(self->scene_us1 == NULL)
	{
		goto failure;
	}

	vkk_vertexBufferInfo_t sky_vbi_array[] =
	{
		// vertex
		{
			.location   = 0,
			.components = 3,
			.format     = VKK_VERTEX_FORMAT_FLOAT,
		},
		// V
		{
			.location   = 1,
			.components = 3,
			.format     = VKK_VERTEX_FORMAT_FLOAT,
		},
	};

	vkk_graphicsPipelineInfo_t sky_flat_gpi =
	{
		.renderer          = rend,
		.pl                = self->scene_pl,
		.vs                = "shaders/sky_vert.spv",
		.fs                = "shaders/sky_flat_frag.spv",
		.vb_count          = 2,
		.vbi               = sky_vbi_array,
		.primitive         = VKK_PRIMITIVE_TRIANGLE_STRIP,
		.primitive_restart = 0,
		.cull_mode         = VKK_CULL_MODE_NONE,
		.depth_test        = 0,
		.depth_write       = 0,
		.blend_mode        = VKK_BLEND_MODE_DISABLED,
	};
	self->sky_flat_gp = vkk_graphicsPipeline_new(engine, &sky_flat_gpi);
	if(self->sky_flat_gp == NULL)
	{
		goto failure;
	}

	vkk_graphicsPipelineInfo_t sky_atmo_gpi =
	{
		.renderer          = rend,
		.pl                = self->scene_pl,
		.vs                = "shaders/sky_vert.spv",
		.fs                = "shaders/sky_atmo_frag.spv",
		.vb_count          = 2,
		.vbi               = sky_vbi_array,
		.primitive         = VKK_PRIMITIVE_TRIANGLE_STRIP,
		.primitive_restart = 0,
		.cull_mode         = VKK_CULL_MODE_NONE,
		.depth_test        = 0,
		.depth_write       = 0,
		.blend_mode        = VKK_BLEND_MODE_DISABLED,
	};
	self->sky_atmo_gp = vkk_graphicsPipeline_new(engine, &sky_atmo_gpi);
	if(self->sky_atmo_gp == NULL)
	{
		goto failure;
	}

	// success
	return self;

	// failure
	failure:
		atmo_renderer_delete(&self);
	return NULL;
}

void atmo_renderer_delete(atmo_renderer_t** _self)
{
	ASSERT(_self);

	atmo_renderer_t* self = *_self;
	if(self)
	{
		vkk_graphicsPipeline_delete(&self->sky_atmo_gp);
		vkk_graphicsPipeline_delete(&self->sky_flat_gp);
		vkk_uniformSet_delete(&self->scene_us1);
		vkk_uniformSet_delete(&self->scene_us0);
		vkk_buffer_delete(&self->scene_ub103_phase_g_mie);
		vkk_buffer_delete(&self->scene_ub102_II4);
		vkk_buffer_delete(&self->scene_ub101_Zenith4);
		vkk_buffer_delete(&self->scene_ub100_Unused);
		vkk_buffer_delete(&self->scene_ub003_P0H);
		vkk_buffer_delete(&self->scene_ub002_L4);
		vkk_buffer_delete(&self->scene_ub001_RaRp);
		vkk_buffer_delete(&self->scene_ub000_mvp);
		vkk_pipelineLayout_delete(&self->scene_pl);
		vkk_uniformSetFactory_delete(&self->scene_usf1);
		vkk_uniformSetFactory_delete(&self->scene_usf0);
		vkk_buffer_delete(&self->vb_V);
		vkk_buffer_delete(&self->vb_vertex);
		FREE(self);
		*_self = NULL;
	}
}

void atmo_renderer_draw(atmo_renderer_t* self,
                        atmo_solver_t* solver,
                        float width, float height)
{
	ASSERT(self);
	ASSERT(solver);

	vkk_renderer_t* rend;
	rend = vkk_engine_defaultRenderer(self->engine);

	atmo_solverParam_t* param = &solver->param;

	float aspect = width/height;
	float fovy   = 90.0f;
	if(aspect < 1.0f)
	{
		fovy /= aspect;
	}

	float    h     = atmo_renderer_getH(self, param);
	float    phi   = atmo_renderer_getPhi(self);
	float    delta = atmo_renderer_getDelta(self);
	float    omega = atmo_renderer_getOmega(self);
	uint32_t k     = atmo_renderer_getK(self);

	cc_vec3f_t eye =
	{
		.z = h + param->Rp,
	};

	cc_vec4f_t P0h =
	{
		.x = eye.x,
		.y = eye.y,
		.z = eye.z,
		.w = h,
	};

	cc_vec3f_t at =
	{
		.x = sin(phi),
		.z = cos(phi),
	};
	cc_vec3f_normalize(&at);

	cc_vec3f_t y =
	{
		.y = 1.0f,
	};

	// dist is used to compute the plane corners
	// a larger distance improves numerical stability
	float dist = 1000.0f;

	cc_vec3f_t up;
	cc_vec3f_cross_copy(&at, &y, &up);
	cc_vec3f_normalize(&up);
	cc_vec3f_muls(&at, dist);
	cc_vec3f_addv(&at, &eye);

	cc_mat4f_t mvm;
	cc_mat4f_lookat(&mvm, 1,
	                eye.x, eye.y, eye.z,
	                at.x,  at.y,  at.z,
	                up.x,  up.y,  up.z);

	// view plane normal
	cc_vec3f_t vpn;
	cc_vec3f_load(&vpn, -mvm.m20, -mvm.m21, -mvm.m22);

	// compute the vectors for the plane
	cc_vec3f_t vx;
	cc_vec3f_t vy;
	cc_vec3f_t vz;
	cc_vec3f_load(&vx, mvm.m00, mvm.m01, mvm.m02);
	cc_vec3f_load(&vy, mvm.m10, mvm.m11, mvm.m12);
	cc_vec3f_muls_copy(&vpn, dist, &vz);
	cc_vec3f_muls(&vy, tanf((fovy/2.0f)*(M_PI/180.0f))*
	              cc_vec3f_mag(&vz));
	cc_vec3f_muls(&vx, aspect*cc_vec3f_mag(&vy));

	// compute the vectors of the plane corners
	cc_vec3f_t v00;
	cc_vec3f_t v01;
	cc_vec3f_t v10;
	cc_vec3f_t v11;
	cc_vec3f_subv_copy(&vz, &vx, &v00);
	cc_vec3f_subv_copy(&v00, &vy, &v10);
	cc_vec3f_addv(&v00, &vy);
	cc_vec3f_addv_copy(&vz, &vx, &v01);
	cc_vec3f_subv_copy(&v01, &vy, &v11);
	cc_vec3f_addv(&v01, &vy);

	// compute direction rays for sky shaders
	cc_vec3f_t V[4];
	cc_vec3f_normalize_copy(&v00, &V[0]);
	cc_vec3f_normalize_copy(&v10, &V[1]);
	cc_vec3f_normalize_copy(&v01, &V[2]);
	cc_vec3f_normalize_copy(&v11, &V[3]);

	vkk_renderer_updateBuffer(rend, self->vb_V,
	                          4*sizeof(cc_vec3f_t), V);

	cc_vec2f_t RaRp = { .x = param->Ra, .y = param->Rp };

	cc_vec4f_t L =
	{
		.x = -sin(delta)*cos(omega),
		.y = -sin(delta)*sin(omega),
		.z = -cos(delta),
	};
	cc_vec4f_normalize(&L);

	cc_vec4f_t II4 =
	{
		.r = param->spectral_irradiance_r*
		     param->spectral_to_rgb_r*
		     powf(2.0f, param->exposure),
		.g = param->spectral_irradiance_g*
		     param->spectral_to_rgb_g*
		     powf(2.0f, param->exposure),
		.b = param->spectral_irradiance_b*
		     param->spectral_to_rgb_b*
		     powf(2.0f, param->exposure),
	};

	cc_mat4f_t mvp;
	cc_mat4f_ortho(&mvp, 1, -1.0f, 1.0f,
	               -1.0f, 1.0f, 0.0f, 2.0f);
	vkk_renderer_updateBuffer(rend, self->scene_ub000_mvp,
	                          sizeof(cc_mat4f_t),
	                          (const void*) &mvp);
	vkk_renderer_updateBuffer(rend, self->scene_ub001_RaRp,
	                          sizeof(cc_vec2f_t),
	                          (const void*) &RaRp);
	vkk_renderer_updateBuffer(rend, self->scene_ub002_L4,
	                          sizeof(cc_vec4f_t),
	                          (const void*) &L);
	vkk_renderer_updateBuffer(rend, self->scene_ub003_P0H,
	                          sizeof(cc_vec4f_t),
	                          (const void*) &P0h);

	vkk_uniformSet_t* us_array[] =
	{
		self->scene_us0,
		self->scene_us1,
	};

	vkk_image_t* image = atmo_solver_image(solver, k);
	if(image)
	{
		cc_vec4f_t Unused = { 0 };
		cc_vec4f_t Zenith4 =
		{
			.z = 1.0f,
		};
		vkk_renderer_updateBuffer(rend, self->scene_ub100_Unused,
		                          sizeof(cc_vec4f_t),
		                          (const void*) &Unused);
		vkk_renderer_updateBuffer(rend, self->scene_ub101_Zenith4,
		                          sizeof(cc_vec4f_t),
		                          (const void*) &Zenith4);
		vkk_renderer_updateBuffer(rend, self->scene_ub102_II4,
		                          sizeof(cc_vec4f_t),
		                          (const void*) &II4);
		vkk_renderer_updateBuffer(rend, self->scene_ub103_phase_g_mie,
		                          sizeof(float),
		                          (const void*) &param->phase_g_mie);
		vkk_uniformAttachment_t ua_array[] =
		{
			// sampler104_fIS
			{
				.binding = 4,
				.type    = VKK_UNIFORM_TYPE_IMAGE_REF,
				.image   = image,
			},
		};
		vkk_renderer_updateUniformSetRefs(rend, self->scene_us1,
		                                  1, ua_array);

		vkk_renderer_bindGraphicsPipeline(rend,
		                                  self->sky_atmo_gp);
		vkk_renderer_bindUniformSets(rend, 2, us_array);
	}
	else
	{
		vkk_renderer_bindGraphicsPipeline(rend,
		                                  self->sky_flat_gp);
		vkk_renderer_bindUniformSets(rend, 1, us_array);
	}

	vkk_buffer_t* vertex_buffers[] =
	{
		self->vb_vertex,
		self->vb_V,
	};
	vkk_renderer_draw(rend, 4, 2, vertex_buffers);
}

int atmo_renderer_event(atmo_renderer_t* self,
                        vkk_platformEvent_t* event)
{
	ASSERT(self);
	ASSERT(event);

	vkk_platformEventKey_t* e = &event->key;
	if((event->type == VKK_PLATFORM_EVENTTYPE_KEY_UP) ||
	   ((event->type == VKK_PLATFORM_EVENTTYPE_KEY_DOWN) &&
	    (event->key.repeat)))
	{
		if(e->keycode == 'i')
		{
			self->ctrl_h = cc_clamp(self->ctrl_h - 0.005f,
			                        0.0f, 1.0f);
			return 1;
		}
		else if(e->keycode == 'o')
		{
			self->ctrl_h = cc_clamp(self->ctrl_h + 0.005f,
			                        0.0f, 1.0f);
			return 1;
		}
		else if(e->keycode == 'j')
		{
			self->ctrl_phi = cc_clamp(self->ctrl_phi + 0.005f,
			                          0.0f, 1.0f);
			return 1;
		}
		else if(e->keycode == 'k')
		{
			self->ctrl_phi = cc_clamp(self->ctrl_phi - 0.005f,
			                          0.0f, 1.0f);
			return 1;
		}
		else if(e->keycode == 'w')
		{
			self->ctrl_delta = cc_clamp(self->ctrl_delta - 0.005f,
			                            0.0f, 1.0f);
			return 1;
		}
		else if(e->keycode == 's')
		{
			self->ctrl_delta = cc_clamp(self->ctrl_delta + 0.005f,
			                            0.0f, 1.0f);
			return 1;
		}
		else if(e->keycode == 'a')
		{
			self->ctrl_omega = self->ctrl_omega + 0.005f;
			while(self->ctrl_omega >= 1.0f)
			{
				self->ctrl_omega -= 1.0f;
			}
			return 1;
		}
		else if(e->keycode == 'd')
		{
			self->ctrl_omega = self->ctrl_omega - 0.005f;
			while(self->ctrl_omega < 0.0f)
			{
				self->ctrl_omega += 1.0f;
			}
			return 1;
		}
		else if(e->keycode == 'r')
		{
			atmo_renderer_resetCtrl(self);
			return 1;
		}
		else if((e->keycode >= '0') && (e->keycode <= '9'))
		{
			self->ctrl_k = e->keycode - '0';
			return 1;
		}
	}

	return 0;
}

float atmo_renderer_getH(atmo_renderer_t* self,
                         atmo_solverParam_t* param)
{
	ASSERT(self);

	// clamp h to avoid near clipping plane
	float Ha = param->Ra - param->Rp;
	return cc_clamp(self->ctrl_h*Ha, 3.0f, Ha - 10.0f);
}

float atmo_renderer_getPhi(atmo_renderer_t* self)
{
	ASSERT(self);

	return cc_deg2rad(180.0f*self->ctrl_phi);
}

float atmo_renderer_getDelta(atmo_renderer_t* self)
{
	ASSERT(self);

	return cc_deg2rad(180.0f*self->ctrl_delta);
}

float atmo_renderer_getOmega(atmo_renderer_t* self)
{
	ASSERT(self);

	return cc_deg2rad(360.0f*self->ctrl_omega);
}

uint32_t atmo_renderer_getK(atmo_renderer_t* self)
{
	ASSERT(self);

	return (self->ctrl_k == 0) ? 10 : self->ctrl_k;
}
