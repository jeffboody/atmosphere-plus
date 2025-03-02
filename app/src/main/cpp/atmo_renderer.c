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
#include "libcc/math/cc_vec2f.h"
#include "libcc/math/cc_vec3f.h"
#include "libcc/cc_log.h"
#include "libcc/cc_memory.h"
#include "libgltf/gltf.h"
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

static int
atmo_renderer_newSphere(atmo_renderer_t* self,
                        gltf_file_t* glb)
{
	ASSERT(self);
	ASSERT(glb);

	gltf_accessorType_e  at1  = GLTF_ACCESSOR_TYPE_SCALAR;
	gltf_accessorType_e  at3  = GLTF_ACCESSOR_TYPE_VEC3;
	gltf_componentType_e ctus = GLTF_COMPONENT_TYPE_UNSIGNED_SHORT;
	gltf_componentType_e ctf  = GLTF_COMPONENT_TYPE_FLOAT;

	gltf_accessor_t* a_vb;
	a_vb = gltf_file_getAccessor(glb, 0);
	if((a_vb == NULL) ||
	   (a_vb->has_bufferView == 0)   ||
	   (a_vb->bufferView     != 0)   ||
	   (a_vb->type           != at3) ||
	   (a_vb->componentType  != ctf))
	{
		LOGE("invalid");
		return 0;
	}

	gltf_accessor_t* a_ib;
	a_ib = gltf_file_getAccessor(glb, 1);
	if((a_ib == NULL) ||
	   (a_ib->has_bufferView == 0)   ||
	   (a_ib->bufferView     != 1)   ||
	   (a_ib->type           != at1) ||
	   (a_ib->componentType  != ctus))
	{
		LOGE("invalid");
		return 0;
	}

	gltf_bufferView_t* bv_vb;
	gltf_bufferView_t* bv_ib;
	bv_vb = gltf_file_getBufferView(glb, 0);
	bv_ib = gltf_file_getBufferView(glb, 1);
	if((bv_vb == NULL) || (bv_ib == NULL))
	{
		LOGE("invalid");
		return 0;
	}

	const char* buf_vb;
	const char* buf_ib;
	buf_vb = gltf_file_getBuffer(glb, bv_vb);
	buf_ib = gltf_file_getBuffer(glb, bv_ib);
	if((buf_vb == NULL) || (buf_ib == NULL))
	{
		LOGE("invalid");
		return 0;
	}

	self->sphere_ic = a_ib->count;
	self->sphere_it = VKK_INDEX_TYPE_USHORT;

	self->sphere_ib = vkk_buffer_new(self->engine,
	                                 VKK_UPDATE_MODE_STATIC,
	                                 VKK_BUFFER_USAGE_INDEX,
	                                 a_ib->count*sizeof(unsigned short),
	                                 buf_ib);
	if(self->sphere_ib == NULL)
	{
		return 0;
	}

	self->sphere_vb = vkk_buffer_new(self->engine,
	                                 VKK_UPDATE_MODE_STATIC,
	                                 VKK_BUFFER_USAGE_VERTEX,
	                                 a_vb->count*sizeof(cc_vec3f_t),
	                                 buf_vb);
	if(self->sphere_vb == NULL)
	{
		goto fail_sphere_vb;
	}

	// success
	return 1;

	// failure
	fail_sphere_vb:
		vkk_buffer_delete(&self->sphere_ib);
	return 0;
}

static int atmo_renderer_importSphere(atmo_renderer_t* self)
{
	ASSERT(self);

	char resource[256];
	snprintf(resource, 256, "%s/resource.bfs",
	         vkk_engine_internalPath(self->engine));

	bfs_file_t* bfs;
	bfs = bfs_file_open(resource, 1, BFS_MODE_RDONLY);
	if(bfs == NULL)
	{
		return 0;
	}

	size_t size = 0;
	void*  data = NULL;
	if(bfs_file_blobGet(bfs, 0, "models/Sphere.glb",
	                    &size, &data) == 0)
	{
		goto fail_blob;
	}

	if(size == 0)
	{
		LOGE("invalid");
		goto fail_blob;
	}

	gltf_file_t* glb;
	glb = gltf_file_openb((char*) data, size,
	                      GLTF_FILEMODE_REFERENCE);
	if(glb == NULL)
	{
		goto fail_glb;
	}

	if(atmo_renderer_newSphere(self, glb) == 0)
	{
		goto fail_parse;
	}

	gltf_file_close(&glb);
	FREE(data);
	bfs_file_close(&bfs);

	// success
	return 1;

	// failure
	fail_parse:
		gltf_file_close(&glb);
	fail_glb:
		FREE(data);
	fail_blob:
		bfs_file_close(&bfs);
	return 0;
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

	if(atmo_renderer_importSphere(self) == 0)
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
	};

	self->scene_usf0 = vkk_uniformSetFactory_new(engine, um, 3,
	                                             scene_ub0_array);
	if(self->scene_usf0 == NULL)
	{
		goto failure;
	}

	vkk_uniformBinding_t scene_ub1_array[] =
	{
		// ub100_P0H
		{
			.binding = 0,
			.type    = VKK_UNIFORM_TYPE_BUFFER,
			.stage   = VKK_STAGE_VSFS,
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

	self->scene_ub100_P0H = vkk_buffer_new(engine, um,
	                                       VKK_BUFFER_USAGE_UNIFORM,
	                                       sizeof(cc_vec4f_t),
	                                       NULL);
	if(self->scene_ub100_P0H == NULL)
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
	};
	self->scene_us0 = vkk_uniformSet_new(engine, 0, 3,
	                                     scene_ua0_array,
	                                     self->scene_usf0);
	if(self->scene_us0 == NULL)
	{
		goto failure;
	}

	vkk_uniformAttachment_t scene_ua1_array[] =
	{
		// ub100_P0H
		{
			.binding = 0,
			.type    = VKK_UNIFORM_TYPE_BUFFER,
			.buffer  = self->scene_ub100_P0H,
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

	vkk_vertexBufferInfo_t planet_vbi_array[] =
	{
		// vertex
		{
			.location   = 0,
			.components = 3,
			.format     = VKK_VERTEX_FORMAT_FLOAT,
		},
	};

	vkk_graphicsPipelineInfo_t planet_gpi =
	{
		.renderer          = rend,
		.pl                = self->scene_pl,
		.vs                = "shaders/planet_vert.spv",
		.fs                = "shaders/planet_frag.spv",
		.vb_count          = 1,
		.vbi               = planet_vbi_array,
		.primitive         = VKK_PRIMITIVE_TRIANGLE_LIST,
		.primitive_restart = 0,
		.cull_mode         = VKK_CULL_MODE_BACK,
		.depth_test        = 1,
		.depth_write       = 1,
		.blend_mode        = VKK_BLEND_MODE_DISABLED,
	};
	self->planet_gp = vkk_graphicsPipeline_new(engine,
	                                           &planet_gpi);
	if(self->planet_gp == NULL)
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
	};

	vkk_graphicsPipelineInfo_t sky_flat_gpi =
	{
		.renderer          = rend,
		.pl                = self->scene_pl,
		.vs                = "shaders/sky_flat_vert.spv",
		.fs                = "shaders/sky_flat_frag.spv",
		.vb_count          = 1,
		.vbi               = sky_vbi_array,
		.primitive         = VKK_PRIMITIVE_TRIANGLE_LIST,
		.primitive_restart = 0,
		.cull_mode         = VKK_CULL_MODE_FRONT,
		.depth_test        = 1,
		.depth_write       = 1,
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
		.vs                = "shaders/sky_atmo_vert.spv",
		.fs                = "shaders/sky_atmo_frag.spv",
		.vb_count          = 1,
		.vbi               = sky_vbi_array,
		.primitive         = VKK_PRIMITIVE_TRIANGLE_LIST,
		.primitive_restart = 0,
		.cull_mode         = VKK_CULL_MODE_FRONT,
		.depth_test        = 1,
		.depth_write       = 1,
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
		vkk_graphicsPipeline_delete(&self->planet_gp);
		vkk_uniformSet_delete(&self->scene_us1);
		vkk_uniformSet_delete(&self->scene_us0);
		vkk_buffer_delete(&self->scene_ub103_phase_g_mie);
		vkk_buffer_delete(&self->scene_ub102_II4);
		vkk_buffer_delete(&self->scene_ub101_Zenith4);
		vkk_buffer_delete(&self->scene_ub100_P0H);
		vkk_buffer_delete(&self->scene_ub002_L4);
		vkk_buffer_delete(&self->scene_ub001_RaRp);
		vkk_buffer_delete(&self->scene_ub000_mvp);
		vkk_pipelineLayout_delete(&self->scene_pl);
		vkk_uniformSetFactory_delete(&self->scene_usf1);
		vkk_uniformSetFactory_delete(&self->scene_usf0);
		vkk_buffer_delete(&self->sphere_vb);
		vkk_buffer_delete(&self->sphere_ib);
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

	cc_vec4f_t Zenith4 =
	{
		.z = 1.0f,
	};

	cc_vec3f_t at =
	{
		.x = sin(phi),
		.z = cos(phi),
	};
	cc_vec3f_normalize(&at);

	cc_vec3f_t vy =
	{
		.y = 1.0f,
	};

	cc_vec3f_t up;
	cc_vec3f_cross_copy(&at, &vy, &up);
	cc_vec3f_normalize(&up);
	cc_vec3f_muls(&at, 1000.0f);
	cc_vec3f_addv(&at, &eye);

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
		     param->exposure*param->spectral_to_rgb_r,
		.g = param->spectral_irradiance_g*
		     param->exposure*param->spectral_to_rgb_g,
		.b = param->spectral_irradiance_b*
		     param->exposure*param->spectral_to_rgb_b,
	};

	cc_mat4f_t mvp;
	cc_mat4f_perspective(&mvp, 1, fovy, aspect,
	                     1.0f, 2.0f*param->Ra);
	cc_mat4f_lookat(&mvp, 0,
	                eye.x, eye.y, eye.z,
	                at.x,  at.y,  at.z,
	                up.x,  up.y,  up.z);
	vkk_renderer_updateBuffer(rend, self->scene_ub000_mvp,
	                          sizeof(cc_mat4f_t),
	                          (const void*) &mvp);
	vkk_renderer_updateBuffer(rend, self->scene_ub001_RaRp,
	                          sizeof(cc_vec2f_t),
	                          (const void*) &RaRp);
	vkk_renderer_updateBuffer(rend, self->scene_ub002_L4,
	                          sizeof(cc_vec4f_t),
	                          (const void*) &L);
	vkk_renderer_bindGraphicsPipeline(rend, self->planet_gp);
	vkk_renderer_bindUniformSets(rend, 1, &self->scene_us0);

	vkk_buffer_t* vertex_buffers[] =
	{
		self->sphere_vb,
	};
	vkk_renderer_drawIndexed(rend, self->sphere_ic, 1,
	                         self->sphere_it,
	                         self->sphere_ib,
	                         vertex_buffers);

	vkk_image_t* image = atmo_solver_image(solver, k);
	if(image)
	{
		vkk_renderer_updateBuffer(rend, self->scene_ub100_P0H,
		                          sizeof(cc_vec4f_t),
		                          (const void*) &P0h);
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

		vkk_uniformSet_t* us_array[] =
		{
			self->scene_us0,
			self->scene_us1,
		};
		vkk_renderer_bindUniformSets(rend, 2, us_array);
		vkk_renderer_bindGraphicsPipeline(rend,
		                                  self->sky_atmo_gp);
	}
	else
	{
		vkk_renderer_bindGraphicsPipeline(rend,
		                                  self->sky_flat_gp);
	}

	vkk_renderer_drawIndexed(rend, self->sphere_ic, 1,
	                         self->sphere_it,
	                         self->sphere_ib,
	                         vertex_buffers);
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
