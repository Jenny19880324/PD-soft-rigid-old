// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "MarqueeGL.h"
#include "bind_vertex_attrib_array.h"
#include "create_shader_program.h"
#include "destroy_shader_program.h"
#include "gl.h"
#include <iostream>

IGL_INLINE void igl::opengl::MarqueeGL::init()
{
	if (is_initialized) 
	{
		return;
	}
	is_initialized = true;
	std::string marquee_vertex_shader_string = 
	R"(#version 150
	in vec3 position;
	
	void main()
	{
		gl_Position = vec4(position, 1.0);
	}
	)";
	
	std::string marquee_fragment_shader_string = 
	R"(#version 150
	out vec4 outColor;
	void main()
	{
		outColor = vec4(0.2, 0.2, 0.2, 1.0);
	})";
	
	init_buffers();
	create_shader_program(
		marquee_vertex_shader_string,
		marquee_fragment_shader_string,
		{},
		shader_marquee);


}


IGL_INLINE void igl::opengl::MarqueeGL::free()
{
	const auto free = [](GLuint & id)
	{
		if (id)
		{
			destroy_shader_program(id);
			id = 0;
		}
	};
	
	if (is_initialized) 
	{
		free(shader_marquee);
		free_buffers();
	}
}


IGL_INLINE void igl::opengl::MarqueeGL::init_buffers()
{
	glGenVertexArrays(1, &vao_marquee);
	glBindVertexArray(vao_marquee);
	glGenBuffers(1, &vbo_V);
}


IGL_INLINE void igl::opengl::MarqueeGL::bind()
{
	glBindVertexArray(vao_marquee);
	glUseProgram(shader_marquee);
	bind_vertex_attrib_array(shader_marquee, "position", vbo_V, V_vbo, dirty);	
}


IGL_INLINE void igl::opengl::MarqueeGL::draw()
{
	glDisable(GL_DEPTH_TEST);

	double lw;
	glGetDoublev(GL_LINE_WIDTH, &lw);
	glLineWidth(1);
	glDrawArrays(GL_LINE_LOOP, 0, 4);
	glLineWidth(lw);
	glEnable(GL_DEPTH_TEST);
}


IGL_INLINE void igl::opengl::MarqueeGL::free_buffers()
{
	if (is_initialized) 
	{
		glDeleteVertexArrays(1, &vao_marquee);
		
		glDeleteBuffers(1, &vbo_V);
	}
}