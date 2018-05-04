// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "FanGL.h"
#include "bind_vertex_attrib_array.h"
#include "create_shader_program.h"
#include "destroy_shader_program.h"
#include "gl.h"
#include <iostream>
#include <cmath>
#include <math.h>
#include <Eigen/geometry>

IGL_INLINE void igl::opengl::FanGL::init()
{
	if (is_initialized) 
	{
		return;
	}
	is_initialized = true;
	std::string marquee_vertex_shader_string = 
	R"(#version 150
	in vec2 position;
	uniform mat3 proj;
	
	void main()
	{
		gl_Position = vec4(vec2(proj * vec3(position, 1.0)), 0.0, 1.0);
	}
	)";
	
	std::string marquee_fragment_shader_string = 
	R"(#version 150
	out vec4 outColor;
	void main()
	{
		outColor = vec4(0.2, 0.2, 0.2, 0.5);
	})";
	
	init_buffers();
	create_shader_program(
		marquee_vertex_shader_string,
		marquee_fragment_shader_string,
		{},
		shader_fan);


}


IGL_INLINE void igl::opengl::FanGL::free()
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
		free(shader_fan);
		free_buffers();
	}
}


IGL_INLINE void igl::opengl::FanGL::init_buffers()
{
	glGenVertexArrays(1, &vao_fan);
	glBindVertexArray(vao_fan);
	glGenBuffers(1, &vbo_V);
}


IGL_INLINE void igl::opengl::FanGL::bind()
{
	glBindVertexArray(vao_fan);
	glUseProgram(shader_fan);

	if (V3.rows() == 2) {
		V_vbo = V3;
	}
	else if (V3.rows() == 3){
		Eigen::Vector2f v0(0., 1.);
		Eigen::Vector2f v1 = V3.row(1) - V3.row(0);
		float r = v1.norm();
		v1.normalize();
		Eigen::Vector2f v2 = V3.row(2) - V3.row(0);
		v2.normalize();

		float rotate_angle = acos(v1.dot(v2));
		float start_angle = acos(v0.dot(v1));
		float end_angle = acos(v0.dot(v2));

		//assert(fabs(start_angle + rotate_angle - end_angle) < 1e-3);
		
		int number_of_triangles = (int)ceil(rotate_angle * 180.f / M_PI);
		V_vbo.resize(number_of_triangles + 2, Eigen::NoChange);
		V_vbo.setZero();
		V_vbo.row(0) = V3.row(0);
		V_vbo.row(1) = V3.row(1);
		for (int i = 0; i < number_of_triangles; i++) {
			V_vbo(i+2, 0) = V_vbo(0, 0) + r * sin(start_angle - i * M_PI / 180.f);
			V_vbo(i+2, 1) = V_vbo(0, 1) + r * cos(start_angle - i * M_PI / 180.f);
		}
	}


	bind_vertex_attrib_array(shader_fan, "position", vbo_V, V_vbo, dirty);	
}


IGL_INLINE void igl::opengl::FanGL::draw()
{
	GLint proji = glGetUniformLocation(shader_fan, "proj");
	glUniformMatrix3fv(proji, 1, GL_FALSE, proj.data());

	double lw;
	glGetDoublev(GL_LINE_WIDTH, &lw);
	glLineWidth(1);
	glDrawArrays(GL_LINE_LOOP, 0, V_vbo.rows());
	glDrawArrays(GL_TRIANGLE_FAN, 0, V_vbo.rows());
	glLineWidth(lw);
}


IGL_INLINE void igl::opengl::FanGL::free_buffers()
{
	if (is_initialized) 
	{
		glDeleteVertexArrays(1, &vao_fan);
		
		glDeleteBuffers(1, &vbo_V);
	}
}