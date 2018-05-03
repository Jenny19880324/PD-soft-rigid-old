// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_OPENGL_FAN_GL_H
#define IGL_OPENGL_FAN_GL_H

// Create a fan from screen points
// compatible format 
//  The class includes a shader and the opengl calls to plot the data

#include <igl/igl_inline.h>
#include <Eigen/Core>


namespace igl
{
  namespace opengl
  {
	class FanGL 
	{
		public:
		typedef unsigned int GLuint;
		
		bool is_initialized = false;
		bool dirty = false;
		GLuint vao_fan;
		GLuint shader_fan;
		GLuint vbo_V; // Vertices of the current fan

		Eigen::Matrix3f proj;
		
		// Three points of a fan
		Eigen::Matrix<float, Eigen::Dynamic, 2, Eigen::RowMajor> V3;
		
		// Interpolated temporary copy of the content of each VBO
		Eigen::Matrix<float, Eigen::Dynamic, 2, Eigen::RowMajor> V_vbo;
		
		// Initialize shaders and buffers
		IGL_INLINE void init();
		
		// Release all resorces
		IGL_INLINE void free();
		
		// Create a new set of OpenGL buffer objects
		IGL_INLINE void init_buffers();
		
		// Bind the undrlying OpenGL buffer objects for subsequent marquee draw calls
		IGL_INLINE void bind();
		
		// Draw the fan
		IGL_INLINE void draw();
	
		// Release the OpenGL buffer objects
		IGL_INLINE void free_buffers();
		
	};
  }
}
#ifndef IGL_STATIC_LIBRARY
#include "FanGL.cpp"
#endif
#endif
