/*
The MIT License (MIT)

Copyright (c) 2024 sigma-axis

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <cstdint>
#include "../buffer_base.hpp"

namespace Calculation::bin2x
{
	Bounds inflate(int src_w, int src_h,
		i16 const* src_buf, bool src_colored, size_t src_stride, i16 thresh,
		i16* dst_buf, bool dst_colored, size_t dst_stride,
		void* heap, int size2_sq);

	template<int denom>
	constexpr int inflate_radius(int numer) { return (numer + (denom >> 1)) / denom; }
	// size = floor((size2_sq^(1/2) + 1)/2), dst_(w/h) = src_(w/h) + 2*size.
	constexpr size_t inflate_heap_size(int dst_w, int dst_h, int size) {
		return sizeof(i32) * (dst_w * dst_h + 2 * size + 3);
	}

	Bounds deflate(int src_w, int src_h,
		i16 const* src_buf, bool src_colored, size_t src_stride, i16 thresh,
		i16* dst_buf, bool dst_colored, size_t dst_stride,
		void* heap, int size2_sq);

	template<int denom>
	constexpr int deflate_radius(int numer) { return numer / denom; }
	// size = floor(size2_sq^(1/2)/2).
	constexpr size_t deflate_heap_size(int src_w, int src_h, int size) {
		return sizeof(i32) * ((src_w - 2 * size) * (src_h + 2) + 2 * size + 4);
	}
}

