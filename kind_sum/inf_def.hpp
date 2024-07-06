/*
The MIT License (MIT)

Copyright (c) 2024 sigma-axis

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <cstdint>
#include <exedit/pixel.hpp>
#include "../buffer_base.hpp"

namespace Calculation::sum
{
	Bounds inflate(int src_w, int src_h,
		i16* src_buf, size_t src_stride,
		i16* dst_buf, bool dst_colored, size_t dst_stride, int a_sum_cap_rate,
		void* heap, int size_sq);
	Bounds inflate(int src_w, int src_h,
		ExEdit::PixelYCA const* src_buf, size_t src_stride,
		i16* dst_buf, bool dst_colored, size_t dst_stride, int a_sum_cap_rate,
		void* heap, int size_sq, void* alpha_space);

	template<int denom>
	constexpr int inflate_radius(int numer) { return numer / denom; }
	// size = floor(size_sq^(1/2)).
	constexpr size_t inflate_heap_size(int dst_w, int dst_h, int size) {
		return sizeof(int8_t) * (dst_w * dst_h) + sizeof(i32) * (2 * size + 1);
	}

	Bounds deflate(int src_w, int src_h,
		i16* src_buf, size_t src_stride,
		i16* dst_buf, bool dst_colored, size_t dst_stride, int a_sum_cap_rate,
		void* heap, int size_sq);
	Bounds deflate(int src_w, int src_h,
		ExEdit::PixelYCA const* src_buf, size_t src_stride,
		i16* dst_buf, bool dst_colored, size_t dst_stride, int a_sum_cap_rate,
		void* heap, int size_sq, void* alpha_space);

	template<int denom>
	constexpr int deflate_radius(int numer) { return std::max(0, numer / denom - 1); }
	// size = floor(size_sq^(1/2)).
	constexpr size_t deflate_heap_size(int src_w, int src_h, int size) {
		return sizeof(int8_t) * (src_w * src_h) + sizeof(i32) * (2 * size + 1);
	}

	constexpr size_t alpha_space_size(int src_w, int src_h) {
		return sizeof(i16) * ((src_w + 1) & (-2)) * src_h;
	}

	size_t constexpr log2_den_cap_rate = 12,
		den_cap_rate = 1 << log2_den_cap_rate;

	inline int a_sum_cap_from_rate(int a_sum_cap_rate, int size_sq) {
		return max_alpha + static_cast<int>((
			(2 * std::sqrt(2 * std::sqrt(size_sq) - 1)) - 1
			) * a_sum_cap_rate * ((1.0 * max_alpha) / den_cap_rate));
	}
}


