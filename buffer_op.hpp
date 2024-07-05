/*
The MIT License (MIT)

Copyright (c) 2024 sigma-axis

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <cstdint>
#include <cmath>

#include <exedit/pixel.hpp>
#include "buffer_base.hpp"


////////////////////////////////
// よくあるバッファ操作の宣言．
////////////////////////////////
namespace Calculation::buff
{
	void copy_alpha(ExEdit::PixelYCA const* src, size_t src_stride, int src_x, int src_y, int src_w, int src_h,
		ExEdit::PixelYCA* dst, size_t dst_stride, int dst_x, int dst_y);
	void copy_alpha(i16 const* a_src, size_t a_stride, int src_x, int src_y, int src_w, int src_h,
		ExEdit::PixelYCA* dst, size_t dst_stride, int dst_x, int dst_y);
	void copy_alpha(ExEdit::PixelYCA const* src, size_t src_stride, int src_x, int src_y, int src_w, int src_h,
		i16* a_dst, size_t a_stride, int dst_x, int dst_y);

	void clear_alpha(ExEdit::PixelYCA* dst, size_t dst_stride, int x, int y, int w, int h);
	void clear_alpha(i16* a_dst, size_t a_stride, int x, int y, int w, int h, bool strict = false);

	void clear_alpha_chrome(ExEdit::PixelYCA* dst, size_t dst_stride, Bounds const& outer, Bounds const& inner);
	void clear_alpha_chrome(i16* a_dst, size_t a_stride, Bounds const& outer, Bounds const& inner);

	void mult_alpha(i16 const* a_src, size_t a_stride, int src_x, int src_y, int src_w, int src_h,
		ExEdit::PixelYCA* dst, size_t dst_stride, int dst_x, int dst_y);
	void mult_alpha(i16 alpha, ExEdit::PixelYCA* dst, size_t dst_stride, int x, int y, int w, int h);
	void mult_alpha(i16 alpha, i16* a_dst, size_t a_stride, int x, int y, int w, int h);

	void binarize(ExEdit::PixelYCA const* src, size_t src_stride, int src_x, int src_y, int src_w, int src_h,
		ExEdit::PixelYCA* dst, size_t dst_stride, int dst_x, int dst_y, i16 thresh);
	void binarize(ExEdit::PixelYCA const* src, size_t src_stride, int src_x, int src_y, int src_w, int src_h,
		i16* a_dst, size_t a_stride, int dst_x, int dst_y, i16 thresh);
	void binarize(i16 const* a_src, size_t a_stride, int src_x, int src_y, int src_w, int src_h,
		ExEdit::PixelYCA* dst, size_t dst_stride, int dst_x, int dst_y, i16 thresh);

	constexpr size_t log2_den_blur_px = 12,
		den_blur_px = 1 << log2_den_blur_px;
	// returns the inflation size of each side.
	void blur_alpha(ExEdit::PixelYCA* dst, size_t dst_stride,
		int dst_x, int dst_y, int dst_w, int dst_h, int blur_px);
	// returns the inflation size of each side.
	void blur_alpha(i16* a_dst, size_t a_stride,
		int dst_x, int dst_y, int dst_w, int dst_h, int blur_px);
	constexpr int blur_displace(int blur_px) {
		return ((blur_px - 1) >> (log2_den_blur_px + 1)) + 1;
	}

	// ripped a piece of code from exedit/pixel.hpp.
	constexpr ExEdit::PixelYC fromRGB(uint8_t r, uint8_t g, uint8_t b) {
		auto r_ = (r << 6) + 18;
		auto g_ = (g << 6) + 18;
		auto b_ = (b << 6) + 18;
		return {
			static_cast<int16_t>(((r_* 4918)>>16)+((g_* 9655)>>16)+((b_* 1875)>>16)-3),
			static_cast<int16_t>(((r_*-2775)>>16)+((g_*-5449)>>16)+((b_* 8224)>>16)+1),
			static_cast<int16_t>(((r_* 8224)>>16)+((g_*-6887)>>16)+((b_*-1337)>>16)+1),
		};
	}

	constexpr ExEdit::PixelYCA* alpha_to_pixel(i16* a_ptr) {
		constexpr auto ofs = offsetof(ExEdit::PixelYCA, a) / sizeof(i16);
		return (ExEdit::PixelYCA*)(a_ptr - ofs);
	}
	constexpr ExEdit::PixelYCA const* alpha_to_pixel(i16 const* a_ptr) {
		constexpr auto ofs = offsetof(ExEdit::PixelYCA, a) / sizeof(i16);
		return (ExEdit::PixelYCA const*)(a_ptr - ofs);
	}
	constexpr i16* pixel_to_alpha(ExEdit::PixelYCA* p_ptr) {
		return &p_ptr->a;
	}
	constexpr i16 const* pixel_to_alpha(ExEdit::PixelYCA const* p_ptr) {
		return &p_ptr->a;
	}
}

