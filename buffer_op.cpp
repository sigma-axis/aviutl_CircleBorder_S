/*
The MIT License (MIT)

Copyright (c) 2024 sigma-axis

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <cstdint>
#include <algorithm>
#include <tuple>
#include <cmath>
#include <bit>

#include "exedit/pixel.hpp"
#include "multi_thread.hpp"
#include "buffer_op.hpp"

using namespace Calculation;


////////////////////////////////
// よくあるバッファ操作の実装．
////////////////////////////////
void buff::copy_alpha(ExEdit::PixelYCA const* src, size_t src_stride, int src_x, int src_y, int src_w, int src_h,
	ExEdit::PixelYCA* dst, size_t dst_stride, int dst_x, int dst_y)
{
	if (src_w <= 0 || src_h <= 0) return;

	src += src_x + src_y * src_stride;
	dst += dst_x + dst_y * dst_stride;
	multi_thread(src_h, [=](int thread_id, int thread_num) {
		for (int y = thread_id * src_h / thread_num, y1 = (thread_id + 1) * src_h / thread_num;
			y < y1; y++) {
			auto src_y = src + y * src_stride;
			auto dst_y = dst + y * dst_stride;
			for (int x = src_w; --x >= 0; src_y++, dst_y++) dst_y->a = src_y->a;
		}
	});
}

void buff::copy_alpha(i16 const* a_src, size_t a_stride, int src_x, int src_y, int src_w, int src_h,
	ExEdit::PixelYCA* dst, size_t dst_stride, int dst_x, int dst_y)
{
	if (src_w <= 0 || src_h <= 0) return;

	a_src += src_x + src_y * a_stride;
	dst += dst_x + dst_y * dst_stride;
	multi_thread(src_h, [=](int thread_id, int thread_num) {
		for (int y = thread_id * src_h / thread_num, y1 = (thread_id + 1) * src_h / thread_num;
			y < y1; y++) {
			auto src_y = a_src + y * a_stride;
			auto dst_y = dst + y * dst_stride;
			for (int x = src_w; --x >= 0; src_y++, dst_y++) dst_y->a = *src_y;
		}
	});
}

void buff::copy_alpha(ExEdit::PixelYCA const* src, size_t src_stride, int src_x, int src_y, int src_w, int src_h,
	i16* a_dst, size_t a_stride, int dst_x, int dst_y)
{
	if (src_w <= 0 || src_h <= 0) return;

	src += src_x + src_y * src_stride;
	a_dst += dst_x + dst_y * a_stride;
	multi_thread(src_h, [=](int thread_id, int thread_num) {
		for (int y = thread_id * src_h / thread_num, y1 = (thread_id + 1) * src_h / thread_num;
			y < y1; y++) {
			auto src_y = src + y * src_stride;
			auto dst_y = a_dst + y * a_stride;
			for (int x = src_w; --x >= 0; src_y++, dst_y++) *dst_y = src_y->a;
		}
	});
}

void buff::clear_alpha(ExEdit::PixelYCA* dst, size_t dst_stride, int x, int y, int w, int h)
{
	if (w <= 0 || h <= 0) return;

	dst += x + y * dst_stride;
	multi_thread(std::min(w, h), [=](int thread_id, int thread_num) {
		for (int y = h * thread_id / thread_num, y1 = h * (thread_id + 1) / thread_num;
			y < y1; y++) {
			auto dst_y = dst + y * dst_stride;
			for (int x = w; --x >= 0; dst_y++) dst_y->a = 0;
		}
	});
}

void buff::clear_alpha(i16* dst, size_t a_stride, int x, int y, int w, int h, bool strict)
{
	if (w <= 0 || h <= 0) return;

	dst += x + y * a_stride;
	if (!strict && a_stride <= static_cast<size_t>(w + (w >> 1)))
		// if the stride is small enough w.r.t. the width, use memset() instead.
		std::memset(dst, 0, sizeof(*dst) * (w + (h - 1) * a_stride));
	else {
		multi_thread(std::min(w, h), [=](int thread_id, int thread_num) {
			for (int y = h * thread_id / thread_num, y1 = h * (thread_id + 1) / thread_num;
				y < y1; y++)
				std::memset(dst + y * a_stride, 0, sizeof(*dst) * w);
		});
	}
}

void buff::clear_alpha_chrome(ExEdit::PixelYCA* dst, size_t dst_stride, Bounds const& outer, Bounds const& inner)
{
	if (outer.T < inner.T)
		clear_alpha(dst, dst_stride, outer.L, outer.T, outer.wd(), inner.T - outer.T);
	if (inner.B < outer.B)
		clear_alpha(dst, dst_stride, outer.L, inner.B, outer.wd(), outer.B - inner.B);
	auto t = std::max(outer.T, inner.T), h = std::min(outer.B, inner.B) - t;
	if (outer.L < inner.L)
		clear_alpha(dst, dst_stride, outer.L, t, inner.L - outer.L, h);
	if (inner.R < outer.R)
		clear_alpha(dst, dst_stride, inner.R, t, outer.R - inner.R, h);
}

void buff::clear_alpha_chrome(i16* a_dst, size_t a_stride, Bounds const& outer, Bounds const& inner)
{
	if (outer.T < inner.T)
		clear_alpha(a_dst, a_stride, outer.L, outer.T, outer.wd(), inner.T - outer.T);
	if (inner.B < outer.B)
		clear_alpha(a_dst, a_stride, outer.L, inner.B, outer.wd(), outer.B - inner.B);
	auto t = std::max(outer.T, inner.T), h = std::min(outer.B, inner.B) - t;
	if (outer.L < inner.L)
		clear_alpha(a_dst, a_stride, outer.L, t, inner.L - outer.L, h, true);
	if (inner.R < outer.R)
		clear_alpha(a_dst, a_stride, inner.R, t, outer.R - inner.R, h, true);
}

void buff::mult_alpha(i16 const* a_src, size_t a_stride, int src_x, int src_y, int src_w, int src_h,
	ExEdit::PixelYCA* dst, size_t dst_stride, int dst_x, int dst_y)
{
	a_src += src_x + src_y * a_stride;
	dst += dst_x + dst_y * dst_stride;
	multi_thread(src_h, [=](int thread_id, int thread_num) {
		for (int y = thread_id * src_h / thread_num, y1 = (thread_id + 1) * src_h / thread_num;
			y < y1; y++) {
			auto src_y = a_src + y * a_stride;
			auto dst_y = dst + y * dst_stride;
			for (int x = src_w; --x >= 0; src_y++, dst_y++)
				dst_y->a = (*src_y * dst_y->a) >> log2_max_alpha;
		}
	});
}

void buff::mult_alpha(i16 alpha, ExEdit::PixelYCA* dst, size_t dst_stride, int x, int y, int w, int h)
{
	if (w <= 0 || h <= 0) return;

	dst += x + y * dst_stride;
	multi_thread(h, [=](int thread_id, int thread_num) {
		for (int y = thread_id * h / thread_num, y1 = (thread_id + 1) * h / thread_num;
			y < y1; y++) {
			auto dst_y = dst + y * dst_stride;
			for (int x = w; --x >= 0; dst_y++)
				dst_y->a = (alpha * dst_y->a) >> log2_max_alpha;
		}
	});
}

void buff::mult_alpha(i16 alpha, i16* a_dst, size_t a_stride, int x, int y, int w, int h)
{
	if (w <= 0 || h <= 0) return;

	a_dst += x + y * a_stride;
	multi_thread(h, [=](int thread_id, int thread_num) {
		for (int y = thread_id * h / thread_num, y1 = (thread_id + 1) * h / thread_num;
			y < y1; y++) {
			auto dst_y = a_dst + y * a_stride;
			for (int x = w; --x >= 0; dst_y++)
				*dst_y = (alpha * (*dst_y)) >> log2_max_alpha;
		}
	});
}

void buff::binarize(ExEdit::PixelYCA const* src, size_t src_stride, int src_x, int src_y, int src_w, int src_h,
	ExEdit::PixelYCA* dst, size_t dst_stride, int dst_x, int dst_y, i16 thresh)
{
	if (src_w <= 0 || src_h <= 0) return;

	multi_thread(src_h, [=](int thread_id, int thread_num) {
		for (int y = src_h * thread_id / thread_num, y1 = src_h * (thread_id + 1) / thread_num;
			y < y1; y++) {
			auto src1 = src + src_x + (src_y + y) * src_stride;
			auto dst1 = dst + dst_x + (dst_y + y) * dst_stride;
			for (int x = src_w; --x >= 0; src1++, dst1++)
				dst1->a = src1->a > thresh ? max_alpha : 0;
		}
	});
}

void buff::binarize(ExEdit::PixelYCA const* src, size_t src_stride, int src_x, int src_y, int src_w, int src_h,
	i16* a_dst, size_t a_stride, int dst_x, int dst_y, i16 thresh)
{
	if (src_w <= 0 || src_h <= 0) return;

	multi_thread(src_h, [=](int thread_id, int thread_num) {
		for (int y = src_h * thread_id / thread_num, y1 = src_h * (thread_id + 1) / thread_num;
			y < y1; y++) {
			auto src1 = src + src_x + (src_y + y) * src_stride;
			auto dst1 = a_dst + dst_x + (dst_y + y) * a_stride;
			for (int x = src_w; --x >= 0; src1++, dst1++)
				*dst1 = src1->a > thresh ? max_alpha : 0;
		}
	});
}

void buff::binarize(i16 const* a_src, size_t a_stride, int src_x, int src_y, int src_w, int src_h,
	ExEdit::PixelYCA* dst, size_t dst_stride, int dst_x, int dst_y, i16 thresh)
{
	if (src_w <= 0 || src_h <= 0) return;

	multi_thread(src_h, [=](int thread_id, int thread_num) {
		for (int y = src_h * thread_id / thread_num, y1 = src_h * (thread_id + 1) / thread_num;
			y < y1; y++) {
			auto src1 = a_src + src_x + (src_y + y) * a_stride;
			auto dst1 = dst + dst_x + (dst_y + y) * dst_stride;
			for (int x = src_w; --x >= 0; src1++, dst1++)
				dst1->a = *src1 > thresh ? max_alpha : 0;
		}
	});
}


template<size_t a_step>
static inline void blur_alpha_core(i16* a_dst, size_t a_stride, int w, int h, int blur_px)
{
	using namespace buff;
	if (w <= 0 || h <= 0 || blur_px <= 0) return;

	int denom_len = std::bit_width(static_cast<uint32_t>(blur_px + den_blur_px)) + log2_max_alpha;
	uint32_t numer = (1ULL << denom_len) / (blur_px + den_blur_px);
	int denom_len2 = std::max<int>(0, denom_len + log2_max_alpha - log2_den_blur_px - 31);
	denom_len -= denom_len2 + log2_den_blur_px;
	// final alpha will be calculated as: (((weighted sum of alpha) >> denom_len2) * numer) >> denom_len.

	int size_i = (blur_px - 1) >> log2_den_blur_px;
	int size_f = (((blur_px - 1) & (den_blur_px - 1)) + 2) >> 1; // from 1 to den_blur_px/2 (inclusive).
	if ((size_i & 1) != 0) size_f += den_blur_px >> 1;
	size_i = (size_i & (-2)) + 1; // length of fully weighted pixels.
	// weight for partially weighted pixels ~ size_f/(den_blur_px).

	auto calc_alpha = [&](uint32_t sum, uint32_t edge_sum) {
		sum += (edge_sum * size_f) >> log2_den_blur_px;
		return ((sum >> denom_len2) * numer) >> denom_len;
	};

	int displace = (size_i + 1) >> 1; // inflated lengths on the four sides.

	// perform vertical convolution.
	auto const D = 2 * displace;
	multi_thread(w, [&, hh = std::min(h, D), H = h - D](int thread_id, int thread_num) {
		for (int x = w * thread_id / thread_num, x1 = w * (thread_id + 1) / thread_num;
			x < x1; x++) {
			auto a_y0 = a_dst + x * a_step + (H + D - 1) * a_stride,
				a_y1 = a_y0 + D * a_stride;
			uint32_t sum = 0;

			for (int y = hh; --y >= 0; a_y1 -= a_stride, a_y0 -= a_stride) {
				*a_y1 = calc_alpha(sum, *a_y0);
				sum += *a_y0;
			}
			if (H >= 0) {
				for (int y = H; --y >= 0; a_y1 -= a_stride, a_y0 -= a_stride) {
					sum -= *a_y1;
					*a_y1 = calc_alpha(sum, *a_y0 + *a_y1);
					sum += *a_y0;
				}
			}
			else {
				auto a = calc_alpha(sum, 0);
				for (int y = -H; --y >= 0; a_y1 -= a_stride) *a_y1 = a;
			}
			for (int y = hh; --y >= 0; a_y1 -= a_stride) {
				sum -= *a_y1;
				*a_y1 = calc_alpha(sum, *a_y1);
			}
		}
	});

	// perform horizontal convolution.
	multi_thread(h + D, [&, ww = std::min(w, D), W = w - D](int thread_id, int thread_num) {
		for (int y = (h + D) * thread_id / thread_num, y1 = (h + D) * (thread_id + 1) / thread_num;
			y < y1; y++) {
			auto a_x0 = a_dst + (W + D - 1) * a_step + y * a_stride,
				a_x1 = a_x0 + D * a_step;
			uint32_t sum = 0;

			for (int x = ww; --x >= 0; a_x1 -= a_step, a_x0 -= a_step) {
				*a_x1 = std::clamp<i16>(calc_alpha(sum, *a_x0), 0, max_alpha);
				sum += *a_x0;
			}
			if (W >= 0) {
				for (int x = W; --x >= 0; a_x1 -= a_step, a_x0 -= a_step) {
					sum -= *a_x1;
					*a_x1 = std::clamp<i16>(calc_alpha(sum, *a_x0 + *a_x1), 0, max_alpha);
					sum += *a_x0;
				}
			}
			else {
				auto a = std::clamp<i16>(calc_alpha(sum, 0), 0, max_alpha);
				for (int x = -W; --x >= 0; a_x1 -= a_step) *a_x1 = a;
			}
			for (int x = ww; --x >= 0; a_x1 -= a_step) {
				sum -= *a_x1;
				*a_x1 = std::clamp<i16>(calc_alpha(sum, *a_x1), 0, max_alpha);
			}
		}
	});

	return;
}

void buff::blur_alpha(ExEdit::PixelYCA* dst, size_t dst_stride,
	int dst_x, int dst_y, int dst_w, int dst_h, int blur_px)
{
	blur_alpha_core<4>(&dst[dst_x + dst_y * dst_stride].a, 4 * dst_stride, dst_w, dst_h, blur_px);
}

void buff::blur_alpha(i16* a_dst, size_t a_stride,
	int dst_x, int dst_y, int dst_w, int dst_h, int blur_px)
{
	blur_alpha_core<1>(a_dst + dst_x + dst_y * a_stride, a_stride, dst_w, dst_h, blur_px);
}


