/*
The MIT License (MIT)

Copyright (c) 2024 sigma-axis

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <cstdint>
#include <tuple>

#include <exedit/pixel.hpp>
#include "../multi_thread.hpp"
#include "../buffer_base.hpp"

namespace Calculation::masking
{
	enum class mask : uint8_t {
		gray,
		zero, // all pixels nearby are fully transparent.
		full, // all pixels nearby are fully opaque.
	};
	static_assert(sizeof(mask) == 1);
}

namespace Calculation::masking::inflation
{
	inline auto mask_v_color(int src_w, int src_h, int size,
		ExEdit::PixelYCA const* src_buf, size_t src_stride,
		mask* mask_buf, size_t mask_stride,
		i16* a_buf, size_t a_stride)
	{
		// has a second task to copy alpha values to a_buf
		// --- those values are referred so many times in later processes
		//     that it seems to be faster if they are placed within a compact space.
		auto bounds = multi_thread(src_w, [&](int thread_id, int thread_num) {
			int left = src_w, right = -1;

			for (int x = thread_id; x < src_w; x += thread_num) {
				auto s_buf_x = src_buf + x;
				auto m_buf_x = mask_buf + x;
				auto d_buf_x = a_buf + x;

				int cnt_o = 0, cnt_i = 2 * size;
				for (int y = src_h; --y >= 0; s_buf_x += src_stride, m_buf_x += mask_stride, d_buf_x += a_stride) {
					*d_buf_x = s_buf_x->a;

					cnt_o--; cnt_i--;
					if (*d_buf_x > 0) cnt_o = 2 * size; else *d_buf_x = 0;
					if (*d_buf_x < max_alpha) cnt_i = 2 * size; else *d_buf_x = max_alpha;
					*m_buf_x = cnt_o < 0 ? mask::zero :
						cnt_i < 0 ? mask::full : mask::gray;
				}
				for (int y = 2 * size; --y >= 0; m_buf_x += mask_stride) {
					cnt_o--;
					*m_buf_x = cnt_o < 0 ? mask::zero : mask::gray;
				}

				if (cnt_o > -src_h) {
					// some (partially) opaque pixels found.
					if (right < 0) left = right = x; else right = x;
				}
			}

			return std::pair{ left, right };
		});

		// aggregate the returned bounds.
		return unite_interval_alt<int>(bounds);
	}
	inline auto mask_v_alpha(int src_w, int src_h, int size,
		i16* a_buf, size_t a_stride,
		mask* mask_buf, size_t mask_stride)
	{
		auto bounds = multi_thread(src_w, [&](int thread_id, int thread_num) {
			int left = src_w, right = -1;

			for (int x = thread_id; x < src_w; x += thread_num) {
				auto a_buf_x = a_buf + x;
				auto m_buf_x = mask_buf + x;

				int cnt_o = 0, cnt_i = 2 * size;
				for (int y = src_h; --y >= 0; a_buf_x += a_stride, m_buf_x += mask_stride) {
					cnt_o--; cnt_i--;
					if (*a_buf_x > 0) cnt_o = 2 * size; else *a_buf_x = 0;
					if (*a_buf_x < max_alpha) cnt_i = 2 * size; else *a_buf_x = max_alpha;
					*m_buf_x = cnt_o < 0 ? mask::zero :
						cnt_i < 0 ? mask::full : mask::gray;
				}
				for (int y = 2 * size; --y >= 0; m_buf_x += mask_stride) {
					cnt_o--;
					*m_buf_x = cnt_o < 0 ? mask::zero : mask::gray;
				}

				if (cnt_o > -src_h) {
					// some (partially) opaque pixels found.
					if (right < 0) left = right = x; else right = x;
				}
			}

			return std::pair{ left, right };
		});

		// aggregate the returned bounds.
		return unite_interval_alt<int>(bounds);
	}

	inline auto mask_h(int src_w, int src_h, int size,
		mask* mask_buf, size_t mask_stride)
	{
		auto dst_h = src_h + 2 * size;
		auto bounds = multi_thread(dst_h, [&](int thread_id, int thread_num) {
			int top = dst_h, bottom = -1;

			for (int y = thread_id; y < dst_h; y += thread_num) {
				auto m_buf_y = mask_buf + y * mask_stride;

				int cnt_o = 0, cnt_i = 2 * size;
				for (int x = src_w; --x >= 0; m_buf_y++) {
					cnt_o--; cnt_i--;
					if (*m_buf_y != mask::zero) cnt_o = 2 * size;
					if (*m_buf_y != mask::full) cnt_i = 2 * size;
					*m_buf_y = cnt_o < 0 ? mask::zero :
						cnt_i < 0 ? mask::full : mask::gray;
				}
				for (int x = 2 * size; --x >= 0; m_buf_y++) {
					cnt_o--;
					*m_buf_y = cnt_o < 0 ? mask::zero : mask::gray;
				}

				if (cnt_o > -src_w) {
					// some (partially) opaque pixels found.
					if (bottom < 0) top = bottom = y; else bottom = y;
				}
			}

			return std::pair{ top, bottom };
		});

		// aggregate the returned bounds.
		return unite_interval_alt<int>(bounds);
	}
}

namespace Calculation::masking::deflation
{
	// diff_size == size_mask - size_canvas, >= 0 (either 0 or 1)
	template<int diff_size>
	inline auto mask_h_color(int src_w, int src_h, int size_mask,
		ExEdit::PixelYCA const* src_buf, size_t src_stride,
		mask* mask_buf, size_t mask_stride,
		i16* a_buf, size_t a_stride)
	{
		// has a second task to copy alpha values to a_buf
		// --- those values are referred so many times in later processes
		//     that it seems to be faster if they are placed within a compact space.

		using Calculation::max_alpha;

		int const inner_w1 = 2 * size_mask - diff_size,
			inner_w2 = src_w - inner_w1;
		auto bounds = multi_thread(src_h, [&](int thread_id, int thread_num) {
			int top = src_h, bottom = -1;

			for (int y = thread_id; y < src_h; y += thread_num) {
				auto s_buf_y = src_buf + y * src_stride;
				auto m_buf_y = mask_buf + y * mask_stride;
				auto d_buf_y = a_buf + y * a_stride;

				if constexpr (diff_size > 0) {
					for (int x = diff_size; --x >= 0;)
						d_buf_y[-1 - x] = 0;
				}

				int cnt_o = 0, cnt_i = 2 * size_mask;
				for (int x = inner_w1; --x >= 0; s_buf_y++, d_buf_y++) {
					*d_buf_y = s_buf_y->a;

					cnt_o--; cnt_i--;
					if (*d_buf_y > 0) cnt_o = 2 * size_mask; else *d_buf_y = 0;
					if (*d_buf_y < max_alpha) cnt_i = 2 * size_mask; else *d_buf_y = max_alpha;
				}
				for (int x = inner_w2; --x >= 0; s_buf_y++, m_buf_y++, d_buf_y++) {
					*d_buf_y = s_buf_y->a;

					cnt_o--; cnt_i--;
					if (*d_buf_y > 0) cnt_o = 2 * size_mask; else *d_buf_y = 0;
					if (*d_buf_y < max_alpha) cnt_i = 2 * size_mask; else *d_buf_y = max_alpha;
					*m_buf_y = cnt_o < 0 ? mask::zero :
						cnt_i < 0 ? mask::full : mask::gray;
				}
				if (cnt_o > -src_w) {
					// some (partially) opaque pixels found.
					if (bottom < 0) top = bottom = y; else bottom = y;
				}

				if constexpr (diff_size > 0) {
					for (int x = diff_size; --x >= 0; m_buf_y++, d_buf_y++) {
						*d_buf_y = 0;

						cnt_o--;
						*m_buf_y = cnt_o < 0 ? mask::zero : mask::gray;
					}
				}
			}

			return std::pair{ top, bottom };
		});

		// aggregate the returned bounds.
		return unite_interval_alt<int>(bounds);
	}
	// diff_size == size_mask - size_canvas, >= 0 (either 0 or 1)
	template<int diff_size>
	inline auto mask_h_alpha(int src_w, int src_h, int size_mask,
		i16* a_buf, size_t a_stride,
		mask* mask_buf, size_t mask_stride)
	{
		using Calculation::max_alpha;

		int const inner_w1 = 2 * size_mask - diff_size,
			inner_w2 = src_w - inner_w1;
		auto bounds = multi_thread(src_h, [&](int thread_id, int thread_num) {
			int top = src_h, bottom = -1;

			for (int y = thread_id; y < src_h; y += thread_num) {
				auto s_buf_y = a_buf + y * a_stride;
				auto m_buf_y = mask_buf + y * mask_stride;

				if constexpr (diff_size > 0) {
					for (int x = diff_size; --x >= 0;)
						s_buf_y[-1 - x] = 0;
				}

				int cnt_o = 0, cnt_i = 2 * size_mask;
				for (int x = inner_w1; --x >= 0; s_buf_y++) {
					cnt_o--; cnt_i--;
					if (*s_buf_y > 0) cnt_o = 2 * size_mask; else *s_buf_y = 0;
					if (*s_buf_y < max_alpha) cnt_i = 2 * size_mask; else *s_buf_y = max_alpha;
				}
				for (int x = inner_w2; --x >= 0; s_buf_y++, m_buf_y++) {
					cnt_o--; cnt_i--;
					if (*s_buf_y > 0) cnt_o = 2 * size_mask; else *s_buf_y = 0;
					if (*s_buf_y < max_alpha) cnt_i = 2 * size_mask; else *s_buf_y = max_alpha;
					*m_buf_y = cnt_o < 0 ? mask::zero :
						cnt_i < 0 ? mask::full : mask::gray;
				}
				if (cnt_o > -src_w) {
					// some (partially) opaque pixels found.
					if (bottom < 0) top = bottom = y; else bottom = y;
				}

				if constexpr (diff_size > 0) {
					for (int x = diff_size; --x >= 0; s_buf_y++, m_buf_y++) {
						*s_buf_y = 0;

						cnt_o--;
						*m_buf_y = cnt_o < 0 ? mask::zero : mask::gray;
					}
				}
			}

			return std::pair{ top, bottom };
		});

		// aggregate the returned bounds.
		return unite_interval_alt<int>(bounds);
	}

	// diff_size == size_mask - size_canvas, >= 0 (either 0 or 1)
	template<int diff_size>
	inline auto mask_v(int src_w, int src_h, int size_mask,
		mask* mask_buf, size_t mask_stride)
	{
		int const dst_w = src_w - 2 * size_mask + 2 * diff_size,
			inner_h1 = 2 * size_mask - diff_size,
			inner_h2 = src_h - inner_h1;
		auto bounds = multi_thread(dst_w, [&](int thread_id, int thread_num) {
			int left = dst_w, right = -1;

			for (int x = thread_id; x < dst_w; x += thread_num) {
				auto m_buf_s = mask_buf + x, m_buf_d = m_buf_s;

				int cnt_o = 0, cnt_i = 2 * size_mask;
				for (int y = inner_h1; --y >= 0; m_buf_s += mask_stride) {
					cnt_o--; cnt_i--;
					if (*m_buf_s != mask::zero) cnt_o = 2 * size_mask;
					if (*m_buf_s != mask::full) cnt_i = 2 * size_mask;
				}
				for (int y = inner_h2; --y >= 0; m_buf_s += mask_stride, m_buf_d += mask_stride) {
					cnt_o--; cnt_i--;
					if (*m_buf_s != mask::zero) cnt_o = 2 * size_mask;
					if (*m_buf_s != mask::full) cnt_i = 2 * size_mask;
					*m_buf_d = cnt_o < 0 ? mask::zero :
						cnt_i < 0 ? mask::full : mask::gray;
				}
				if (cnt_o > -src_h) {
					// some (partially) opaque pixels found.
					if (right < 0) left = right = x; else right = x;
				}

				if constexpr (diff_size > 0) {
					for (int y = diff_size; --y >= 0; m_buf_d += mask_stride) {
						cnt_o--;
						*m_buf_d = cnt_o < 0 ? mask::zero : mask::gray;
					}
				}
			}

			return std::pair{ left, right };
		});

		// aggregate the returned bounds.
		return unite_interval_alt<int>(bounds);
	}
}
