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

#include "../multi_thread.hpp"
#include "../arithmetics.hpp"
#include "inf_def.hpp"

using namespace Calculation;

template<size_t a_step>
static inline void pass1(int src_w, int src_h, int size,
	i16 const* a_buf, size_t a_stride, i16 thresh,
	i32* med_buf, size_t med_stride)
{
	multi_thread(src_h, [=](int thread_id, int thread_num) {
		int dst_w = src_w - 2 * size;

		for (int y = thread_id; y < src_h; y += thread_num) {
			auto a_buf_y = a_buf + size * a_step + y * a_stride;
			auto m_buf_y = med_buf + y * med_stride;

			// left -> right
			int count = 0;
			auto a_buf_x = a_buf_y - a_step;
			for (int x = size; --x >= 0; a_buf_x -= a_step) {
				if (*a_buf_x <= thresh) break;
				count++;
			}
			a_buf_x = a_buf_y;
			auto m_buf_x = m_buf_y;
			for (int x = dst_w; --x >= 0; a_buf_x += a_step, m_buf_x++) {
				count++;
				if (*a_buf_x <= thresh) count = 0;
				*m_buf_x = count;
			}

			// left <- right
			count = 0;
			for (int x = size; --x >= 0; a_buf_x += a_step) {
				if (*a_buf_x <= thresh) break;
				count++;
			}
			a_buf_x = a_buf_y + (dst_w - 1) * a_step;
			m_buf_x--;
			for (int x = dst_w; --x >= 0; a_buf_x -= a_step, m_buf_x--) {
				count++;
				if (*a_buf_x <= thresh) count = 0;
				*m_buf_x = std::min(*m_buf_x, count);
			}
		}
	});
}

template<size_t a_step>
static inline void pass2(int src_w, int src_h, int size,
	i32 const* med_buf, size_t med_stride,
	i16* a_buf, size_t a_stride, i32 const* arc, i32* cnt_buf)
{
	auto dst_w = src_w - 2 * size, dst_h = src_h - 2 * size;
	multi_thread(dst_w, [=](int thread_id, int thread_num) {
		int x0 = dst_w * thread_id / thread_num, x1 = dst_w * (thread_id + 1) / thread_num;
		auto m_buf_x0 = med_buf + x0;
		auto a_buf_x0 = a_buf + x0 * a_step;
		auto const count_x0 = cnt_buf + x0;

		auto set_count = [&](int val) {
			auto count = count_x0;
			for (int x = x1 - x0; --x >= 0; count++) *count = { val };
		};

		// top -> bottom
		set_count(size);
		for (int y = size; --y >= 0; m_buf_x0 += med_stride) {
			auto count = count_x0;
			auto m_buf_y = m_buf_x0;
			for (int x = x1 - x0; --x >= 0; count++, m_buf_y++) {
				--*count;
				if (*m_buf_y <= size) *count = std::max(*count, arc[*m_buf_y]);
			}
		}
		for (int y = dst_h; --y >= 0; m_buf_x0 += med_stride, a_buf_x0 += a_stride) {
			auto count = count_x0;
			auto m_buf_y = m_buf_x0;
			auto a_buf_y = a_buf_x0;

			for (int x = x1 - x0; --x >= 0; count++, m_buf_y++, a_buf_y += a_step) {
				--*count;
				if (*m_buf_y <= size) *count = std::max(*count, arc[*m_buf_y]);
				*a_buf_y = *count < 0 ? max_alpha : 0;
			}
		}

		// top <- bottom
		set_count(size);
		m_buf_x0 += size * med_stride; m_buf_x0 -= med_stride; a_buf_x0 -= a_stride;
		for (int y = size; --y >= 0; m_buf_x0 -= med_stride) {
			auto count = count_x0;
			auto m_buf_y = m_buf_x0;
			for (int x = x1 - x0; --x >= 0; count++, m_buf_y++) {
				--*count;
				if (*m_buf_y <= size) *count = std::max(*count, arc[*m_buf_y]);
			}
		}
		for (int y = dst_h; --y >= 0; m_buf_x0 -= med_stride, a_buf_x0 -= a_stride) {
			auto count = count_x0;
			auto m_buf_y = m_buf_x0;
			auto a_buf_y = a_buf_x0;

			for (int x = x1 - x0; --x >= 0; count++, m_buf_y++, a_buf_y += a_step) {
				--*count;
				if (*m_buf_y <= size) *count = std::max(*count, arc[*m_buf_y]);
				if (*count >= 0) *a_buf_y = 0; // std::min(*a_buf_y, 0)
			}
		}
	});
}

Bounds bin::deflate(int src_w, int src_h,
	i16 const* src_buf, bool src_colored, size_t src_stride, i16 thresh,
	i16* dst_buf, bool dst_colored, size_t dst_stride,
	void* heap, int size_sq)
{
	auto* const arc = reinterpret_cast<i32*>(heap);
	int size = arith::arc::quarter(size_sq, arc);
	int const dst_w = src_w - 2 * size;
	auto* cnt_buf = arc + size + 1;
	auto* const med_buf = cnt_buf + dst_w;

	(src_colored ? pass1<4> : pass1<1>)
		(src_w, src_h, size, src_buf, src_stride, thresh, med_buf, dst_w);

	(dst_colored ? pass2<4> : pass2<1>)
		(src_w, src_h, size, med_buf, dst_w, dst_buf, dst_stride, arc, cnt_buf);

	return { 0, 0, src_w - 2 * size, src_h - 2 * size };
}
