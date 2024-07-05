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

enum class flg : i16 {
	left = 0,
	equiv = 1,
	right = 2,
};
struct med_data {
	i16 d = 0;
	flg f = flg::equiv;

	med_data(int d, flg f) : d{ static_cast<i16>(d) }, f{ f } {}
	// assuming that left -> right goes before right -> left.
	void push(int d2) {
		if (d >= d2) {
			f = d > d2 ? flg::right : flg::equiv;
			d = d2;
		}
	}
};
static_assert(sizeof(med_data) == sizeof(i32));

template<size_t a_step>
static inline void pass1(int src_w, int src_h, int size, int size1,
	i16 const* a_buf, size_t a_stride, i16 thresh,
	med_data* med_buf, size_t med_stride)
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
				*m_buf_x = { count, flg::left };
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
				m_buf_x->push(count);
			}
		}
	});
}

template<size_t a_step>
static inline void pass2(int src_w, int src_h, int size, int size1,
	med_data const* med_buf, size_t med_stride,
	i16* a_buf, size_t a_stride, i32 const* arc)
{
	struct fill_count {
		int l, r;
		fill_count(int d) : l{ d }, r{ d } {}

		void operator--() { l -= 2; r -= 2; }
		void operator--(int) { --(*this); }

		void update(med_data const& m, int32_t const* arc) {
			int d1 = arc[2 * m.d], d2 = arc[2 * m.d + 1];
			switch (m.f) {
			case flg::left: break;
			case flg::right: std::swap(d1, d2); break;
			case flg::equiv: default: d2 = d1; break;
			}
			l = std::max(l, d1); r = std::max(r, d2);
		}

		i16 alpha() const {
			return (-std::clamp(l, -2, 0) - std::clamp(r, -2, 0)) << (log2_max_alpha - 2);
		}
	};

	auto dst_w = src_w - 2 * size, dst_h = src_h - 2 * size;
	multi_thread(dst_w, [=](int thread_id, int thread_num) {
		for (int x = thread_id; x < dst_w; x += thread_num) {
			auto m_buf_x = med_buf + x;
			auto a_buf_x = a_buf + x * a_step;

			// top -> bottom
			fill_count count{ arc[0] };
			auto m_buf_y = m_buf_x;
			auto a_buf_y = a_buf_x;
			for (int y = size; --y >= 0; m_buf_y += med_stride) {
				count--;
				if (m_buf_y->d <= size1) count.update(*m_buf_y, arc);
			}
			for (int y = dst_h; --y >= 0; m_buf_y += med_stride, a_buf_y += a_stride) {
				count--;
				if (m_buf_y->d <= size1) count.update(*m_buf_y, arc);
				*a_buf_y = count.alpha();
			}

			// top <- bottom
			count = { arc[0] };
			m_buf_y += size * med_stride; m_buf_y -= med_stride; a_buf_y -= a_stride;
			for (int y = size; --y >= 0; m_buf_y -= med_stride) {
				count--;
				if (m_buf_y->d <= size1) count.update(*m_buf_y, arc);
			}
			for (int y = dst_h; --y >= 0; m_buf_y -= med_stride, a_buf_y -= a_stride) {
				count--;
				if (m_buf_y->d <= size1) count.update(*m_buf_y, arc);
				*a_buf_y = std::min(*a_buf_y, count.alpha());
			}
		}
	});
}


Bounds bin2x::deflate(int src_w, int src_h,
	i16 const* src_buf, bool src_colored, size_t src_stride, i16 thresh,
	i16* dst_buf, bool dst_colored, size_t dst_stride,
	void* heap, int size2_sq)
{
	auto* const arc = reinterpret_cast<i32*>(heap);
	arc[0] = arith::arc::quarter(size2_sq, &arc[1]);
	arc[arc[0] + 2] = -2;
	int const size = arc[0] >> 1,
		size1 = (arc[0] + 1) >> 1; // the length upto which searching should reach.
	auto* const med_buf = reinterpret_cast<med_data*>(arc + (arc[0] + 3));

	(src_colored ? pass1<4> : pass1<1>)
		(src_w, src_h, size, size1, src_buf, src_stride, thresh, med_buf, src_w);

	(dst_colored ? pass2<4> : pass2<1>)
		(src_w, src_h, size, size1, med_buf, src_w, dst_buf, dst_stride, arc);

	return { 0, 0, src_w - 2 * size, src_h - 2 * size };
}
