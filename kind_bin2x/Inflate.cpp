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
	upper = 0,
	equiv = 1,
	lower = 2,
};
struct med_data {
	i16 d = 0;
	flg f = flg::equiv;

	med_data(int d, flg f) : d{ static_cast<i16>(d) }, f{ f } {}
	// assuming that top -> bottom goes before bottom -> top.
	void push(int d2) {
		if (d >= d2) {
			f = d > d2 ? flg::lower : flg::equiv;
			d = d2;
		}
	}
};
static_assert(sizeof(med_data) == sizeof(i32));

template<size_t a_step>
static inline auto pass1(int src_w, int src_h, int size,
	i16 const* a_buf, size_t a_stride, i16 thresh,
	med_data* med_buf, size_t med_stride)
{
	auto const bounds = multi_thread(src_w, [=](int thread_id, int thread_num) {
		int x0 = src_w * thread_id / thread_num, x1 = src_w * (thread_id + 1) / thread_num;
		auto a_buf_x0 = a_buf + x0 * a_step;
		auto m_buf_x0 = med_buf + x0;
		auto const count_x0 = reinterpret_cast<i32*>(m_buf_x0);

		auto set_count = [&](i32 val) {
			auto count = count_x0;
			for (int x = x1 - x0; --x >= 0; count++) *count = val;
		};

		// top -> bottom
		set_count(size);
		m_buf_x0 += size * med_stride;
		for (int y = src_h; --y >= 0; a_buf_x0 += a_stride, m_buf_x0 += med_stride) {
			auto count = count_x0;
			auto a_buf_y = a_buf_x0;
			auto m_buf_y = m_buf_x0;

			for (int x = x1 - x0; --x >= 0; count++, a_buf_y += a_step, m_buf_y++) {
				++*count;
				if (*a_buf_y > thresh) *count = 0;
				*m_buf_y = { *count, flg::upper };
			}
		}
		for (int y = size; --y >= 0; m_buf_x0 += med_stride) {
			auto count = count_x0;
			auto m_buf_y = m_buf_x0;

			for (int x = x1 - x0; --x >= 0; count++, m_buf_y++) {
				++*count;
				*m_buf_y = { *count, flg::upper };
			}
		}

		int left = src_w, right = -1;
		{
			int x = x0;
			for (auto count = count_x0; x < x1; x++, count++) {
				if (*count >= src_h + 2 * size) continue;
				if (right < 0) left = right = x; else right = x;
			}
		}
		if (left > right) {
			// no opaque pixels were found.
			m_buf_x0 -= med_stride; m_buf_x0 -= (size + src_h) * med_stride;
		}
		else {
			// opaque pixels exist.
			set_count(size);
			a_buf_x0 -= a_stride; m_buf_x0 -= med_stride; m_buf_x0 -= size * med_stride;

			// top <- bottom
			for (int y = src_h; --y >= 0; a_buf_x0 -= a_stride, m_buf_x0 -= med_stride) {
				auto count = count_x0;
				auto a_buf_y = a_buf_x0;
				auto m_buf_y = m_buf_x0;

				for (int x = x1 - x0; --x >= 0; count++, a_buf_y += a_step, m_buf_y++) {
					++*count;
					if (*a_buf_y > thresh) *count = 0;
					m_buf_y->push(*count);
				}
			}
		}
		for (int y = size; --y >= 0; m_buf_x0 -= med_stride) {
			auto count = count_x0;
			auto m_buf_y = m_buf_x0;

			for (int x = x1 - x0; --x >= 0; count++, m_buf_y++) {
				++*count;
				*m_buf_y = { *count, flg::lower };
			}
		}

		return std::pair{ left, right };
	});

	// aggregate the returned bounds.
	return unite_interval_alt<int>(bounds);
}

template<size_t a_step>
static inline auto pass2(int src_w, int src_h, int size,
	med_data const* med_buf, size_t med_stride,
	i16* a_buf, size_t a_stride, i32 const* arc)
{
	struct fill_count {
		int u, l;
		fill_count(int d) : u{ d }, l{ d } {}

		void operator--() { u -= 2; l -= 2; }
		void operator--(int) { --(*this); }

		void update(med_data m, int32_t const* arc) {
			int d1 = arc[2 * m.d], d2 = arc[2 * m.d + 1];
			switch (m.f) {
			case flg::upper: break;
			case flg::lower: std::swap(d1, d2); break;
			case flg::equiv: default: d2 = d1; break;
			}
			u = std::max(u, d1); l = std::max(l, d2);
		}

		i16 alpha() const {
			return (4 + std::clamp(u, -2, 0) + std::clamp(l, -2, 0)) << (log2_max_alpha - 2);
		}
	};

	int dst_h = src_h + 2 * size;
	auto const bounds = multi_thread(dst_h, [=](int thread_id, int thread_num) {
		int top = dst_h, bottom = -1;

		for (int y = thread_id; y < dst_h; y += thread_num) {
			auto m_buf_y = med_buf + y * med_stride;
			auto a_buf_y = a_buf + y * a_stride;

			// left -> right
			fill_count count{ 0 };
			auto m_buf_x = m_buf_y;
			auto a_buf_x = a_buf_y + size * a_step;
			for (int x = src_w; --x >= 0; m_buf_x++, a_buf_x += a_step) {
				count--;
				if (m_buf_x->d <= size) count.update(*m_buf_x, arc);
				*a_buf_x = count.alpha();
			}
			for (int x = size; --x >= 0; a_buf_x += a_step) {
				count--;
				*a_buf_x = count.alpha();
			}

			if (count.u <= -2 * (src_w + size) && count.l <= -2 * (src_w + size)) {
				// no opaque pixel were written.
				a_buf_x = a_buf_y + size * a_step; a_buf_x -= a_step;
			}
			else {
				// some opaque pixels were written. update the bounds.
				if (bottom < 0) top = bottom = y; else bottom = y;

				// left <- right
				count = { 0 };
				m_buf_x--; a_buf_x -= a_step; a_buf_x -= size * a_step;
				for (int x = src_w; --x >= 0; m_buf_x--, a_buf_x -= a_step) {
					count--;
					if (m_buf_x->d <= size) count.update(*m_buf_x, arc);
					*a_buf_x = std::max(*a_buf_x, count.alpha());
				}
			}
			for (int x = size; --x >= 0; a_buf_x -= a_step) {
				count--;
				*a_buf_x = count.alpha();
			}
		}

		return std::pair{ top, bottom };
	});

	// aggregate the returned bounds.
	return unite_interval_alt<int>(bounds);
}


Bounds bin2x::inflate(int src_w, int src_h,
	i16 const* src_buf, bool src_colored, size_t src_stride, i16 thresh,
	i16* dst_buf, bool dst_colored, size_t dst_stride,
	void* heap, int size2_sq)
{
	auto* const arc = reinterpret_cast<i32*>(heap);
	arc[0] = arith::arc::quarter(size2_sq, &arc[1]);
	arc[arc[0] + 2] = -2;
	int const size = (arc[0] + 1) >> 1;

	auto* med_buf = reinterpret_cast<med_data*>(arc + (arc[0] + 3));
	size_t med_stride = src_w + 2 * size;

	auto [left, right] = (src_colored ? pass1<4> : pass1<1>)
		(src_w, src_h, size, src_buf, src_stride, thresh, med_buf, med_stride);

	if (left >= right) return { 0,0,0,0 };

	med_buf += left; dst_buf += left * (dst_colored ? 4 : 1);
	src_w = right - left;

	auto [top, bottom] = (dst_colored ? pass2<4> : pass2<1>)
		(src_w, src_h, size, med_buf, med_stride, dst_buf, dst_stride, arc);

	right += 2 * size;
	return { left, top, right, bottom };
}
