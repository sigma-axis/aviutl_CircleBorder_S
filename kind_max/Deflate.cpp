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
#include <numeric>

#include "../multi_thread.hpp"
#include "../arithmetics.hpp"
#include "inf_def.hpp"
#include "masking.hpp"

using namespace Calculation;
using mask = masking::mask;

template<size_t src_step, size_t a_step>
static inline void find_min(int src_w, int src_h, int size,
	i16 const* src_buf, size_t src_stride,
	mask const* mask_buf, size_t mask_stride,
	i16* a_buf, size_t a_stride, i32 const* arc)
{
	// arc[i]: i ranges from -size to size.

	int const dst_w = src_w - 2 * size, dst_h = src_h - 2 * size,
		disk_area = 1 + 4 * (size + std::accumulate(arc + 1, arc + size + 1, 0));
#pragma warning(suppress : 6262) // allocating > 16 KiB on stack.
	multi_thread(dst_h, [&](int thread_id, int thread_num)
	{
		// the buckets that count pixels at each alpha value (except alpha == full).
		uint32_t bucket[max_alpha + 1]{}; // bucket[max_alpha] is simply ignored.
		int curr_min = max_alpha;
		auto add = [&](int alpha) {
			alpha = std::max(alpha, 0);
			bucket[alpha]++;
			curr_min = std::min(curr_min, alpha);
		};
		auto pop = [&](int alpha) {
			alpha = std::max(alpha, 0);
			bucket[alpha]--;
		};
		auto update_min = [&] {
			while (curr_min < max_alpha && bucket[curr_min] == 0) curr_min++;
		};

		int const y0 = dst_h * thread_id / thread_num, y1 = dst_h * (thread_id + 1) / thread_num;

		// first state of buckets.
		switch (mask_buf[y0 * mask_stride]) {
		case mask::zero: bucket[0] = disk_area - 2 * size - 1; curr_min = 0; break;
		case mask::full: curr_min = max_alpha; break;
		case mask::gray:
		default:
			for (int dy = -size; dy <= size; dy++) {
				int secant = arc[dy];
				for (int dx = -secant; dx <= secant; dx++)
					add(src_buf[(dx + size) * src_step + (y0 + dy + size) * src_stride]);
			}
			break;
		}

		bool l2r = true;
		auto s_buf_pt = &src_buf[size * src_step + (y0 + size) * src_stride];
		auto m_buf_pt = mask_buf + y0 * mask_stride;
		auto a_buf_pt = a_buf + y0 * a_stride;
		for (int y = y0; /*y < y1*/; y++, l2r ^= true,
			s_buf_pt += src_stride, m_buf_pt += mask_stride, a_buf_pt += a_stride) {
			// aggregate the points on the "incoming arc".
			if (y > y0) {
				switch (*m_buf_pt) {
				case mask::zero: curr_min = 0; break;
				case mask::full: break; // ignore full alpha.
				case mask::gray:
				default:
					for (int dx = -size; dx <= size; dx++)
						add(s_buf_pt[dx * src_step + arc[dx] * src_stride]);
					break;
				}
			}

			if (l2r) {
				for (int x = 0; x < dst_w; x++,
					s_buf_pt += src_step, m_buf_pt++, a_buf_pt += a_step) {
					switch (*m_buf_pt) {
					case mask::zero: *a_buf_pt = curr_min = 0; continue;
					case mask::full: *a_buf_pt = curr_min = max_alpha; continue;
					}

					// aggregate the points on the "incoming arc".
					if (x > 0) {
						for (int dy = -size; dy <= size; dy++)
							add(s_buf_pt[+arc[dy] * src_step + dy * src_stride]);
					}

					// write the alpha value.
					update_min();
					*a_buf_pt = curr_min;

					// aggregate the points on the "outgoing arc".
					if (x < dst_w - 1) {
						for (int dy = -size; dy <= size; dy++)
							pop(s_buf_pt[-arc[dy] * src_step + dy * src_stride]);
					}
				}
				s_buf_pt -= src_step; m_buf_pt--; a_buf_pt -= a_step;
			}
			else {
				for (int x = dst_w - 1; x >= 0; x--,
					s_buf_pt -= src_step, m_buf_pt--, a_buf_pt -= a_step) {
					switch (*m_buf_pt) {
					case mask::zero: *a_buf_pt = curr_min = 0; continue;
					case mask::full: *a_buf_pt = curr_min = max_alpha; continue;
					}

					// aggregate the points on the "incoming arc".
					if (x < dst_w - 1) {
						for (int dy = -size; dy <= size; dy++)
							add(s_buf_pt[-arc[dy] * src_step + dy * src_stride]);
					}

					// write the alpha value.
					update_min();
					*a_buf_pt = curr_min;

					// aggregate the points on the "outgoing arc".
					if (x > 0) {
						for (int dy = -size; dy <= size; dy++)
							pop(s_buf_pt[+arc[dy] * src_step + dy * src_stride]);
					}
				}
				s_buf_pt += src_step; m_buf_pt++; a_buf_pt += a_step;
			}

			// aggregate the points on the "outgoing arc".
			if (y < y1 - 1) {
				switch (*m_buf_pt) {
				case mask::zero: break;
				case mask::full: break; // ignore full alpha.
				case mask::gray:
				default:
					for (int dx = -size; dx <= size; dx++)
						pop(s_buf_pt[dx * src_step - arc[dx] * src_stride]);
					break;
				}
			}
			else break;
		}
	});
}


inline static Bounds deflate_common(auto&& alloc_and_mask_h,
	int src_w, int src_h,
	i16* dst_buf, bool dst_colored, size_t dst_stride,
	void* heap, int size_sq)
{
	using namespace max;
	using namespace masking::deflation;

	auto* const arc = reinterpret_cast<i32*>(heap);
	int const size = arith::arc::half(size_sq, arc);

	auto* mask_buf = reinterpret_cast<mask*>(arc + (2 * size + 1));
	size_t mask_stride = (src_w - 2 * size + 3) & (-4);

	auto [src_buf, src_stride, top, bottom] = alloc_and_mask_h(size, mask_buf, mask_stride);
	if (top >= bottom - 2 * size) return { 0,0,0,0 };

	mask_buf += top * mask_stride;
	src_buf += top * src_stride;
	dst_buf += top * dst_stride;
	src_h = bottom - top;

	auto [left, right] = mask_v<0>(src_w, src_h, size, mask_buf, mask_stride);
	if (left >= right) return { 0,0,0,0 };

	bottom -= 2 * size;
	mask_buf += left;
	src_buf += left;
	dst_buf += left * (dst_colored ? 4 : 1);
	src_w = right - left + 2 * size;

	(dst_colored ? find_min<1, 4> : find_min<1, 1>)
		(src_w, src_h, size, src_buf, src_stride,
			mask_buf, mask_stride, dst_buf, dst_stride, arc + size);

	return { left, top, right, bottom };
}

Bounds max::deflate(int src_w, int src_h,
	i16* src_buf, size_t src_stride,
	i16* dst_buf, bool dst_colored, size_t dst_stride,
	void* heap, int size_sq)
{
	using namespace masking::deflation;
	return deflate_common([&](int size, mask* mask_buf, size_t mask_stride) {
		auto [top, bottom] = mask_h_alpha<0>(src_w, src_h, size, src_buf, src_stride, mask_buf, mask_stride);
		return std::tuple{ src_buf, src_stride, top, bottom };
	}, src_w, src_h, dst_buf, dst_colored, dst_stride, heap, size_sq);
}

Bounds max::deflate(int src_w, int src_h,
	ExEdit::PixelYCA const* src_buf, size_t src_stride,
	i16* dst_buf, bool dst_colored, size_t dst_stride,
	void* heap, int size_sq, void* alpha_space)
{
	using namespace masking::deflation;
	return deflate_common([&](int size, mask* mask_buf, size_t mask_stride) {
		i16* med_buf = reinterpret_cast<i16*>(alpha_space);
		size_t med_stride = (src_w + 1) & (-2);
		auto [top, bottom] = mask_h_color<0>(src_w, src_h, size, src_buf, src_stride, mask_buf, mask_stride,
			med_buf, med_stride);
		return std::tuple{ med_buf, med_stride, top, bottom };
	}, src_w, src_h, dst_buf, dst_colored, dst_stride, heap, size_sq);
}

