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
#include <bit>

#include "../multi_thread.hpp"
#include "../arithmetics.hpp"
#include "../buffer_op.hpp"
#include "inf_def.hpp"
#include "../kind_max/masking.hpp"

using namespace Calculation;
using mask = masking::mask;

template<size_t src_step, size_t a_step>
static inline void take_inv_sum(int src_w, int src_h, int size_canvas, int size_disk,
	i16 const* src_buf, size_t src_stride,
	mask const* mask_buf, size_t mask_stride,
	i16* a_buf, size_t a_stride, int a_sum_cap, i32 const* arc)
{
	// assumably, size_canvas = max(0, size_disk-1).
	// arc[i]: i ranges from -size_disk to size_disk.

	int64_t const max_sum_alpha = static_cast<int64_t>(max_alpha)
		* (1 + 4 * (size_disk + std::accumulate(arc + 1, arc + size_disk + 1, 0))),
		sum_full_val = max_sum_alpha - (1 + 2 * size_disk) * max_alpha;

	int const denom_bits = [](int N) {
		if (N <= 15) return 0;
		return N - 16;
	}(std::bit_width(static_cast<uint32_t>(a_sum_cap)) /* at most 21. */);
	int const numer = ((1ULL << 31) << denom_bits) / a_sum_cap; // numer * (a_sum_cap>>denom_bits) ~ 2^31.
	auto alpha_from_sum = [=](int64_t const& sum) -> i16 {
		auto inv = max_sum_alpha - sum;
		if (inv >= a_sum_cap) return 0;
		return max_alpha - static_cast<i16>(((static_cast<uint32_t>(
			inv) >> denom_bits) * numer + ((1 << 19) - 1)) >> 19);
	};

	int const dst_w = src_w - 2 * size_canvas, dst_h = src_h - 2 * size_canvas;
	multi_thread(dst_h, [&](int thread_id, int thread_num)
	{
		int const y0 = dst_h * thread_id / thread_num, y1 = dst_h * (thread_id + 1) / thread_num;

		int64_t sum_alpha = 0;
		// first state of sum_alpha.
		switch (mask_buf[y0 * mask_stride]) {
		case mask::zero: sum_alpha = 0; break;
		case mask::full: sum_alpha = max_sum_alpha; break;
		case mask::gray:
		default:
			for (int dy = -size_disk; dy <= size_disk; dy++) {
				int secant = arc[dy];
				int_fast32_t diff = 0;
				for (int dx = -secant; dx <= secant; dx++)
					diff += src_buf[(dx + size_canvas) * src_step + (y0 + dy + size_canvas) * src_stride];
				sum_alpha += diff;
			}
			break;
		}

		bool l2r = true;
		auto s_buf_pt = &src_buf[size_canvas * src_step + (y0 + size_canvas) * src_stride];
		auto m_buf_pt = mask_buf + y0 * mask_stride;
		auto a_buf_pt = a_buf + y0 * a_stride;
		for (int y = y0; y < y1; y++, l2r ^= true,
			s_buf_pt += src_stride, m_buf_pt += mask_stride, a_buf_pt += a_stride) {
			// aggregate the points on the "incoming" and "outgoing" arcs.
			if (y > y0) {
				switch (*m_buf_pt) {
				case mask::zero: sum_alpha = 0; break;
				case mask::full: sum_alpha = max_sum_alpha; break;
				case mask::gray:
				default:
					auto s_buf_dx = s_buf_pt - size_disk * src_step;
					int_fast32_t diff = 0;
					for (int dx = -size_disk; dx <= size_disk; dx++, s_buf_dx += src_step)
						diff += -s_buf_dx[-(arc[dx] + 1) * src_stride] + s_buf_dx[+arc[dx] * src_stride];
					sum_alpha += diff;
					break;
				}
			}

			if (l2r) {
				for (int x = 0; x < dst_w; x++,
					s_buf_pt += src_step, m_buf_pt++, a_buf_pt += a_step) {
					if (x > 0) {
						switch (*m_buf_pt) {
						case mask::zero: *a_buf_pt = 0; sum_alpha = 0; continue;
						case mask::full: *a_buf_pt = max_alpha; sum_alpha = max_sum_alpha; continue;
						}

						// aggregate the points on the "incoming" and "outgoing" arcs.
						auto s_buf_dy = s_buf_pt - size_disk * src_stride;
						int_fast32_t diff = 0;
						for (int dy = -size_disk; dy <= size_disk; dy++, s_buf_dy += src_stride)
							diff += -s_buf_dy[-(arc[dy] + 1) * src_step] + s_buf_dy[+arc[dy] * src_step];
						sum_alpha += diff;
					}

					// write the alpha value.
					*a_buf_pt = alpha_from_sum(sum_alpha);
				}
				s_buf_pt -= src_step; m_buf_pt--; a_buf_pt -= a_step;
			}
			else {
				for (int x = dst_w - 1; x >= 0; x--,
					s_buf_pt -= src_step, m_buf_pt--, a_buf_pt -= a_step) {
					if (x < dst_w - 1) {
						switch (*m_buf_pt) {
						case mask::zero: *a_buf_pt = 0; sum_alpha = 0; continue;
						case mask::full: *a_buf_pt = max_alpha; sum_alpha = max_sum_alpha; continue;
						}

						// aggregate the points on the "incoming" and "outgoing" arcs.
						auto s_buf_dy = s_buf_pt - size_disk * src_stride;
						int_fast32_t diff = 0;
						for (int dy = -size_disk; dy <= size_disk; dy++, s_buf_dy += src_stride)
							diff += s_buf_dy[-arc[dy] * src_step] - s_buf_dy[+(arc[dy] + 1) * src_step];
						sum_alpha += diff;
					}

					// write the alpha value.
					*a_buf_pt = alpha_from_sum(sum_alpha);
				}
				s_buf_pt += src_step; m_buf_pt++; a_buf_pt += a_step;
			}
		}
	});
}


inline static Bounds deflate_common(auto&& alloc_and_mask_h,
	int src_w, int src_h,
	i16* dst_buf, bool dst_colored, size_t dst_stride, int a_sum_cap_rate,
	void* heap, int size_sq)
{
	using namespace sum;
	using namespace masking::deflation;

	auto* const arc = reinterpret_cast<i32*>(heap);
	int const size_disk = arith::arc::half(size_sq, arc),
		size = std::max(size_disk - 1, 0);

	auto* mask_buf = reinterpret_cast<mask*>(arc + (2 * size_disk + 1));
	size_t mask_stride = (src_w - 2 * size + 3) & (-4);

	auto [src_buf, src_stride, top, bottom] = alloc_and_mask_h(size, size_disk, mask_buf, mask_stride);
	if (top >= bottom - 2 * size) return { 0,0,0,0 };

	mask_buf += top * mask_stride;
	src_buf += top * src_stride;
	dst_buf += top * dst_stride;
	src_h = bottom - top;

	auto [left, right] = (size < size_disk ? mask_v<1> : mask_v<0>)(src_w, src_h, size_disk, mask_buf, mask_stride);
	if (left >= right) return { 0,0,0,0 };

	bottom -= 2 * size;
	mask_buf += left;
	src_buf += left;
	dst_buf += left * (dst_colored ? 4 : 1);
	src_w = right - left + 2 * size;

	// clear the 1-dot chrome of outside the source image,
	// as searching may reach that area otherwise the code would be too complicated.
	if (size < size_disk) {
		// assuming difference is known to be 1.
		constexpr int diff = 1;

		// left and right margins are already cleared during the call to mask_h_****<>();
		auto len = right - left + 2 * size_disk;
		if (top - diff < 0)
			std::memset(src_buf + (left - diff) + (top - diff) * src_stride, 0, sizeof(*src_buf) * len);
		if (bottom + 2 * size_disk - diff > src_h)
			std::memset(src_buf + (left - diff) + src_h * src_stride, 0, sizeof(*src_buf) * len);
	}

	(dst_colored ? take_inv_sum<1, 4> : take_inv_sum<1, 1>)
		(src_w, src_h, size, size_disk, src_buf, src_stride,
			mask_buf, mask_stride, dst_buf, dst_stride,
			a_sum_cap_from_rate(a_sum_cap_rate, size_sq), arc + size_disk);

	return { left, top, right, bottom };
}

Bounds sum::deflate(int src_w, int src_h,
	i16* src_buf, size_t src_stride,
	i16* dst_buf, bool dst_colored, size_t dst_stride, int a_sum_cap_rate,
	void* heap, int size_sq)
{
	using namespace masking::deflation;
	return deflate_common([&](int size, int size_disk, mask* mask_buf, size_t mask_stride) {
		auto [top, bottom] = (size < size_disk ? mask_h_alpha<1> : mask_h_alpha<0>)
			(src_w, src_h, size_disk, src_buf, src_stride, mask_buf, mask_stride);
		return std::tuple{ src_buf, src_stride, top, bottom };
	}, src_w, src_h, dst_buf, dst_colored, dst_stride, a_sum_cap_rate, heap, size_sq);
}

Bounds sum::deflate(int src_w, int src_h,
	ExEdit::PixelYCA const* src_buf, size_t src_stride,
	i16* dst_buf, bool dst_colored, size_t dst_stride, int a_sum_cap_rate,
	void* heap, int size_sq, void* alpha_space)
{
	using namespace masking::deflation;
	return deflate_common([&](int size, int size_disk, mask* mask_buf, size_t mask_stride) {
		size_t med_stride = (src_w + 3) & (-2);
		i16* med_buf = reinterpret_cast<i16*>(alpha_space) + (1 + med_stride);

		auto [top, bottom] = (size < size_disk ? mask_h_color<1> : mask_h_color<0>)
			(src_w, src_h, size_disk, src_buf, src_stride, mask_buf, mask_stride, med_buf, med_stride);
		return std::tuple{ med_buf, med_stride, top, bottom };
	}, src_w, src_h, dst_buf, dst_colored, dst_stride, a_sum_cap_rate, heap, size_sq);
}

