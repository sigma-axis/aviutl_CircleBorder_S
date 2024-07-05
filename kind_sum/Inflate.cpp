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
#include "inf_def.hpp"
#include "../kind_max/masking.hpp"

using namespace Calculation;
using mask = masking::mask;

template<size_t src_step, size_t a_step>
static inline void take_sum(int src_w, int src_h, int size,
	i16 const* src_buf, size_t src_stride,
	mask const* mask_buf, size_t mask_stride,
	i16* a_buf, size_t a_stride, int a_sum_cap, i32 const* arc)
{
	// arc[i]: i ranges from -size to size.

	int64_t const max_sum_alpha = static_cast<int64_t>(max_alpha)
		* (1 + 4 * (size + std::accumulate(arc + 1, arc + size + 1, 0))),
		sum_full_val = max_sum_alpha - (1 + 2 * size) * max_alpha;

	int const denom_bits = [](int N) {
		if (N <= 15) return 0;
		return N - 16;
	}(std::bit_width(static_cast<uint32_t>(a_sum_cap)) /* at most 21. */);
	int const numer = ((1ULL << 31) << denom_bits) / a_sum_cap; // numer * (a_sum_cap>>denom_bits) ~ 2^31.
	auto alpha_from_sum = [=](int64_t const& sum) -> i16 {
		if (sum >= a_sum_cap) return max_alpha;
		return static_cast<i16>((static_cast<uint32_t>(
			sum >> denom_bits) * numer + ((1 << 19) - 1)) >> 19);
	};

	int const dst_w = src_w + 2 * size, dst_h = src_h + 2 * size;
	multi_thread(dst_h, [&](int thread_id, int thread_num)
	{
		for (int y = thread_id; y < dst_h; y += thread_num) {
			auto s_buf_pt = src_buf - size * src_step + (y - size) * src_stride;
			auto m_buf_pt = mask_buf + y * mask_stride;
			auto a_buf_pt = a_buf + y * a_stride;

			int64_t sum_alpha = 0;
			int const dy_min = std::max(-size, size - y), dy_max = std::min(size, dst_h - size - 1 - y);
			for (int x = 0; x < dst_w; x++,
				s_buf_pt += src_step, m_buf_pt++, a_buf_pt += a_step) {

				// hints by masking.
				switch (*m_buf_pt) {
				case mask::zero: *a_buf_pt = 0; sum_alpha = 0; continue;
				case mask::full: *a_buf_pt = max_alpha; sum_alpha = sum_full_val; continue;
				}

				// aggregate the points on the "incoming arc".
				if (x < dst_w - size) {
					int dy0 = dy_min, dy1 = dy_max;
					if (x < size) {
						auto c = arc[size - x];
						dy0 = std::max(-c, dy_min);
						dy1 = std::min(+c, dy_max);
					}
					int_fast32_t diff = 0;
					if (x < dst_w - 2 * size) {
						for (int dy = dy0; dy <= dy1; dy++)
							diff += s_buf_pt[arc[dy] * src_step + dy * src_stride];
					}
					else {
						auto c = arc[dst_w - size - x] + 1;
						for (int dy = dy0, c1 = std::min(-c, dy1); dy <= c1; dy++)
							diff += s_buf_pt[arc[dy] * src_step + dy * src_stride];
						for (int dy = std::max(+c, dy0); dy <= dy1; dy++)
							diff += s_buf_pt[arc[dy] * src_step + dy * src_stride];
					}
					sum_alpha += diff;
				}

				// write the alpha value.
				*a_buf_pt = alpha_from_sum(sum_alpha);

				// aggregate the points on the "outgoing arc".
				if (x >= size) {
					int dy0 = dy_min, dy1 = dy_max;
					if (x >= dst_w - size) {
						auto c = arc[x - dst_w + size + 1];
						dy0 = std::max(-c, dy_min);
						dy1 = std::min(+c, dy_max);
					}
					int_fast32_t diff = 0;
					if (x >= 2 * size) {
						for (int dy = dy0; dy <= dy1; dy++)
							diff += s_buf_pt[-arc[dy] * src_step + dy * src_stride];
					}
					else {
						auto c = arc[x - size + 1] + 1;
						for (int dy = dy0, c1 = std::min(-c, dy1); dy <= c1; dy++)
							diff += s_buf_pt[-arc[dy] * src_step + dy * src_stride];
						for (int dy = std::max(+c, dy0); dy <= dy1; dy++)
							diff += s_buf_pt[-arc[dy] * src_step + dy * src_stride];
					}
					sum_alpha -= diff;
				}
			}
		}
	});
}

Bounds sum::inflate(int src_w, int src_h,
	i16* src_buf, bool src_colored, size_t src_stride,
	i16* dst_buf, bool dst_colored, size_t dst_stride, int a_sum_cap_rate,
	void* heap, int size_sq)
{
	using namespace masking::inflation;

	auto* const arc = reinterpret_cast<i32*>(heap);
	int const size = arith::arc::half(size_sq, arc);

	auto* mask_buf = reinterpret_cast<mask*>(arc + (2 * size + 1));
	size_t mask_stride = (src_w + 2 * size + 3) & (-4);

	auto [left, right] = (src_colored ? mask_v<4> : mask_v<1>)
		(src_w, src_h, size, src_buf, src_stride, mask_buf, mask_stride);
	if (left >= right) return { 0,0,0,0 };

	mask_buf += left;
	src_buf += left * (src_colored ? 4 : 1);
	dst_buf += left * (dst_colored ? 4 : 1);
	src_w = right - left;

	auto [top, bottom] = mask_h(src_w, src_h, size, mask_buf, mask_stride);

	right += 2 * size;
	mask_buf += top * mask_stride;
	src_buf += top * src_stride;
	dst_buf += top * dst_stride;
	src_h = bottom - top - 2 * size;

	if (size <= 0) std::unreachable();
	int a_sum_cap = max_alpha + static_cast<int>((
		(2 * std::sqrt(2 * std::sqrt(size_sq) - 1)) - 1
		) * a_sum_cap_rate * ((1.0 * max_alpha) / den_cap_rate));

	(src_colored ?
		dst_colored ? take_sum<4, 4> : take_sum<4, 1> :
		dst_colored ? take_sum<1, 4> : take_sum<1, 1>)
		(src_w, src_h, size, src_buf, src_stride,
			mask_buf, mask_stride, dst_buf, dst_stride, a_sum_cap, arc + size);

	return { left, top, right, bottom };
}

