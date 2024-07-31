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
#include "../kind_max/masking.hpp"

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
	multi_thread(dst_h, [&](int thread_id, int thread_num)
	{
		for (int y = thread_id; y < dst_h; y += thread_num) {
			auto s_buf_pt = src_buf + size * src_step + (y + size) * src_stride;
			auto m_buf_pt = mask_buf + y * mask_stride;
			auto a_buf_pt = a_buf + y * a_stride;

			int curr_min = max_alpha, curr_min_dur = -1;
			for (int x = 0; x < dst_w; x++,
				s_buf_pt += src_step, m_buf_pt++, a_buf_pt += a_step) {

				// hints by masking.
				switch (*m_buf_pt) {
				case mask::zero:
					*a_buf_pt = curr_min = 0;
					curr_min_dur = 2 * size;
					continue;
				case mask::full:
					*a_buf_pt = curr_min = max_alpha;
					curr_min_dur = 2 * size;
					continue;
				}

				if (--curr_min_dur >= 0) {
					// search the points on the "incoming arc".
					auto examine = [&](int dy) {
						int dx = arc[dy];
						int a = s_buf_pt[+dx * src_step + dy * src_stride];
						if (a <= curr_min) {
							int dur = 2 * dx;
							if (a < curr_min || dur > curr_min_dur) {
								curr_min = a; curr_min_dur = dur;
								if (a <= 0) return true;
							}
						}
						return false;
					};

					for (int dy = 0; dy >= -size; dy--) {
						if (examine(dy)) goto search_end1;
					}
					for (int dy = 1; dy <= size; dy++) {
						if (examine(dy)) goto search_end1;
					}

				search_end1:
					*a_buf_pt = curr_min;
				}
				else {
					// search the entire disc.
					int expiring_min = max_alpha; curr_min = max_alpha;
					if constexpr (true) {
						for (int dy = -size; dy <= size; dy++) {
							int const secant = arc[dy];
							for (int dx = -secant; dx <= secant; dx++) {
								int a = s_buf_pt[dx * src_step + dy * src_stride];
								if (a <= curr_min) {
									int dur = secant + dx;
									if (a < curr_min || dur > curr_min_dur) {
										if (dur > 0) {
											curr_min = a; curr_min_dur = dur;
											if (a <= 0) goto search_end2;
										}
										else if (a < expiring_min) expiring_min = a;
									}
								}
							}
						}
					}
					else {
						for (int dx = size; dx >= -size; dx--) {
							int const secant = arc[(-size - dx) >> 1],
								dur = dx + size;
							for (int dy = -secant; dy <= secant; dy++) {
								int a = s_buf_pt[(dx + size - arc[dy]) * src_step + dy * src_stride];
								if (a <= curr_min) {
									if (a < curr_min || dur > curr_min_dur) {
										if (dur > 0) {
											curr_min = a; curr_min_dur = dur;
											if (a <= 0) goto search_end2;
										}
										else if (a < expiring_min) expiring_min = a;
									}
								}
							}
						}
					}

				search_end2:
					*a_buf_pt = std::min(curr_min, expiring_min);
				}
			}
		}
	});
}

template<size_t src_step, size_t a_step>
static inline void find_min_2(int src_w, int src_h, int size,
	i16 const* src_buf, size_t src_stride,
	mask const* mask_buf, size_t mask_stride,
	i16* a_buf, size_t a_stride, i32 const* arc)
{
	// arc[i]: i ranges from -size to size.

	int const dst_w = src_w - 2 * size, dst_h = src_h - 2 * size,
		disk_area = 1 + 4 * (size + std::accumulate(arc + 1, arc + size + 1, 0));
	multi_thread(dst_w, [&](int thread_id, int thread_num)
	{
		for (int x = thread_id; x < dst_w; x += thread_num) {
			auto s_buf_pt = src_buf + (x + size) * src_step + size * src_stride;
			auto m_buf_pt = mask_buf + x;
			auto a_buf_pt = a_buf + x * a_step;

			int curr_min = max_alpha, curr_min_dur = -1;
			for (int y = 0; y < dst_h; y++,
				s_buf_pt += src_stride, m_buf_pt += mask_stride, a_buf_pt += a_stride) {

				// hints by masking.
				switch (*m_buf_pt) {
				case mask::zero:
					*a_buf_pt = curr_min = 0;
					curr_min_dur = 2 * size;
					continue;
				case mask::full:
					*a_buf_pt = curr_min = max_alpha;
					curr_min_dur = 2 * size;
					continue;
				}

				if (--curr_min_dur >= 0) {
					// search the points on the "incoming arc".
					auto examine = [&](int dx) {
						int dy = arc[dx];
						int a = s_buf_pt[+dx * src_step + dy * src_stride];
						if (a <= curr_min) {
							int dur = 2 * dy;
							if (a < curr_min || dur > curr_min_dur) {
								curr_min = a; curr_min_dur = dur;
								if (a <= 0) return true;
							}
						}
						return false;
					};

					for (int dx = 0; dx >= -size; dx--) {
						if (examine(dx)) goto search_end1;
					}
					for (int dx = 1; dx <= size; dx++) {
						if (examine(dx)) goto search_end1;
					}

				search_end1:
					*a_buf_pt = curr_min;
				}
				else {
					// search the entire disc.
					int expiring_min = max_alpha; curr_min = max_alpha;
					if constexpr (true) {
						for (int dx = -size; dx <= size; dx++) {
							int const secant = arc[dx];
							for (int dy = -secant; dy <= secant; dy++) {
								int a = s_buf_pt[dx * src_step + dy * src_stride];
								if (a <= curr_min) {
									int dur = secant + dy;
									if (a < curr_min || dur > curr_min_dur) {
										if (dur > 0) {
											curr_min = a; curr_min_dur = dur;
											if (a <= 0) goto search_end2;
										}
										else if (a < expiring_min) expiring_min = a;
									}
								}
							}
						}
					}
					else if constexpr (false) {
						for (int dy = -size; dy <= size; dy++) {
							int const secant = arc[dy];
							for (int dx = -secant; dx <= secant; dx++) {
								int a = s_buf_pt[dx * src_step + dy * src_stride];
								if (a <= curr_min) {
									int dur = arc[dx] + dy;
									if (a < curr_min || dur > curr_min_dur) {
										if (dur > 0) {
											curr_min = a; curr_min_dur = dur;
											if (a <= 0) goto search_end2;
										}
										else if (a < expiring_min) expiring_min = a;
									}
								}
							}
						}
					}
					else {
						for (int dy = size; dy >= -size; dy--) {
							int const secant = arc[(-size - dy) >> 1],
								dur = dy + size;
							for (int dx = -secant; dx <= secant; dx++) {
								int a = s_buf_pt[dx * src_step + (dy + size - arc[dx]) * src_stride];
								if (a <= curr_min) {
									if (a < curr_min || dur > curr_min_dur) {
										if (dur > 0) {
											curr_min = a; curr_min_dur = dur;
											if (a <= 0) goto search_end2;
										}
										else if (a < expiring_min) expiring_min = a;
									}
								}
							}
						}
					}

				search_end2:
					*a_buf_pt = std::min(curr_min, expiring_min);
				}
			}
		}
	});
}

inline static Bounds deflate_common(auto&& alloc_and_mask_h,
	int src_w, int src_h,
	i16* dst_buf, bool dst_colored, size_t dst_stride,
	void* heap, int size_sq)
{
	using namespace masking::deflation;

	auto* const arc = reinterpret_cast<i32*>(heap);
	int const size = arith::arc::half(size_sq, arc);

	auto* mask_heap = reinterpret_cast<i32*>(arc + (2 * size + 1));

	auto* mask_buf = reinterpret_cast<mask*>(mask_heap + 2 * (src_w - 2 * size));
	size_t mask_stride = (src_w - 2 * size + 3) & (-4);

	auto [src_buf, src_stride, top, bottom] = alloc_and_mask_h(size, mask_buf, mask_stride);
	if (top >= bottom - 2 * size) return { 0,0,0,0 };

	mask_buf += top * mask_stride;
	src_buf += top * src_stride;
	dst_buf += top * dst_stride;
	src_h = bottom - top;

	auto [left, right] = mask_v<0>(src_w, src_h, size, mask_buf, mask_stride, mask_heap);
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

Bounds max_fast::deflate(int src_w, int src_h,
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

Bounds max_fast::deflate(int src_w, int src_h,
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

