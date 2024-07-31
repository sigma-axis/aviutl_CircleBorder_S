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

#include "../multi_thread.hpp"
#include "../arithmetics.hpp"
#include "inf_def.hpp"
#include "../kind_max/masking.hpp"

using namespace Calculation;
using mask = masking::mask;

template<size_t src_step, size_t a_step>
static inline void find_max(int src_w, int src_h, int size,
	i16 const* src_buf, size_t src_stride,
	mask const* mask_buf, size_t mask_stride,
	i16* a_buf, size_t a_stride, i32 const* arc)
{
	// arc[i]: i ranges from -size to size.

	int dst_w = src_w + 2 * size, dst_h = src_h + 2 * size;
	multi_thread(dst_h, [&](int thread_id, int thread_num)
	{
		for (int y = thread_id; y < dst_h; y += thread_num) {
			auto s_buf_pt = src_buf - size * src_step + (y - size) * src_stride;
			auto m_buf_pt = mask_buf + y * mask_stride;
			auto a_buf_pt = a_buf + y * a_stride;

			int curr_max = 0, curr_max_dur = -1;
			int const dy_min = std::max(-size, size - y), dy_max = std::min(size, dst_h - size - 1 - y);
			for (int x = 0; x < dst_w; x++,
				s_buf_pt += src_step, m_buf_pt++, a_buf_pt += a_step) {

				// hints by masking.
				switch (*m_buf_pt) {
				case mask::zero:
					*a_buf_pt = curr_max = 0;
					curr_max_dur = 2 * size;
					if (x >= dst_w - size) {
						for (int rest = dst_w - x; --rest >= 0; a_buf_pt += a_step) *a_buf_pt = 0;
						x = dst_w; // to break the for-loop, not the switch-block.
					}
					continue;
				case mask::full:
					*a_buf_pt = curr_max = max_alpha;
					curr_max_dur = 2 * size;
					continue;
				}

				if (--curr_max_dur >= 0) {
					// search the points on the "incoming arc".
					if (x >= dst_w - size) {
						if (curr_max > 0) {
							*a_buf_pt = curr_max;
							continue;
						}
						for (int rest = dst_w - x; --rest >= 0; a_buf_pt += a_step) *a_buf_pt = 0;
						break;
					}

					int dy0 = dy_min, dy1 = dy_max, c0 = 0, c1 = 1;
					if (x < size) {
						auto c = arc[size - x];
						dy0 = std::max(-c, dy0);
						dy1 = std::min(+c, dy1);
					}
					if (x >= dst_w - 2 * size) {
						c1 = arc[dst_w - size - x] + 1; c0 = -c1;
					}
					c0 = std::min(c0, dy1); c1 = std::max(c1, dy0);

					auto examine = [&](int dy) {
						int dx = arc[dy];
						int a = s_buf_pt[dx * src_step + dy * src_stride];
						if (a >= curr_max) {
							int dur = 2 * dx;
							if (a > curr_max || dur > curr_max_dur) {
								curr_max = a; curr_max_dur = dur;
								if (a >= max_alpha) return true;
							}
						}
						return false;
					};

					for (int dy = c0; dy >= dy0; dy--) {
						if (examine(dy)) goto search_end1;
					}
					for (int dy = c1; dy <= dy1; dy++) {
						if (examine(dy)) goto search_end1;
					}

				search_end1:
					*a_buf_pt = curr_max;
				}
				else {
					// search the entire disc.
					int expiring_max = 0; curr_max = 0;
					if constexpr (true) {
						for (int dy = dy_min; dy <= dy_max; dy++) {
							int const secant = arc[dy],
								dx0 = std::max(-secant, size - x),
								dx1 = std::min(+secant, dst_w - size - x - 1);
							for (int dx = dx0; dx <= dx1; dx++) {
								int a = s_buf_pt[dx * src_step + dy * src_stride];
								if (a >= curr_max) {
									int dur = secant + dx;
									if (a > curr_max || dur > curr_max_dur) {
										if (dur > 0) {
											curr_max = a; curr_max_dur = dur;
											if (a >= max_alpha) goto search_end2;
										}
										else if (a > expiring_max) expiring_max = a;
									}
								}
							}
						}
					}
					else {
						for (int dx = std::min(size, dst_w - size - x - 1); dx >= -size; dx--) {
							if (x + ((dx + size) >> 1) < size) break;

							int const secant = arc[(-size - dx) >> 1];
							int dy0 = std::max(-secant, dy_min),
								dy1 = std::min(+secant, dy_max),
								c0 = 0, c1 = 1;
							if (x + dx >= dst_w - 2 * size) {
								auto c = arc[dst_w - 2 * size - x - dx - 1];
								dy0 = std::max(-c, dy0);
								dy1 = std::min(+c, dy1);
							}

							if (x + dx < size) {
								c1 = arc[x + dx + 1] + 1; c0 = -c1;
							}
							c0 = std::min(c0, dy1); c1 = std::max(c1, dy0);

							auto examine = [&, dur = dx + size, s_buf = s_buf_pt + (dx + size) * src_step](int dy) {
								int a = s_buf[-arc[dy] * src_step + dy * src_stride];
								if (a >= curr_max) {
									if (a > curr_max || dur > curr_max_dur) {
										if (dur > 0) {
											curr_max = a; curr_max_dur = dur;
											if (a >= max_alpha) return true;
										}
										else if (a > expiring_max) expiring_max = a;
									}
								}
								return false;
							};
							for (int dy = dy0; dy <= c0; dy++) {
								if (examine(dy)) goto search_end2;
							}
							for (int dy = c1; dy <= dy1; dy++) {
								if (examine(dy)) goto search_end2;
							}
						}
					}

				search_end2:
					*a_buf_pt = std::max(curr_max, expiring_max);
				}
			}
		}
	});
}
template<size_t src_step, size_t a_step>
static inline void find_max_2(int src_w, int src_h, int size,
	i16 const* src_buf, size_t src_stride,
	mask const* mask_buf, size_t mask_stride,
	i16* a_buf, size_t a_stride, i32 const* arc)
{
	// arc[i]: i ranges from -size to size.

	int dst_w = src_w + 2 * size, dst_h = src_h + 2 * size;
	multi_thread(dst_w, [&](int thread_id, int thread_num)
	{
		for (int x = thread_id; x < dst_w; x += thread_num) {
			auto s_buf_pt = src_buf + (x - size) * src_step - size * src_stride;
			auto m_buf_pt = mask_buf + x;
			auto a_buf_pt = a_buf + x * a_step;

			int curr_max = 0, curr_max_dur = -1;
			int const dx_min = std::max(-size, size - x), dx_max = std::min(size, dst_w - size - 1 - x);
			for (int y = 0; y < dst_h; y++,
				s_buf_pt += src_stride, m_buf_pt += mask_stride, a_buf_pt += a_stride) {

				// hints by masking.
				switch (*m_buf_pt) {
				case mask::zero:
					*a_buf_pt = curr_max = 0;
					curr_max_dur = 2 * size;
					if (y >= dst_h - size) {
						for (int rest = dst_h - y; --rest >= 0; a_buf_pt += a_stride) *a_buf_pt = 0;
						y = dst_h; // to break the for-loop, not the switch-block.
					}
					continue;
				case mask::full:
					*a_buf_pt = curr_max = max_alpha;
					curr_max_dur = 2 * size;
					continue;
				}

				if (--curr_max_dur >= 0) {
					// search the points on the "incoming arc".
					if (y >= dst_h - size) {
						if (curr_max > 0) {
							*a_buf_pt = curr_max;
							continue;
						}
						for (int rest = dst_h - y; --rest >= 0; a_buf_pt += a_stride) *a_buf_pt = 0;
						break;
					}

					int dx0 = dx_min, dx1 = dx_max, c0 = 0, c1 = 1;
					if (y < size) {
						auto c = arc[size - y];
						dx0 = std::max(-c, dx0);
						dx1 = std::min(+c, dx1);
					}
					if (y >= dst_h - 2 * size) {
						c1 = arc[dst_h - size - y] + 1; c0 = -c1;
					}
					c0 = std::min(c0, dx1); c1 = std::max(c1, dx0);

					auto examine = [&](int dx) {
						int dy = arc[dx];
						int a = s_buf_pt[dx * src_step + dy * src_stride];
						if (a >= curr_max) {
							int dur = 2 * dy;
							if (a > curr_max || dur > curr_max_dur) {
								curr_max = a; curr_max_dur = dur;
								if (a >= max_alpha) return true;
							}
						}
						return false;
					};

					for (int dx = c0; dx >= dx0; dx--) {
						if (examine(dx)) goto search_end1;
					}
					for (int dx = c1; dx <= dx1; dx++) {
						if (examine(dx)) goto search_end1;
					}

				search_end1:
					*a_buf_pt = curr_max;
				}
				else {
					// search the entire disc.
					int expiring_max = 0; curr_max = 0;
					if constexpr (true) {
						for (int dx = dx_min; dx <= dx_max; dx++) {
							int const secant = arc[dx],
								dy0 = std::max(-secant, size - y),
								dy1 = std::min(+secant, dst_h - size - y - 1);
							for (int dy = dy0; dy <= dy1; dy++) {
								int a = s_buf_pt[dx * src_step + dy * src_stride];
								if (a >= curr_max) {
									int dur = secant + dy;
									if (a > curr_max || dur > curr_max_dur) {
										if (dur > 0) {
											curr_max = a; curr_max_dur = dur;
											if (a >= max_alpha) goto search_end2;
										}
										else if (a > expiring_max) expiring_max = a;
									}
								}
							}
						}
					}
					else if constexpr (false) {
						int const
							dy0 = std::max(-size, size - y),
							dy1 = std::min(size, dst_h - size - y - 1);
						for (int dy = dy1; dy >= dy0; dy--) {
							int const secant = arc[dy],
								dx0 = std::max(-secant, dx_min),
								dx1 = std::min(+secant, dx_max);
							for (int dx = dx0; dx < dx1; dx++) {
								int a = s_buf_pt[dx * src_step + dy * src_stride];
								if (a >= curr_max) {
									int dur = arc[dx] + dy;
									if (a > curr_max || dur > curr_max_dur) {
										if (dur > 0) {
											curr_max = a; curr_max_dur = dur;
											if (a >= max_alpha) goto search_end2;
										}
										else if (a > expiring_max) expiring_max = a;
									}
								}
							}
						}
					}
					else {
						for (int dy = std::min(size, dst_h - size - y - 1); dy >= -size; dy--) {
							if (y + ((dy + size) >> 1) < size) break;

							int const secant = arc[(-size - dy) >> 1];
							int dx0 = std::max(-secant, dx_min),
								dx1 = std::min(+secant, dx_max),
								c0 = 0, c1 = 1;
							if (y + dy >= dst_h - 2 * size) {
								auto c = arc[dst_h - 2 * size - y - dy - 1];
								dx0 = std::max(-c, dx0);
								dx1 = std::min(+c, dx1);
							}

							if (y + dy < size) {
								c1 = arc[y + dy + 1] + 1; c0 = -c1;
							}
							c0 = std::min(c0, dx1); c1 = std::max(c1, dx0);

							auto examine = [&, dur = dy + size, s_buf = s_buf_pt + (dy + size) * src_stride](int dx) {
								int a = s_buf[dx * src_step - arc[dx] * src_stride];
								if (a >= curr_max) {
									if (a > curr_max || dur > curr_max_dur) {
										if (dur > 0) {
											curr_max = a; curr_max_dur = dur;
											if (a >= max_alpha) return true;
										}
										else if (a > expiring_max) expiring_max = a;
									}
								}
								return false;
							};
							for (int dx = dx0; dx <= c0; dx++) {
								if (examine(dx)) goto search_end2;
							}
							for (int dx = c1; dx <= dx1; dx++) {
								if (examine(dx)) goto search_end2;
							}
						}
					}

				search_end2:
					*a_buf_pt = std::max(curr_max, expiring_max);
				}
			}
		}
	});
}


inline static Bounds inflate_common(auto&& alloc_and_mask_v,
	int src_w, int src_h,
	i16* dst_buf, bool dst_colored, size_t dst_stride,
	void* heap, int size_sq)
{
	using namespace masking::inflation;

	auto* const arc = reinterpret_cast<i32*>(heap);
	int const size = arith::arc::half(size_sq, arc);

	auto* const mask_heap = reinterpret_cast<i32*>(arc + 2 * size + 1);

	auto* mask_buf = reinterpret_cast<mask*>(mask_heap + 2 * (src_w + 2 * size));
	size_t const mask_stride = (src_w + 2 * size + 3) & (-4);

	auto [src_buf, src_stride, left, right] = alloc_and_mask_v(size, mask_buf, mask_stride, mask_heap);
	if (left >= right) return { 0,0,0,0 };

	mask_buf += left;
	src_buf += left;
	dst_buf += left * (dst_colored ? 4 : 1);
	src_w = right - left;

	auto [top, bottom] = mask_h(src_w, src_h, size, mask_buf, mask_stride);

	right += 2 * size;
	mask_buf += top * mask_stride;
	src_buf += top * src_stride;
	dst_buf += top * dst_stride;
	src_h = bottom - top - 2 * size;

	(dst_colored ? find_max_2<1, 4> : find_max_2<1, 1>)
		(src_w, src_h, size, src_buf, src_stride,
			mask_buf, mask_stride, dst_buf, dst_stride, arc + size);

	return { left, top, right, bottom };
}

Bounds max_fast::inflate(int src_w, int src_h,
	i16* src_buf, size_t src_stride,
	i16* dst_buf, bool dst_colored, size_t dst_stride,
	void* heap, int size_sq)
{
	using namespace masking::inflation;
	return inflate_common([&](int size, mask* mask_buf, size_t mask_stride, void* mask_heap) {
		auto [left, right] = mask_v_alpha(src_w, src_h, size, src_buf, src_stride, mask_buf, mask_stride, mask_heap);
		return std::tuple{ src_buf, src_stride, left, right };
	}, src_w, src_h, dst_buf, dst_colored, dst_stride, heap, size_sq);
}

Bounds max_fast::inflate(int src_w, int src_h,
	ExEdit::PixelYCA const* src_buf, size_t src_stride,
	i16* dst_buf, bool dst_colored, size_t dst_stride,
	void* heap, int size_sq, void* alpha_space)
{
	using namespace masking::inflation;
	return inflate_common([&](int size, mask* mask_buf, size_t mask_stride, void* mask_heap) {
		i16* med_buf = reinterpret_cast<i16*>(alpha_space);
		size_t med_stride = (src_w + 1) & (-2);
		auto [left, right] = mask_v_color(src_w, src_h, size, src_buf, src_stride, mask_buf, mask_stride, mask_heap,
			med_buf, med_stride);

		return std::tuple{ med_buf, med_stride, left, right };
	}, src_w, src_h, dst_buf, dst_colored, dst_stride, heap, size_sq);
}
