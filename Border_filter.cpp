/*
The MIT License (MIT)

Copyright (c) 2024 sigma-axis

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <cstdint>
#include <algorithm>

#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
using byte = uint8_t;
#include <exedit.hpp>

#include "multi_thread.hpp"
#include "buffer_op.hpp"
#include "tiled_image.hpp"

#include "kind_bin/inf_def.hpp"
#include "kind_bin2x/inf_def.hpp"
#include "kind_max/inf_def.hpp"
#include "kind_sum/inf_def.hpp"

#include "kind_max_fast/inf_def.hpp"

#include "filter_defl.hpp"
#include "Border.hpp"

using namespace Filter::Border;
using namespace Calculation;
namespace Border_filter::params
{
	using namespace impl;
	namespace exp
	{
		static constexpr int
			den_size		= track_den[idx_track::size		],
			den_neg_size	= track_den[idx_track::neg_size	],
			den_transp		= track_den[idx_track::transp	],
			den_f_transp	= track_den[idx_track::f_transp	],
			den_blur		= track_den[idx_track::blur		],
			den_param_a		= track_den[idx_track::param_a	],
			den_img_x		= track_den[idx_track::img_x	],
			den_img_y		= track_den[idx_track::img_y	],

			min_size		= track_min[idx_track::size		],
			min_neg_size	= track_min[idx_track::neg_size	],
			min_transp		= track_min[idx_track::transp	],
			min_f_transp	= track_min[idx_track::f_transp	],
			min_blur		= track_min[idx_track::blur		],
			min_param_a		= track_min[idx_track::param_a	],
			min_img_x		= track_min[idx_track::img_x	],
			min_img_y		= track_min[idx_track::img_y	],

			max_size		= track_max[idx_track::size		],
			max_neg_size	= track_max[idx_track::neg_size	],
			max_transp		= track_max[idx_track::transp	],
			max_f_transp	= track_max[idx_track::f_transp	],
			max_blur		= track_max[idx_track::blur		],
			max_param_a		= track_max[idx_track::param_a	],
			max_img_x		= track_max[idx_track::img_x	],
			max_img_y		= track_max[idx_track::img_y	];
	
		static_assert(den_size == den_neg_size);

		using defl_base = Filter::Common::defl_base<den_size>;
	}
}
using namespace Border_filter::params::exp;


////////////////////////////////
// フィルタ処理．
////////////////////////////////
struct infl_base {
protected:
	struct process_spec {
		int sum_displace, neg_displace;
		bool do_infl, do_defl,
			allows_buffer_overlap; // whether intermediate result can overlap final result.
	};

	static size_t obj_mem_max() {
		return sizeof(ExEdit::PixelYCA) * ((exedit.yca_max_w + 8) * (exedit.yca_max_h + 4) - 4);
	}
	static std::pair<int, int> max_size_cand(size_t pix_max) {
		int w, h = w = static_cast<int>(std::sqrt(pix_max)) & (-4);
		if (w < exedit.yca_max_w) {
			w = (exedit.yca_max_w + 3) & (-4);
			h = (pix_max / w + 3) & (-4);
		}
		else if (h < exedit.yca_max_h) {
			h = (exedit.yca_max_h + 3) & (-4);
			w = (pix_max / h + 3) & (-4);
		}
		return { w, h };
	}
	virtual std::pair<int, int> max_size() const = 0;
	virtual process_spec tell_spec(int sum_size, int neg_size) const = 0;

	// copies alpha values from `efpip->obj_edit` to `efpip->obj_temp`,
	// in a manner that is contiguous with non-zero-sized filter.
	virtual void zero_op(int param_a, ExEdit::PixelYCA const* src_buf, size_t src_stride,
		int src_w, int src_h, i16* dst_buf, bool dst_colored, size_t dst_stride) const
	{
		if (dst_colored)
			buff::copy_alpha(src_buf, src_stride, 0, 0, src_w, src_h,
				buff::alpha_to_pixel(dst_buf), dst_stride / 4, 0, 0);
		else buff::copy_alpha(src_buf, src_stride, 0, 0, src_w, src_h,
			dst_buf, dst_stride, 0, 0);
	}
	virtual Bounds inflate_1(int sum_size_raw, int param_a, ExEdit::PixelYCA* src_buf, size_t stride,
		int src_w, int src_h, ExEdit::PixelYCA* dst_buf, void* heap) const = 0;
	virtual Bounds inflate_2(int sum_size_raw, int param_a, ExEdit::PixelYCA* src_buf, size_t src_stride,
		int src_w, int src_h, i16* dst_buf, size_t dst_stride, void* heap, void* alpha_space) const = 0;
	virtual Bounds deflate_2(int neg_size_raw, int param_a, i16* src_buf, size_t src_stride,
		int src_w, int src_h, ExEdit::PixelYCA* dst_buf, size_t dst_stride, void* heap) const = 0;

private:
	struct sizing {
		int sum_size_raw, sum_displace,
			neg_size_raw, neg_displace,
			blur_size_raw, blur_displace;
		bool do_infl, do_defl,
			allows_buffer_overlap;
		bool invalid;
	};
	sizing measure(int size, int neg_size, int blur_px, int src_w, int src_h) const
	{
		constexpr sizing invalid_size{ .invalid = true };
		if (size <= 0 && neg_size <= 0) return invalid_size;

		auto [mem_max_w, mem_max_h] = max_size();
		// to fit with the intermediate buffers.
		while (true) {
			blur_px = std::clamp(blur_px, 0, size);
			auto sum_size = size - (blur_px >> 1) + neg_size;
			auto spec = tell_spec(sum_size, neg_size);

			// calculate the inflated size.
			auto w = src_w + 2 * spec.sum_displace, h = src_h + 2 * spec.sum_displace;

			// check if it exceeds the limit.
			auto diff = std::max(w - mem_max_w, h - mem_max_h);
			if (diff <= 0) {
				// calculate the final size.
				int blur_displace = buff::blur_displace((blur_px * buff::den_blur_px) / den_size);
				w +=2 * (-spec.neg_displace + blur_displace),
				h +=2 * (-spec.neg_displace + blur_displace);

				// check if it exceeds the limit.
				diff = std::max(w - exedit.yca_max_w, h - exedit.yca_max_h);
				if (diff <= 0) return {
					.sum_size_raw = sum_size, .sum_displace = spec.sum_displace,
					.neg_size_raw = neg_size, .neg_displace = spec.neg_displace,
					.blur_size_raw = blur_px, .blur_displace = blur_displace,
					.do_infl = spec.do_infl, .do_defl = spec.do_defl,
					.allows_buffer_overlap = spec.allows_buffer_overlap,
				}; // it's OK.
			}

			diff = (diff + 1) >> 1;

			// either lack of memory or excess of final size.
			// try reducing the deflation size first.
			if (neg_size >= diff * den_size) neg_size -= diff * den_size;
			else if (neg_size > 0) neg_size = 0;

			// then the inflation size.
			else if (size >= diff * den_size) size -= diff * den_size;
			else if (size > 0) size = 0;
			else return invalid_size; // nothing that I could do.
		}
	}

public:
	int measure_displace(int size, int neg_size, int blur_px, int src_w, int src_h) const
	{
		auto sz = measure(size, neg_size, blur_px, src_w, src_h);
		return sz.sum_displace - sz.neg_displace + sz.blur_displace;
	}

	struct infl_result {
		int displace;
		bool is_empty; // entire image is found transparent.
		bool invalid;
	};
	infl_result operator()(int size, int neg_size, int blur_px, int param_a, ExEdit::FilterProcInfo* efpip) const
	{
		constexpr infl_result invalid{ .invalid = true };
		// calculate sizing values.
		int const src_w = efpip->obj_w, src_h = efpip->obj_h;
		auto const sz = measure(size, neg_size, blur_px, src_w, src_h);
		if (sz.invalid) return invalid;

		// manipulate the final size so chatterings wouldn't occur.
		int const displace = (sz.sum_size_raw - sz.neg_size_raw + (sz.blur_size_raw >> 1) + (den_size >> 1)) / den_size,
			diff_displace = displace - (sz.sum_displace - sz.neg_displace + sz.blur_displace),
			diff_disp_cnt = diff_displace * (1 + efpip->obj_line);

		Bounds bd{ 0, 0, src_w, src_h };
		if (!sz.do_infl && !sz.do_defl) zero_op(param_a, efpip->obj_edit, efpip->obj_line,
			bd.wd(), bd.ht(), &efpip->obj_temp->a, true, 4 * efpip->obj_line);
		else {
			if (sz.do_defl) {
				// allocate memory layout.
				size_t med_stride = (bd.wd() + 2 * sz.sum_displace + 1) & (-2);
				i16* med_buffer; void* heap;
				if (sz.allows_buffer_overlap) {
					med_buffer = reinterpret_cast<i16*>(efpip->obj_temp);
					heap = *exedit.memory_ptr;
				}
				else {
					med_buffer = reinterpret_cast<i16*>(*exedit.memory_ptr);
					heap = med_buffer + med_stride * (bd.ht() + 2 * sz.sum_displace);
				}

				// then process by two passes.
				if (sz.do_infl) {
					bd = inflate_2(sz.sum_size_raw, param_a,
						efpip->obj_edit, efpip->obj_line, bd.wd(), bd.ht(),
						med_buffer, med_stride, heap, efpip->obj_temp);
					if (bd.is_empty()) return {
						.displace = displace,
						.is_empty = true,
					};
				}
				else zero_op(param_a, efpip->obj_edit, efpip->obj_line, bd.wd(), bd.ht(),
					med_buffer, false, med_stride);
				bd = deflate_2(sz.neg_size_raw, param_a,
					med_buffer + bd.L + bd.T * med_stride, med_stride, bd.wd(), bd.ht(),
					&efpip->obj_temp[bd.L + bd.T * efpip->obj_line + diff_disp_cnt], efpip->obj_line, heap)
					.move(bd.L + diff_displace, bd.T + diff_displace);
			}
			else {
				// process by one pass.
				bd = inflate_1(sz.sum_size_raw, param_a,
					efpip->obj_edit, efpip->obj_line, bd.wd(), bd.ht(),
					&efpip->obj_temp[diff_disp_cnt], *exedit.memory_ptr)
					.move(diff_displace, diff_displace);
			}
			if (bd.is_empty()) return {
				.displace = displace,
				.is_empty = true,
			};

			// apply blur.
			if (sz.blur_size_raw > 0) {
				buff::blur_alpha(efpip->obj_temp, efpip->obj_line,
					bd.L, bd.T, bd.wd(), bd.ht(), (sz.blur_size_raw * buff::den_blur_px) / den_size,
					*exedit.memory_ptr);

				bd = bd.inflate_br(2 * sz.blur_displace);
			}

			// clear the four sides of margins if present.
			int dst_w = src_w + 2 * displace, dst_h = src_h + 2 * displace;
			buff::clear_alpha_chrome(efpip->obj_temp, efpip->obj_line,
				{ 0, 0, dst_w, dst_h }, bd);
		}

		return { .displace = displace };
	}
};

// base class for "bin" and "bin2x".
struct infl_bin_base : infl_base {
protected:
	static constexpr i16 to_thresh(int param_a) {
		return (max_alpha - 1) * param_a / max_param_a;
	}
	void zero_op(int param_a, ExEdit::PixelYCA const* src_buf, size_t src_stride,
		int src_w, int src_h, i16* dst_buf, bool dst_colored, size_t dst_stride) const override
	{
		if (dst_colored)
			buff::binarize(src_buf, src_stride, 0, 0, src_w, src_h,
				buff::alpha_to_pixel(dst_buf), dst_stride / 4, 0, 0, to_thresh(param_a));
		else buff::binarize(src_buf, src_stride, 0, 0, src_w, src_h,
			dst_buf, dst_stride, 0, 0, to_thresh(param_a));
	}
	Bounds inflate_2(int sum_size_raw, int param_a, ExEdit::PixelYCA* src_buf, size_t src_stride,
		int src_w, int src_h, i16* dst_buf, size_t dst_stride, void* heap, void* alpha_space) const override
	{
		return bin::inflate(src_w, src_h,
			&src_buf->a, true, 4 * src_stride, to_thresh(param_a),
			dst_buf, false, dst_stride,
			heap, (sum_size_raw * sum_size_raw) / (den_size * den_size));
	}
};

// algorithm "bin".
constexpr struct : infl_bin_base {
private:
	static inline constinit int mem_max_w = 0, mem_max_h = 0;
	static void init_mem_max()
	{
		size_t const mem_max = obj_mem_max();
		std::tie(mem_max_w, mem_max_h) = max_size_cand(mem_max / sizeof(i32));

		// trim by 4 dots until it fits within available space.
		while (sizeof(i16) * ((mem_max_w + 1) & (-2)) * mem_max_h > mem_max ||
			bin::inflate_heap_size(mem_max_w, mem_max_h, std::min(mem_max_w, mem_max_h) >> 1) > mem_max)
			mem_max_w -= 4, mem_max_h -= 4;
	}

protected:
	std::pair<int, int> max_size() const override
	{
		if (mem_max_w == 0) init_mem_max();
		return { mem_max_w, mem_max_h };
	}
	process_spec tell_spec(int sum_size, int neg_size) const override
	{
		int sum_displace = bin::inflate_radius<den_size>(sum_size),
			neg_displace = bin::deflate_radius<den_size>(neg_size);
		return {
			.sum_displace = sum_displace,
			.neg_displace = neg_displace,
			.do_infl = sum_displace > 0,
			.do_defl = neg_displace > 0,
			.allows_buffer_overlap = true,
		};
	}
	Bounds inflate_1(int sum_size_raw, int param_a, ExEdit::PixelYCA* src_buf, size_t stride,
		int src_w, int src_h, ExEdit::PixelYCA* dst_buf, void* heap) const override
	{
		return bin::inflate(src_w, src_h,
			&src_buf->a, true, 4 * stride, to_thresh(param_a),
			&dst_buf->a, true, 4 * stride,
			heap, (sum_size_raw * sum_size_raw) / (den_size * den_size));
	}
	Bounds deflate_2(int neg_size_raw, int param_a, i16* src_buf, size_t src_stride,
		int src_w, int src_h, ExEdit::PixelYCA* dst_buf, size_t dst_stride, void* heap) const override
	{
		return bin::deflate(src_w, src_h,
			src_buf, false, src_stride, to_thresh(param_a),
			&dst_buf->a, true, 4 * dst_stride,
			heap, (neg_size_raw * neg_size_raw) / (den_size * den_size));
	}
} infl_bin{};

// algorithm "bin2x".
constexpr struct : infl_bin_base {
private:
	static inline constinit int mem_max_w = 0, mem_max_h = 0;
	static void init_mem_max()
	{
		size_t const mem_max = obj_mem_max();
		std::tie(mem_max_w, mem_max_h) = max_size_cand(mem_max / sizeof(i32));

		// trim by 4 dots until it fits within available space.
		while (sizeof(i16) * ((mem_max_w + 1) & (-2)) * mem_max_h > mem_max ||
			bin2x::inflate_heap_size(mem_max_w, mem_max_h, std::min(mem_max_w, mem_max_h) >> 1) + sizeof(i32) > mem_max)
			mem_max_w -= 4, mem_max_h -= 4;
	}

protected:
	std::pair<int, int> max_size() const override
	{
		if (mem_max_w == 0) init_mem_max();
		return { mem_max_w, mem_max_h };
	}
	process_spec tell_spec(int sum_size, int neg_size) const override
	{
		bool do_defl = neg_size >= den_size / 2;
		int sum_displace = (do_defl ? bin::inflate_radius<den_size> : bin2x::inflate_radius<den_size>)(sum_size),
			neg_displace = bin2x::deflate_radius<den_size>(neg_size);
		return {
			.sum_displace = sum_displace,
			.neg_displace = neg_displace,
			.do_infl = sum_displace > 0,
			.do_defl = do_defl,
			.allows_buffer_overlap = true,
		};
	}
	Bounds inflate_1(int sum_size_raw, int param_a, ExEdit::PixelYCA* src_buf, size_t stride,
		int src_w, int src_h, ExEdit::PixelYCA* dst_buf, void* heap) const override
	{
		return bin2x::inflate(src_w, src_h,
			&src_buf->a, true, 4 * stride, to_thresh(param_a),
			&dst_buf->a, true, 4 * stride,
			heap, (4 * sum_size_raw * sum_size_raw) / (den_size * den_size));
	}
	Bounds deflate_2(int neg_size_raw, int param_a, i16* src_buf, size_t src_stride,
		int src_w, int src_h, ExEdit::PixelYCA* dst_buf, size_t dst_stride, void* heap) const override
	{
		return bin2x::deflate(src_w, src_h,
			src_buf, false, src_stride, to_thresh(param_a),
			&dst_buf->a, true, 4 * dst_stride,
			heap, (4 * neg_size_raw * neg_size_raw) / (den_size * den_size));
	}
} infl_bin2x{};

// algorithm "max".
constexpr struct : infl_base {
private:
	static inline constinit int mem_max_w = 0, mem_max_h = 0;
	static void init_mem_max()
	{
		size_t const mem_max = obj_mem_max();
		std::tie(mem_max_w, mem_max_h) = max_size_cand(mem_max / (sizeof(i16) + sizeof(uint8_t)));

		// trim by 4 dots until it fits within available space.
		while (sizeof(i16) * ((mem_max_w + 1) & (-2)) * mem_max_h
			+ max::inflate_heap_size(mem_max_w, mem_max_h, std::min(mem_max_w, mem_max_h) >> 1) > mem_max)
			mem_max_w -= 4, mem_max_h -= 4;
	}

protected:
	std::pair<int, int> max_size() const override
	{
		if (mem_max_w == 0) init_mem_max();
		return { mem_max_w, mem_max_h };
	}
	process_spec tell_spec(int sum_size, int neg_size) const override
	{
		int sum_displace = max::inflate_radius<den_size>(sum_size),
			neg_displace = max::deflate_radius<den_size>(neg_size);
		return {
			.sum_displace = sum_displace,
			.neg_displace = neg_displace,
			.do_infl = sum_displace > 0,
			.do_defl = neg_displace > 0,
			.allows_buffer_overlap = false,
		};
	}
	Bounds inflate_1(int sum_size_raw, int param_a, ExEdit::PixelYCA* src_buf, size_t stride,
		int src_w, int src_h, ExEdit::PixelYCA* dst_buf, void* heap) const override
	{
		return max::inflate(src_w, src_h, src_buf, stride,
			&dst_buf->a, true, 4 * stride,
			reinterpret_cast<void*>(reinterpret_cast<uintptr_t>(heap) + max::alpha_space_size(src_w, src_h)),
			(sum_size_raw * sum_size_raw) / (den_size * den_size), heap);
	}
	Bounds inflate_2(int sum_size_raw, int param_a, ExEdit::PixelYCA* src_buf, size_t src_stride,
		int src_w, int src_h, i16* dst_buf, size_t dst_stride, void* heap, void* alpha_space) const override
	{
		return max::inflate(src_w, src_h, src_buf, src_stride,
			dst_buf, false, dst_stride,
			heap, (sum_size_raw * sum_size_raw) / (den_size * den_size), alpha_space);
	}
	Bounds deflate_2(int neg_size_raw, int param_a, i16* src_buf, size_t src_stride,
		int src_w, int src_h, ExEdit::PixelYCA* dst_buf, size_t dst_stride, void* heap) const override
	{
		return max::deflate(src_w, src_h, src_buf, src_stride,
			&dst_buf->a, true, 4 * dst_stride,
			heap, (neg_size_raw * neg_size_raw) / (den_size * den_size));
	}
} infl_max;

// algorithm "max_fast".
constexpr struct : infl_base {
private:
	static inline constinit int mem_max_w = 0, mem_max_h = 0;
	static void init_mem_max()
	{
		size_t const mem_max = obj_mem_max();
		std::tie(mem_max_w, mem_max_h) = max_size_cand(mem_max / (sizeof(i16) + sizeof(uint8_t)));

		// trim by 4 dots until it fits within available space.
		while (sizeof(i16) * ((mem_max_w + 1) & (-2)) * mem_max_h
			+ max_fast::inflate_heap_size(mem_max_w, mem_max_h, std::min(mem_max_w, mem_max_h) >> 1) > mem_max)
			mem_max_w -= 4, mem_max_h -= 4;
	}

protected:
	std::pair<int, int> max_size() const override
	{
		if (mem_max_w == 0) init_mem_max();
		return { mem_max_w, mem_max_h };
	}
	process_spec tell_spec(int sum_size, int neg_size) const override
	{
		int sum_displace = max_fast::inflate_radius<den_size>(sum_size),
			neg_displace = max_fast::deflate_radius<den_size>(neg_size);
		return {
			.sum_displace = sum_displace,
			.neg_displace = neg_displace,
			.do_infl = sum_displace > 0,
			.do_defl = neg_displace > 0,
			.allows_buffer_overlap = false,
		};
	}
	Bounds inflate_1(int sum_size_raw, int param_a, ExEdit::PixelYCA* src_buf, size_t stride,
		int src_w, int src_h, ExEdit::PixelYCA* dst_buf, void* heap) const override
	{
		return max_fast::inflate(src_w, src_h, src_buf, stride,
			&dst_buf->a, true, 4 * stride,
			reinterpret_cast<void*>(reinterpret_cast<uintptr_t>(heap) + max_fast::alpha_space_size(src_w, src_h)),
			(sum_size_raw * sum_size_raw) / (den_size * den_size), heap);
	}
	Bounds inflate_2(int sum_size_raw, int param_a, ExEdit::PixelYCA* src_buf, size_t src_stride,
		int src_w, int src_h, i16* dst_buf, size_t dst_stride, void* heap, void* alpha_space) const override
	{
		return max_fast::inflate(src_w, src_h, src_buf, src_stride,
			dst_buf, false, dst_stride,
			heap, (sum_size_raw * sum_size_raw) / (den_size * den_size), alpha_space);
	}
	Bounds deflate_2(int neg_size_raw, int param_a, i16* src_buf, size_t src_stride,
		int src_w, int src_h, ExEdit::PixelYCA* dst_buf, size_t dst_stride, void* heap) const override
	{
		return max_fast::deflate(src_w, src_h, src_buf, src_stride,
			&dst_buf->a, true, 4 * dst_stride,
			heap, (neg_size_raw * neg_size_raw) / (den_size * den_size));
	}
} infl_max_fast;

// algorithm "sum".
constexpr struct : infl_base {
private:
	static inline constinit int mem_max_w = 0, mem_max_h = 0;
	static void init_mem_max()
	{
		size_t const mem_max = obj_mem_max();
		std::tie(mem_max_w, mem_max_h) = max_size_cand(mem_max / (sizeof(i16) + sizeof(uint8_t)));

		// trim by 4 dots until it fits within available space.
		while (sizeof(i16) * ((mem_max_w + 1) & (-2)) * mem_max_h
			+ sum::inflate_heap_size(mem_max_w, mem_max_h, std::min(mem_max_w, mem_max_h) >> 1) > mem_max)
			mem_max_w -= 4, mem_max_h -= 4;
	}

protected:
	std::pair<int, int> max_size() const override
	{
		if (mem_max_w == 0) init_mem_max();
		return { mem_max_w, mem_max_h };
	}
	process_spec tell_spec(int sum_size, int neg_size) const override
	{
		int sum_displace = sum::inflate_radius<den_size>(sum_size),
			neg_displace = sum::deflate_radius<den_size>(neg_size);
		return {
			.sum_displace = sum_displace + (neg_size >= den_size ? 1 : 0),
			.neg_displace = neg_displace + (neg_size >= den_size ? 1 : 0),
			.do_infl = sum_displace > 0,
			.do_defl = neg_size >= den_size,
			.allows_buffer_overlap = false,
		};
	}
	Bounds inflate_1(int sum_size_raw, int param_a, ExEdit::PixelYCA* src_buf, size_t stride,
		int src_w, int src_h, ExEdit::PixelYCA* dst_buf, void* heap) const override
	{
		return sum::inflate(src_w, src_h, src_buf, stride,
			&dst_buf->a, true, 4 * stride, (param_a * sum::den_cap_rate) / max_param_a,
			reinterpret_cast<void*>(reinterpret_cast<uintptr_t>(heap) + sum::alpha_space_size(src_w, src_h)),
			(sum_size_raw * sum_size_raw) / (den_size * den_size), heap);
	}
	Bounds inflate_2(int sum_size_raw, int param_a, ExEdit::PixelYCA* src_buf, size_t src_stride,
		int src_w, int src_h, i16* dst_buf, size_t dst_stride, void* heap, void* alpha_space) const override
	{
		return sum::inflate(src_w, src_h, src_buf, src_stride,
			dst_buf + (1 + dst_stride), false, dst_stride, (param_a * sum::den_cap_rate) / max_param_a,
			heap, (sum_size_raw * sum_size_raw) / (den_size * den_size), alpha_space)
			.inflate_br(2);
	}
	Bounds deflate_2(int neg_size_raw, int param_a, i16* src_buf, size_t src_stride,
		int src_w, int src_h, ExEdit::PixelYCA* dst_buf, size_t dst_stride, void* heap) const override
	{
		return sum::deflate(src_w - 2, src_h - 2, src_buf + (1 + src_stride), src_stride,
			&dst_buf->a, true, 4 * dst_stride, (param_a * sum::den_cap_rate) / max_param_a,
			heap, (neg_size_raw * neg_size_raw) / (den_size * den_size));
	}
} infl_sum;

static constexpr infl_base const& choose_infl(Filter::Algorithm algorithm) {
	switch (algorithm) {
		using algo = Filter::Algorithm;
	case algo::bin: return infl_bin;
	default:
	case algo::bin2x: return infl_bin2x;
	case algo::max: return infl_max;
	case algo::max_fast: return infl_max_fast;
	case algo::sum: return infl_sum;
	}
}

constexpr Filter::Common::defl_bin<den_size, max_param_a> defl_bin{};
constexpr Filter::Common::defl_bin2x<den_size, max_param_a> defl_bin2x{};
constexpr Filter::Common::defl_max<den_size> defl_max{};
constexpr Filter::Common::defl_max_fast<den_size> defl_max_fast{};
constexpr Filter::Common::defl_sum<den_size, max_param_a> defl_sum{};

static constexpr defl_base const& choose_defl(Filter::Algorithm algorithm) {
	switch (algorithm) {
		using algo = Filter::Algorithm;
	case algo::bin: return defl_bin;
	default:
	case algo::bin2x: return defl_bin2x;
	case algo::max: return defl_max;
	case algo::max_fast: return defl_max_fast;
	case algo::sum: return defl_sum;
	}
}

static inline void expand_foursides(int displace, int f_alpha, ExEdit::Filter* efp, ExEdit::FilterProcInfo* efpip)
{
	if (displace <= 0) return;

	efp->exfunc->fill(efpip->obj_temp,
		0, 0, efpip->obj_w + 2 * displace, efpip->obj_h + 2 * displace,
		0, 0, 0, 0, 0x02);
	if (f_alpha > 0) {
		efp->exfunc->bufcpy(efpip->obj_temp, displace, displace,
			efpip->obj_edit, 0, 0, efpip->obj_w, efpip->obj_h, max_alpha - f_alpha,
			0x12000003
			| (f_alpha < max_alpha ? 0 : 1 << 24) // flag to ignore alpha multiplication.
		);
	}
	efpip->obj_w += 2 * displace; efpip->obj_h += 2 * displace;
	std::swap(efpip->obj_temp, efpip->obj_edit);
}


////////////////////////////////
// フィルタ処理のエントリポイント．
////////////////////////////////
BOOL impl::func_proc(ExEdit::Filter* efp, ExEdit::FilterProcInfo* efpip)
{
	int const src_w = efpip->obj_w, src_h = efpip->obj_h;
	if (src_w <= 0 || src_h <= 0) return TRUE;

	int const
		size		= std::clamp(efp->track[idx_track::size		], min_size		, max_size		),
		neg_size	= std::clamp(efp->track[idx_track::neg_size	], min_neg_size	, max_neg_size	),
		transp		= std::clamp(efp->track[idx_track::transp	], min_transp	, max_transp	),
		f_transp	= std::clamp(efp->track[idx_track::f_transp	], min_f_transp	, max_f_transp	),
		blur		= std::clamp(efp->track[idx_track::blur		], min_blur		, max_blur		),
		param_a		= std::clamp(efp->track[idx_track::param_a	], min_param_a	, max_param_a	),
		img_x		= std::clamp(efp->track[idx_track::img_x	], min_img_x	, max_img_x		),
		img_y		= std::clamp(efp->track[idx_track::img_y	], min_img_y	, max_img_y		);
	auto* const exdata = reinterpret_cast<Exdata*>(efp->exdata_ptr);

	int const
		alpha = std::clamp<int>(max_alpha * (max_transp - transp) / max_transp, 0, max_alpha),
		f_alpha = std::clamp<int>(max_alpha * (max_f_transp - f_transp) / max_f_transp, 0, max_alpha),
		blur_px = blur == 0 || arith::abs(size) <= den_size ? 0 :
			(std::max(arith::abs(size) - den_size - 1, 0) * (blur - 1)) / (max_blur - 1) + 1;

	// "fake" size by 0.4 so integral-sized shape will more likely look smooth.
	int const lifted_size = size + (size >= 0 ? +1 : -1) * ((den_size >> 1) - 1);

	// handle trivial cases.
	if (size == 0 || alpha <= 0) {
		if (size > 0) {
			if (int displace = choose_infl(exdata->algorithm)
				.measure_displace(lifted_size, neg_size, blur_px, src_w, src_h);
				displace > 0) {
				expand_foursides(displace, f_alpha, efp, efpip);
				return TRUE;
			}
		}
		if (f_alpha >= max_alpha); // nothing to do.
		else if (f_alpha <= 0)
			buff::clear_alpha(efpip->obj_edit, efpip->obj_line, 0, 0, src_w, src_h);
		else
			buff::mult_alpha(f_alpha, efpip->obj_edit, efpip->obj_line, 0, 0, src_w, src_h);
		return TRUE;
	}

	// general cases.
	if (lifted_size > 0) {
		auto result = choose_infl(exdata->algorithm)(lifted_size, neg_size, blur_px, param_a, efpip);
		if (result.invalid) return TRUE;

		int const
			dst_w = efpip->obj_w += 2 * result.displace,
			dst_h = efpip->obj_h += 2 * result.displace;
		std::swap(efpip->obj_temp, efpip->obj_edit);

		if (result.is_empty) {
			// the entire image is transparent.
			if (result.displace > 0) {
				efp->exfunc->fill(efpip->obj_temp,
					src_w, 0, dst_w - src_w, src_h,
					0, 0, 0, 0, 0x02);
				efp->exfunc->fill(efpip->obj_temp,
					0, src_h, dst_w, dst_h - src_h,
					0, 0, 0, 0, 0x02);
			}
			return TRUE;
		}

		if (tiled_image const img{ exdata->file, img_x, img_y, result.displace, efp, *exedit.memory_ptr }) {
			// image seems to have been successfully loaded.
			// fill with the pattern image.
			auto blend = [&, img_stride = efpip->obj_line](ExEdit::PixelYCA const& infl, ExEdit::PixelYCA const& orig, int px_x, int px_y) noexcept -> ExEdit::PixelYCA {
				auto& col = img[px_x + px_y * img_stride];
				int a = ((alpha * ((infl.a * col.a) >> log2_max_alpha)) >> log2_max_alpha) - orig.a;
				int A = static_cast<uint32_t>(f_alpha * orig.a) >> log2_max_alpha; // making sure A >= 0.
				if (a > 0) return {
					.y  = static_cast<i16>((a * col.y + A * orig.y) / (a + A)),
					.cb = static_cast<i16>((a * col.cb + A * orig.cb) / (a + A)),
					.cr = static_cast<i16>((a * col.cr + A * orig.cr) / (a + A)),
					.a  = static_cast<i16>(a + A),
				};
				else return {
					.y = orig.y, .cb = orig.cb, .cr = orig.cr,
					.a = static_cast<i16>(A),
				};
			};
			auto paint = [&, img_stride = efpip->obj_line](ExEdit::PixelYCA const& infl, int px_x, int px_y) noexcept -> ExEdit::PixelYCA {
				auto& col = img[px_x + px_y * img_stride];
				return {
					.y = col.y, .cb = col.cb, .cr = col.cr,
					.a = static_cast<i16>((alpha * ((infl.a * col.a) >> log2_max_alpha)) >> log2_max_alpha),
				};
			};

			multi_thread(dst_h, [&](int thread_id, int thread_num) {
				for (int y = thread_id; y < dst_h; y += thread_num) {
					int i_x = img.ox, i_y = (y + img.oy) % img.h;
					auto incr_x = [&] { i_x++; if (i_x >= img.w) i_x -= img.w; };

					auto* dst = efpip->obj_edit + y * efpip->obj_line;
					if (y < result.displace || y >= dst_h - result.displace) {
						for (int x = dst_w; --x >= 0; dst++, incr_x()) *dst = paint(*dst, i_x, i_y);
					}
					else {
						auto* src = efpip->obj_temp + (y - result.displace) * efpip->obj_line;
						for (int x = result.displace; --x >= 0; dst++, incr_x()) *dst = paint(*dst, i_x, i_y);
						for (int x = dst_w - 2 * result.displace; --x >= 0; src++, dst++, incr_x()) *dst = blend(*dst, *src, i_x, i_y);
						for (int x = result.displace; --x >= 0; dst++, incr_x()) *dst = paint(*dst, i_x, i_y);
					}
				}
			});
		}
		else {
			// fill with specified color.
			auto const col = buff::fromRGB(exdata->color.r, exdata->color.g, exdata->color.b);
			auto blend = [&](ExEdit::PixelYCA const& infl, ExEdit::PixelYCA const& orig) noexcept -> ExEdit::PixelYCA {
				int a = ((alpha * infl.a) >> log2_max_alpha) - orig.a;
				int A = static_cast<uint32_t>(f_alpha * orig.a) >> log2_max_alpha; // making sure A >= 0.
				if (a > 0) return {
					.y  = static_cast<i16>((a * col.y  + A * orig.y ) / (a + A)),
					.cb = static_cast<i16>((a * col.cb + A * orig.cb) / (a + A)),
					.cr = static_cast<i16>((a * col.cr + A * orig.cr) / (a + A)),
					.a  = static_cast<i16>(a + A),
				};
				else return {
					.y = orig.y, .cb = orig.cb, .cr = orig.cr,
					.a = static_cast<i16>(A),
				};
			};
			auto paint = [&](ExEdit::PixelYCA const& infl) noexcept -> ExEdit::PixelYCA {
				return {
					.y = col.y, .cb = col.cb, .cr = col.cr,
					.a = static_cast<i16>((alpha * infl.a) >> log2_max_alpha),
				};
			};

			multi_thread(dst_h, [&](int thread_id, int thread_num) {
				for (int y = thread_id; y < dst_h; y+= thread_num) {
					auto* dst = efpip->obj_edit + y * efpip->obj_line;
					if (y < result.displace || y >= dst_h - result.displace) {
						for (int x = dst_w; --x >= 0; dst++) *dst = paint(*dst);
					}
					else {
						auto* src = efpip->obj_temp + (y - result.displace) * efpip->obj_line;
						for (int x = result.displace; --x >= 0; dst++) *dst = paint(*dst);
						for (int x = dst_w - 2 * result.displace; --x >= 0; src++, dst++) *dst = blend(*dst, *src);
						for (int x = result.displace; --x >= 0; dst++) *dst = paint(*dst);
					}
				}
			});
		}
	}
	else {
		auto result = choose_defl(exdata->algorithm)(
			-lifted_size, neg_size, blur_px, param_a, false, false, efpip);
		if (result.invalid) return TRUE;

		if (tiled_image const img{ exdata->file, img_x, img_y, -result.displace, efp, *exedit.memory_ptr }) {
			// image seems to have been successfully loaded.
			// fill with the pattern image.
			auto blend = [&, img_stride = efpip->obj_line](i16 defl, ExEdit::PixelYCA const& orig, int px_x, int px_y) noexcept -> ExEdit::PixelYCA {
				auto& col = img[px_x + px_y * img_stride];
				int a = (alpha * (((max_alpha - defl) * col.a) >> log2_max_alpha)) >> log2_max_alpha,
					A = static_cast<uint32_t>(orig.a * (max_alpha - a)) >> log2_max_alpha; // making sure A >= 0.
				a = static_cast<uint32_t>(a * orig.a) >> log2_max_alpha; // making sure a >= 0.
				A = (f_alpha * A) >> log2_max_alpha;
				if (a > 0) return {
					.y  = static_cast<i16>((a * col.y  + A * orig.y ) / (a + A)),
					.cb = static_cast<i16>((a * col.cb + A * orig.cb) / (a + A)),
					.cr = static_cast<i16>((a * col.cr + A * orig.cr) / (a + A)),
					.a  = static_cast<i16>(a + A),
				};
				else return {
					.y = orig.y, .cb = orig.cb, .cr = orig.cr,
					.a = static_cast<i16>(A),
				};
			};
			auto paint = [&](ExEdit::PixelYCA const& orig, int px_x, int px_y) noexcept -> ExEdit::PixelYCA {
				return blend(0, orig, px_x, px_y);
			};

			if (result.is_empty) {
				// the entire image is re-colored.
				multi_thread(src_h, [&](int thread_id, int thread_num) {
					int const y0 = src_h * thread_id / thread_num, y1 = src_h * (thread_id + 1) / thread_num;
					int i_y = (y0 + img.oy) % img.h;
					auto incr_y = [&] { i_y++; if (i_y >= img.w) i_y -= img.h; };
					auto* dst_y = efpip->obj_edit + y0 * efpip->obj_line;
					for (int y = y1 - y0; --y >= 0; dst_y += efpip->obj_line, incr_y()) {
						int i_x = img.ox;
						auto incr_x = [&] { i_x++; if (i_x >= img.w) i_x -= img.w; };

						auto* dst = dst_y;
						for (int x = src_w; --x >= 0; dst++, incr_x()) *dst = paint(*dst, i_x, i_y);
					}
				});
			}
			else {
				multi_thread(src_h, [&, in_w = src_w - 2 * result.displace](int thread_id, int thread_num) {
					for (int y = thread_id; y < src_h; y += thread_num) {
						int i_x = img.ox, i_y = (y + img.oy) % img.h;
						auto incr_x = [&] { i_x++; if (i_x >= img.w) i_x -= img.w; };

						auto* dst = efpip->obj_edit + y * efpip->obj_line;
						if (y < result.displace || y >= src_h - result.displace) {
							for (int x = src_w; --x >= 0; dst++, incr_x()) *dst = paint(*dst, i_x, i_y);
						}
						else {
							auto* src = reinterpret_cast<i16*>(efpip->obj_temp) + (y - result.displace) * result.a_stride;
							for (int x = result.displace; --x >= 0; dst++, incr_x()) *dst = paint(*dst, i_x, i_y);
							for (int x = in_w; --x >= 0; src++, dst++, incr_x()) *dst = blend(*src, *dst, i_x, i_y);
							for (int x = result.displace; --x >= 0; dst++, incr_x()) *dst = paint(*dst, i_x, i_y);
						}
					}
				});
			}
		}
		else {
			// fill with specified color.
			auto col = buff::fromRGB(exdata->color.r, exdata->color.g, exdata->color.b);
			auto blend = [&](i16 defl, ExEdit::PixelYCA const& orig) noexcept -> ExEdit::PixelYCA {
				int a = (alpha * (max_alpha - defl)) >> log2_max_alpha,
					A = static_cast<uint32_t>(orig.a * (max_alpha - a)) >> log2_max_alpha; // making sure A >= 0.
				a = static_cast<uint32_t>(a * orig.a) >> log2_max_alpha; // making sure a >= 0.
				A = (f_alpha * A) >> log2_max_alpha;
				if (a > 0) return {
					.y  = static_cast<i16>((a * col.y  + A * orig.y ) / (a + A)),
					.cb = static_cast<i16>((a * col.cb + A * orig.cb) / (a + A)),
					.cr = static_cast<i16>((a * col.cr + A * orig.cr) / (a + A)),
					.a  = static_cast<i16>(a + A),
				};
				else return {
					.y = orig.y, .cb = orig.cb, .cr = orig.cr,
					.a = static_cast<i16>(A),
				};
			};
			auto paint = [&](ExEdit::PixelYCA const& orig) noexcept -> ExEdit::PixelYCA {
				return blend(0, orig);
			};

			if (result.is_empty) {
				// the entire image is re-colored.
				multi_thread(src_h, [&](int thread_id, int thread_num) {
					int const y0 = src_h * thread_id / thread_num, y1 = src_h * (thread_id + 1) / thread_num;
					auto* dst_y = efpip->obj_edit + y0 * efpip->obj_line;
					for (int y = y1 - y0; --y >= 0; dst_y += efpip->obj_line) {
						auto* dst = dst_y;
						for (int x = src_w; --x >= 0; dst++) *dst = paint(*dst);
					}
				});
			}
			else {
				multi_thread(src_h, [&, in_w = src_w - 2 * result.displace](int thread_id, int thread_num) {
					for (int y = thread_id; y < src_h; y += thread_num) {
						auto* dst = efpip->obj_edit + y * efpip->obj_line;
						if (y < result.displace || y >= src_h - result.displace) {
							for (int x = src_w; --x >= 0; dst++) *dst = paint(*dst);
						}
						else {
							auto* src = reinterpret_cast<i16*>(efpip->obj_temp) + (y - result.displace) * result.a_stride;
							for (int x = result.displace; --x >= 0; dst++) *dst = paint(*dst);
							for (int x = in_w; --x >= 0; src++, dst++) *dst = blend(*src, *dst);
							for (int x = result.displace; --x >= 0; dst++) *dst = paint(*dst);
						}
					}
				});
			}
		}
	}

	return TRUE;
}

