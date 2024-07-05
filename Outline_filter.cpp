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

#include "arithmetics.hpp"
#include "multi_thread.hpp"
#include "buffer_op.hpp"
#include "tiled_image.hpp"

#include "kind_bin/inf_def.hpp"
#include "kind_bin2x/inf_def.hpp"
#include "kind_max/inf_def.hpp"
#include "kind_sum/inf_def.hpp"

#include "Outline.hpp"

using namespace Filter::Outline;
using namespace Calculation;
namespace Outline_filter::params
{
	using namespace impl;
	namespace exp
	{
		static constexpr int
			den_distance	= track_den[idx_track::distance	],
			den_thickness	= track_den[idx_track::thickness],
			den_pos_rad		= track_den[idx_track::pos_rad	],
			den_neg_rad		= track_den[idx_track::neg_rad	],
			den_blur		= track_den[idx_track::blur		],
			den_param_a		= track_den[idx_track::param_a	],
			den_img_x		= track_den[idx_track::img_x	],
			den_img_y		= track_den[idx_track::img_y	],

			min_distance	= track_min[idx_track::distance	],
			min_thickness	= track_min[idx_track::thickness],
			min_pos_rad		= track_min[idx_track::pos_rad	],
			min_neg_rad		= track_min[idx_track::neg_rad	],
			min_blur		= track_min[idx_track::blur		],
			min_param_a		= track_min[idx_track::param_a	],
			min_img_x		= track_min[idx_track::img_x	],
			min_img_y		= track_min[idx_track::img_y	],

			max_distance	= track_max[idx_track::distance	],
			max_thickness	= track_max[idx_track::thickness],
			max_pos_rad		= track_max[idx_track::pos_rad	],
			max_neg_rad		= track_max[idx_track::neg_rad	],
			max_blur		= track_max[idx_track::blur		],
			max_param_a		= track_max[idx_track::param_a	],
			max_img_x		= track_max[idx_track::img_x	],
			max_img_y		= track_max[idx_track::img_y	];
	
		static_assert(
			den_distance == den_thickness &&
			den_distance == den_pos_rad &&
			den_distance == den_neg_rad);

		using impl::FilterOrder;
	}
}
using namespace Outline_filter::params::exp;


////////////////////////////////
// フィルタ処理．
////////////////////////////////
struct outline_base {
protected:
	struct process_spec {
		int displace;
		bool valid;
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
	virtual std::tuple<int, int> max_size() const = 0;
	virtual process_spec tell_spec(int size_raw, bool is_final) const = 0;

	virtual Bounds inflate_med(int size_raw, int param_a, i16* src_buf, bool src_colored, size_t src_stride,
		int src_w, int src_h, i16* dst_buf, size_t dst_stride, void* heap) const = 0;
	virtual Bounds deflate_med(int size_raw, int param_a, i16* src_buf, bool src_colored, size_t src_stride,
		int src_w, int src_h, i16* dst_buf, size_t dst_stride, void* heap) const = 0;
	virtual Bounds inflate(int size_raw, int param_a, i16* src_buf, bool src_colored, size_t src_stride,
		int src_w, int src_h, i16* dst_buf, size_t dst_stride, void* heap) const {
		return inflate_med(size_raw, param_a, src_buf, src_colored, src_stride,
			src_w, src_h, dst_buf, dst_stride, heap);
	}
	virtual Bounds deflate(int size_raw, int param_a, i16* src_buf, bool src_colored, size_t src_stride,
		int src_w, int src_h, i16* dst_buf, size_t dst_stride, void* heap) const {
		return deflate_med(size_raw, param_a, src_buf, src_colored, src_stride,
			src_w, src_h, dst_buf, dst_stride, heap);
	}

private:
	struct sizing {
		int pass1_infl[3], pass1_displace[3], pass1_cnt,
			pass2_infl, pass2_displace,
			blur_size_raw, blur_displace,
			final_displace;
		bool has_hole, is_empty, zero_sized;
	};
	sizing measure(int distance_raw, int pos_rad_raw, int neg_rad_raw, int thick_raw, int blur_px,
		int src_w, int src_h, FilterOrder order) const
	{
		constexpr sizing zero_sized = { .zero_sized = true };

		auto const [max_w, max_h] = max_size();
		int const
			max_displace = std::min(max_w - src_w, max_h - src_h) >> 1,
			min_displace = -(std::min(src_w, src_h) >> 1);
		int const max_final_displace = std::min(exedit.yca_max_w - src_w, exedit.yca_max_h - src_h) >> 1;
		while (true) {
			sizing ret{};
			blur_px = std::clamp(blur_px, 0, std::max((arith::abs(thick_raw) >> 1) - den_distance, 0));
			static_assert(((-min_thickness) >> 1) - den_distance > max_blur);

			int const
				adj_dist_raw = distance_raw + (thick_raw >= 0 ? +1 : -1) * (blur_px >> 1),
				adj_thick_raw = thick_raw - (thick_raw >= 0 ? +1 : -1) * blur_px;

			int diff;
			// about the primary outline.
			{
				int min_diff_size, max_diff_size;
				if (adj_dist_raw >= 0) {
					min_diff_size = -std::max(pos_rad_raw - adj_dist_raw, 0);
					max_diff_size = neg_rad_raw + adj_dist_raw;
				}
				else {
					min_diff_size = -pos_rad_raw + adj_dist_raw;
					max_diff_size = std::max(neg_rad_raw + adj_dist_raw, 0);
				}
				int diff_sizes[3] = { min_diff_size, max_diff_size, adj_dist_raw };
				if (order == FilterOrder::defl_once || (order != FilterOrder::infl_once && adj_dist_raw < 0))
					std::swap(diff_sizes[0], diff_sizes[1]);

				// remove unnecessary inflate-deflate cycles.
				if (diff_sizes[1] == 0) diff_sizes[1] = diff_sizes[0];
				if (diff_sizes[0] == diff_sizes[2]) diff_sizes[0] = 0;

				// differentiate.
				for (int i = 3; --i > 0; ) diff_sizes[i] -= diff_sizes[i - 1];

				// precalculate sizes for each step.
				process_spec specs[3]; bool is_final = true;
				for (int i = 3; --i >= 0; ) {
					specs[i] = tell_spec(diff_sizes[i], is_final);
					is_final &= !specs[i].valid;
				}

				// extract valid operations, write them to the returned array.
				auto infl_1_push_alt = [&](int size_raw, int displace) {
					ret.pass1_infl[ret.pass1_cnt] = size_raw;
					ret.pass1_displace[ret.pass1_cnt] = displace;
					ret.final_displace += displace;
					ret.pass1_cnt++;
					return size_raw < 0 ?
						std::min(ret.final_displace - min_displace, 0) :
						std::max(ret.final_displace - max_displace, 0);
				};
				for (int i = 0; i < 3; i++) {
					if (specs[i].valid) {
						diff = infl_1_push_alt(diff_sizes[i], specs[i].displace);
						if (diff > 0) goto fail1;
						else if (diff < 0)
							// results in completely transparent image.
							// still continues calculating the final size.
							ret.is_empty = true;
					}
				}
			}

			// about line width.
			if (thick_raw <= min_thickness) ret.has_hole = false;
			else {
				auto spec = tell_spec(adj_thick_raw, true);
				ret.pass2_infl = adj_thick_raw;
				ret.pass2_displace = spec.displace;

				if (spec.valid) {
					if (ret.pass2_infl < 0)
						ret.has_hole = ret.final_displace + ret.pass2_displace >= min_displace;
					else {
						ret.has_hole = true;
						ret.final_displace += ret.pass2_displace;
						if (diff = ret.final_displace - max_displace;
							diff > 0) goto fail2;
					}
				}
				else ret.is_empty = true; // the line is zero-width, thus transparent.
			}

			// about blurring.
			ret.blur_size_raw = blur_px;
			ret.blur_displace = buff::blur_displace((ret.blur_size_raw * buff::den_blur_px) / den_blur);
			ret.final_displace += ret.blur_displace;
			if (diff = ret.final_displace - max_final_displace;
				diff > 0) goto fail1;

			// determine whether it's zero-sized.
			if (ret.final_displace < min_displace) return zero_sized;

			// it's OK.
			return ret;

		fail1:
			// trim the primary size.
			distance_raw -= diff * den_distance;
			continue;

		fail2:
			// trim the secondary size.
			thick_raw -= diff * den_distance;
			continue;
		}
	}

public:
	struct outline_result {
		int displace, stride;
		Bounds bd;
		bool is_empty, zero_sized;
	};

	outline_result operator()(int distance_raw, int pos_rad_raw, int neg_rad_raw, int thick_raw, int blur_px, int param_a,
		FilterOrder order, ExEdit::FilterProcInfo* efpip) const
	{
		constexpr outline_result zero_sized{ .zero_sized = true };
		int const src_w = efpip->obj_w, src_h = efpip->obj_h;

		auto sz = measure(distance_raw, pos_rad_raw, neg_rad_raw, thick_raw, blur_px, src_w, src_h, order);
		if (sz.zero_sized) return zero_sized; // zero-sized.
		if (sz.is_empty) {
			// valid size, filled with transparent pixels.
			return {
				.displace = sz.final_displace,
				.is_empty = true
			};
		}

		outline_result ret{
			.displace = sz.final_displace,
			.stride = (src_w + 2 * sz.final_displace + 1) & (-2),
			.bd = { 0, 0, src_w, src_h },
		};

		size_t stride = 4 * efpip->obj_line;
		// first pass to create a curve at the specified distance and curvatures.
		if (sz.pass1_cnt == 0) {
			stride = (efpip->obj_w + 1) & (-2);
			buff::copy_alpha(efpip->obj_edit, efpip->obj_line, 0, 0, src_w, src_h,
				reinterpret_cast<i16*>(efpip->obj_temp), stride, 0, 0);
			std::swap(efpip->obj_edit, efpip->obj_temp);
		}
		else {
			for (int i = 0; i < sz.pass1_cnt; i++) {
				bool dst_final = i == sz.pass1_cnt - 1;
				size_t src_stride = std::exchange(stride,
					dst_final && (!sz.has_hole || sz.pass2_infl < 0) ?
					ret.stride : ((ret.bd.wd() + 2 * sz.pass1_displace[i] + 1) & (-2)));
				infdef(sz.pass1_infl[i], sz.pass1_displace[i], param_a,
					i == 0, src_stride, dst_final, stride, ret.bd, efpip);
				if (ret.bd.is_empty()) return {
					.displace = sz.final_displace,
					.is_empty = true,
				};
			}
		}

		// carving the hole by the second pass.
		if (sz.has_hole) {
			Bounds bd2 = ret.bd;
			size_t stride2 = sz.pass2_infl > 0 ?
				ret.stride : (ret.bd.wd() + 2 * sz.pass2_displace + 1) & (-2);
			infdef(sz.pass2_infl, sz.pass2_displace, param_a, false, stride, true, stride2, bd2, efpip);
			if (bd2.is_empty()) {
				if (sz.pass2_infl > 0) return {
					.displace = sz.final_displace,
					.is_empty = true,
				};
			}
			else {
				// there is a hole.
				// move the larger image to obj_temp.
				if (sz.pass2_infl > 0) {
					std::swap(efpip->obj_edit, efpip->obj_temp);
					std::swap(ret.bd, bd2);
					std::swap(stride, stride2);
				}

				// carve the image.
				multi_thread(bd2.ht(), [&](int thread_id, int thread_num) {
					int y0 = bd2.T + bd2.ht() * thread_id / thread_num,
						y1 = bd2.T + bd2.ht() * (thread_id + 1) / thread_num;
					auto src = reinterpret_cast<i16*>(efpip->obj_edit)
							+ bd2.L + y0 * stride2,
						dst = reinterpret_cast<i16*>(efpip->obj_temp)
							+ bd2.L + y0 * stride
							+ arith::abs(sz.pass2_displace) * (1 + stride);

					for (int y = y1 - y0; --y >= 0; src += stride2, dst += stride) {
						auto src_y = src, dst_y = dst;
						for (int x = bd2.wd(); --x >= 0; src_y++, dst_y++)
							*dst_y = std::max(*dst_y - *src_y, 0);
					}
				});
			}
			std::swap(efpip->obj_edit, efpip->obj_temp);
		}

		// apply blur.
		if (sz.blur_size_raw > 0) {
			buff::blur_alpha(reinterpret_cast<i16*>(efpip->obj_edit), ret.stride,
				ret.bd.L, ret.bd.T, ret.bd.wd(), ret.bd.ht(), (sz.blur_size_raw * buff::den_blur_px) / den_blur);
			ret.bd = ret.bd.inflate_br(2 * sz.blur_displace, 2 * sz.blur_displace);
		}

		return ret;
	}

private:
	void infdef(int size, int displace, int param_a, bool src_colored, size_t src_stride, bool dst_final, size_t dst_stride,
		Bounds& bd, ExEdit::FilterProcInfo* efpip) const
	{
		auto const
			src_buf = (src_colored ? &efpip->obj_edit->a : reinterpret_cast<i16*>(efpip->obj_edit))
				+ bd.L * (src_colored ? 4 : 1) + bd.T * src_stride,
			dst_buf = reinterpret_cast<i16*>(efpip->obj_temp)
				+ bd.L + bd.T * dst_stride;

		bd = (this->*(size > 0 ?
			dst_final ? &outline_base::inflate : &outline_base::inflate_med :
			dst_final ? &outline_base::deflate : &outline_base::deflate_med))(
				size, param_a, src_buf, src_colored, src_stride,
				bd.wd(), bd.ht(), dst_buf, dst_stride, *exedit.memory_ptr)
			.move(bd.L, bd.T);

		std::swap(efpip->obj_edit, efpip->obj_temp);
	}
};

struct outline_bin_base : outline_base {
protected:
	static constexpr i16 to_thresh(int param_a) {
		return (max_alpha - 1) * param_a / max_param_a;
	}
	process_spec tell_spec(int size_raw, bool is_final) const override {
		int const displace = size_raw > 0 ?
			+bin::inflate_radius<den_distance>(+size_raw) :
			-bin::deflate_radius<den_distance>(-size_raw);
		return {
			.displace = displace,
			.valid = displace != 0,
		};
	}
	Bounds inflate_med(int size_raw, int param_a, i16* src_buf, bool src_colored, size_t src_stride,
		int src_w, int src_h, i16* dst_buf, size_t dst_stride, void* heap) const {
		return bin::inflate(src_w, src_h, src_buf,
			src_colored, src_stride, to_thresh(param_a),
			dst_buf, false, dst_stride, *exedit.memory_ptr, (size_raw * size_raw) / (den_distance * den_distance));
	}
	Bounds deflate_med(int size_raw, int param_a, i16* src_buf, bool src_colored, size_t src_stride,
		int src_w, int src_h, i16* dst_buf, size_t dst_stride, void* heap) const {
		return bin::deflate(src_w, src_h, src_buf,
			src_colored, src_stride, to_thresh(param_a),
			dst_buf, false, dst_stride, *exedit.memory_ptr, (size_raw * size_raw) / (den_distance * den_distance));
	}
};

// algorithm "bin".
constexpr struct : outline_bin_base {
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
	std::tuple<int, int> max_size() const override {
		if (mem_max_w == 0) init_mem_max();
		return { mem_max_w, mem_max_h };
	}
} outline_bin{};

// algorithm "bin2x".
constexpr struct : outline_bin_base {
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
	std::tuple<int, int> max_size() const override {
		if (mem_max_w == 0) init_mem_max();
		return { mem_max_w, mem_max_h };
	}
	process_spec tell_spec(int size_raw, bool is_final) const override {
		if (is_final) {
			int const displace = size_raw > 0 ?
				+bin2x::inflate_radius<den_distance>(+size_raw) :
				-bin2x::deflate_radius<den_distance>(-size_raw);
			return {
				.displace = displace,
				.valid = size_raw >= (den_distance >> 1) || size_raw <= -(den_distance >> 1),
			};
		}
		else return outline_bin_base::tell_spec(size_raw, is_final);
	}
	Bounds inflate(int size_raw, int param_a, i16* src_buf, bool src_colored, size_t src_stride,
		int src_w, int src_h, i16* dst_buf, size_t dst_stride, void* heap) const override {
		return bin2x::inflate(src_w, src_h, src_buf,
			src_colored, src_stride, to_thresh(param_a),
			dst_buf, false, dst_stride, *exedit.memory_ptr, (4 * size_raw * size_raw) / (den_distance * den_distance));
	}
	Bounds deflate(int size_raw, int param_a, i16* src_buf, bool src_colored, size_t src_stride,
		int src_w, int src_h, i16* dst_buf, size_t dst_stride, void* heap) const override {
		return bin2x::deflate(src_w, src_h, src_buf,
			src_colored, src_stride, to_thresh(param_a),
			dst_buf, false, dst_stride, *exedit.memory_ptr, (4 * size_raw * size_raw) / (den_distance * den_distance));
	}
} outline_bin2x{};

// algorithm "max".
constexpr struct : outline_base {
private:
	static inline constinit int mem_max_w = 0, mem_max_h = 0;
	static void init_mem_max()
	{
		size_t const mem_max = obj_mem_max();
		std::tie(mem_max_w, mem_max_h) = max_size_cand(mem_max / sizeof(i16));

		// trim by 4 dots until it fits within available space.
		while (sizeof(i16) * ((mem_max_w + 1) & (-2)) * mem_max_h > mem_max ||
			max::inflate_heap_size(mem_max_w, mem_max_h, std::min(mem_max_w, mem_max_h) >> 1) > mem_max)
			mem_max_w -= 4, mem_max_h -= 4;
	}

protected:
	std::tuple<int, int> max_size() const override {
		if (mem_max_w == 0) init_mem_max();
		return { mem_max_w, mem_max_h };
	}
	process_spec tell_spec(int size_raw, bool is_final) const override {
		int const displace = size_raw > 0 ?
			+max::inflate_radius<den_distance>(+size_raw) :
			-max::deflate_radius<den_distance>(-size_raw);
		return {
			.displace = displace,
			.valid = displace != 0,
		};
	}
	Bounds inflate_med(int size_raw, int param_a, i16* src_buf, bool src_colored, size_t src_stride,
		int src_w, int src_h, i16* dst_buf, size_t dst_stride, void* heap) const {
		return max::inflate(src_w, src_h, src_buf,
			src_colored, src_stride,
			dst_buf, false, dst_stride, *exedit.memory_ptr, (size_raw * size_raw) / (den_distance * den_distance));
	}
	Bounds deflate_med(int size_raw, int param_a, i16* src_buf, bool src_colored, size_t src_stride,
		int src_w, int src_h, i16* dst_buf, size_t dst_stride, void* heap) const {
		return max::deflate(src_w, src_h, src_buf,
			src_colored, src_stride,
			dst_buf, false, dst_stride, *exedit.memory_ptr, (size_raw * size_raw) / (den_distance * den_distance));
	}
} outline_max{};

// algorithm "sum".
constexpr struct : outline_base {
	constexpr static int to_cap_rate(int param_a) {
		return sum::den_cap_rate * param_a / max_param_a;
	}

private:
	static inline constinit int mem_max_w = 0, mem_max_h = 0;
	static void init_mem_max()
	{
		size_t const mem_max = obj_mem_max();
		std::tie(mem_max_w, mem_max_h) = max_size_cand(mem_max / sizeof(i16));

		// trim by 4 dots until it fits within available space.
		while (sizeof(i16) * ((mem_max_w + 1) & (-2)) * mem_max_h > mem_max ||
			sum::inflate_heap_size(mem_max_w, mem_max_h, std::min(mem_max_w, mem_max_h) >> 1) > mem_max)
			mem_max_w -= 4, mem_max_h -= 4;
	}

protected:
	std::tuple<int, int> max_size() const override {
		if (mem_max_w == 0) init_mem_max();
		return { mem_max_w, mem_max_h };
	}
	process_spec tell_spec(int size_raw, bool is_final) const override {
		int const displace = size_raw > 0 ?
			+sum::inflate_radius<den_distance>(+size_raw) :
			-sum::deflate_radius<den_distance>(-size_raw);
		return {
			.displace = displace,
			.valid = size_raw >= den_distance || size_raw <= -den_distance,
		};
	}
	Bounds inflate_med(int size_raw, int param_a, i16* src_buf, bool src_colored, size_t src_stride,
		int src_w, int src_h, i16* dst_buf, size_t dst_stride, void* heap) const {
		return sum::inflate(src_w, src_h, src_buf,
			src_colored, src_stride,
			dst_buf, false, dst_stride, to_cap_rate(param_a),
			*exedit.memory_ptr, (size_raw * size_raw) / (den_distance * den_distance));
	}
	Bounds deflate_med(int size_raw, int param_a, i16* src_buf, bool src_colored, size_t src_stride,
		int src_w, int src_h, i16* dst_buf, size_t dst_stride, void* heap) const {
		return sum::deflate(src_w, src_h, src_buf,
			src_colored, src_stride,
			dst_buf, false, dst_stride, to_cap_rate(param_a),
			*exedit.memory_ptr, (size_raw * size_raw) / (den_distance * den_distance));
	}
} outline_sum{};

constexpr outline_base const& choose_outline(Filter::Algorithm algorithm) {
	switch (algorithm) {
		using algo = Filter::Algorithm;
	case algo::bin: return outline_bin;
	default:
	case algo::bin2x: return outline_bin2x;
	case algo::max: return outline_max;
	case algo::sum: return outline_sum;
	}
}


////////////////////////////////
// フィルタ処理のエントリポイント．
////////////////////////////////
BOOL impl::func_proc(ExEdit::Filter* efp, ExEdit::FilterProcInfo* efpip)
{
	int const src_w = efpip->obj_w, src_h = efpip->obj_h;
	if (src_w <= 0 || src_h <= 0) return TRUE;

	int const
		distance	= std::clamp(efp->track[idx_track::distance	], min_distance	, max_distance	),
		thickness	= std::clamp(efp->track[idx_track::thickness], min_thickness, max_thickness	),
		pos_rad		= std::clamp(efp->track[idx_track::pos_rad	], min_pos_rad	, max_pos_rad	),
		neg_rad		= std::clamp(efp->track[idx_track::neg_rad	], min_neg_rad	, max_neg_rad	),
		blur_px		= std::clamp(efp->track[idx_track::blur		], min_blur		, max_blur		),
		param_a		= std::clamp(efp->track[idx_track::param_a	], min_param_a	, max_param_a	),
		img_x		= std::clamp(efp->track[idx_track::img_x	], min_img_x	, max_img_x		),
		img_y		= std::clamp(efp->track[idx_track::img_y	], min_img_y	, max_img_y		);
	auto* const exdata = reinterpret_cast<Exdata*>(efp->exdata_ptr);

	// handle trivial cases.
	if (std::min(src_w, src_h) <= 2 * ((-distance - std::max(thickness, 0)) / den_distance)) {
		// should turn empty.
		efpip->obj_w = efpip->obj_h = 0;
		return FALSE; // invalidate subsequent filters.
	}

	auto result = choose_outline(exdata->algorithm)(
		distance, pos_rad, neg_rad, thickness,
		blur_px, param_a, exdata->order, efpip);
	if (result.zero_sized) {
		// should turn empty.
		efpip->obj_w = efpip->obj_h = 0;
		return FALSE; // invalidate subsequent filters.
	}

	int const
		displace = arith::floor_div(distance + (den_distance - 1), den_distance)
			+ (thickness >= 0 ? (thickness + (den_distance - 1)) / den_distance : 0),
		dst_w = efpip->obj_w += 2 * displace,
		dst_h = efpip->obj_h += 2 * displace;

	if (result.is_empty) {
		// resized but completely transparent.
		buff::clear_alpha(efpip->obj_edit, efpip->obj_line, 0, 0, dst_w, dst_h);
		return TRUE;
	}

	int const diff_displace = displace - result.displace,
		T = result.bd.T + diff_displace, B = result.bd.B + diff_displace,
		L = std::max(result.bd.L + diff_displace, 0), R = std::min(result.bd.R + diff_displace, dst_w),
		r = dst_w - R, in_w = R - L;
	i16 const* const src0 = reinterpret_cast<i16*>(efpip->obj_edit)
		+ L - diff_displace * (1 + result.stride);
	if (tiled_image const img{ exdata->file, img_x, img_y, displace, efp, *exedit.memory_ptr }) {
		// image seems to have been successfully loaded.
		// fill with the pattern image.
		auto paint = [&, img_stride = efpip->obj_line](i16 src_a, int px_x, int px_y) noexcept -> ExEdit::PixelYCA {
			ExEdit::PixelYCA col = img[px_x + px_y * img_stride];
			col.a = static_cast<i16>((col.a * src_a) >> log2_max_alpha);
			return col;
		};

		multi_thread(dst_h, [&](int thread_id, int thread_num) {
			for (int y = thread_id; y < dst_h; y += thread_num) {
				auto* dst = efpip->obj_temp + y * efpip->obj_line;
				if (y < T || y >= B) {
					for (int x = dst_w; --x >= 0; dst++) dst->a = 0;
				}
				else {
					int i_x = (img.ox + L) % img.w, i_y = (y + img.oy) % img.h;
					auto incr_x = [&] { i_x++; if (i_x >= img.w) i_x -= img.w; };

					auto* src = src0 + y * result.stride;
					for (int x = L; --x >= 0; dst++) dst->a = 0;
					for (int x = in_w; --x >= 0; dst++, src++, incr_x()) *dst = paint(*src, i_x, i_y);
					for (int x = r; --x >= 0; dst++) dst->a = 0;
				}
			}
		});
	}
	else {
		// fill with specified color.
		auto const col = buff::fromRGB(exdata->color.r, exdata->color.g, exdata->color.b);
		auto paint = [&](i16 src_a) noexcept -> ExEdit::PixelYCA {
			return { .y = col.y, .cb = col.cb, .cr = col.cr, .a = src_a };
		};

		multi_thread(/*dst_h*/true, [&](int thread_id, int thread_num) {
			for (int y = thread_id; y < dst_h; y += thread_num) {
				auto* dst = efpip->obj_temp + y * efpip->obj_line;
				if (y < T || y >= B) {
					for (int x = dst_w; --x >= 0; dst++) dst->a = 0;
				}
				else {
					auto* src = src0 + y * result.stride;
					for (int x = L; --x >= 0; dst++) dst->a = 0;
					for (int x = in_w; --x >= 0; dst++, src++) *dst = paint(*src);
					for (int x = r; --x >= 0; dst++) dst->a = 0;
				}
			}
		});
	}
	std::swap(efpip->obj_edit, efpip->obj_temp);

	return TRUE;
}

