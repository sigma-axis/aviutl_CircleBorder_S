/*
The MIT License (MIT)

Copyright (c) 2024 sigma-axis

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <cstdint>
#include <algorithm>
//#include <chrono>
//
//struct stopwatch {
//	using clock = std::chrono::high_resolution_clock;
//	clock::time_point start_at;
//	stopwatch() : start_at{ clock::now() } {}
//	void lap(const char* mes) const {
//		auto dur = clock::now() - start_at;
//		printf_s("%s: %lld\n", mes, dur.count());
//	}
//	void resume() { start_at = clock::now(); }
//};

#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
using byte = uint8_t;
#include <exedit/FilterProcInfo.hpp>

#include "arithmetics.hpp"
#include "multi_thread.hpp"
#include "buffer_op.hpp"

#include "kind_bin/inf_def.hpp"
#include "kind_bin2x/inf_def.hpp"
#include "kind_max/inf_def.hpp"
#include "kind_sum/inf_def.hpp"

#include "CircleBorder_S.hpp"

using i16 = int16_t;
using i32 = int32_t;


////////////////////////////////
// 内側縁取り・角丸め共通処理．
////////////////////////////////
namespace Filter::Common::impl
{
	using namespace Calculation;
	namespace exp {}
}
namespace Filter::Common { using namespace impl::exp; }
namespace Filter::Common::impl::exp
{
	template<size_t den_radius>
	struct defl_base {
	protected:
		struct process_spec {
			int sum_displace, neg_displace;
			bool do_defl, do_infl,
				allows_buffer_overlap; // whether intermediate result can overlap final result.
		};

		virtual process_spec tell_spec(int sum_size, int neg_size) const = 0;

		// copy alpha values from `efpip->obj_edit` to `efpip->obj_temp`, contiguous with non-zero.
		// returns the stride of the destination buffer (if not dst_colored).
		virtual void zero_op(int param_a, ExEdit::PixelYCA const* src_buf, size_t src_stride,
			int src_w, int src_h, i16* dst_buf, bool dst_colored, size_t dst_stride) const
		{
			if (dst_colored)
				buff::copy_alpha(src_buf, src_stride, 0, 0, src_w, src_h,
					buff::alpha_to_pixel(dst_buf), dst_stride / 4, 0, 0);
			else buff::copy_alpha(src_buf, src_stride, 0, 0, src_w, src_h,
				dst_buf, dst_stride, 0, 0);
		}
		virtual Bounds deflate_1(int sum_size_raw, int param_a, ExEdit::PixelYCA* src_buf, size_t src_stride,
			int src_w, int src_h, i16* dst_buf, bool dst_colored, size_t dst_stride, void* heap) const = 0;
		virtual Bounds deflate_2(int sum_size_raw, int param_a, ExEdit::PixelYCA* src_buf, size_t src_stride,
			int src_w, int src_h, i16* dst_buf, size_t dst_stride, void* heap, void* alpha_space) const = 0;
		virtual Bounds inflate_2(int neg_size_raw, int param_a, i16* src_buf, size_t src_stride,
			int src_w, int src_h, i16* dst_buf, bool dst_colored, size_t dst_stride, void* heap) const = 0;

	private:
		struct sizing {
			int sum_size_raw, sum_displace,
				neg_size_raw, neg_displace,
				blur_size_raw, blur_displace;
			bool do_defl, do_infl,
				allows_buffer_overlap;
			bool invalid;
		};
		sizing measure(int size, int neg_size, int blur_px) const
		{
			constexpr sizing invalid_size{ .invalid = true };
			if (size <= 0 && neg_size <= 0) return invalid_size;

			blur_px = std::clamp(blur_px, 0, size);
			auto sum_size = size - (blur_px >> 1) + neg_size;
			auto spec = tell_spec(sum_size, neg_size);
			auto blur_displace = buff::blur_displace((blur_px * buff::den_blur_px) / den_radius);
			return {
				.sum_size_raw = sum_size, .sum_displace = spec.sum_displace,
				.neg_size_raw = neg_size, .neg_displace = spec.neg_displace,
				.blur_size_raw = blur_px, .blur_displace = blur_displace,
				.do_defl = spec.do_defl, .do_infl = spec.do_infl,
				.allows_buffer_overlap = spec.allows_buffer_overlap,
			};
		}

	public:
		struct defl_result {
			int displace,
				a_stride; // will be 4*efpip->obj_line if `colored` is true.
			bool is_empty, // entire image turned transparent.
				colored;
			bool invalid;
		};
		defl_result operator()(int size, int neg_size, int blur_px, int param_a,
			bool dst_colored, bool tamely_diplace, ExEdit::FilterProcInfo* efpip) const
		{
			constexpr defl_result invalid{ .invalid = true };
			// calculate sizing values.
			auto const sz = measure(size, neg_size, blur_px);
			if (sz.invalid) return invalid;

			int const displace = sz.sum_displace - sz.neg_displace - sz.blur_displace,
				result_displace = tamely_diplace ?
					std::max(sz.sum_size_raw + (sz.blur_size_raw >> 1) - sz.neg_size_raw - sz.blur_size_raw, 0) / den_radius :
					std::max(displace, 0),
				diff_displace = displace - result_displace;

			Bounds bd{ 0, 0, efpip->obj_w, efpip->obj_h };
			if (!sz.do_defl && !sz.do_infl) {
				int const a_stride = dst_colored ? 4 * efpip->obj_line : (bd.wd() + 1) & (-2);
				zero_op(param_a, efpip->obj_edit, efpip->obj_line, bd.wd(), bd.ht(),
					dst_colored ? &efpip->obj_temp->a : reinterpret_cast<i16*>(efpip->obj_temp), dst_colored, a_stride);
				return { .displace = 0, .a_stride = a_stride, .colored = dst_colored };
			}
			else if (2 * sz.sum_displace >= std::min(efpip->obj_w, efpip->obj_h))
				return { .is_empty = true };
			else {
				// adjust the destination to handle with negative deflation.
				int const dst_stride = dst_colored ? 4 * efpip->obj_line : 
					((efpip->obj_w - 2 * displace + 1) & (-2)),
					dst_step = dst_colored ? 4 : 1;
				i16* const dst_buf = (dst_colored ?
					&efpip->obj_temp->a : reinterpret_cast<i16*>(efpip->obj_temp))
						+ diff_displace * (dst_step + dst_stride);

				if (sz.do_infl) {
					// allocate memory layout.
					size_t const med_stride = (bd.wd() - 2 * sz.sum_displace + 1) & (-2);
					i16* med_buffer; void* heap; void* alpha_space;
					if (sz.allows_buffer_overlap) {
						med_buffer = reinterpret_cast<i16*>(efpip->obj_temp);
						heap = *exedit.memory_ptr;
						alpha_space = nullptr;
					}
					else {
						med_buffer = reinterpret_cast<i16*>(*exedit.memory_ptr);
						heap = med_buffer + med_stride * (bd.ht() - 2 * sz.sum_displace);
						alpha_space = efpip->obj_temp;
					}

					// then process by two passes.
					if (sz.do_defl) {
						bd = deflate_2(sz.sum_size_raw, param_a,
							efpip->obj_edit, efpip->obj_line, bd.wd(), bd.ht(),
							med_buffer, med_stride, heap, alpha_space);
						if (bd.is_empty()) return {
							.displace = result_displace,
							.is_empty = true,
						};
					}
					else zero_op(param_a, efpip->obj_edit, efpip->obj_line, bd.wd(), bd.ht(),
						med_buffer, false, med_stride);
					bd = inflate_2(sz.neg_size_raw, param_a,
						med_buffer + bd.L + bd.T * med_stride, med_stride, bd.wd(), bd.ht(),
						dst_buf + bd.L * dst_step + bd.T * dst_stride, dst_colored, dst_stride, heap)
						.move(bd.L + diff_displace, bd.T + diff_displace);
				}
				else {
					// process by one pass.
					bd = deflate_1(sz.sum_size_raw, param_a,
						efpip->obj_edit, efpip->obj_line, bd.wd(), bd.ht(),
						dst_buf, dst_colored, dst_stride, *exedit.memory_ptr)
						.move(diff_displace, diff_displace);
				}
				if (bd.is_empty()) return {
					.displace = result_displace,
					.is_empty = true,
				};

				// apply blur.
				if (sz.blur_size_raw > 0) {
					if (dst_colored) buff::blur_alpha(efpip->obj_temp, efpip->obj_line,
						bd.L, bd.T, bd.wd(), bd.ht(), (blur_px * buff::den_blur_px) / den_radius);
					else buff::blur_alpha(reinterpret_cast<i16*>(efpip->obj_temp), dst_stride,
						bd.L, bd.T, bd.wd(), bd.ht(), (blur_px * buff::den_blur_px) / den_radius);

					bd = bd.inflate_br(2 * sz.blur_displace);
				}

				// clear the four sides of margins if present.
				int const dst_w = efpip->obj_w - 2 * result_displace, dst_h = efpip->obj_h - 2 * result_displace;
				if (dst_colored) buff::clear_alpha_chrome(efpip->obj_temp, efpip->obj_line,
					{ 0, 0, dst_w, dst_h }, bd);
				else buff::clear_alpha_chrome(reinterpret_cast<i16*>(efpip->obj_temp), dst_stride,
					{ 0, 0, dst_w, dst_h }, bd);

				return { .displace = result_displace, .a_stride = dst_stride, .colored = dst_colored };
			}
		}
	};

	// base class for "bin" and "bin2x".
	template<size_t den_radius, size_t max_param_a>
	struct defl_bin_base : defl_base<den_radius> {
	protected:
		using defl_base<den_radius>::process_spec;
		static constexpr i16 to_thresh(int param_a) {
			return static_cast<i16>((max_alpha - 1) * param_a / max_param_a);
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
		Bounds deflate_2(int sum_size_raw, int param_a, ExEdit::PixelYCA* src_buf, size_t src_stride,
			int src_w, int src_h, i16* dst_buf, size_t dst_stride, void* heap, void* alpha_space) const override
		{
			return bin::deflate(src_w, src_h,
				&src_buf->a, true, 4 * src_stride, to_thresh(param_a),
				dst_buf, false, dst_stride,
				heap, (sum_size_raw * sum_size_raw) / (den_radius * den_radius));
		}
	};

	// algorithm "bin".
	template<size_t den_radius, size_t max_param_a>
	struct defl_bin : defl_bin_base<den_radius, max_param_a> {
	protected:
		using defl_bin_base<den_radius, max_param_a>::process_spec;
		using defl_bin_base<den_radius, max_param_a>::to_thresh;

		process_spec tell_spec(int sum_size, int neg_size) const override
		{
			int sum_displace = bin::deflate_radius<den_radius>(sum_size),
				neg_displace = bin::inflate_radius<den_radius>(neg_size);
			return {
				.sum_displace = sum_displace,
				.neg_displace = neg_displace,
				.do_defl = sum_displace > 0,
				.do_infl = neg_displace > 0,
				.allows_buffer_overlap = true,
			};
		}

		Bounds deflate_1(int sum_size_raw, int param_a, ExEdit::PixelYCA* src_buf, size_t src_stride,
			int src_w, int src_h, i16* dst_buf, bool dst_colored, size_t dst_stride, void* heap) const override
		{
			return bin::deflate(src_w, src_h,
				&src_buf->a, true, 4 * src_stride, to_thresh(param_a),
				dst_buf, dst_colored, dst_stride,
				heap, (sum_size_raw * sum_size_raw) / (den_radius * den_radius));
		}
		Bounds inflate_2(int neg_size_raw, int param_a, i16* src_buf, size_t src_stride,
			int src_w, int src_h, i16* dst_buf, bool dst_colored, size_t dst_stride, void* heap) const override
		{
			return bin::inflate(src_w, src_h,
				src_buf, false, src_stride, 0,
				dst_buf, dst_colored, dst_stride,
				heap, (neg_size_raw * neg_size_raw) / (den_radius * den_radius));
		}
	};

	// algorithm "bin2x".
	template<size_t den_radius, size_t max_param_a>
	struct defl_bin2x : defl_bin_base<den_radius, max_param_a> {
	protected:
		using defl_bin_base<den_radius, max_param_a>::process_spec;
		using defl_bin_base<den_radius, max_param_a>::to_thresh;

		process_spec tell_spec(int sum_size, int neg_size) const override
		{
			bool do_infl = neg_size >= den_radius / 2;
			int sum_displace = (do_infl ? bin::deflate_radius<den_radius> : bin2x::deflate_radius<den_radius>)(sum_size),
				neg_displace = bin2x::inflate_radius<den_radius>(neg_size);
			return {
				.sum_displace = sum_displace,
				.neg_displace = neg_displace,
				.do_defl = (sum_size >= static_cast<int>(do_infl ? den_radius : den_radius / 2)),
				.do_infl = do_infl,
				.allows_buffer_overlap = true,
			};
		}

		Bounds deflate_1(int sum_size_raw, int param_a, ExEdit::PixelYCA* src_buf, size_t src_stride,
			int src_w, int src_h, i16* dst_buf, bool dst_colored, size_t dst_stride, void* heap) const override
		{
			return bin2x::deflate(src_w, src_h,
				&src_buf->a, true, 4 * src_stride, to_thresh(param_a),
				dst_buf, dst_colored, dst_stride,
				heap, (4 * sum_size_raw * sum_size_raw) / (den_radius * den_radius));
		}
		Bounds inflate_2(int neg_size_raw, int param_a, i16* src_buf, size_t src_stride,
			int src_w, int src_h, i16* dst_buf, bool dst_colored, size_t dst_stride, void* heap) const override
		{
			return bin2x::inflate(src_w, src_h,
				src_buf, false, src_stride, 0,
				dst_buf, dst_colored, dst_stride,
				heap, (4 * neg_size_raw * neg_size_raw) / (den_radius * den_radius));
		}
	};

	// algorithm "max".
	template<size_t den_radius>
	struct defl_max : defl_base<den_radius> {
	protected:
		using defl_base<den_radius>::process_spec;

		process_spec tell_spec(int sum_size, int neg_size) const override
		{
			int sum_displace = max::deflate_radius<den_radius>(sum_size),
				neg_displace = max::inflate_radius<den_radius>(neg_size);
			return {
				.sum_displace = sum_displace,
				.neg_displace = neg_displace,
				.do_defl = sum_displace > 0,
				.do_infl = neg_displace > 0,
				.allows_buffer_overlap = false,
			};
		}

		Bounds deflate_1(int sum_size_raw, int param_a, ExEdit::PixelYCA* src_buf, size_t src_stride,
			int src_w, int src_h, i16* dst_buf, bool dst_colored, size_t dst_stride, void* heap) const
		{
			return max::deflate(src_w, src_h,
				src_buf, src_stride,
				dst_buf, dst_colored, dst_stride,
				reinterpret_cast<void*>(reinterpret_cast<uintptr_t>(heap) + max::alpha_space_size(src_w, src_h)),
				(sum_size_raw * sum_size_raw) / (den_radius * den_radius), heap);
		}
		Bounds deflate_2(int sum_size_raw, int param_a, ExEdit::PixelYCA* src_buf, size_t src_stride,
			int src_w, int src_h, i16* dst_buf, size_t dst_stride, void* heap, void* alpha_space) const
		{
			return max::deflate(src_w, src_h, src_buf, src_stride,
				dst_buf, false, dst_stride,
				heap, (sum_size_raw * sum_size_raw) / (den_radius * den_radius), alpha_space);
		}
		Bounds inflate_2(int neg_size_raw, int param_a, i16* src_buf, size_t src_stride,
			int src_w, int src_h, i16* dst_buf, bool dst_colored, size_t dst_stride, void* heap) const
		{
			return max::inflate(src_w, src_h, src_buf, src_stride,
				dst_buf, dst_colored, dst_stride,
				heap, (neg_size_raw * neg_size_raw) / (den_radius * den_radius));
		}
	};

	// algorithm "sum".
	template<size_t den_radius, size_t max_param_a>
	struct defl_sum : defl_base<den_radius> {
		constexpr static int to_cap_rate(int param_a) {
			return sum::den_cap_rate * param_a / max_param_a;
		}

	protected:
		using defl_base<den_radius>::process_spec;

		process_spec tell_spec(int sum_size, int neg_size) const override
		{
			int sum_displace = sum::deflate_radius<den_radius>(sum_size),
				neg_displace = sum::inflate_radius<den_radius>(neg_size);
			return {
				.sum_displace = sum_displace,
				.neg_displace = neg_displace,
				.do_defl = sum_size >= den_radius,
				.do_infl = neg_size >= den_radius,
				.allows_buffer_overlap = false,
			};
		}

		Bounds deflate_1(int sum_size_raw, int param_a, ExEdit::PixelYCA* src_buf, size_t src_stride,
			int src_w, int src_h, i16* dst_buf, bool dst_colored, size_t dst_stride, void* heap) const
		{
			return sum::deflate(src_w, src_h, src_buf, src_stride,
				dst_buf, dst_colored, dst_stride, to_cap_rate(param_a),
				reinterpret_cast<void*>(reinterpret_cast<uintptr_t>(heap) + max::alpha_space_size(src_w, src_h)),
				(sum_size_raw * sum_size_raw) / (den_radius * den_radius), heap);
		}
		Bounds deflate_2(int sum_size_raw, int param_a, ExEdit::PixelYCA* src_buf, size_t src_stride,
			int src_w, int src_h, i16* dst_buf, size_t dst_stride, void* heap, void* alpha_space) const
		{
			return sum::deflate(src_w, src_h, src_buf, src_stride,
				dst_buf, false, dst_stride, to_cap_rate(param_a),
				heap, (sum_size_raw * sum_size_raw) / (den_radius * den_radius), alpha_space);
		}
		Bounds inflate_2(int neg_size_raw, int param_a, i16* src_buf, size_t src_stride,
			int src_w, int src_h, i16* dst_buf, bool dst_colored, size_t dst_stride, void* heap) const
		{
			return sum::inflate(src_w, src_h, src_buf, src_stride,
				dst_buf, dst_colored, dst_stride, to_cap_rate(param_a),
				heap, (neg_size_raw * neg_size_raw) / (den_radius * den_radius));
		}
	};
}

