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

#include "kind_bin/inf_def.hpp"
#include "kind_bin2x/inf_def.hpp"
#include "kind_max/inf_def.hpp"
#include "kind_sum/inf_def.hpp"

#include "filter_defl.hpp"
#include "Rounding.hpp"

using namespace Filter::Rounding;
using namespace Calculation;
namespace Rounding_filter::params
{
	using namespace impl;
	namespace exp
	{
		static constexpr int
			den_radius	= track_den[idx_track::radius	],
			den_shrink	= track_den[idx_track::shrink	],
			den_transp	= track_den[idx_track::transp	],
			den_blur	= track_den[idx_track::blur		],
			den_param_a	= track_den[idx_track::param_a	],

			min_radius	= track_min[idx_track::radius	],
			min_shrink	= track_min[idx_track::shrink	],
			min_transp	= track_min[idx_track::transp	],
			min_blur	= track_min[idx_track::blur		],
			min_param_a	= track_min[idx_track::param_a	],

			max_radius	= track_max[idx_track::radius	],
			max_shrink	= track_max[idx_track::shrink	],
			max_transp	= track_max[idx_track::transp	],
			max_blur	= track_max[idx_track::blur		],
			max_param_a	= track_max[idx_track::param_a	];
	
		static_assert(den_radius == den_shrink);

		using defl_base = Filter::Common::defl_base<den_radius>;
	}
}
using namespace Rounding_filter::params::exp;


////////////////////////////////
// フィルタ処理の補助クラス．
////////////////////////////////
constexpr Filter::Common::defl_bin<den_radius, max_param_a> defl_bin{};
constexpr Filter::Common::defl_bin2x<den_radius, max_param_a> defl_bin2x{};
constexpr Filter::Common::defl_max<den_radius> defl_max{};
constexpr Filter::Common::defl_max_fast<den_radius> defl_max_fast{};
constexpr Filter::Common::defl_sum<den_radius, max_param_a> defl_sum{};

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


////////////////////////////////
// フィルタ処理のエントリポイント．
////////////////////////////////
BOOL impl::func_proc(ExEdit::Filter* efp, ExEdit::FilterProcInfo* efpip)
{
	int const src_w = efpip->obj_w, src_h = efpip->obj_h;
	if (src_w <= 0 || src_h <= 0) return TRUE;

	int const
		radius	= std::clamp(efp->track[idx_track::radius	], min_radius	, max_radius	),
		shrink	= std::clamp(efp->track[idx_track::shrink	], min_shrink	, max_shrink	),
		transp	= std::clamp(efp->track[idx_track::transp	], min_transp	, max_transp	),
		blur_px	= std::clamp(efp->track[idx_track::blur		], min_blur		, max_blur		),
		param_a	= std::clamp(efp->track[idx_track::param_a	], min_param_a	, max_param_a	);
	bool const crop = efp->check[idx_check::crop] != check_data::unchecked;
	auto const algorithm = reinterpret_cast<Exdata*>(efp->exdata_ptr)->algorithm;

	int const alpha = std::clamp<int>(max_alpha * (max_transp - transp) / max_transp, 0, max_alpha);

	// "fake" radius by 0.4 so integral-radius will more likely look smooth.
	int const lifted_radius = radius == 0 ? 0 : radius + ((den_radius >> 1) - 1);

	// handle trivial cases.
	if (!crop && alpha >= max_alpha) return TRUE; // entire effect is disabled.
	if (lifted_radius <= 0 && shrink <= 0 && blur_px <= 0) return TRUE; // none of the values are effective.

	// create the shape of alpha values onto efpip->obj_temp.
	auto result = choose_defl(algorithm)(shrink + blur_px, lifted_radius, blur_px, param_a, crop, crop, efpip);
	if (result.invalid) return TRUE;

	int const dst_w = src_w - 2 * result.displace, dst_h = src_h - 2 * result.displace;
	if (crop) {
		// final size shrinks.
		efpip->obj_w = dst_w; efpip->obj_h = dst_h;

		if (result.is_empty) {
			// simply clear the current image.
			if (dst_w <= 0 || dst_h <= 0) {
				efpip->obj_w = efpip->obj_h = 0;
				return FALSE; // invalidate subsequent filters.
			}
			buff::clear_alpha(efpip->obj_edit, efpip->obj_line, 0, 0, dst_w, dst_h);
			return TRUE;
		}

		// place color onto that shape.
		constexpr auto combine = [](i16 defl, ExEdit::PixelYCA const& src) noexcept -> ExEdit::PixelYCA {
			return {
				.y = src.y, .cb = src.cb, .cr = src.cr,
				.a = std::min(defl, src.a)
			};
		};
		multi_thread(dst_h, [&](int thread_id, int thread_num) {
			int y0 = dst_h * thread_id / thread_num,
				y1 = dst_h * (thread_id + 1) / thread_num;
			auto src = efpip->obj_edit + result.displace + (result.displace + y0) * efpip->obj_line,
				dst = efpip->obj_temp + y0 * efpip->obj_line;
			for (int y = y1 - y0; --y >= 0; src += efpip->obj_line, dst += efpip->obj_line) {
				auto src_y = src, dst_y = dst;
				for (int x = dst_w; --x >= 0; src_y++, dst_y++)
					*dst_y = combine(dst_y->a, *src_y);
			}
		});

		std::swap(efpip->obj_edit, efpip->obj_temp);
		return TRUE;
	}
	else {
		// final size keeps unchanged.

		if (result.is_empty) {
			// simply apply alpha value to the current image.
			buff::mult_alpha(alpha, efpip->obj_edit, efpip->obj_line, 0, 0, src_w, src_h);
			return TRUE;
		}

		// place that alpha values onto the current image.
		auto combine = [&](i16 defl, ExEdit::PixelYCA& dst) noexcept {
			auto diff_a = dst.a - defl;
			if (diff_a <= 0) return;
			dst.a = static_cast<i16>(defl + ((diff_a * alpha) >> log2_max_alpha));
		};
		auto decay = [&](ExEdit::PixelYCA& dst) noexcept {
			dst.a = static_cast<i16>((dst.a * alpha) >> log2_max_alpha);
		};
		multi_thread(src_h, [&](int thread_id, int thread_num) {
			for (int y = thread_id; y < src_h; y += thread_num) {
				auto dst = efpip->obj_edit + y * efpip->obj_line;
				if (y < result.displace || y >= dst_h + result.displace) {
					for (int x = src_w; --x >= 0; dst++) decay(*dst);
				}
				else {
					auto src = reinterpret_cast<i16*>(efpip->obj_temp) + (y - result.displace) * result.a_stride;
					for (int x = result.displace; --x >= 0; dst++) decay(*dst);
					for (int x = dst_w; --x >= 0; dst++, src++) combine(*src, *dst);
					for (int x = result.displace; --x >= 0; dst++) decay(*dst);
				}
			}
		});

		return TRUE;
	}
}

