/*
The MIT License (MIT)

Copyright (c) 2024 sigma-axis

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <cstdint>
#include <numeric>

#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>

using byte = uint8_t;
#include <exedit.hpp>

#include "CircleBorder_S.hpp"


////////////////////////////////
// 仕様書．
////////////////////////////////
namespace Filter::Border
{
	namespace impl
	{
		FILTER_INFO("縁取りσ");

		// trackbars.
		constexpr const char* track_names[]
			= { "サイズ", "凹半径", "透明度", "内透明度", "ぼかし", "param_a", "画像X", "画像Y" };
		constexpr int32_t
			track_den[]			= {    10,   10,   10,   10,   10,   10,     1,     1 },
			track_min[]			= { -5000,    0,    0,    0,    0,    0, -4000, -4000 },
			track_min_drag[]	= { -2000,    0,    0,    0,    0,    0, -1000, -1000 },
			track_default[]		= {    50,    0,    0,    0,    0,  500,     0,     0 },
			track_max_drag[]	= {  2000, 2000, 1000, 1000, 1000, 1000,  1000,  1000 },
			track_max[]			= {  5000, 5000, 1000, 1000, 1000, 1000,  4000,  4000 };
		constexpr int track_link[] = { 0, 0, 0, 0, 0, 0, 1, -1, };
		namespace idx_track
		{
			enum id : int {
				size,
				neg_size,
				transp,
				f_transp,
				blur,
				param_a,
				img_x,
				img_y,
			};
			constexpr int count_entries = std::size(track_names);
		};

		static_assert(
			std::size(track_names) == std::size(track_den) &&
			std::size(track_names) == std::size(track_min) &&
			std::size(track_names) == std::size(track_min_drag) &&
			std::size(track_names) == std::size(track_default) &&
			std::size(track_names) == std::size(track_max_drag) &&
			std::size(track_names) == std::size(track_max) &&
			std::size(track_names) == std::size(track_link));

		// checks.
		constexpr const char* check_names[]
			= { gui::algorithm_names, "縁色の設定", "パターン画像ファイル" };
		constexpr int32_t
			check_default[] = { check_data::dropdown, check_data::button, check_data::button };
		namespace idx_check {
			enum id : int {
				algorithm,
				color,
				file,
			};
			constexpr int count_entries = std::size(check_names);
		};

		static_assert(std::size(check_names) == std::size(check_default));

		// exdata.
		constexpr ExEdit::ExdataUse exdata_use[] =
		{
			{ .type = ExEdit::ExdataUse::Type::Number, .size = sizeof(Algorithm), .name = "kind" },
			{ .type = ExEdit::ExdataUse::Type::Binary, .size = 3, .name = "color" },
			{ .type = ExEdit::ExdataUse::Type::Padding, .size = 1, .name = nullptr },
			{ .type = ExEdit::ExdataUse::Type::String, .size = 256, .name = "file" },
		};
		namespace idx_data
		{
			namespace _impl
			{
				static consteval size_t idx(auto name) {
					auto ret = std::find_if(std::begin(exdata_use), std::end(exdata_use),
						[name](auto d) { return d.name != nullptr && std::string_view{ d.name } == name; }) - exdata_use;
					if (ret < std::size(exdata_use)) return ret;
					std::unreachable();
				}
			}
			enum id : int {
				algorithm = _impl::idx("kind"),
				color = _impl::idx("color"),
				file = _impl::idx("file"),
			};
			constexpr int count_entries = 3;
		}
		//#pragma pack(push, 1)
		struct Exdata {
			Algorithm algorithm = Algorithm::bin2x;
			ExEdit::Exdata::ExdataColor color{ .r = 0, .g = 0, .b = 0 };
			char file[exdata_use[idx_data::file].size]{};
		};

		//#pragma pack(pop)
		constexpr static Exdata exdata_def = {};

		static_assert(sizeof(Exdata) == std::accumulate(
			std::begin(exdata_use), std::end(exdata_use), size_t{ 0 }, [](auto v, auto d) { return v + d.size; }));

		// callbacks.
		BOOL func_proc(ExEdit::Filter* efp, ExEdit::FilterProcInfo* efpip);
		BOOL func_WndProc(HWND hwnd, UINT message, WPARAM wparam, LPARAM lparam, AviUtl::EditHandle* editp, ExEdit::Filter* efp);
		int32_t func_window_init(HINSTANCE hinstance, HWND hwnd, int y, int base_id, int sw_param, ExEdit::Filter* efp);
		static inline BOOL func_init(ExEdit::Filter* efp) { exedit.init(efp->exedit_fp); return TRUE; }
	}

	inline constinit ExEdit::Filter filter = {
		.flag = ExEdit::Filter::Flag::Effect,
		.name = const_cast<char*>(impl::filter_name),
		.track_n = std::size(impl::track_names),
		.track_name = const_cast<char**>(impl::track_names),
		.track_default = const_cast<int*>(impl::track_default),
		.track_s = const_cast<int*>(impl::track_min),
		.track_e = const_cast<int*>(impl::track_max),
		.check_n = std::size(impl::check_names),
		.check_name = const_cast<char**>(impl::check_names),
		.check_default = const_cast<int*>(impl::check_default),
		.func_proc = &impl::func_proc,
		.func_init = &impl::func_init,
		.func_WndProc = &impl::func_WndProc,
		.exdata_size = sizeof(impl::exdata_def),
		.information = const_cast<char*>(impl::info),
		.func_window_init = &impl::func_window_init,
		.exdata_def = const_cast<impl::Exdata*>(&impl::exdata_def),
		.exdata_use = impl::exdata_use,
		.track_scale = const_cast<int*>(impl::track_den),
		.track_link = const_cast<int*>(impl::track_link),
		.track_drag_min = const_cast<int*>(impl::track_min_drag),
		.track_drag_max = const_cast<int*>(impl::track_max_drag),
	};
}

