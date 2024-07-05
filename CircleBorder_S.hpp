/*
The MIT License (MIT)

Copyright (c) 2024 sigma-axis

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <cstdint>
#include <algorithm>
#include <numeric>

#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <commdlg.h>

using byte = uint8_t;
#include <exedit.hpp>


////////////////////////////////
// 主要情報源の変数アドレス．
////////////////////////////////
inline constinit struct ExEdit092 {
	void init(AviUtl::FilterPlugin* efp)
	{
		if (fp == nullptr)
			init_pointers(efp);
	}
	AviUtl::FilterPlugin* fp = nullptr;

	HWND*	hwnd_setting_dlg;	// 0x1539c8
	AviUtl::EditHandle** editp;	// 0x1a532c
	void**	memory_ptr;			// 0x1a5328 // 少なくとも最大画像サイズと同サイズは保証されるっぽい．

	int32_t	yca_max_w;			// 0x196748
	int32_t	yca_max_h;			// 0x1920e0

	void(*update_any_exdata)(ExEdit::ObjectFilterIndex processing, const char* exdata_use_name);	// 0x04a7e0
	void(*nextundo)();			// 0x08d150

private:
	void init_pointers(AviUtl::FilterPlugin* efp);
} exedit{};


////////////////////////////////
// フィルタ共通．
////////////////////////////////
namespace Filter
{
	struct check_data {
		enum def : int32_t {
			unchecked = 0,
			checked = 1,
			button = -1,
			dropdown = -2,
		};
	};

	enum class Algorithm : uint32_t {
		bin = 0,
		bin2x = 1,
		sum = 2,
		max = 3,
	};
	constexpr int algorithm_count = 4;

	namespace gui
	{
		constexpr auto algorithm_names = "2値化\0002値化倍精度\0総和\0最大値\0";
		constexpr auto algorithm_caption = L"方式",
			param_a_name = L"αしきい値", param_a_name_alt = L"基準α和",
			invalid_name = L"----";
		constexpr auto choose_param_a_name(Algorithm algorithm) {
			switch (algorithm) {
				using enum Algorithm;
			case bin: case bin2x: return param_a_name;
			case sum: return param_a_name_alt;
			default: return invalid_name;
			}
		}
		constexpr auto invalid_name_narrow = "----";

		constexpr auto color_format = L"RGB ( %d , %d , %d )";
		constexpr size_t size_col_fmt = std::wstring_view{ color_format }.size() + 1
			+ 3 * (std::size(L"255")- std::size(L"%d"));
	}
}


////////////////////////////////
// バージョン情報．
////////////////////////////////
#define PLUGIN_VERSION	"v1.00"
#define PLUGIN_AUTHOR	"sigma-axis"
#define FILTER_INFO_FMT(name, ver, author)	(name##" "##ver##" by "##author)
#define FILTER_INFO(name)	constexpr char filter_name[] = name, info[] = FILTER_INFO_FMT(name, PLUGIN_VERSION, PLUGIN_AUTHOR)

