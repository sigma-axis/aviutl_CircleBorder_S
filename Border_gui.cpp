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
#include <Windows.h>
#include <commdlg.h>

using byte = uint8_t;
#include <exedit/Filter.hpp>

#include "Border.hpp"
#include "relative_path.hpp"

using namespace Filter;
using namespace Border;


////////////////////////////////
// ウィンドウ状態の保守．
////////////////////////////////
/*
efp->exfunc->get_hwnd(efp->processing, i, j):
	i = 0:		j 番目のスライダーの中央ボタン．
	i = 1:		j 番目のスライダーの左トラックバー．
	i = 2:		j 番目のスライダーの右トラックバー．
	i = 3:		j 番目のチェック枠のチェックボックス．
	i = 4:		j 番目のチェック枠のボタン．
	i = 5, 7:	j 番目のチェック枠の右にある static (テキスト).
	i = 6:		j 番目のチェック枠のコンボボックス．
	otherwise -> nullptr.
*/

static inline void update_extendedfilter_wnd(ExEdit::Filter* efp)
{
	using namespace Border::impl;
	auto* exdata = reinterpret_cast<Exdata*>(efp->exdata_ptr);

	// choose button text from "αしきい値", "基準α和", or "----".
	::SetWindowTextW(efp->exfunc->get_hwnd(efp->processing, 0, idx_track::param_a),
		gui::choose_param_a_name(exdata->algorithm));

	// select item in the combo box for algorithm.
	if (int idx = static_cast<int>(exdata->algorithm); 0 <= idx || idx < algorithm_count) {
		::SendMessageW(efp->exfunc->get_hwnd(efp->processing, 6, idx_check::algorithm), CB_SETCURSEL, idx, 0);
	}

	// whether the background pattern image is specified, or single color.
	wchar_t col_fmt[gui::size_col_fmt] = L"";
	auto file = "", img_x = track_names[idx_track::img_x], img_y = track_names[idx_track::img_y];
	if (exdata->file[0] == '\0') {
		// single color.
		::swprintf_s(col_fmt, gui::color_format,
			exdata->color.r, exdata->color.g, exdata->color.b);
		img_x = img_y = gui::invalid_name_narrow;
	}
	else {
		// background image.
		std::end(exdata->file)[-1] = '\0';
		file = relative_path::ptr_file_name(exdata->file);
	}

	// choose button text from "画像X/Y" or "----".
	::SetWindowTextA(efp->exfunc->get_hwnd(efp->processing, 0, idx_track::img_x), img_x);
	::SetWindowTextA(efp->exfunc->get_hwnd(efp->processing, 0, idx_track::img_y), img_y);

	// set label text next to the combo box and buttons.
	::SetWindowTextW(efp->exfunc->get_hwnd(efp->processing, 7, idx_check::algorithm), gui::algorithm_caption);
	::SetWindowTextW(efp->exfunc->get_hwnd(efp->processing, 5, idx_check::color), col_fmt);
	::SetWindowTextA(efp->exfunc->get_hwnd(efp->processing, 5, idx_check::file), file);
}

BOOL impl::func_WndProc(HWND, UINT message, WPARAM wparam, LPARAM lparam, AviUtl::EditHandle*, ExEdit::Filter* efp)
{
	// return TRUE if the image needs be re-rendered.

	if (message != ExEdit::ExtendedFilter::Message::WM_EXTENDEDFILTER_COMMAND) return FALSE;

	auto* exdata = reinterpret_cast<Exdata*>(efp->exdata_ptr);
	auto chk = static_cast<idx_check::id>(wparam >> 16);
	auto cmd = wparam & 0xffff;

	switch (cmd) {
		using namespace ExEdit::ExtendedFilter::CommandId;
	case EXTENDEDFILTER_SELECT_DROPDOWN:
		if (chk == idx_check::algorithm) {
			auto const alg = static_cast<Algorithm>(std::clamp(static_cast<int>(lparam), 0, algorithm_count - 1));
			if (exdata->algorithm != alg) {
				efp->exfunc->set_undo(efp->processing, 0);
				exdata->algorithm = alg;
				exedit.update_any_exdata(efp->processing, exdata_use[idx_data::algorithm].name);

				update_extendedfilter_wnd(efp);
			}
			return TRUE;
		}
		break;
	case EXTENDEDFILTER_PUSH_BUTTON:
		switch (chk) {
		case idx_check::color:
		{
			efp->exfunc->set_undo(efp->processing, 0);
			char const heading = std::exchange(exdata->file[0], '\0');
			if (efp->exfunc->x6c(efp, &exdata->color, 0x002)) { // color_dialog
				std::memset(exdata->file, 0, sizeof(exdata->file));
				exedit.update_any_exdata(efp->processing, exdata_use[idx_data::color].name);
				exedit.update_any_exdata(efp->processing, exdata_use[idx_data::file].name);

				update_extendedfilter_wnd(efp);
			}
			else exdata->file[0] = heading;
			return TRUE;
		}
		case idx_check::file:
		{
			decltype(exdata->file) file{};
			::strcpy_s(file, exdata->file);
			OPENFILENAMEA ofn{
				.lStructSize = sizeof(ofn),
				.hwndOwner = *exedit.hwnd_setting_dlg,
				.hInstance = nullptr,
				.lpstrFilter = efp->exfunc->get_loadable_image_extension(),
				.lpstrFile = file,
				.nMaxFile = std::size(file),
			};
			if (::GetOpenFileNameA(&ofn)) {
				efp->exfunc->set_undo(efp->processing, 0);
				std::memset(exdata->file, 0, sizeof(exdata->file));
				::strcpy_s(exdata->file, relative_path::relative{ file }.str_rel.c_str());
				exedit.update_any_exdata(efp->processing, exdata_use[idx_data::file].name);

				update_extendedfilter_wnd(efp);
			}
			return TRUE;
		}
		}
		break;
	case EXTENDEDFILTER_D_AND_D:
	{
		auto file = reinterpret_cast<char const*>(lparam);
		if (file == nullptr) break;

		efp->exfunc->set_undo(efp->processing, 0);
		std::memset(exdata->file, 0, sizeof(exdata->file));
		::strcpy_s(exdata->file, relative_path::relative{ file }.str_rel.c_str());
		exedit.update_any_exdata(efp->processing, exdata_use[idx_data::file].name);

		update_extendedfilter_wnd(efp);
		return TRUE;
	}
	}
	return FALSE;
}

int impl::func_window_init(HINSTANCE hinstance, HWND hwnd, int y, int base_id, int sw_param, ExEdit::Filter* efp)
{
	if (sw_param != 0) update_extendedfilter_wnd(efp);
	return 0;
}

