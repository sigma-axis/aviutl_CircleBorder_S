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

#include "Rounding.hpp"

using namespace Filter;
using namespace Rounding;


////////////////////////////////
// ウィンドウ状態の保守．
////////////////////////////////
static inline void update_extendedfilter_wnd(ExEdit::Filter* efp)
{
	using namespace Rounding::impl;
	auto* exdata = reinterpret_cast<Exdata*>(efp->exdata_ptr);

	// choose button text from "透明度", or "----".
	::SetWindowTextA(efp->exfunc->get_hwnd(efp->processing, 0, idx_track::transp),
		efp->check[idx_check::crop] != check_data::unchecked ?
		gui::invalid_name_narrow : track_names[idx_track::transp]);

	// choose button text from "αしきい値", "基準α和", or "----".
	::SetWindowTextW(efp->exfunc->get_hwnd(efp->processing, 0, idx_track::param_a), [&] {
		switch (exdata->algorithm) {
			using enum Filter::Algorithm;
		case bin: case bin2x: return gui::param_a_name;
		case sum: return gui::param_a_name_alt;
		default: return gui::invalid_name;
		}
	}());

	// select item in the combo box for algorithm.
	if (int idx = static_cast<int>(exdata->algorithm); 0 <= idx || idx < algorithm_count) {
		::SendMessageW(efp->exfunc->get_hwnd(efp->processing, 6, idx_check::algorithm), CB_SETCURSEL, idx, 0);
	}

	// set label text next to the combo box.
	::SetWindowTextW(efp->exfunc->get_hwnd(efp->processing, 7, idx_check::algorithm), gui::algorithm_caption);
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

	case EXTENDEDFILTER_UPDATE_CHECK:
		if (chk == idx_check::crop) {
			update_extendedfilter_wnd(efp);
			return TRUE;
		}
		break;
	}
	return FALSE;
}

int impl::func_window_init(HINSTANCE hinstance, HWND hwnd, int y, int base_id, int sw_param, ExEdit::Filter* efp)
{
	if (sw_param != 0) update_extendedfilter_wnd(efp);
	return 0;
}

