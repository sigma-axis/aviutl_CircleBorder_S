/*
The MIT License (MIT)

Copyright (c) 2024 sigma-axis

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <cstdint>

#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>

#include "multi_thread.hpp"
#include "Border.hpp"
#include "Rounding.hpp"
#include "Outline.hpp"


////////////////////////////////
// 変数アドレス初期化．
////////////////////////////////
void ExEdit092::init_pointers(AviUtl::FilterPlugin* efp)
{
	fp = efp;

	auto pick_addr = [exedit_base = reinterpret_cast<int32_t>(efp->dll_hinst), this]
		<class T>(T& target, ptrdiff_t offset) { target = reinterpret_cast<T>(exedit_base + offset); };
	auto pick_val  = [exedit_base = reinterpret_cast<int32_t>(efp->dll_hinst), this]
		<class T>(T& target, ptrdiff_t offset) { target = *reinterpret_cast<T*>(exedit_base + offset); };

	pick_addr(hwnd_setting_dlg,		0x1539c8);
	pick_addr(editp,				0x1a532c);
	pick_addr(memory_ptr,			0x1a5328);

	pick_val (yca_max_w,			0x196748);
	pick_val (yca_max_h,			0x1920e0);

	pick_addr(update_any_exdata,	0x04a7e0);
	pick_addr(nextundo,				0x08d150);

	// make ready the multi-thread function feature.
	constexpr uint32_t ofs_num_threads_address = 0x086384;
	auto aviutl_base = reinterpret_cast<uint32_t>(fp->hinst_parent);
	multi_thread.init(fp->exfunc->exec_multi_thread_func,
		reinterpret_cast<int32_t*>(aviutl_base + ofs_num_threads_address));
}


////////////////////////////////
// DLL 初期化．
////////////////////////////////
BOOL WINAPI DllMain(HINSTANCE hinst, DWORD fdwReason, LPVOID lpvReserved)
{
	switch (fdwReason) {
	case DLL_PROCESS_ATTACH:
		::DisableThreadLibraryCalls(hinst);
		break;
	}
	return TRUE;
}


////////////////////////////////
// エントリポイント．
////////////////////////////////
EXTERN_C __declspec(dllexport) ExEdit::Filter* const* __stdcall GetFilterTableList() {
	constexpr static ExEdit::Filter* filter_list[] = {
		&Filter::Border::filter,
		&Filter::Rounding::filter,
		&Filter::Outline::filter,
		nullptr,
	};

	return filter_list;
}

