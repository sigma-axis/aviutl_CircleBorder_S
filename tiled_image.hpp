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
using byte = uint8_t;
#include <exedit/Filter.hpp>
#include <exedit/Exfunc.hpp>

#include "relative_path.hpp"


////////////////////////////////
// パターン画像の取得．
////////////////////////////////
struct tiled_image {
	int w = 0, h = 0, ox = 0, oy = 0;
	ExEdit::PixelYCA* buff = nullptr;

	operator bool() const { return buff != nullptr; }
	tiled_image(char const* path, int img_x, int img_y, int displace, ExEdit::Filter* efp, void* buffer)
	{
		if (path == nullptr || path[0] == '\0') return;

		buff = reinterpret_cast<ExEdit::PixelYCA*>(buffer);
		if (efp->exfunc->load_image(buff, relative_path::absolute{ path }.abs_path.data(), &w, &h, 0, 0) == 0) {
			buff = nullptr;
			return;
		}

		ox = (-img_x - displace) % w;
		oy = (-img_y - displace) % h;
		if (ox < 0) ox += w; if (oy < 0) oy += h;
	}
	auto& operator[](int idx) const { return buff[idx]; }
};

