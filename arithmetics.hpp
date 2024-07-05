/*
The MIT License (MIT)

Copyright (c) 2024 sigma-axis

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <cstdint>
#include <cmath>


////////////////////////////////
// 単純計算補助．
////////////////////////////////
namespace Calculation::arith
{
	constexpr auto square(auto x) { return x * x; }
	constexpr auto square_sum(auto... x) { return (square(x) + ...); }
	constexpr auto floor_div(auto dividend, auto divisor) {
		auto d = std::div(dividend, divisor);
		if (d.rem < 0) d.quot -= divisor >= 0 ? +1 : -1;
		return d.quot;
	}
	constexpr auto abs(auto x) { return x < 0 ? -x : x; }
	namespace arc
	{
		inline size_t half(size_t size_sq, int32_t* arc)
		{
			int sz = static_cast<int>(std::sqrt(size_sq));
			arc[sz] = sz;
			for (int i = 1, K = size_sq - 1, d = 3; K >= 0; i++, K -= d, d += 2)
				arc[sz + i] = arc[sz - i] = static_cast<int32_t>(std::sqrt(K));

			return sz;
		}
		inline size_t quarter(size_t size_sq, int32_t* arc)
		{
			for (int i = 0, K = size_sq, d = 1; K >= 0; i++, K -= d, d += 2)
				arc[i] = static_cast<int32_t>(std::sqrt(K));
			return arc[0];
		}
	}
}
