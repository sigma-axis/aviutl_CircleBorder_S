/*
The MIT License (MIT)

Copyright (c) 2024 sigma-axis

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <cstdint>
#include <limits>
#include <concepts>
#include <tuple>


////////////////////////////////
// バッファ操作の基礎ヘッダ．
////////////////////////////////
namespace Calculation
{
	using i16 = int16_t;
	using i32 = int32_t;

	constexpr size_t log2_max_alpha = 12;
	constexpr i16 max_alpha = 1 << log2_max_alpha;

	struct Bounds {
		int L, T; // inclusive.
		int R, B; // exclusive.

		constexpr int wd() const { return R - L; }
		constexpr int ht() const { return B - T; }
		constexpr bool is_empty() const { return wd() <= 0 || ht() <= 0; }

		[[nodiscard]] constexpr Bounds move(int x, int y) const { return { L + x, T + y, R + x, B + y }; }
		[[nodiscard]] constexpr Bounds inflate(int x, int y) const { return { L - x, T - y, R + x, B + y }; }
		[[nodiscard]] constexpr Bounds inflate(int len) const { return inflate(len, len); }
		[[nodiscard]] constexpr Bounds inflate_tl(int x, int y) const { return { L - x, T - y, R, B }; }
		[[nodiscard]] constexpr Bounds inflate_tl(int len) const { return inflate_tl(len, len); }
		[[nodiscard]] constexpr Bounds inflate_br(int x, int y) const { return { L, T, R + x, B + y }; }
		[[nodiscard]] constexpr Bounds inflate_br(int len) const { return inflate_br(len, len); }

		// for the use of std::tie().
		constexpr operator std::tuple<int&, int&, int&, int&>() { return { L, T, R, B }; }
	};

	// each element in `range` is a left-closed, right-open interval.
	// returned interval is left-closed, right-open.
	template<std::totally_ordered T>
	constexpr auto unite_interval(auto&& range) {
		T lbd = std::numeric_limits<T>::max(), ubd = std::numeric_limits<T>::min();
		for (auto& [l, u] : range) {
			if (l < u) {
				lbd = std::min(lbd, l);
				ubd = std::max(ubd, u);
			}
		}
		return std::pair{ lbd, ubd };
	}
	// each element in `range` is a both-end-closed interval.
	// returned interval is left-closed, right-open.
	template<std::totally_ordered T>
	constexpr auto unite_interval_alt(auto&& range) {
		T lbd = std::numeric_limits<T>::max(), ubd = std::numeric_limits<T>::min();
		for (auto& [l, u] : range) {
			if (l <= u) {
				lbd = std::min(lbd, l);
				ubd = std::max(ubd, u);
			}
		}
		return std::pair{ lbd, ubd + 1 };
	}
}

