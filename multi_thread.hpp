/*
The MIT License (MIT)

Copyright (c) 2024 sigma-axis

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <cstdint>
#include <tuple>
#include <vector>
#include <thread>
#include <mutex>


////////////////////////////////
// AviUtl のマルチスレッド関数のラッパー．
////////////////////////////////
inline constinit struct MultiThread {
	auto operator()(int num_parallel, auto&&... args, auto&& func) const {
		return (*this)(num_parallel < num_threads(), args..., func);
	}
	auto operator()(bool single_thread, auto&&... args, auto&& func) const
	{
		using RetT = std::invoke_result_t<decltype(func), int, int, decltype(args)...>;
		if (single_thread) {
			if constexpr (std::is_void_v<RetT>)
				return func(0, 1, args...);
			else return std::vector<RetT>{ func(0, 1, args...) };
		}

		auto cxt = std::tuple{ &func, &args... };
		constexpr auto invoke = [](auto& cxt, auto... params) {
			return [&]<size_t... I>(std::index_sequence<I...>) {
				return (*std::get<0>(cxt))(params..., *std::get<1 + I>(cxt)...);
			}(std::make_index_sequence<sizeof...(args)>{});
		};

		if constexpr (std::is_void_v<RetT>) {
			exec_multi_thread_func([](int thread_id, int thread_num, void* param1, void*) {
				invoke(*reinterpret_cast<decltype(cxt)*>(param1), thread_id, thread_num);
			}, &cxt, nullptr);
		}
		else {
			std::vector<RetT> ret(num_threads());

			exec_multi_thread_func([](int thread_id, int thread_num, void* param1, void* param2) {
				// assign the return value to a std::vector<>.
				(*reinterpret_cast<decltype(ret)*>(param2))[thread_id]
					= invoke(*reinterpret_cast<decltype(cxt)*>(param1), thread_id, thread_num);
			}, &cxt, &ret);

			return ret;
		}
	}

	int32_t num_threads() const {
		return *ptr_num_threads != 0 ? *ptr_num_threads : def_num_threads;
	}

private:
	//decltype(AviUtl::ExFunc::exec_multi_thread_func) exec_multi_thread_func = nullptr;
	int32_t (*exec_multi_thread_func)(void(*func)(int thread_id, int thread_num, void* param1, void* param2), void* param1, void* param2) = nullptr;
	int32_t* ptr_num_threads = nullptr; // 0x086384
	int32_t def_num_threads = 0;

	friend struct ExEdit092;
	void init(decltype(exec_multi_thread_func) mt_func, int32_t* num_threads) {
		if (def_num_threads > 0) return;

		exec_multi_thread_func = mt_func;
		ptr_num_threads = num_threads;
		def_num_threads = std::thread::hardware_concurrency();
	}
} multi_thread{};

