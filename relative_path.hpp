/*
The MIT License (MIT)

Copyright (c) 2024 sigma-axis

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <cstdint>
#include <string>


////////////////////////////////
// ファイルパスの相対化．
////////////////////////////////
namespace relative_path
{
	namespace header
	{
		using namespace std::string_view_literals;
		constexpr auto aviutl = "<exe>"sv, project = "<aup>"sv;
	}

	struct relative {
		std::string_view abs_path;
		std::string str_rel;
		explicit relative(std::string_view const& src);
	};
	struct absolute {
		std::string_view str_rel;
		std::string abs_path;
		explicit absolute(std::string_view const& src);
	};

	size_t pos_file_name(std::string_view const& file);
	inline auto* ptr_file_name(auto* file) { return file + pos_file_name(file); }
}

