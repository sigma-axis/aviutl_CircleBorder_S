/*
The MIT License (MIT)

Copyright (c) 2024 sigma-axis

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <cstdint>
#include <string>
#include <vector>

#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>

#include "relative_path.hpp"
#include "CircleBorder_S.hpp"

using namespace relative_path;


////////////////////////////////
// ファイルパスの相対化．
////////////////////////////////
static std::string path_aviutl_dir{};
static std::vector<char> path_aviutl_dir_lower{};
static inline void init() {
	if (!path_aviutl_dir.empty()) return;
	path_aviutl_dir.resize_and_overwrite(MAX_PATH - 1,
		[](auto p, auto c) { return ::GetModuleFileNameA(nullptr, p, c + 1); });
	path_aviutl_dir.erase(pos_file_name(path_aviutl_dir));
	path_aviutl_dir.shrink_to_fit();

	path_aviutl_dir_lower = { path_aviutl_dir.begin(), path_aviutl_dir.end() };
	::CharLowerBuffA(path_aviutl_dir_lower.data(), path_aviutl_dir_lower.size());
}

// 絶対 → 相対．
relative::relative(std::string_view const& src) : abs_path{ src }, str_rel{}
{
	init();
	size_t src_pos = 0;
	std::vector<char> src_buf{ src.begin(), src.end() };
	::CharLowerBuffA(src_buf.data(), src_buf.size());

	// compare with the project file path.
	if (size_t len = *exedit.editp == nullptr ? 0 : pos_file_name((*exedit.editp)->project_filename);
		0 < len && len <= src.size()) {
		std::vector<char> cmp_buf{ (*exedit.editp)->project_filename, (*exedit.editp)->project_filename + len };
		::CharLowerBuffA(cmp_buf.data(), cmp_buf.size());
		if (std::memcmp(src_buf.data(), cmp_buf.data(), sizeof(char) * cmp_buf.size()) == 0) {
			src_pos = len;
			str_rel = header::project;
			goto copy_rest;
		}
	}

	// compare with the path to the parent directory of `aviutl.exe`.
	if (std::memcmp(src_buf.data(), path_aviutl_dir_lower.data(), sizeof(char) * path_aviutl_dir_lower.size()) == 0) {
		src_pos = path_aviutl_dir.size();
		str_rel = header::aviutl;
		goto copy_rest;
	}

copy_rest:
	// simply copy otherwise.
	str_rel += src.substr(src_pos);
}

// 相対 → 絶対．
absolute::absolute(std::string_view const& src) : str_rel{ src }, abs_path{}
{
	init();
	size_t src_pos = 0;

	// see if it starts with "<exe>", which represents the parent directory of `aviutl.exe`.
	if (std::memcmp(src.data(), header::aviutl.data(), sizeof(char) * header::aviutl.size()) == 0) {
		src_pos = header::aviutl.size();
		abs_path = path_aviutl_dir;
	}

	// then see "<aup>", which represents the parent directory of the current project file.
	else if (std::memcmp(src.data(), header::project.data(), sizeof(char) * header::project.size()) == 0) {
		src_pos = header::project.size();
		if (size_t len; *exedit.editp != nullptr &&
			(len = pos_file_name((*exedit.editp)->project_filename)) > 0)
			abs_path = { (*exedit.editp)->project_filename, len };
		else
			// as a fallback, use the path to the parent directory of `aviutl.exe`.
			abs_path = path_aviutl_dir;
	}

	// copy the rest.
	abs_path += src.substr(src_pos);
}

size_t relative_path::pos_file_name(std::string_view const& file)
{
	constexpr std::string_view delimiters = ">\\/";
	auto name = &*file.end(), st = file.data();
	while (name != st) {
		auto p = ::CharPrevA(st, name);
		if (!delimiters.contains(*p)) name = p;
		else break;
	}
	return name - st;
}

