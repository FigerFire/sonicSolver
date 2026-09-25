/// @file SF_casePath.cpp
/// @brief case 字典解析器的内部无求解逻辑辅助实现。

#include "private/SF_casePath.h"

#include "private/SF_foamParser.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sys/stat.h>

namespace SF::IOPrivate {

std::string stripTrailingSlash(std::string path) {
    while (path.size() > 1 && (path.back() == '/' || path.back() == '\\')) {
        path.pop_back();
    }
    return path;
}

bool isDirectoryPath(const std::string& path) {
    struct stat st {};
    return !path.empty() && stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}

std::string parentPath(const std::string& path) {
    std::string p = stripTrailingSlash(path);
    const size_t pos = p.find_last_of("/\\");
    if (pos == std::string::npos) return ".";
    if (pos == 0) return p.substr(0, 1);
    return p.substr(0, pos);
}

std::string basePathName(const std::string& path) {
    std::string p = stripTrailingSlash(path);
    const size_t pos = p.find_last_of("/\\");
    return (pos == std::string::npos) ? p : p.substr(pos + 1);
}

std::string basePathNameNoExtension(const std::string& path) {
    std::string name = basePathName(path);
    const size_t dot = name.find_last_of('.');
    if (dot != std::string::npos) name = name.substr(0, dot);
    return name;
}

std::string upperToken(std::string s) {
    s.erase(std::remove(s.begin(), s.end(), '"'), s.end());
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return (char)std::toupper(c); });
    return s;
}

std::string trimToken(std::string s) {
    auto notSpace = [](unsigned char c) { return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), notSpace));
    s.erase(std::find_if(s.rbegin(), s.rend(), notSpace).base(), s.end());
    return s;
}

std::string stripTokenQuotes(std::string s) {
    s = trimToken(s);
    while (!s.empty() && (s.back() == ',' || s.back() == ';')) s.pop_back();
    if (s.size() >= 2 &&
        ((s.front() == '"' && s.back() == '"') ||
         (s.front() == '\'' && s.back() == '\''))) {
        s = s.substr(1, s.size() - 2);
    }
    return trimToken(s);
}

bool fileExists(const std::string& path) {
    std::ifstream in(path);
    return (bool)in;
}

std::string resolveCaseFile(const std::string& caseDir, const std::string& fileName) {
    if (fileName.empty()) return fileName;
    if (fileName[0] == '/' || fileName[0] == '\\') {
        return fileName;
    }
    return caseDir + "/" + fileName;
}

bool isPlainRelativeFileName(const std::string& path) {
    return !path.empty()
        && path[0] != '/'
        && path[0] != '\\'
        && path.rfind("..", 0) != 0
        && path.find('/') == std::string::npos
        && path.find('\\') == std::string::npos;
}

[[noreturn]] void fatalILWConfig(const std::string& message) {
    std::cerr << "[SF FATAL] " << message << std::endl;
    std::exit(1);
}

[[noreturn]] void fatalIBMConfig(const std::string& message) {
    std::cerr << "[SF FATAL] " << message << std::endl;
    std::exit(1);
}

[[noreturn]] void fatalCaseConfig(const std::string& message) {
    std::cerr << "[SF FATAL] " << message << std::endl;
    std::exit(1);
}

} // namespace SF::IOPrivate
