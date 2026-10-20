#pragma once

/// @file SF_casePath.h
/// @brief Case path resolution and fatal configuration diagnostics.

#include <string>

namespace SF::IOPrivate {

std::string stripTrailingSlash(std::string path);
bool isDirectoryPath(const std::string& path);
std::string parentPath(const std::string& path);
std::string basePathName(const std::string& path);
std::string basePathNameNoExtension(const std::string& path);
std::string upperToken(std::string value);
std::string stripTokenQuotes(std::string value);
bool fileExists(const std::string& path);
std::string resolveCaseFile(const std::string& caseDir,
                            const std::string& fileName);
bool isPlainRelativeFileName(const std::string& path);
[[noreturn]] void fatalILWConfig(const std::string& message);
[[noreturn]] void fatalIBMConfig(const std::string& message);
[[noreturn]] void fatalCaseConfig(const std::string& message);

} // namespace SF::IOPrivate

