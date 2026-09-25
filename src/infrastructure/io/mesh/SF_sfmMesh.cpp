/// @file SF_sfmMesh.cpp
/// @brief SFM/生成网格文件的 IO 实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/
/*--------------Sonic Fluid-------------------*/

#include "SF_sfmMesh.h"
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace SF {

// ============================================================
//  内部工具函数
// ============================================================

namespace {

std::string upperSection(std::string section) {
    std::transform(section.begin(), section.end(), section.begin(),
                   [](unsigned char c) { return (char)std::toupper(c); });
    return section;
}

bool isReservedSection(const std::string& section) {
    const std::string upper = upperSection(section);
    return upper == "INFORMATION" ||
           upper == "POINT" ||
           upper == "CELL" ||
           upper == "CELLS" ||
           upper == "FACE" ||
           upper == "FACES" ||
           upper == "PATCH" ||
           upper == "PATCHES" ||
           upper == "INTERFACE";
}

} // namespace

std::string SFMMeshIO::resolvePath(const std::string& filePath) const {
    if (filePath.empty()) return filePath;
    // 绝对路径或已含..的路径直接返回
    if (filePath[0] == '/' || filePath[0] == '\\'
        || filePath.rfind("..", 0) == 0) {
        return filePath;
    }
    return caseDir_ + "/" + filePath;
}

void SFMMeshIO::cleanLine(std::string& line) {
    // 去除注释 (以 '#' 开头但非 section 标记的行, 和行内 '#' 注释)
    size_t comment = line.find("//");
    if (comment != std::string::npos) line = line.substr(0, comment);
    // 去除首尾空白
    while (!line.empty() && std::isspace(static_cast<unsigned char>(line.front())))
        line.erase(line.begin());
    while (!line.empty() && std::isspace(static_cast<unsigned char>(line.back())))
        line.pop_back();
}

// ============================================================
//  核心读取实现
// ============================================================

bool SFMMeshIO::readSFMInternal(const std::string& filePath,
                                bool readPoints, bool readSets,
                                SFMMeshBlock& outBlock) const {
    std::string absPath = resolvePath(filePath);
    std::ifstream in(absPath);
    if (!in.is_open()) {
        std::cerr << "[SFM] 无法打开文件: " << absPath << std::endl;
        return false;
    }

    outBlock.sourceFile = filePath;
    outBlock.hasPoints = false;
    outBlock.hasSets = false;

    std::string line, section;
    bool inPoints = false;
    int pointCount = 0;
    int expectedPoints = 0;

    while (std::getline(in, line)) {
        cleanLine(line);
        if (line.empty()) continue;

        // 检测section标记
        if (line[0] == '#') {
            section = line.substr(1);
            cleanLine(section);

            if (upperSection(section) == "POINT") {
                inPoints = readPoints;
                if (readPoints) {
                    expectedPoints = (outBlock.nx + 1) * (outBlock.ny + 1) * (outBlock.nz + 1);
                    outBlock.x.reserve(expectedPoints);
                    outBlock.y.reserve(expectedPoints);
                    outBlock.z.reserve(expectedPoints);
                    outBlock.hasPoints = true;
                }
                continue;
            }
            if (upperSection(section) == "INFORMATION") {
                inPoints = false;
                continue; // 下一行读取维度
            }

            if (isReservedSection(section)) {
                inPoints = false;
                continue;
            }

            // 其他#SetName标记
            if (readSets && !section.empty()) {
                // 这是一个集合段；空集合也要保留，表示该patch存在但本块无点。
                outBlock.pointSets[section];
                outBlock.hasSets = true;
                inPoints = false;
            }
            continue;
        }

        // 读取Information行 (维度)
        if (upperSection(section) == "INFORMATION" && outBlock.nx == 0) {
            std::istringstream dimIn(line);
            if (dimIn >> outBlock.nx >> outBlock.ny >> outBlock.nz >> outBlock.ng) {
                continue;
            }
        }

        // 读取坐标点
        if (inPoints && readPoints) {
            std::istringstream ptIn(line);
            double px, py, pz;
            if (ptIn >> px >> py >> pz) {
                outBlock.x.push_back(px);
                outBlock.y.push_back(py);
                outBlock.z.push_back(pz);
                ++pointCount;
            }
        }

        if (readSets && !inPoints && !section.empty()
            && !isReservedSection(section)) {
            std::istringstream ss(line);
            int idx;
            while (ss >> idx) outBlock.pointSets[section].push_back(idx);
        }
    }

    in.close();

    // 校验
    if (readPoints && outBlock.hasPoints && !outBlock.isValid()) {
        std::cerr << "[SFM] 坐标点数不匹配: 期望" << expectedPoints
                  << ", 实际" << pointCount
                  << " (" << absPath << ")" << std::endl;
        return false;
    }

    return readPoints ? outBlock.hasPoints : outBlock.hasSets || outBlock.nx > 0;
}

// ============================================================
//  公开读接口
// ============================================================

bool SFMMeshIO::readSingleBlock(const std::string& filePath,
                                SFMMeshBlock& outBlock) const {
    return readSFMInternal(filePath, true, true, outBlock);
}

bool SFMMeshIO::readPointsOnly(const std::string& filePath,
                               SFMMeshBlock& outBlock) const {
    return readSFMInternal(filePath, true, false, outBlock);
}

bool SFMMeshIO::readSetsOnly(const std::string& filePath,
                             std::map<std::string, std::vector<int>>& outSets) const {
    SFMMeshBlock temp;
    if (!readSFMInternal(filePath, false, true, temp)) return false;
    outSets = std::move(temp.pointSets);
    return !outSets.empty();
}

bool SFMMeshIO::hasSetsSection(const std::string& filePath) const {
    std::string absPath = resolvePath(filePath);
    std::ifstream in(absPath);
    if (!in.is_open()) return false;

    std::string line;
    while (std::getline(in, line)) {
        cleanLine(line);
        if (!line.empty() && line[0] == '#') {
            std::string sec = line.substr(1);
            cleanLine(sec);
            if (!isReservedSection(sec) && !sec.empty()) {
                return true;
            }
        }
    }
    return false;
}

bool SFMMeshIO::readDimensions(const std::string& filePath,
                               int& nx, int& ny, int& nz, int& ng) const {
    std::string absPath = resolvePath(filePath);
    std::ifstream in(absPath);
    if (!in.is_open()) return false;

    std::string line;
    while (std::getline(in, line)) {
        cleanLine(line);
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        std::istringstream ss(line);
        if (ss >> nx >> ny >> nz >> ng) return true;
        break;
    }
    return false;
}

// ============================================================
//  公开写接口
// ============================================================

bool SFMMeshIO::writeSingleBlock(const std::string& filePath,
                                 const SFMMeshBlock& block,
                                 bool writeSets,
                                 bool useGlobalIndex) const {
    std::string absPath = resolvePath(filePath);
    std::ofstream out(absPath);
    if (!out.is_open()) {
        std::cerr << "[SFM] 无法创建文件: " << absPath << std::endl;
        return false;
    }

    // #Information
    out << "#Information\n"
        << block.nx << " " << block.ny << " " << block.nz << " " << block.ng
        << "\n\n";

    // #Point
    out << "#Point\n";
    int np = static_cast<int>(block.x.size());
    out << std::scientific << std::setprecision(12);
    for (int i = 0; i < np; ++i) {
        out << block.x[i] << " " << block.y[i] << " " << block.z[i] << "\n";
    }
    out << "\n";

    // #SetName (可选)
    if (writeSets && !block.pointSets.empty()) {
        for (const auto& [name, indices] : block.pointSets) {
            out << "#" << name << "\n";
            // 每行最多输出20个索引
            int count = 0;
            for (int idx : indices) {
                out << idx;
                ++count;
                if (count % 20 == 0 || count == static_cast<int>(indices.size()))
                    out << "\n";
                else
                    out << " ";
            }
        }
        out << "\n";
    }

    out.close();
    return true;
}

bool SFMMeshIO::writeSetsFile(const std::string& filePath,
                              const std::map<std::string, std::vector<int>>& sets,
                              int nx, int ny, int nz, int ng) const {
    std::string absPath = resolvePath(filePath);
    std::ofstream out(absPath);
    if (!out.is_open()) {
        std::cerr << "[SFM] 无法创建sets文件: " << absPath << std::endl;
        return false;
    }

    out << "#Information\n"
        << nx << " " << ny << " " << nz << " " << ng << "\n\n";

    for (const auto& [name, indices] : sets) {
        out << "#" << name << "\n";
        int count = 0;
        for (int idx : indices) {
            out << idx;
            ++count;
            if (count % 20 == 0 || count == static_cast<int>(indices.size()))
                out << "\n";
            else
                out << " ";
        }
    }
    out << "\n";
    out.close();
    return true;
}

bool SFMMeshIO::writeMultiBlockLocal(
    const std::string& basePath,
    const std::vector<SFMMeshBlock>& blocks,
    const std::vector<std::map<std::string, std::vector<int>>>& localSets) const
{
    if (blocks.size() != localSets.size()) {
        std::cerr << "[SFM] 块数(" << blocks.size()
                  << ")与sets数(" << localSets.size() << ")不匹配" << std::endl;
        return false;
    }

    for (size_t b = 0; b < blocks.size(); ++b) {
        SFMMeshBlock block = blocks[b];
        block.pointSets = localSets[b];

        char suffix[32];
        std::snprintf(suffix, sizeof(suffix), "_block%03zu.sfm", b);
        std::string fileName = basePath + suffix;

        if (!writeSingleBlock(fileName, block, true, false)) return false;
        std::cout << "[SFM] 写出: " << fileName << std::endl;
    }

    return true;
}

bool SFMMeshIO::writeMultiBlockGlobal(
    const std::string& basePath,
    const std::vector<SFMMeshBlock>& blocks,
    const std::map<std::string, std::vector<int>>& globalSets,
    int globalNx, int globalNy, int globalNz, int ng) const
{
    // 先写出各块(不含sets段)
    for (size_t b = 0; b < blocks.size(); ++b) {
        char suffix[32];
        std::snprintf(suffix, sizeof(suffix), "_block%03zu.sfm", b);
        std::string fileName = basePath + suffix;

        if (!writeSingleBlock(fileName, blocks[b], false, false)) return false;
        std::cout << "[SFM] 写出: " << fileName << std::endl;
    }

    // 写出共享sets文件
    std::string setsFilePath = basePath + "_sets.sfm";
    if (!writeSetsFile(setsFilePath, globalSets, globalNx, globalNy, globalNz, ng))
        return false;
    std::cout << "[SFM] 写出共享sets: " << setsFilePath << std::endl;

    return true;
}

} // namespace SF
