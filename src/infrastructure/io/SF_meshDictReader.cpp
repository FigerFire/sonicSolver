/// @file SF_meshDictReader.cpp
/// @brief case 配置、场或 VTK 结果的基础设施 IO 实现。

/*--------------Sonic Fluid-------------------*/
/*----------IO private implementation---------*/

#include "mesh/SF_meshParameters.h"
#include "core/interfaces/SF_log.h"
#include "private/SF_casePath.h"
#include "private/SF_foamParser.h"

#include "SF_meshGen.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace SF {

using namespace IOPrivate;

bool MeshIO::readParameters(const std::string& path, MeshParameters& out) {
    std::ifstream file(path);
    if (!file.is_open()) {
        broadcast("Error: ", "Cannot open " + path);
        return false;
    }
    const std::string text((std::istreambuf_iterator<char>(file)),
                           std::istreambuf_iterator<char>());
    const std::string number = R"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)";
    const std::regex foamFile(R"(\bFoamFile\b)");
    const std::regex blockMeshObject(R"(\bobject\s+blockMeshDict\s*;)");
    if (!std::regex_search(text, foamFile) || !std::regex_search(text, blockMeshObject)) {
        broadcast("Error: ", "createMesh accepts only an OpenFOAM blockMeshDict: " + path);
        return false;
    }

    auto listBody = [](const std::string& source, const std::string& name,
                       std::string* body) {
        const std::regex key("\\b" + name + R"(\s*\()");
        std::smatch match;
        if (!std::regex_search(source, match, key)) return false;
        const size_t open = (size_t)match.position() + match.length() - 1;
        int depth = 0;
        for (size_t i = open; i < source.size(); ++i) {
            if (source[i] == '(') ++depth;
            else if (source[i] == ')' && --depth == 0) {
                *body = source.substr(open + 1, i - open - 1);
                return true;
            }
        }
        return false;
    };
    auto parseDoubles = [](const std::string& source) {
        std::vector<double> values;
        const std::regex value(R"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)");
        for (std::sregex_iterator it(source.begin(), source.end(), value), end;
             it != end; ++it) values.push_back(std::stod(it->str()));
        return values;
    };
    auto parseIntegers = [](const std::string& source) {
        std::vector<int> values;
        const std::regex value(R"([-+]?\d+)");
        for (std::sregex_iterator it(source.begin(), source.end(), value), end;
             it != end; ++it) values.push_back(std::stoi(it->str()));
        return values;
    };

    out = MeshParameters{};
    std::smatch scale;
    if (std::regex_search(text, scale,
                          std::regex("\\bconvertToMeters\\s+(" + number + R"()\s*;)"))) {
        out.scale = std::stod(scale[1].str());
    }
    if (!std::isfinite(out.scale) || out.scale <= 0.0) {
        broadcast("Error: ", "blockMeshDict convertToMeters must be finite and positive.");
        return false;
    }

    // ghost 容量属于网格输入，不绑定某个离散格式。Runtime 仍由当前
    // numerics contract 请求实际 halo depth；网格只负责提供足够容量。
    std::smatch ghostLayers;
    if (std::regex_search(text, ghostLayers,
                          std::regex(R"(\bnGhost\s+([-+]?\d+)\s*;)"))) {
        out.nGhost = std::stoi(ghostLayers[1].str());
    }
    if (out.nGhost <= 0) {
        broadcast("Error: ", "blockMeshDict nGhost must be a positive integer.");
        return false;
    }

    std::string verticesBody;
    if (!listBody(text, "vertices", &verticesBody)) {
        broadcast("Error: ", "blockMeshDict is missing vertices.");
        return false;
    }
    const std::regex vertex("\\(\\s*(" + number + ")\\s+(" + number
                            + ")\\s+(" + number + R"()\s*\))");
    for (std::sregex_iterator it(verticesBody.begin(), verticesBody.end(), vertex), end;
         it != end; ++it) {
        out.vertices.emplace_back(std::stod((*it)[1].str()), std::stod((*it)[2].str()),
                                  std::stod((*it)[3].str()));
    }

    std::string blocksBody;
    if (!listBody(text, "blocks", &blocksBody)) {
        broadcast("Error: ", "blockMeshDict is missing blocks.");
        return false;
    }
    const std::regex unsupportedBlock("\\b(wedge|prism)\\s*\\(",
                                      std::regex::icase);
    if (std::regex_search(blocksBody, unsupportedBlock)) {
        broadcast("Error: ", "createMesh currently supports only hex blocks; wedge/prism are not implemented.");
        return false;
    }
    const std::regex block("\\bhex\\s*\\(([^)]*)\\)\\s*\\(([^)]*)\\)\\s*"
                           "simpleGrading\\s*\\(([^)]*)\\)", std::regex::icase);
    for (std::sregex_iterator it(blocksBody.begin(), blocksBody.end(), block), end;
         it != end; ++it) {
        const std::vector<int> vertices = parseIntegers((*it)[1].str());
        const std::vector<int> cells = parseIntegers((*it)[2].str());
        const std::vector<double> grading = parseDoubles((*it)[3].str());
        if (vertices.size() != 8 || cells.size() != 3 || grading.size() != 3) {
            broadcast("Error: ", "Invalid hex block definition in blockMeshDict.");
            return false;
        }
        MeshBlock meshBlock{};
        for (int i = 0; i < 8; ++i) meshBlock.verts[i] = vertices[(size_t)i];
        for (int i = 0; i < 3; ++i) {
            meshBlock.cells[i] = cells[(size_t)i];
            meshBlock.grade[i] = grading[(size_t)i];
        }
        out.blocks.push_back(meshBlock);
    }

    std::string edgesBody;
    if (!listBody(text, "edges", &edgesBody)) {
        broadcast("Error: ", "blockMeshDict is missing edges.");
        return false;
    }
    const std::regex arc("\\barc\\s+([-+]?\\d+)\\s+([-+]?\\d+)\\s*\\(\\s*("
                         + number + ")\\s+(" + number + ")\\s+(" + number
                         + R"()\s*\))", std::regex::icase);
    for (std::sregex_iterator it(edgesBody.begin(), edgesBody.end(), arc), end;
         it != end; ++it) {
        MeshEdge edge;
        edge.type = EdgeType::ARC;
        edge.v0 = std::stoi((*it)[1].str());
        edge.v1 = std::stoi((*it)[2].str());
        edge.arcP = Vector3(std::stod((*it)[3].str()), std::stod((*it)[4].str()),
                            std::stod((*it)[5].str()));
        out.edges.push_back(edge);
    }

    std::string boundaryBody;
    if (!listBody(text, "boundary", &boundaryBody)) {
        broadcast("Error: ", "blockMeshDict is missing boundary.");
        return false;
    }
    const std::regex patch(R"(([A-Za-z_][A-Za-z0-9_.-]*)\s*\{([^{}]*)\})");
    const std::regex face(R"(\(\s*([-+]?\d+)\s+([-+]?\d+)\s+([-+]?\d+)\s+([-+]?\d+)\s*\))");
    for (std::sregex_iterator it(boundaryBody.begin(), boundaryBody.end(), patch), end;
         it != end; ++it) {
        std::string facesBody;
        if (!listBody((*it)[2].str(), "faces", &facesBody)) continue;
        MeshFacePatch meshPatch;
        meshPatch.name = (*it)[1].str();
        for (std::sregex_iterator fit(facesBody.begin(), facesBody.end(), face), fend;
             fit != fend; ++fit) {
            meshPatch.faces.push_back({std::stoi((*fit)[1].str()), std::stoi((*fit)[2].str()),
                                       std::stoi((*fit)[3].str()), std::stoi((*fit)[4].str())});
        }
        if (!meshPatch.faces.empty()) out.facePatches.push_back(std::move(meshPatch));
    }

    if (out.vertices.empty() || out.blocks.empty()) {
        broadcast("Error: ", "blockMeshDict must contain at least one vertex and one hex block.");
        return false;
    }
    const auto validVertex = [&out](int index) {
        return index >= 0 && index < (int)out.vertices.size();
    };
    for (const MeshBlock& meshBlock : out.blocks) {
        for (int vertexIndex : meshBlock.verts) {
            if (!validVertex(vertexIndex)) {
                broadcast("Error: ", "blockMeshDict block references an invalid vertex index.");
                return false;
            }
        }
        for (int direction = 0; direction < 3; ++direction) {
            if (meshBlock.cells[direction] <= 0 || !std::isfinite(meshBlock.grade[direction])
                || meshBlock.grade[direction] <= 0.0) {
                broadcast("Error: ", "blockMeshDict cells and simpleGrading must be positive.");
                return false;
            }
        }
    }
    for (const MeshEdge& edge : out.edges) {
        if (!validVertex(edge.v0) || !validVertex(edge.v1)) {
            broadcast("Error: ", "blockMeshDict arc references an invalid vertex index.");
            return false;
        }
    }
    for (const MeshFacePatch& meshPatch : out.facePatches) {
        for (const auto& faceVertices : meshPatch.faces) {
            for (int vertexIndex : faceVertices) {
                if (!validVertex(vertexIndex)) {
                    broadcast("Error: ", "blockMeshDict boundary face references an invalid vertex index.");
                    return false;
                }
            }
        }
    }
    broadcast("Parsed blockMeshDict: ", std::to_string(out.vertices.size()) + " verts, "
             + std::to_string(out.blocks.size()) + " hex blocks, scale="
             + std::to_string(out.scale) + ", nGhost=" + std::to_string(out.nGhost));
    return true;
}

} // namespace SF
