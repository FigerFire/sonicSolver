/// @file SF_documents.cpp
/// @brief native case 的磁盘查询委托。native 语义不再经过任何中间文档表示。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_compatibility.h"
#include "private/SF_casePath.h"

namespace SF {

// native case 的语义全部来自 Model::Description 的 typed sections；
// 这里只保留真正的磁盘查询，用于 mesh / geometry / result 路径校验。
bool CaseAdapter::fileExists(const std::string& path) const {
    return IOPrivate::fileExists(path);
}

bool CaseAdapter::isDirectoryPath(const std::string& path) const {
    return IOPrivate::isDirectoryPath(path);
}

} // namespace SF
