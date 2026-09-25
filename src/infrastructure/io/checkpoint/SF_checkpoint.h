#pragma once
#include "infrastructure/io/field/SF_field.h"
namespace SF::CheckpointIO {
/// @brief checkpoint 仅消费 restart 字段，不包含任何物理量白名单。
inline void write(const Model::FieldRegistry& registry,const std::string& path) { FieldIO::write(registry,path,true); }
inline void read(Model::FieldRegistry& registry,const std::string& path) { FieldIO::read(registry,path,true); }
}
