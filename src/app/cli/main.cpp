/// @file main.cpp
/// @brief 命令行入口、参数解析与应用启动职责。

#include "SF_cli.h"

int main(int argc, char* argv[])
{
    return SF::CLI::run(argc, argv);
}
