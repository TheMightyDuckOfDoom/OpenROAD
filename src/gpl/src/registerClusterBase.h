// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#pragma once

#include <memory>

namespace cts {
class TritonCTS;
}

namespace utl {
class Logger;
}

namespace gpl {

class NesterovBaseCommon;
class GNet;

class RegisterClusterBase
{
 public:
  RegisterClusterBase();
  RegisterClusterBase(std::shared_ptr<NesterovBaseCommon> nbc,
             cts::TritonCTS* cts,
             utl::Logger* log);

  void executeRegisterClustering();

 private:
  utl::Logger* log_ = nullptr;
  std::shared_ptr<NesterovBaseCommon> nbc_;
  cts::TritonCTS* cts_ = nullptr;
};

}  // namespace gpl
