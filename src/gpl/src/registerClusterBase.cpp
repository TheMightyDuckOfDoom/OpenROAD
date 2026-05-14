// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "registerClusterBase.h"

#include <memory>

#include "cts/TritonCTS.h"
#include "nesterovBase.h"
#include "odb/db.h"
#include "utl/Logger.h"

namespace gpl {

using utl::GPL;

// RegisterClusterBase
RegisterClusterBase::RegisterClusterBase() = default;

RegisterClusterBase::RegisterClusterBase(std::shared_ptr<NesterovBaseCommon> nbc,
              cts::TritonCTS* cts,
                                         utl::Logger* log)
    : RegisterClusterBase()
{
  nbc_ = std::move(nbc);
  cts_ = cts;
  log_ = log;
}

void RegisterClusterBase::executeRegisterClustering()
{
  log_->info(GPL, 9999, "---- Execute Register Clustering.");
  if (cts_ == nullptr) {
    log_->error(GPL, 9998, "TritonCTS instance is null. Cannot execute register clustering.");
    return;
  }
  cts_->runTritonCts();
  nbc_->fixPointers();
}

}  // namespace gpl
