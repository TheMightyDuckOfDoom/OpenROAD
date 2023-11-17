///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (c) 2023, Tobias Senti
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include "rop/Rop.hh"

#include <iostream>

#include "odb/db.h"
#include "utl/Logger.h"

namespace rop {

using utl::ROP;

Rop::Rop()
{
}

void Rop::init(odb::dbDatabase* db, utl::Logger* logger)
{
  db_ = db;
  block_ = nullptr;
  logger_ = logger;
}

double Rop::dbuToMicrons(int64_t dbu)
{
  // Check if block_ is valid 
  if (block_ == nullptr) {
    logger_->error(ROP, 2, "dbuToMicrons: block_ is nullptr.");
    return 0.0;
  }
  return (double) dbu / block_->getDbUnitsPerMicron();
}

void Rop::optimize_net_routing(odb::dbNet* net)
{
  // Get associated wire
  odb::dbWire* wire = net->getWire();

  // Check if wire exists
  if (wire == nullptr) {
    logger_->warn(
      ROP, 3, "Net {} does not have detailed route.", net->getName());
    return;
  }

  // Get dbBlock -> Needed in dbuToMicrons
  block_ = db_->getChip()->getBlock();

  // Print net name
  logger_->info(ROP, 4, "Net: {}", net->getName());

  // Print starting wire length
  int64_t starting_wl = wire->getLength();
  logger_->info(ROP, 5, "\tStarting Length: {:.2f}um", dbuToMicrons(starting_wl));
}

}  // namespace rop
