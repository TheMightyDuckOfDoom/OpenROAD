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
#include "odb/dbShape.h"
#include "utl/Logger.h"

namespace rop {

using utl::ROP;
using namespace odb;

Rop::Rop()
{
}

void Rop::init(dbDatabase* db, utl::Logger* logger)
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

void Rop::build_maps(dbWire* wire, segment_map_t &segment_map, 
  via_map_t &via_map)
{
  // Loop over all shapes
  dbWireShapeItr shapes;
  dbShape s;
  for(shapes.begin(wire); shapes.next(s);) {
    // Check if shape is a via or a segment
    if(s.isVia()) {
      // Via
      Point pos = s.getViaXY();
      logger_->info(ROP, 6, "Via x: {} y: {}", pos.x(), pos.y());
      // Add via to via_map
      via_map.insert(std::pair<Point, dbShape>(s.getViaXY(), s));
    } else {
      // Segment
      Point min = Point(s.xMin(), s.yMin());
      Point max = Point(s.xMax(), s.yMax());

      logger_->info(ROP, 7, "Segment x: {} y: {} to x2: {} y2: {}", min.x(), min.y(), max.x(), max.y());

      int pitch = s.getTechLayer()->getMinWidth();

      logger_->info(ROP, 8, "Pitch: {}", pitch);

      Point start, end;

      // Determine Direction
      if (s.getDX() == pitch) {
        // Segment is Vertical
        logger_->info(ROP, 9, "Segment is vertical");

        // Center X
        int x = (min.x() + max.x()) / 2;
        start.setX(x);
        end.setX(x);

        // Y
        // Add pitch to minY
        // Subtract pitch from maxY
        start.setY(min.y() + pitch / 2);
        end.setY(max.y() - pitch / 2);
      } else if (s.getDY() == pitch) {
        // Segment is Horizontal
        logger_->info(ROP, 10, "Segment is horizontal");

        // Center Y
        int y = (min.y() + max.y()) / 2;
        start.setY(y);
        end.setY(y);

        // X
        // Add pitch to minX
        // Subtract pitch from maxX
        start.setX(min.x() + pitch / 2);
        end.setX(max.x() - pitch / 2);
      } else {
        // Cannot determine direction
        logger_->info(ROP, 11, "Cannot determine segment direction");
        continue;
      }

      logger_->info(ROP, 12, "Inserted segment x: {} y: {} to x2: {} y2: {}", start.x(), start.y(), end.x(), end.y());

      segment_map.insert(std::pair<Point, odb::dbShape>(start, s));
      segment_map.insert(std::pair<Point, odb::dbShape>(end, s));
    }
  }
} 

void Rop::optimize_net_routing(dbNet* net)
{
  // Get associated wire
  dbWire* wire = net->getWire();

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

  // Create Via and Segment Maps
  segment_map_t segment_map;
  via_map_t via_map;
  build_maps(wire, segment_map, via_map);
}

}  // namespace rop
