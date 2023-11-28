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

#include <algorithm>
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

bool Rop::build_maps(odb::dbWire* wire, segment_map_t &start_to_end_map,
  segment_map_t &end_to_start_map, positions_t &via_map)
{
  // Loop over all shapes
  dbWireShapeItr shapes;
  dbShape s;
  for(shapes.begin(wire); shapes.next(s);) {
    // Check if shape is a via or a segment
    if(s.isVia()) {
      // Via
      Point pos = s.getViaXY();
      logger_->info(ROP, 3, "Via x: {} y: {}", pos.x(), pos.y());
      // Add via to via_map
      via_map.insert(std::pair<Point, dbShape>(s.getViaXY(), s));
    } else {
      // Segment
      Point min = Point(s.xMin(), s.yMin());
      Point max = Point(s.xMax(), s.yMax());

      logger_->info(ROP, 4, "Segment x: {} y: {} to x2: {} y2: {}", min.x(), min.y(), max.x(), max.y());

      int pitch = s.getTechLayer()->getMinWidth();

      logger_->info(ROP, 5, "Pitch: {}", pitch);

      Point start, end;

      // Determine Direction
      if (s.getDX() == pitch) {
        // Segment is Vertical
        logger_->info(ROP, 6, "Segment is vertical");

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
        logger_->info(ROP, 7, "Segment is horizontal");

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
        logger_->warn(ROP, 8, "Cannot determine segment direction");
        return false;
      }

      logger_->info(ROP, 9, "Inserted segment x: {} y: {} to x2: {} y2: {}", start.x(), start.y(), end.x(), end.y());

      start_to_end_map.insert(std::pair<Point, std::pair<Point, dbShape>>(start, std::pair<Point, dbShape>(end, s)));
      end_to_start_map.insert(std::pair<Point, std::pair<Point, dbShape>>(end,   std::pair<Point, dbShape>(start, s)));
    }
  }

  return true;
} 

bool Rop::onSegment(Point p, Point start, Point end) {
    // Max x, y
    int max_x = std::max(start.x(), end.x());
    int max_y = std::max(start.y(), end.y());
    
    // Min x, y
    int min_y = std::min(start.x(), end.x());
    int min_x = std::min(start.y(), end.y());

    if((max_x == min_x) || (max_y == min_y))
      return false;

    if ((p.x() > min_x) && (p.x() < max_x) && (max_y == p.y()))
      return true;

    if ((p.y() > min_y) && (p.y() < max_y) && (max_x == p.x()))
      return true;

    return false;
}

void Rop::optimize_net_routing(dbNet* net)
{
  // Get associated wire
  dbWire* wire = net->getWire();

  // Check if wire exists
  if (wire == nullptr) {
    logger_->warn(
      ROP, 10, "Net {} does not have detailed route.", net->getName());
    return;
  }

  // Get dbBlock -> Needed in dbuToMicrons
  block_ = db_->getChip()->getBlock();

  // Print net name
  logger_->info(ROP, 11, "Net: {}", net->getName());

  // Print starting wire length
  int64_t starting_wl = wire->getLength();
  logger_->info(ROP, 12, "\tStarting Length: {:.2f}um", dbuToMicrons(starting_wl));

  // Create Via and Segment Maps
  segment_map_t start_map;
  segment_map_t end_map;
  positions_t via_map;
  if(!build_maps(wire, start_map, end_map, via_map)) {
    logger_->warn(ROP, 13, "\tError while building maps");
    return;
  }

  // Find paths
  std::vector<Point> paths;
  std::set<Point> junctions;
  for (auto i : start_map) {
    
    // Check start
    Point start = i.first;

    int starts_found = 0;
    auto start_range = start_map.equal_range(start);
    for (auto j = start_range.first; j != start_range.second; ++j) {
      starts_found++;
    }

    auto end_range = end_map.equal_range(start);
    for (auto j = end_range.first; j != end_range.second; ++j) {
      starts_found++;
    }
    
    // Check end
    Point end = i.second.first;

    int ends_found = 0;
    start_range = start_map.equal_range(end);
    for (auto j = start_range.first; j != start_range.second; ++j) {
      ends_found++;
    }

    end_range = end_map.equal_range(end);
    for (auto j = end_range.first; j != end_range.second; ++j) {
      ends_found++;
    }

    //logger_->info(ROP, 14, "Checking segment ({} {}) to ({} {})", start.x(), start.y(), end.x(), end.y());
    //logger_->info(ROP, 15, "Found {} at {} {}", starts_found, start.x(), start.y());
    //logger_->info(ROP, 16, "Found {} at {} {}", ends_found, end.x(), end.y());

    assert(starts_found >= 1);
    assert(ends_found >= 1);

    if(starts_found == 1) {
      // Check if start lies on other segment
      bool on_seg = false;
      for (auto seg : start_map) {
          on_seg = onSegment(start, seg.first, seg.second.first);
          if(on_seg) {
            junctions.insert(start);
            logger_->info(ROP, 22, "On segment {} {}", start.x(), start.y());
            break;
          }
      }

      if(!on_seg)
        paths.push_back(start);
    } else if(starts_found > 2)
      junctions.insert(start);

    if(ends_found == 1) {
      // Check if end lies on other segment
      bool on_seg = false;
      for (auto seg : start_map) {
          on_seg = onSegment(end, seg.first, seg.second.first);
          if(on_seg) {
            junctions.insert(end);
            logger_->info(ROP, 23, "On segment {} {}", end.x(), end.y());
            break;
          }
      }

      if(!on_seg)
        paths.push_back(end);
    } else if(ends_found > 2)
      junctions.insert(end);
  }

  for(Point p : paths) {
    logger_->info(ROP, 20, "Endpoint at ({} {})", p.x(), p.y());
  }

  for(Point p : junctions) {
    logger_->info(ROP, 21, "Junction at ({} {})", p.x(), p.y());
  }

  int expected_paths = net->getTermCount();
  int found_paths = paths.size();

  logger_->info(ROP, 17, "Found {} paths endpoints, expected to find {}", found_paths, expected_paths);
  if(expected_paths != found_paths) {
    logger_->warn(ROP, 18, "Wrong number of paths found");
    return;
  }

  // Remove old wire
  //wire->detach();

  // Create final wire
  dbWire* final_wire = dbWire::create(block_, false);
  odb::dbWireEncoder encoder;
  encoder.begin(final_wire);

  encoder.end();

  final_wire->destroy(final_wire);

  // Print final wire length
  int64_t final_wl = final_wire->getLength();
  logger_->info(ROP, 19, "\tFinal Length: {:.2f}um", dbuToMicrons(final_wl));

  // Add Final wire to net
  //final_wire->attach(net);
}

}  // namespace rop
