// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2020-2025, The OpenROAD Authors

#include "ant/AntennaChecker.hh"

#include <omp.h>
#include <tcl.h>

#include <algorithm>
#include <boost/pending/disjoint_sets.hpp>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <queue>
#include <set>
#include <utility>
#include <vector>

#include "Polygon.hh"
#include "odb/db.h"
#include "odb/dbShape.h"
#include "odb/dbTypes.h"
#include "utl/Logger.h"

namespace ant {

using utl::ANT;

// Abbreviations Index:
//   `PAR`: Partial Area Ratio
//   `CAR`: Cumulative Area Ratio
//   `Area`: Gate Area
//   `S. Area`: Side Diffusion Area
//   `C. Area`: Cumulative Gate Area
//   `C. S. Area`: Cumulative Side (Diffusion) Area

struct AntennaModel
{
  odb::dbTechLayer* layer;

  double metal_factor;
  double diff_metal_factor;

  double cut_factor;
  double diff_cut_factor;

  double side_metal_factor;
  double diff_side_metal_factor;

  double minus_diff_factor;
  double plus_diff_factor;
  double diff_metal_reduce_factor;
};

extern "C" {
extern int Ant_Init(Tcl_Interp* interp);
}

AntennaChecker::AntennaChecker() = default;
AntennaChecker::~AntennaChecker() = default;

void AntennaChecker::init(odb::dbDatabase* db,
                          GlobalRouteSource* global_route_source,
                          utl::Logger* logger)
{
  db_ = db;
  global_route_source_ = global_route_source;
  logger_ = logger;
}

void AntennaChecker::initAntennaRules()
{
  block_ = db_->getChip()->getBlock();
  odb::dbTech* tech = db_->getTech();
  // initialize nets_to_report_ with all nets to avoid issues with
  // multithreading
  if (net_to_report_.empty()) {
    for (odb::dbNet* net : block_->getNets()) {
      if (!net->isSpecial()) {
        net_to_report_[net];
      }
    }
  }

  if (!layer_info_.empty()) {
    return;
  }

  for (odb::dbTechLayer* tech_layer : tech->getLayers()) {
    double metal_factor = 1.0;
    double diff_metal_factor = 1.0;

    double cut_factor = 1.0;
    double diff_cut_factor = 1.0;

    double side_metal_factor = 1.0;
    double diff_side_metal_factor = 1.0;

    double minus_diff_factor = 0.0;
    double plus_diff_factor = 0.0;
    double diff_metal_reduce_factor = 1.0;

    if (tech_layer->hasDefaultAntennaRule()) {
      const odb::dbTechLayerAntennaRule* antenna_rule
          = tech_layer->getDefaultAntennaRule();

      if (antenna_rule->isAreaFactorDiffUseOnly()) {
        diff_metal_factor = antenna_rule->getAreaFactor();

        diff_cut_factor = antenna_rule->getAreaFactor();
      } else {
        metal_factor = antenna_rule->getAreaFactor();
        diff_metal_factor = antenna_rule->getAreaFactor();

        cut_factor = antenna_rule->getAreaFactor();
        diff_cut_factor = antenna_rule->getAreaFactor();
      }
      if (antenna_rule->isSideAreaFactorDiffUseOnly()) {
        diff_side_metal_factor = antenna_rule->getSideAreaFactor();
      } else {
        side_metal_factor = antenna_rule->getSideAreaFactor();
        diff_side_metal_factor = antenna_rule->getSideAreaFactor();
      }

      minus_diff_factor = antenna_rule->getAreaMinusDiffFactor();
      plus_diff_factor = antenna_rule->getGatePlusDiffFactor();

      const double PSR_ratio = antenna_rule->getPSR();
      const odb::dbTechLayerAntennaRule::pwl_pair diffPSR
          = antenna_rule->getDiffPSR();

      uint wire_thickness_dbu = 0;
      tech_layer->getThickness(wire_thickness_dbu);

      const odb::dbTechLayerType layerType = tech_layer->getType();

      // If there is a SIDE area antenna rule, then make sure thickness exists.
      if ((PSR_ratio != 0 || !diffPSR.indices.empty())
          && layerType == odb::dbTechLayerType::ROUTING
          && wire_thickness_dbu == 0) {
        logger_->warn(ANT,
                      13,
                      "No THICKNESS is provided for layer {}.  Checks on this "
                      "layer will not be correct.",
                      tech_layer->getConstName());
      }
    }

    AntennaModel layer_antenna = {tech_layer,
                                  metal_factor,
                                  diff_metal_factor,
                                  cut_factor,
                                  diff_cut_factor,
                                  side_metal_factor,
                                  diff_side_metal_factor,
                                  minus_diff_factor,
                                  plus_diff_factor,
                                  diff_metal_reduce_factor};
    layer_info_[tech_layer] = layer_antenna;
  }
}

double AntennaChecker::gateArea(odb::dbMTerm* mterm)
{
  double max_gate_area = 0;
  if (mterm->hasDefaultAntennaModel()) {
    odb::dbTechAntennaPinModel* pin_model = mterm->getDefaultAntennaModel();
    std::vector<std::pair<double, odb::dbTechLayer*>> gate_areas;
    pin_model->getGateArea(gate_areas);

    for (const auto& [gate_area, layer] : gate_areas) {
      max_gate_area = std::max(max_gate_area, gate_area);
    }
  }
  return max_gate_area;
}

double AntennaChecker::getPwlFactor(
    odb::dbTechLayerAntennaRule::pwl_pair pwl_info,
    double ref_value,
    double default_value)
{
  if (!pwl_info.indices.empty()) {
    if (pwl_info.indices.size() == 1) {
      return pwl_info.ratios[0];
    }
    double pwl_info_index1 = pwl_info.indices[0];
    double pwl_info_ratio1 = pwl_info.ratios[0];
    double slope = 1.0;
    for (int i = 0; i < pwl_info.indices.size(); i++) {
      double pwl_info_index2 = pwl_info.indices[i];
      double pwl_info_ratio2 = pwl_info.ratios[i];
      slope = (pwl_info_ratio2 - pwl_info_ratio1)
              / (pwl_info_index2 - pwl_info_index1);

      if (ref_value >= pwl_info_index1 && ref_value < pwl_info_index2) {
        return pwl_info_ratio1 + (ref_value - pwl_info_index1) * slope;
      }
      pwl_info_index1 = pwl_info_index2;
      pwl_info_ratio1 = pwl_info_ratio2;
    }
    return pwl_info_ratio1 + (ref_value - pwl_info_index1) * slope;
  }
  return default_value;
}

void AntennaChecker::saveGates(odb::dbNet* db_net,
                               LayerToGraphNodes& node_by_layer_map,
                               const int node_count)
{
  std::map<PinType, std::vector<int>, PinTypeCmp> pin_nbrs;
  std::vector<int> ids;
  // iterate all instance pins
  for (odb::dbITerm* iterm : db_net->getITerms()) {
    odb::dbMTerm* mterm = iterm->getMTerm();
    std::string pin_name = fmt::format("  {}/{} ({})",
                                       iterm->getInst()->getConstName(),
                                       mterm->getConstName(),
                                       mterm->getMaster()->getConstName());
    PinType pin = PinType(std::move(pin_name), iterm);
    odb::dbInst* inst = iterm->getInst();
    const odb::dbTransform transform = inst->getTransform();
    for (odb::dbMPin* mterm : mterm->getMPins()) {
      for (odb::dbBox* box : mterm->getGeometry()) {
        odb::dbTechLayer* tech_layer = box->getTechLayer();
        if (tech_layer->getType() != odb::dbTechLayerType::ROUTING) {
          continue;
        }
        // get lower and upper layer
        odb::dbTechLayer* upper_layer = tech_layer->getUpperLayer();
        odb::dbTechLayer* lower_layer = tech_layer->getLowerLayer();

        odb::Rect pin_rect = box->getBox();
        transform.apply(pin_rect);
        // convert rect -> polygon
        Polygon pin_pol = rectToPolygon(pin_rect);
        // if has wire on same layer connect to pin
        ids = findNodesWithIntersection(node_by_layer_map[tech_layer], pin_pol);
        for (const int& index : ids) {
          pin_nbrs[pin].push_back(node_by_layer_map[tech_layer][index]->id);
        }
        // if has via on upper layer connected to pin
        if (upper_layer) {
          ids = findNodesWithIntersection(node_by_layer_map[upper_layer],
                                          pin_pol);
          for (const int& index : ids) {
            pin_nbrs[pin].push_back(node_by_layer_map[upper_layer][index]->id);
          }
        }
        // if has via on lower layer connected to pin
        if (lower_layer) {
          ids = findNodesWithIntersection(node_by_layer_map[lower_layer],
                                          pin_pol);
          for (const int& index : ids) {
            pin_nbrs[pin].push_back(node_by_layer_map[lower_layer][index]->id);
          }
        }
      }
    }
  }
  // run DSU from min_layer to max_layer
  std::vector<int> dsu_parent(node_count);
  std::vector<int> dsu_size(node_count);
  for (int i = 0; i < node_count; i++) {
    dsu_size[i] = 1;
    dsu_parent[i] = i;
  }

  boost::disjoint_sets<int*, int*> dsu(&dsu_size[0], &dsu_parent[0]);

  odb::dbTech* tech = db_->getTech();
  odb::dbTechLayer* iter = tech->findRoutingLayer(1);
  odb::dbTechLayer* lower_layer;
  while (iter) {
    // iterate each node of this layer to union set
    for (auto& node_it : node_by_layer_map[iter]) {
      int id_u = node_it->id;
      // if has lower layer
      lower_layer = iter->getLowerLayer();
      if (lower_layer) {
        // get lower neighbors and union
        for (const int& lower_it : node_it->low_adj) {
          int id_v = node_by_layer_map[lower_layer][lower_it]->id;
          // if they are on different sets then union
          if (dsu.find_set(id_u) != dsu.find_set(id_v)) {
            dsu.union_set(id_u, id_v);
          }
        }
      }
    }
    for (auto& node_it : node_by_layer_map[iter]) {
      int id_u = node_it->id;
      // check gates in same set (first Nodes x gates)
      for (const auto& gate_it : pin_nbrs) {
        for (const int& nbr_id : gate_it.second) {
          if (dsu.find_set(id_u) == dsu.find_set(nbr_id)) {
            node_it->gates.insert(gate_it.first);
            break;
          }
        }
      }
    }
    iter = iter->getUpperLayer();
  }
}

bool AntennaChecker::isValidGate(odb::dbMTerm* mterm)
{
  if (mterm->getIoType() != odb::dbIoType::INPUT) {
    return false;
  }
  if (gateArea(mterm) <= 0.0) {
    if (diffArea(mterm) <= 0.0) {
      // Antenna Diodes have a diffusion area at their inputs
      debugPrint(
          logger_,
          ANT,
          "check_gates",
          1,
          "Gate Area and Diffusion Area for Input {}:{} is 0.0 or not defined. "
          "This pin will be ignored during antenna checking.",
          mterm->getMaster()->getConstName(),
          mterm->getConstName());
    }
    return false;
  }
  return true;
}

void AntennaChecker::calculateWirePar(NodeInfo& info)
{
  // get info from layer map
  odb::dbTechLayer* tech_layer = info.layer;
  const double diff_metal_factor = layer_info_[tech_layer].diff_metal_factor;
  const double diff_side_metal_factor
      = layer_info_[tech_layer].diff_side_metal_factor;
  const double minus_diff_factor = layer_info_[tech_layer].minus_diff_factor;
  const double plus_diff_factor = layer_info_[tech_layer].plus_diff_factor;

  const double metal_factor = layer_info_[tech_layer].metal_factor;
  const double side_metal_factor = layer_info_[tech_layer].side_metal_factor;

  double diff_metal_reduce_factor = 1.0;
  if (tech_layer->hasDefaultAntennaRule()) {
    const odb::dbTechLayerAntennaRule* antenna_rule
        = tech_layer->getDefaultAntennaRule();
    diff_metal_reduce_factor = getPwlFactor(
        antenna_rule->getAreaDiffReduce(), info.iterm_diff_area, 1.0);
  }

  info.PAR = 0.0;
  info.PSR = 0.0;
  info.diff_PAR = 0.0;
  info.diff_PSR = 0.0;

  if (info.iterm_gate_area != 0) {
    if (info.iterm_diff_area != 0) {
      // Calculate PAR
      info.PAR = (diff_metal_factor * info.area) / info.iterm_gate_area;
      info.PSR
          = (diff_side_metal_factor * info.side_area) / info.iterm_gate_area;

      // Calculate PSR
      info.diff_PAR
          = (diff_metal_factor * info.area * diff_metal_reduce_factor
             - minus_diff_factor * info.iterm_diff_area)
            / (info.iterm_gate_area + plus_diff_factor * info.iterm_diff_area);
      info.diff_PSR
          = (diff_side_metal_factor * info.side_area * diff_metal_reduce_factor
             - minus_diff_factor * info.iterm_diff_area)
            / (info.iterm_gate_area + plus_diff_factor * info.iterm_diff_area);
    } else {
      // Calculate PAR
      info.PAR = (metal_factor * info.area) / info.iterm_gate_area;
      info.PSR = (side_metal_factor * info.side_area) / info.iterm_gate_area;

      // Calculate PSR
      info.diff_PAR = (metal_factor * info.area * diff_metal_reduce_factor)
                      / info.iterm_gate_area;
      info.diff_PSR
          = (side_metal_factor * info.side_area * diff_metal_reduce_factor)
            / info.iterm_gate_area;
    }
  }
}

void AntennaChecker::calculateViaPar(NodeInfo& info)
{
  // get info from layer map
  odb::dbTechLayer* tech_layer = info.layer;
  const double diff_cut_factor = layer_info_[tech_layer].diff_cut_factor;
  const double minus_diff_factor = layer_info_[tech_layer].minus_diff_factor;
  const double plus_diff_factor = layer_info_[tech_layer].plus_diff_factor;
  const double cut_factor = layer_info_[tech_layer].cut_factor;

  double diff_metal_reduce_factor = 1.0;
  if (tech_layer->hasDefaultAntennaRule()) {
    const odb::dbTechLayerAntennaRule* antenna_rule
        = tech_layer->getDefaultAntennaRule();
    diff_metal_reduce_factor = getPwlFactor(
        antenna_rule->getAreaDiffReduce(), info.iterm_diff_area, 1.0);
  }

  info.PAR = 0.0;
  info.PSR = 0.0;
  info.diff_PAR = 0.0;
  info.diff_PSR = 0.0;

  if (info.iterm_gate_area != 0) {
    if (info.iterm_diff_area != 0) {
      // Calculate PAR
      info.PAR = (diff_cut_factor * info.area) / info.iterm_gate_area;
      // Calculate diff_PAR
      info.diff_PAR
          = (diff_cut_factor * info.area * diff_metal_reduce_factor
             - minus_diff_factor * info.iterm_diff_area)
            / (info.iterm_gate_area + plus_diff_factor * info.iterm_diff_area);
    } else {
      // Calculate PAR
      info.PAR = (cut_factor * info.area) / info.iterm_gate_area;
      // Calculate diff_PAR
      info.diff_PAR = (cut_factor * info.area * diff_metal_reduce_factor)
                      / info.iterm_gate_area;
    }
  }
}

void AntennaChecker::calculateAreas(const LayerToGraphNodes& node_by_layer_map,
                                    std::vector<NodeInfo>& node_info_list)
{
  for (const auto& it : node_by_layer_map) {
    for (const auto& node_it : it.second) {
      NodeInfo info;

      // Add the gates
      int gates_count = 0;
      std::vector<odb::dbITerm*> iterms;
      for (const auto& gate : node_it->gates) {
        if (!gate.isITerm) {
          continue;
        }

        if (isValidGate(gate.iterm->getMTerm())) {
          info.iterms.push_back(gate.iterm);
        }
        info.iterm_gate_area += gateArea(gate.iterm->getMTerm());
        info.iterm_diff_area += diffArea(gate.iterm->getMTerm());
        gates_count++;
      }
      if (gates_count == 0) {
        continue;
      }

      // Set the area
      double area = gtl::area(node_it->pol);
      // convert from dbu^2 to microns^2
      area = block_->dbuToMicrons(area);
      area = block_->dbuToMicrons(area);
      info.area = area;

      // Calculate the side area of the node
      if (it.first->getRoutingLevel() != 0) {
        // Calculate side area of wire
        uint wire_thickness_dbu = 0;
        it.first->getThickness(wire_thickness_dbu);
        double wire_thickness = block_->dbuToMicrons(wire_thickness_dbu);
        info.side_area = block_->dbuToMicrons(gtl::perimeter(node_it->pol)
                                              * wire_thickness);
      }

      // Set the layer
      info.layer = it.first;

      // If we have the same layer and share atleast one gate, then we can merge
      // it
      bool can_merge = false;
      for (NodeInfo& other_info : node_info_list) {
        // Check if the layer is the same
        if (other_info.layer != info.layer) {
          continue;
        }

        // Check if we share a gate
        for (odb::dbITerm* gate : info.iterms) {
          if (std::find(
                  other_info.iterms.begin(), other_info.iterms.end(), gate)
              != other_info.iterms.end()) {
            can_merge = true;
            break;
          }
        }

        // Merge it with the other node
        if (can_merge) {
          other_info += info;
          break;
        }
      }

      // Add node to list if we have not merged it with another node
      if (!can_merge) {
        node_info_list.push_back(std::move(info));
      }
    }
  }
}

// calculate PAR and PSR of wires and vias
void AntennaChecker::calculatePAR(std::vector<NodeInfo>& node_info_list)
{
  for (NodeInfo& node_info : node_info_list) {
    if (node_info.layer->getRoutingLevel() == 0) {
      calculateViaPar(node_info);
    } else {
      calculateWirePar(node_info);
    }
  }
}

// calculate CAR and CSR of wires and vias
void AntennaChecker::calculateCAR(std::vector<NodeInfo>& node_info_list)
{
  odb::dbTech* tech = db_->getTech();
  assert(tech != nullptr);
  odb::dbTechLayer* iter_layer = tech->findRoutingLayer(1);
  assert(iter_layer != nullptr);

  // Helper map
  std::map<odb::dbTechLayer*, std::vector<NodeInfo*>> layer_to_node_info;
  for (NodeInfo& node_info : node_info_list) {
    if (layer_to_node_info.find(node_info.layer) == layer_to_node_info.end()) {
      layer_to_node_info[node_info.layer] = std::vector<NodeInfo*>();
    }
    layer_to_node_info[node_info.layer].push_back(&node_info);
  }

  // Loop from lowest layer to highest layer
  while (iter_layer) {
    // lower Layer
    odb::dbTechLayer* lower_layer = iter_layer->getLowerLayer();

    // Loop over all nodes in this layer, add the PAR of the lower connected
    // nodes
    for (NodeInfo* node_info : layer_to_node_info[iter_layer]) {
      // CAR(m) = PAR(m) + CAR(m-1)
      node_info->CAR = node_info->PAR;
      node_info->CSR = node_info->PSR;
      node_info->diff_CAR = node_info->diff_PAR;
      node_info->diff_CSR = node_info->diff_PSR;

      for (NodeInfo* lower_node_info : layer_to_node_info[lower_layer]) {
        // Check if the node is connected to the lower layer -> share atleast
        // one gate
        // TODO: This might be wrong
        bool connected = false;
        for (odb::dbITerm* gate : node_info->iterms) {
          if (std::find(lower_node_info->iterms.begin(),
                        lower_node_info->iterms.end(),
                        gate)
              != lower_node_info->iterms.end()) {
            connected = true;
            break;
          }
        }

        // More conservative is to always add the full area of the lower layer
        // bool connected = true;

        if (connected) {
          // Add the CAR of the lower layer to the current layer
          node_info->CAR += lower_node_info->CAR;
          node_info->CSR += lower_node_info->CSR;
          node_info->diff_CAR += lower_node_info->diff_CAR;
          node_info->diff_CSR += lower_node_info->diff_CSR;
        }
      }
    }

    // Get next layer, nullptr if we are at the top layer
    iter_layer = iter_layer->getUpperLayer();
  }
}

bool AntennaChecker::checkPAR(odb::dbNet* db_net,
                              NodeInfo& info,
                              const float ratio_margin,
                              bool verbose,
                              bool report,
                              ViolationReport& net_report)
{
  // get rules
  const odb::dbTechLayer* tech_layer = info.layer;
  const odb::dbTechLayerAntennaRule* antenna_rule
      = tech_layer->getDefaultAntennaRule();
  double PAR_ratio = antenna_rule->getPAR();
  odb::dbTechLayerAntennaRule::pwl_pair diffPAR = antenna_rule->getDiffPAR();
  double diff_PAR_PWL_ratio = getPwlFactor(diffPAR, info.iterm_diff_area, 0.0);
  bool violation = false;
  bool diff_violation = false;

  // apply ratio_margin
  PAR_ratio *= (1.0 - ratio_margin / 100.0);
  diff_PAR_PWL_ratio *= (1.0 - ratio_margin / 100.0);

  // check PAR and diff_PAR
  if (PAR_ratio != 0) {
    // Is only a violation if it is not connected to a diffusion
    if (info.iterm_diff_area == 0) {
      violation = info.PAR > PAR_ratio;
    }
  }
  if (report) {
    std::string par_report = fmt::format(
        "      Partial area ratio: {:7.2f}\n      Required "
        "ratio: "
        "{:7.2f} "
        "(Gate area){}{}",
        info.PAR,
        PAR_ratio,
        violation ? " (VIOLATED)" : "",
        info.iterm_diff_area != 0 ? " (ignored as connected Diffusion area > 0)"
                                  : "");
    net_report.report += par_report + "\n";
  }
  if (diff_PAR_PWL_ratio != 0) {
    // This is always checked, not matter the diffusion area
    diff_violation = info.diff_PAR > diff_PAR_PWL_ratio;
  }
  if (report) {
    std::string diff_par_report = fmt::format(
        "      Partial diffusion area ratio: {:7.2f}\n      "
        "Required "
        "ratio: "
        "{:7.2f} "
        "(Gate area){}",
        info.diff_PAR,
        diff_PAR_PWL_ratio,
        diff_violation ? " (VIOLATED)" : "");
    net_report.report += diff_par_report + "\n";
  }
  return violation || diff_violation;
}

bool AntennaChecker::checkPSR(odb::dbNet* db_net,
                              NodeInfo& info,
                              const float ratio_margin,
                              bool verbose,
                              bool report,
                              ViolationReport& net_report)
{
  // get rules
  const odb::dbTechLayer* tech_layer = info.layer;
  const odb::dbTechLayerAntennaRule* antenna_rule
      = tech_layer->getDefaultAntennaRule();
  double PSR_ratio = antenna_rule->getPSR();
  const odb::dbTechLayerAntennaRule::pwl_pair diffPSR
      = antenna_rule->getDiffPSR();
  double diff_PSR_PWL_ratio = getPwlFactor(diffPSR, info.iterm_diff_area, 0.0);
  bool violation = false;
  bool diff_violation = false;

  // apply ratio_margin
  PSR_ratio *= (1.0 - ratio_margin / 100.0);
  diff_PSR_PWL_ratio *= (1.0 - ratio_margin / 100.0);

  // check PSR and diff_PSR
  if (PSR_ratio != 0) {
    // Is only a violation if it is not connected to a diffusion
    if (info.iterm_diff_area == 0) {
      violation = info.PSR > PSR_ratio;
    }
  }
  if (report) {
    std::string psr_report = fmt::format(
        "      Partial area ratio: {:7.2f}\n      Required ratio: "
        "{:7.2f} "
        "(Side area){}{}",
        info.PSR,
        PSR_ratio,
        violation ? " (VIOLATED)" : "",
        info.iterm_diff_area != 0 ? " (ignored as connected Diffusion area > 0)"
                                  : "");
    net_report.report += psr_report + "\n";
  }
  if (diff_PSR_PWL_ratio != 0) {
    // This is always checked, not matter the diffusion area
    diff_violation = info.diff_PSR > diff_PSR_PWL_ratio;
  }
  if (report) {
    std::string diff_psr_report = fmt::format(
        "      Partial diffusion area ratio: {:7.2f}\n      Required "
        "ratio: "
        "{:7.2f} "
        "(Side area){}",
        info.diff_PSR,
        diff_PSR_PWL_ratio,
        diff_violation ? " (VIOLATED)" : "");
    net_report.report += diff_psr_report + "\n";
  }
  return violation || diff_violation;
}

bool AntennaChecker::checkCAR(odb::dbNet* db_net,
                              const NodeInfo& info,
                              bool verbose,
                              bool report,
                              ViolationReport& net_report)
{
  // get rules
  const odb::dbTechLayer* tech_layer = info.layer;
  const odb::dbTechLayerAntennaRule* antenna_rule
      = tech_layer->getDefaultAntennaRule();
  const double CAR_ratio = antenna_rule->getCAR();
  const odb::dbTechLayerAntennaRule::pwl_pair diffCAR
      = antenna_rule->getDiffCAR();
  const double diff_CAR_PWL_ratio
      = getPwlFactor(diffCAR, info.iterm_diff_area, 0);
  bool violation = false;
  bool diff_violation = false;

  // check CAR and diff_CAR
  if (CAR_ratio != 0) {
    // Is only a violation if it is not connected to a diffusion
    if (info.iterm_diff_area == 0) {
      violation = info.CAR > CAR_ratio;
    }
  }
  if (report) {
    std::string car_report = fmt::format(
        "      Cumulative area ratio: {:7.2f}\n      Required ratio: "
        "{:7.2f} "
        "(Cumulative area){}{}",
        info.CAR,
        CAR_ratio,
        violation ? " (VIOLATED)" : "",
        info.iterm_diff_area != 0 ? " (ignored as connected Diffusion area > 0)"
                                  : "");
    net_report.report += car_report + "\n";
  }
  if (diff_CAR_PWL_ratio != 0) {
    // This is always checked, not matter the diffusion area
    diff_violation = info.diff_CAR > diff_CAR_PWL_ratio;
  }
  if (report) {
    std::string diff_car_report = fmt::format(
        "      Cumulative diffusion area ratio: {:7.2f}\n      Required "
        "ratio: "
        "{:7.2f} "
        "(Cumulative area){}",
        info.diff_CAR,
        diff_CAR_PWL_ratio,
        diff_violation ? " (VIOLATED)" : "");
    net_report.report += diff_car_report + "\n";
  }
  return violation || diff_violation;
}

bool AntennaChecker::checkCSR(odb::dbNet* db_net,
                              const NodeInfo& info,
                              bool verbose,
                              bool report,
                              ViolationReport& net_report)
{
  // get rules
  const odb::dbTechLayer* tech_layer = info.layer;
  const odb::dbTechLayerAntennaRule* antenna_rule
      = tech_layer->getDefaultAntennaRule();
  const double CSR_ratio = antenna_rule->getCSR();
  const odb::dbTechLayerAntennaRule::pwl_pair diffCSR
      = antenna_rule->getDiffCSR();
  const double diff_CSR_PWL_ratio
      = getPwlFactor(diffCSR, info.iterm_diff_area, 0);
  bool violation = false;
  bool diff_violation = false;

  // check CSR and diff_CSR
  if (CSR_ratio != 0) {
    // Is only a violation if it is not connected to a diffusion
    if (info.iterm_diff_area == 0) {
      violation = info.CSR > CSR_ratio;
    }
  }
  if (report) {
    std::string csr_report = fmt::format(
        "      Cumulative area ratio: {:7.2f}\n      Required ratio: "
        "{:7.2f} "
        "(Cumulative side area){}{}",
        info.CSR,
        CSR_ratio,
        violation ? " (VIOLATED)" : "",
        info.iterm_diff_area != 0 ? " (ignored as connected Diffusion area > 0)"
                                  : "");
    net_report.report += csr_report + "\n";
  }
  if (diff_CSR_PWL_ratio != 0) {
    // This is always checked, not matter the diffusion area
    diff_violation = info.diff_CSR > diff_CSR_PWL_ratio;
  }
  if (report) {
    std::string diff_csr_report = fmt::format(
        "      Cumulative diffusion area ratio: {:7.2f}\n      Required "
        "ratio: "
        "{:7.2f} "
        "(Cumulative side area) {}",
        info.diff_CSR,
        diff_CSR_PWL_ratio,
        diff_violation ? "(VIOLATED)" : "");
    net_report.report += diff_csr_report + "\n";
  }
  return violation || diff_violation;
}

bool AntennaChecker::checkRatioViolations(odb::dbNet* db_net,
                                          NodeInfo& node_info,
                                          const float ratio_margin,
                                          bool verbose,
                                          bool report,
                                          ViolationReport& net_report)
{
  bool node_has_violation
      = checkPAR(db_net, node_info, ratio_margin, verbose, report, net_report)
        || checkCAR(db_net, node_info, verbose, report, net_report);
  if (node_info.layer->getRoutingLevel() != 0) {
    bool psr_violation = checkPSR(
        db_net, node_info, ratio_margin, verbose, report, net_report);
    bool csr_violation
        = checkCSR(db_net, node_info, verbose, report, net_report);
    node_has_violation = node_has_violation || psr_violation || csr_violation;
  }

  return node_has_violation;
}

void AntennaChecker::writeReport(std::ofstream& report_file, bool verbose)
{
  std::lock_guard<std::mutex> lock(map_mutex_);
  for (const auto& [net, violation_report] : net_to_report_) {
    if (verbose || violation_report.violated) {
      report_file << violation_report.report;
    }
  }
}

void AntennaChecker::printReport(odb::dbNet* db_net)
{
  if (db_net) {
    logger_->report("{}", net_to_report_[db_net].report);
  } else {
    std::lock_guard<std::mutex> lock(map_mutex_);
    for (const auto& [net, violation_report] : net_to_report_) {
      if (violation_report.violated) {
        logger_->report("{}", violation_report.report);
      }
    }
  }
}

int AntennaChecker::checkGates(odb::dbNet* db_net,
                               bool verbose,
                               bool save_report,
                               odb::dbMTerm* diode_mterm,
                               float ratio_margin,
                               std::vector<NodeInfo>& node_info_list,
                               Violations& antenna_violations)
{
  int pin_violation_count = 0;

  std::map<odb::dbITerm*, std::vector<NodeInfo*>> gate_to_nodes;
  std::map<NodeInfo*, ViolationReport> node_reports;
  std::map<odb::dbTechLayer*, std::vector<NodeInfo*>>
      layer_to_nodes_with_violations;

  // Check all nodes for violations
  size_t violation_count = 0;
  for (NodeInfo& node_info : node_info_list) {
    std::string layer_name
        = fmt::format("    Layer: {}", node_info.layer->getConstName());
    ViolationReport node_report;
    node_report.report += layer_name + "\n";
    node_report.violated = checkRatioViolations(
        db_net, node_info, ratio_margin, verbose, true, node_report);

    // Store node for violation fixing
    if (node_report.violated) {
      if (layer_to_nodes_with_violations.find(node_info.layer)
          == layer_to_nodes_with_violations.end()) {
        layer_to_nodes_with_violations[node_info.layer]
            = std::vector<NodeInfo*>();
      }
      layer_to_nodes_with_violations[node_info.layer].push_back(&node_info);
      violation_count++;

      debugPrint(logger_,
                 ANT,
                 "check_gates",
                 1,
                 "Net {} has antenna violation on layer {}: {}",
                 db_net->getConstName(),
                 node_info.layer->getConstName(),
                 node_report.report);
    }

    // Save report of the node
    node_reports[&node_info] = node_report;

    // Add node info to the gates
    for (odb::dbITerm* iterm : node_info.iterms) {
      if (gate_to_nodes.find(iterm) == gate_to_nodes.end()) {
        gate_to_nodes[iterm] = std::vector<NodeInfo*>();
      }
      gate_to_nodes[iterm].push_back(&node_info);
    }
  }

  // Build the report for all gates with violations
  ViolationReport net_report;
  std::string net_name = fmt::format("Net: {}", db_net->getConstName());
  net_report.report += net_name + "\n";

  // Create reports for each gate
  for (auto& [gate, node_infos] : gate_to_nodes) {
    odb::dbMTerm* mterm = gate->getMTerm();
    std::string pin_name = fmt::format("  Pin:   {}/{} ({})",
                                       gate->getInst()->getConstName(),
                                       mterm->getConstName(),
                                       mterm->getMaster()->getConstName());
    net_report.report += pin_name + "\n";

    // Loop over all nodes for this gate
    bool pin_has_violation = false;
    for (NodeInfo* node_info : node_infos) {
      if (node_info->layer->hasDefaultAntennaRule()) {
        net_report.report += node_reports[node_info].report;

        net_report.report += "\n";
        if (node_reports[node_info].violated) {
          pin_has_violation = true;
          net_report.violated = true;
        }
      }
    }
    if (pin_has_violation) {
      pin_violation_count++;
    }
    net_report.report += "\n";
  }

  // Write report on map
  if (save_report) {
    std::lock_guard<std::mutex> lock(map_mutex_);
    net_to_report_.at(db_net) = net_report;
  }

  if (pin_violation_count > 0) {
    debugPrint(logger_,
               ANT,
               "check_gates",
               1,
               "Net {} has {} pins with antenna violations",
               db_net->getConstName(),
               pin_violation_count);
  }

  // Pass 1: Loop over all violations, determine the diode count

  // How many diodes have we added to the gate
  std::map<odb::dbITerm*, int> gate_to_diode_count;

  if (diode_mterm) {
    // Diffusion are of the diode
    const double diode_diff_area = diffArea(diode_mterm);

    // Start at lowest layer with fixing violation
    odb::dbTechLayer* it_layer = db_->getTech()->findRoutingLayer(1);
    assert(it_layer != nullptr);
    while (it_layer) {
      for (NodeInfo* violation_info :
           layer_to_nodes_with_violations[it_layer]) {
        assert(it_layer == violation_info->layer);

        // If no gates in this violation, skip
        if (violation_info->iterms.empty()) {
          logger_->warn(
              ANT,
              304,
              "No gates found in violation info for net {} on layer {}",
              db_net->getConstName(),
              it_layer->getConstName());

          continue;
        }

        // Step 1: Check if there is a violation
        bool violated = checkRatioViolations(
            db_net, *violation_info, ratio_margin, false, false, net_report);

        // How many new diodes do we need to add for this violation?
        int additional_diode_count = 0;
        NodeInfo previous_info = *violation_info;
        {
          // Step 2: Loop until we no longer have a violation or we exceed the
          // maximum
          // TODO: We should be able to calculate the number of diodes needed
          int gate_idx = 0;
          while (violated) {
            // Increment the diode count for gate at gate_idx
            additional_diode_count++;
            odb::dbITerm* gate = violation_info->iterms[gate_idx];

            // Increment the diode count for this gate
            if (gate_to_diode_count.find(gate) != gate_to_diode_count.end()) {
              // Gate is already in the map, update the count
              gate_to_diode_count[gate] += 1;
            } else {
              // New gate, add it to the map
              gate_to_diode_count[gate] = 1;
            }

            // Update all difusion areas that contain this gate, including this
            for (NodeInfo* node_info_to_add : gate_to_nodes[gate]) {
              node_info_to_add->iterm_diff_area += diode_diff_area;
            }

            // Make sure that we actually have increased the diffusion area
            assert(violation_info->iterm_diff_area
                   > previous_info.iterm_diff_area);

            previous_info = *violation_info;

            // Increment the gate_idx, wrap around to first gate if needed
            gate_idx = (gate_idx + 1) >= (int) violation_info->iterms.size()
                           ? 0
                           : gate_idx + 1;

            debugPrint(logger_,
                       ANT,
                       "check_gates",
                       2,
                       "  checking with {} additional diodes to fix "
                       "on net {}",
                       additional_diode_count,
                       db_net->getConstName());

            // Recalculate the wire parameters
            calculatePAR(node_info_list);
            calculateCAR(node_info_list);

            violated = checkRatioViolations(db_net,
                                            *violation_info,
                                            ratio_margin,
                                            false,
                                            false,
                                            net_report);

            // Check if we have exceeded the maximum number of diodes
            if (additional_diode_count
                >= max_diode_count_per_gate
                       * (int) violation_info->iterms.size()) {
              logger_->warn(
                  ANT,
                  305,
                  "Net {} requires more than {} diodes to "
                  "repair violations, more than {} for each of the {} gates.",
                  db_net->getConstName(),
                  additional_diode_count,
                  max_diode_count_per_gate,
                  violation_info->iterms.size());
              break;
            }
          }
          debugPrint(
              logger_,
              ANT,
              "check_gates",
              1,
              "  {} additional diodes needed to fix violation on layer {} "
              "on net {}",
              additional_diode_count,
              it_layer->getConstName(),
              db_net->getConstName());
        }
      }

      // Get the next layer
      it_layer = it_layer->getUpperLayer();
    }
  }

  // Pass 2: Distribute the diode count of the gates to the violations
  // Greedily assign diodes to the violation
  for (const auto& [layer, nodes_with_violations] :
       layer_to_nodes_with_violations) {
    for (NodeInfo* violation_info : nodes_with_violations) {
      // Get the diode count for this violation
      int diode_count = 0;

      // Check all gates inside the violation for their diode count
      for (const auto& gate : violation_info->iterms) {
        if (gate_to_diode_count.find(gate) != gate_to_diode_count.end()) {
          // Gate is in the map, update the count
          diode_count += gate_to_diode_count[gate];

          // Mark the diodes as already added
          gate_to_diode_count[gate] = 0;
        }
      }

      // Add the violation to the list
      std::vector<odb::dbITerm*> gates;
      gates.reserve(violation_info->iterms.size());
      for (const auto& gate : violation_info->iterms) {
        gates.push_back(gate);
      }
      if (diode_count > 0) {
        // Print the violation info
        debugPrint(
            logger_,
            ANT,
            "check_gates",
            1,
            "  {} diodes needed to fix violation #{} on layer {} on net {}",
            diode_count,
            antenna_violations.size(),
            layer->getConstName(),
            db_net->getConstName());
        // Make sure that we have at least one gate
        assert(!gates.empty());
      }
      antenna_violations.push_back(
          {layer->getRoutingLevel(), std::move(gates), diode_count});
    }
  }

  // Make sure that all diodes have been added to violations
  for (const auto& [gate, diode_count] : gate_to_diode_count) {
    if (diode_count > 0) {
      logger_->warn(ANT,
                    306,
                    "Gate {} has {} diodes that were not added to a violation.",
                    gate->getInst()->getConstName(),
                    diode_count);
    }
    assert(diode_count == 0);
  }

  // Make sure that all violations have been added to the list
  assert(antenna_violations.size() == violation_count);

  return pin_violation_count;
}

void AntennaChecker::buildLayerMaps(odb::dbNet* db_net,
                                    LayerToGraphNodes& node_by_layer_map)
{
  odb::dbWire* wires = db_net->getWire();

  std::map<odb::dbTechLayer*, PolygonSet> set_by_layer;

  wiresToPolygonSetMap(wires, set_by_layer);
  avoidPinIntersection(db_net, set_by_layer);

  int node_count = 0;
  for (const auto& layer_it : set_by_layer) {
    for (const auto& pol_it : layer_it.second) {
      bool isVia = layer_it.first->getRoutingLevel() == 0;
      node_by_layer_map[layer_it.first].push_back(
          std::make_unique<GraphNode>(node_count, isVia, pol_it));
      node_count++;
    }
  }

  // set connections between Polygons ( wire -> via -> wire)
  std::vector<int> upper_index, lower_index;
  for (const auto& layer_it : set_by_layer) {
    // iterate only via layers
    if (layer_it.first->getRoutingLevel() == 0) {
      int via_index = 0;
      for (const auto& via_it : layer_it.second) {
        lower_index = findNodesWithIntersection(
            node_by_layer_map[layer_it.first->getLowerLayer()], via_it);
        upper_index = findNodesWithIntersection(
            node_by_layer_map[layer_it.first->getUpperLayer()], via_it);

        if (upper_index.size() <= 2) {
          // connect upper -> via
          for (int& up_index : upper_index) {
            node_by_layer_map[layer_it.first->getUpperLayer()][up_index]
                ->low_adj.push_back(via_index);
          }
        } else if (upper_index.size() > 2) {
          std::string log_error = fmt::format(
              "ERROR: net {} has via on {} conect with multiple wires on "
              "layer "
              "{} \n",
              db_net->getConstName(),
              layer_it.first->getName(),
              layer_it.first->getUpperLayer()->getName());
          logger_->report("{}", log_error);
        }
        if (lower_index.size() == 1) {
          // connect via -> lower
          for (int& low_index : lower_index) {
            node_by_layer_map[layer_it.first][via_index]->low_adj.push_back(
                low_index);
          }
        } else if (lower_index.size() > 2) {
          std::string log_error = fmt::format(
              "ERROR: net {} has via on {} conect with multiple wires on "
              "layer "
              "{} \n",
              db_net->getConstName(),
              layer_it.first->getName(),
              layer_it.first->getLowerLayer()->getName());
          logger_->report("{}", log_error);
        }
        via_index++;
      }
    }
  }
  saveGates(db_net, node_by_layer_map, node_count);
}

int AntennaChecker::checkNet(odb::dbNet* db_net,
                             bool verbose,
                             bool save_report,
                             odb::dbMTerm* diode_mterm,
                             float ratio_margin,
                             Violations& antenna_violations)
{
  odb::dbWire* wire = db_net->getWire();
  int pin_violations = 0;
  if (wire) {
    LayerToGraphNodes node_by_layer_map;
    buildLayerMaps(db_net, node_by_layer_map);

    std::vector<NodeInfo> node_info_list;
    calculateAreas(node_by_layer_map, node_info_list);

    calculatePAR(node_info_list);
    calculateCAR(node_info_list);

    pin_violations = checkGates(db_net,
                                verbose,
                                save_report,
                                diode_mterm,
                                ratio_margin,
                                node_info_list,
                                antenna_violations);
  }
  return pin_violations;
}

Violations AntennaChecker::getAntennaViolations(odb::dbNet* net,
                                                odb::dbMTerm* diode_mterm,
                                                float ratio_margin)
{
  Violations antenna_violations;
  if (net->isSpecial()) {
    return antenna_violations;
  }

  checkNet(net, false, false, diode_mterm, ratio_margin, antenna_violations);

  return antenna_violations;
}

int AntennaChecker::checkAntennas(odb::dbNet* net,
                                  const int num_threads,
                                  bool verbose)
{
  {
    std::lock_guard<std::mutex> lock(map_mutex_);
    net_to_report_.clear();
  }
  initAntennaRules();

  std::ofstream report_file;
  if (!report_file_name_.empty()) {
    report_file.open(report_file_name_, std::ofstream::out);
  }

  bool drt_routes = haveRoutedNets();
  bool grt_routes = false;
  if (!drt_routes) {
    grt_routes = global_route_source_->haveRoutes();
  }
  bool use_grt_routes = (grt_routes && !drt_routes);
  if (!grt_routes && !drt_routes) {
    logger_->error(ANT,
                   8,
                   "No detailed or global routing found. Run global_route or "
                   "detailed_route first.");
  }

  if (use_grt_routes) {
    global_route_source_->makeNetWires();
  }

  int net_violation_count = 0;
  int pin_violation_count = 0;

  if (net) {
    Violations antenna_violations;
    if (!net->isSpecial()) {
      pin_violation_count
          += checkNet(net, verbose, true, nullptr, 0, antenna_violations);
      if (pin_violation_count > 0) {
        net_violation_count++;
      }
    } else {
      logger_->error(
          ANT, 14, "Skipped net {} because it is special.", net->getName());
    }
  } else {
    nets_.clear();
    for (odb::dbNet* net : block_->getNets()) {
      if (!net->isSpecial()) {
        nets_.push_back(net);
      }
    }
    omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < nets_.size(); i++) {
      odb::dbNet* net = nets_[i];
      Violations antenna_violations;
      int pin_viol_count
          = checkNet(net, verbose, true, nullptr, 0, antenna_violations);
      if (pin_viol_count > 0) {
        std::lock_guard<std::mutex> lock(map_mutex_);
        net_violation_count++;
        pin_violation_count += pin_viol_count;
      }
    }
  }

  if (verbose) {
    printReport(net);
  }

  logger_->info(ANT, 2, "Found {} net violations.", net_violation_count);
  logger_->metric("antenna__violating__nets", net_violation_count);
  logger_->info(ANT, 1, "Found {} pin violations.", pin_violation_count);
  logger_->metric("antenna__violating__pins", pin_violation_count);

  if (!report_file_name_.empty()) {
    writeReport(report_file, verbose);
    report_file.close();
  }

  if (use_grt_routes) {
    global_route_source_->destroyNetWires();
  }

  net_violation_count_ = net_violation_count;
  return net_violation_count;
}

int AntennaChecker::antennaViolationCount() const
{
  return net_violation_count_;
}

bool AntennaChecker::haveRoutedNets()
{
  for (odb::dbNet* net : block_->getNets()) {
    if (!net->isSpecial() && net->getWireType() == odb::dbWireType::ROUTED
        && net->getWire()) {
      return true;
    }
  }
  return false;
}

double AntennaChecker::diffArea(odb::dbMTerm* mterm)
{
  double max_diff_area = 0.0;
  std::vector<std::pair<double, odb::dbTechLayer*>> diff_areas;
  mterm->getDiffArea(diff_areas);
  for (const auto& [area, layer] : diff_areas) {
    max_diff_area = std::max(max_diff_area, area);
  }
  return max_diff_area;
}

void AntennaChecker::setReportFileName(const char* file_name)
{
  report_file_name_ = file_name;
}

}  // namespace ant
