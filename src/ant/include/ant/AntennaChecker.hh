// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2020-2025, The OpenROAD Authors

#pragma once

#include <map>
#include <mutex>
#include <queue>
#include <set>

#include "odb/db.h"
#include "odb/dbWireGraph.h"
#include "utl/Logger.h"

namespace grt {
class GlobalRouter;
}

namespace ant {

using utl::Logger;

struct PARinfo;
struct ARinfo;
struct AntennaModel;

///////////////////////////////////////
struct GraphNode;

struct NodeInfo
{
  double area;
  double side_area;
  double iterm_gate_area;
  double iterm_diff_area;

  double PAR;
  double PSR;
  double CAR;
  double CSR;

  std::vector<odb::dbITerm*> iterms;
  std::vector<odb::dbITerm*> diff_iterms;
  std::vector<double> diff_areas;

  odb::dbTechLayer* layer;

  NodeInfo& operator+=(const NodeInfo& a)
  {
    assert(layer == a.layer);
    // Can only merge if we have the same gates
    for (const auto& iterm : a.iterms) {
      if(std::find(iterms.begin(), iterms.end(), iterm) == iterms.end()) {
        assert(false);
      }
    }

    // Diffusion terminals can be different, make them shared
    assert(diff_iterms.size() == diff_areas.size());
    assert(a.diff_iterms.size() == a.diff_areas.size());
    size_t a_idx = 0;
    for (const auto &a_diff_iterm : a.diff_iterms) {
      if(std::find(diff_iterms.begin(), diff_iterms.end(), a_diff_iterm)
             != diff_iterms.end()) {
        diff_iterms.push_back(a_diff_iterm);
        const double a_diff_area = a.diff_areas[a_idx];
        diff_areas.push_back(a_diff_area);
        iterm_diff_area += a_diff_area;
      }
      a_idx++;
    }

    area += a.area;
    side_area += a.side_area;

    return *this;
  }
  NodeInfo()
  {
    area = 0.0;
    side_area = 0.0;
    iterm_gate_area = 0.0;
    iterm_diff_area = 0.0;

    PAR = 0.0;
    PSR = 0.0;
    CAR = 0.0;
    CSR = 0.0;
  }
};

struct ViolationReport
{
  bool violated;
  std::string report;
  ViolationReport() { violated = false; }
};

class GlobalRouteSource
{
 public:
  virtual ~GlobalRouteSource() = default;

  virtual bool haveRoutes() = 0;
  virtual void makeNetWires() = 0;
  virtual void destroyNetWires() = 0;
};

struct Violation
{
  int routing_level;
  std::vector<odb::dbITerm*> gates;
  int diode_count;
};

using GraphNodes = std::vector<std::unique_ptr<GraphNode>>;
using LayerToGraphNodes = std::map<odb::dbTechLayer*, GraphNodes>;
using Violations = std::vector<Violation>;
using NodeInfoList = std::vector<NodeInfo>;
using NodeInfoPtrList = std::vector<NodeInfo*>;
using GatesToNodes = std::map<odb::dbITerm*, NodeInfoPtrList>;
using LayerToNodes = std::map<odb::dbTechLayer*, NodeInfoPtrList>;
using NodesToReport = std::map<NodeInfo*, ViolationReport>;
using GateToDiodeCount = std::map<odb::dbITerm*, int>;

class AntennaChecker
{
 public:
  AntennaChecker();
  ~AntennaChecker();

  void init(odb::dbDatabase* db,
            GlobalRouteSource* global_route_source,
            utl::Logger* logger);

  // net nullptr -> check all nets
  int checkAntennas(odb::dbNet* net = nullptr,
                    int num_threads = 1,
                    bool verbose = false);
  int antennaViolationCount() const;
  Violations getAntennaViolations(odb::dbNet* net,
                                  odb::dbMTerm* diode_mterm,
                                  float ratio_margin);
  void initAntennaRules();
  void setReportFileName(const char* file_name);

 private:
  bool haveRoutedNets();
  double getPwlFactor(odb::dbTechLayerAntennaRule::pwl_pair pwl_info,
                      double ref_val,
                      double def);
  double diffArea(odb::dbMTerm* mterm);
  double gateArea(odb::dbMTerm* mterm);
  std::vector<std::pair<double, std::vector<odb::dbITerm*>>> parMaxWireLength(
      odb::dbNet* net,
      int layer);
  std::vector<std::pair<double, std::vector<odb::dbITerm*>>>
  getViolatedWireLength(odb::dbNet* net, int routing_level);
  bool isValidGate(odb::dbMTerm* mterm);
  void buildLayerMaps(odb::dbNet* net, LayerToGraphNodes& node_by_layer_map);
  int checkNet(odb::dbNet* net,
               bool verbose,
               bool save_report,
               odb::dbMTerm* diode_mterm,
               float ratio_margin,
               Violations& antenna_violations);
  void saveGates(odb::dbNet* db_net,
                 LayerToGraphNodes& node_by_layer_map,
                 int node_count);
  void calculateAreas(const LayerToGraphNodes& node_by_layer_map,
                      NodeInfoList& node_info_list);
  void calculatePAR(NodeInfoList& node_info_list);
  void calculateCAR(LayerToNodes& layer_to_node_info);
  bool checkRatioViolations(odb::dbNet* db_net,
                            NodeInfo& node_info,
                            float ratio_margin,
                            bool verbose,
                            bool report,
                            ViolationReport& net_report);
  void writeReport(std::ofstream& report_file, bool verbose);
  void printReport(odb::dbNet* db_net);
  int checkGates(odb::dbNet* db_net,
                 bool verbose,
                 bool save_report,
                 odb::dbMTerm* diode_mterm,
                 float ratio_margin,
                 NodeInfoList& node_info_list,
                 LayerToNodes& layer_to_nodes,
                 Violations& antenna_violations);
  void calculateNodePar(NodeInfo& info);
  bool checkPAR(odb::dbNet* db_net,
                NodeInfo& info,
                float ratio_margin,
                bool verbose,
                bool report,
                ViolationReport& net_report);
  bool checkPSR(odb::dbNet* db_net,
                NodeInfo& info,
                float ratio_margin,
                bool verbose,
                bool report,
                ViolationReport& net_report);
  bool checkCAR(odb::dbNet* db_net,
                const NodeInfo& info,
                bool verbose,
                bool report,
                ViolationReport& net_report);
  bool checkCSR(odb::dbNet* db_net,
                const NodeInfo& info,
                bool verbose,
                bool report,
                ViolationReport& net_report);

  odb::dbDatabase* db_{nullptr};
  odb::dbBlock* block_{nullptr};
  GlobalRouteSource* global_route_source_{nullptr};
  utl::Logger* logger_{nullptr};
  std::map<odb::dbTechLayer*, AntennaModel> layer_info_;
  int net_violation_count_{0};
  std::string report_file_name_;
  std::vector<odb::dbNet*> nets_;
  std::map<odb::dbNet*, ViolationReport> net_to_report_;
  std::mutex map_mutex_;
  // consts
  static constexpr int max_diode_count_per_gate = 10;
};

}  // namespace ant
