from openroad import Design, Tech
import helpers
import dpl_aux

tech = Tech()
tech.readLef("Nangate45_data/Nangate45.lef")
design = helpers.make_design(tech)
design.readDef("one_site_gap_disallow.def")

dpl_aux.detailed_placement(design)
design.getOpendp().checkPlacement(False)

def_file = helpers.make_result_file("one_site_gap_disallow.def")
design.writeDef(def_file)
helpers.diff_files(def_file, "one_site_gap_disallow.defok")
