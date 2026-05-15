source helpers.tcl
set test_name large02
read_liberty ./library/nangate45/NangateOpenCellLibrary_typical.lib

read_lef ./nangate45.lef
read_def ./$test_name.def

create_clock -name core_clock -period 2 clk

set_wire_rc -signal -layer metal3
set_wire_rc -clock -layer metal5

set_cts_config \
    -sink_clustering_size 10 \
    -root_buf BUF_X16 \
    -buf_list BUF_X16
cts::set_dummy_load false

global_placement -timing_driven -register_clustering
gui::show

estimate_parasitics -placement
report_worst_slack

set def_file [make_result_file $test_name.def]
write_def $def_file
diff_file $def_file $test_name.defok
source report_hpwl.tcl
