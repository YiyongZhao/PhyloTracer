from phylotracer.GD_Visualizer import get_present_gd_types


def test_get_present_gd_types_excludes_absent_complex():
    count_dic = {
        "S1": {"AABB": 3, "AABX": 1},
        "S2": {"AXBB": 2},
    }
    assert get_present_gd_types(count_dic) == ["AABB", "AXBB", "AABX"]


def test_get_present_gd_types_keeps_complex_when_present():
    count_dic = {
        "S1": {"AABB": 3},
        "S2": {"Complex": 1},
    }
    assert get_present_gd_types(count_dic) == ["AABB", "Complex"]
