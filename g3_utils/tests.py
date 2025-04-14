import pytest
import g3_utils as ut

valid_test_cases = [
    (0,             1, "roach1_0000"),
    (1,             1, "roach1_0001"),
    (213.0,         1, "roach1_0213"),
    (213.9,         1, "roach1_0213"),
    (1000,          1, "roach1_1000"),
    ("roach1_0001", 1, "roach1_0001"),
    ("roach3_9999", 3, "roach3_9999"),
    ("roach5_1231", 5, "roach5_1231"),
    ("roach10_1231", 10, "roach10_1231"),
    ("0001",        1, "roach1_0001"),
    ("01",          1, "roach1_0001"),
    ("1",           1, "roach1_0001"),
    ("1.0",         1, "roach1_0001")
]
@pytest.mark.parametrize("val, roach_id, expected", valid_test_cases)
def test_kid_string_valid(val, roach_id, expected):
    result = ut.kid_string(val, roach_id)
    assert result == expected


invalid_test_cases = [
    (55555,          1, ValueError),
    ("roach01_0001", 1, ValueError),
    ("roach3_6",     3, ValueError),
    ("roach3_00010", 3, ValueError),
    ("roach3_0001",  5, ValueError),
    ("_0101",        1, ValueError),
    (lambda a: True, 1, TypeError)
]
@pytest.mark.parametrize("val, roach_id, expected_exception", invalid_test_cases)
def test_kid_string_invalid(val, roach_id, expected_exception):
    with pytest.raises(expected_exception):
        ut.kid_string(val, roach_id)