import pytest
from minda.ensemble import _infer_strands, parse_bnd_alt


class TestParseBndAlt:
    """Test BND ALT string parsing"""

    def test_forward_bracket(self):
        """Test parsing forward bracket notation"""
        assert parse_bnd_alt("N[chr2:100[") == "["

    def test_reverse_bracket(self):
        """Test parsing reverse bracket notation"""
        assert parse_bnd_alt("N]chr2:100]") == "]"

    def test_bracket_at_start(self):
        """Test parsing when bracket is at start"""
        assert parse_bnd_alt("[chr2:100[N") == "["
        assert parse_bnd_alt("]chr2:100]N") == "]"

    def test_invalid_alt_raises_assertion(self):
        """Test that invalid ALT strings raise assertion error"""
        with pytest.raises(AssertionError):
            parse_bnd_alt("<INV>")
        with pytest.raises(AssertionError):
            parse_bnd_alt("N")


class TestInferStrands:
    """Test strand inference from ALT field"""

    # DEL tests
    def test_del_simple_notation(self):
        """DEL with simple notation should be +-"""
        assert _infer_strands("DEL", "N") == "+-"

    def test_del_bnd_notation(self):
        """DEL with BND notation should be +-"""
        assert _infer_strands("DEL", "N[chr2:100[") == "+-"
    # INS tests
    def test_ins_simple_notation(self):
        """INS with simple notation should be +-"""
        assert _infer_strands("INS", "N") == "+-"

    def test_ins_bnd_notation(self):
        """INS with BND notation should be +-"""
        assert _infer_strands("INS", "N]chr2:100]") == "+-"

    # INV tests - various BND notations
    def test_inv_reverse_reverse_bracket_before(self):
        """INV: [chr2:100[N should give -- (reverse-reverse)"""
        assert _infer_strands("INV", "[chr2:100[N") == "--"

    def test_inv_forward_forward_bracket_after(self):
        """INV: N]chr2:100] should give ++ (forward-forward)"""
        assert _infer_strands("INV", "N]chr2:100]") == "++"

    # DUP tests
    def test_dup_simple_notation(self):
        """DUP with simple notation should be -+"""
        assert _infer_strands("DUP", "<DUP>") == "-+"

    def test_dup_bnd_notation(self):
        """DUP with BND notation should be -+"""
        assert _infer_strands("DUP", "N[chr2:100[") == "-+"

    # BND tests - all four orientations
    def test_bnd_forward_reverse(self):
        """BND: N[chr2:100[ should give +- (forward-reverse)"""
        assert _infer_strands("BND", "N[chr2:100[") == "+-"

    def test_bnd_forward_forward(self):
        """BND: N]chr2:100] should give ++ (forward-forward)"""
        assert _infer_strands("BND", "N]chr2:100]") == "++"

    def test_bnd_reverse_reverse(self):
        """BND: [chr2:100[N should give -- (reverse-reverse)"""
        assert _infer_strands("BND", "[chr2:100[N") == "--"

    def test_bnd_reverse_forward(self):
        """BND: ]chr2:100]N should give -+ (reverse-forward)"""
        assert _infer_strands("BND", "]chr2:100]N") == "-+"


class TestStrandInferenceValidation:
    """Test input validation"""

    def test_invalid_svtype_raises_error(self):
        """Test that invalid SVTYPE raises AssertionError"""
        with pytest.raises(AssertionError, match="invalid svtype"):
            _infer_strands("INVALID", "N")

        with pytest.raises(AssertionError, match="invalid svtype"):
            _infer_strands("del", "N")  # lowercase

        with pytest.raises(AssertionError, match="invalid svtype"):
            _infer_strands("", "N")
