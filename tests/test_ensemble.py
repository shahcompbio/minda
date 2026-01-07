import pytest
from minda.ensemble import _infer_strands, _get_strands_from_alt, _get_strands_from_info


class TestGetStrandsFromInfo:
    """Test STRANDS parsing from INFO field"""

    def test_strands_in_info(self):
        """Test parsing STRANDS from INFO field"""
        assert _get_strands_from_info("SVTYPE=DEL;STRANDS=+-;SVLEN=100") == "+-"
        assert _get_strands_from_info("SVTYPE=BND;STRANDS=++;END=500") == "++"
        assert _get_strands_from_info("STRANDS=--") == "--"
        assert _get_strands_from_info("SVLEN=200;STRANDS=-+;CHR2=chr2") == "-+"

    def test_strands_not_in_info(self):
        """Test when STRANDS is not in INFO field"""
        assert _get_strands_from_info("SVTYPE=DEL;SVLEN=100") is None
        assert _get_strands_from_info("SVTYPE=INS") is None
        assert _get_strands_from_info("") is None

    def test_strands_at_different_positions(self):
        """Test STRANDS at beginning, middle, and end of INFO"""
        assert _get_strands_from_info("STRANDS=+-;SVTYPE=DEL;SVLEN=100") == "+-"
        assert _get_strands_from_info("SVTYPE=DEL;STRANDS=++;SVLEN=100") == "++"
        assert _get_strands_from_info("SVTYPE=DEL;SVLEN=100;STRANDS=--") == "--"


class TestGetStrandsFromAlt:
    """Test strand inference from BND ALT field"""

    def test_forward_bracket(self):
        """Test parsing forward bracket notation: N[chr2:100[ gives +-"""
        assert _get_strands_from_alt("N[chr2:100[") == "+-"

    def test_reverse_bracket(self):
        """Test parsing reverse bracket notation: N]chr2:100] gives ++"""
        assert _get_strands_from_alt("N]chr2:100]") == "++"

    def test_bracket_at_start_forward(self):
        """Test parsing when forward bracket is at start: [chr2:100[N gives --"""
        assert _get_strands_from_alt("[chr2:100[N") == "--"

    def test_bracket_at_start_reverse(self):
        """Test parsing when reverse bracket is at start: ]chr2:100]N gives -+"""
        assert _get_strands_from_alt("]chr2:100]N") == "-+"

    def test_invalid_alt_returns_dots(self):
        """Test that invalid ALT strings return '..'"""
        assert _get_strands_from_alt("<INV>") == ".."
        assert _get_strands_from_alt("N") == ".."
        assert _get_strands_from_alt("") == ".."

    def test_mismatched_brackets_returns_dots(self):
        """Test that mismatched brackets return '..'"""
        assert _get_strands_from_alt("N[chr2:100]") == ".."
        assert _get_strands_from_alt("N]chr2:100[") == ".."


class TestInferStrands:
    """Test strand inference with fallback logic"""

    # DEL tests
    def test_del_simple_notation(self):
        """DEL with simple notation should be +-"""
        assert _infer_strands("DEL", "N", "SVTYPE=DEL;SVLEN=100") == "+-"

    def test_del_symbolic_alt(self):
        """DEL with symbolic ALT should be +-"""
        assert _infer_strands("DEL", "<DEL>", "SVTYPE=DEL;SVLEN=100") == "+-"

    # INS tests
    def test_ins_simple_notation(self):
        """INS with simple notation should be +-"""
        assert _infer_strands("INS", "N", "SVTYPE=INS;SVLEN=100") == "+-"

    def test_ins_symbolic_alt(self):
        """INS with symbolic ALT should be +-"""
        assert _infer_strands("INS", "<INS>", "SVTYPE=INS;SVLEN=100") == "+-"

    # DUP tests
    def test_dup_simple_notation(self):
        """DUP with simple notation should be -+"""
        assert _infer_strands("DUP", "<DUP>", "SVTYPE=DUP;SVLEN=100") == "-+"

    def test_dup_alt_contains_dup(self):
        """DUP in ALT string should be -+"""
        assert _infer_strands("DUP", "N<DUP>", "SVTYPE=DUP;SVLEN=100") == "-+"

    # BND tests - all four orientations inferred from ALT
    def test_bnd_forward_reverse(self):
        """BND: N[chr2:100[ should give +- (forward-reverse)"""
        assert _infer_strands("BND", "N[chr2:100[", "SVTYPE=BND;CHR2=chr2") == "+-"

    def test_bnd_forward_forward(self):
        """BND: N]chr2:100] should give ++ (forward-forward)"""
        assert _infer_strands("BND", "N]chr2:100]", "SVTYPE=BND;CHR2=chr2") == "++"

    def test_bnd_reverse_reverse(self):
        """BND: [chr2:100[N should give -- (reverse-reverse)"""
        assert _infer_strands("BND", "[chr2:100[N", "SVTYPE=BND;CHR2=chr2") == "--"

    def test_bnd_reverse_forward(self):
        """BND: ]chr2:100]N should give -+ (reverse-forward)"""
        assert _infer_strands("BND", "]chr2:100]N", "SVTYPE=BND;CHR2=chr2") == "-+"

    # INV tests - various BND notations
    def test_inv_reverse_reverse_bracket_before(self):
        """INV: [chr2:100[N should give -- (reverse-reverse)"""
        assert _infer_strands("INV", "[chr2:100[N", "SVTYPE=INV;SVLEN=100") == "--"

    def test_inv_forward_forward_bracket_after(self):
        """INV: N]chr2:100] should give ++ (forward-forward)"""
        assert _infer_strands("INV", "N]chr2:100]", "SVTYPE=INV;SVLEN=100") == "++"


class TestStrandInferenceFallback:
    """Test fallback behavior when INFO field has valid STRANDS"""

    def test_info_field_takes_precedence(self):
        """When INFO has valid STRANDS, use it for non-DEL/INS/DUP types"""
        # For BND/INV, INFO field takes precedence over ALT parsing
        # Even though BND ALT says +-, INFO says --
        assert _infer_strands("BND", "N[chr2:100[", "SVTYPE=BND;STRANDS=--") == "--"
        # Even though INV ALT would give ++, INFO says +-
        assert _infer_strands("INV", "N]chr2:100]", "SVTYPE=INV;STRANDS=+-") == "+-"

    def test_fallback_to_alt_when_info_invalid(self):
        """When INFO has invalid STRANDS, fall back to ALT parsing"""
        # INFO has invalid strand value, should parse ALT
        assert _infer_strands("BND", "N]chr2:100]", "SVTYPE=BND;STRANDS=invalid") == "++"
        assert _infer_strands("BND", "[chr2:100[N", "SVTYPE=BND;STRANDS=X") == "--"

    def test_fallback_to_alt_when_no_info_strands(self):
        """When INFO has no STRANDS, fall back to ALT parsing"""
        assert _infer_strands("BND", "N[chr2:100[", "SVTYPE=BND;SVLEN=100") == "+-"
        assert _infer_strands("BND", "]chr2:100]N", "SVTYPE=BND") == "-+"

    def test_returns_dots_when_unparseable(self):
        """When everything fails, return '..'"""
        # Unknown SVTYPE, no INFO strands, unparseable ALT
        assert _infer_strands("UNKNOWN", "<UNKNOWN>", "SVTYPE=UNKNOWN") == ".."
        assert _infer_strands("CTX", "N", "SVTYPE=CTX;SVLEN=100") == ".."
        # Invalid ALT and no INFO strands
        assert _infer_strands("BND", "malformed", "SVTYPE=BND") == ".."
