"""Regression tests for the offline DTU/PBR builders."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


TOOLS = Path(__file__).resolve().parents[1] / "tools"


def _load_tool(name: str):
    spec = importlib.util.spec_from_file_location(name, TOOLS / f"{name}.py")
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load tool {name}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class DTUNameParserTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.builder = _load_tool("build_netmhc_pseudo_reference")

    def test_variable_width_compact_hla_fields_use_support_universe(self):
        known = {
            "HLA-DRB1*13:100",
            "HLA-DRB1*13:101",
            "HLA-DRB5*01:08N",
        }
        self.assertEqual(
            self.builder.parse_official_ii_name("DRB1_13100", known)[2],
            "HLA-DRB1*13:100",
        )
        self.assertEqual(
            self.builder.parse_official_ii_name("DRB1_13101", known)[2],
            "HLA-DRB1*13:101",
        )
        self.assertEqual(
            self.builder.parse_official_ii_name("DRB5_0108N", known)[2],
            "HLA-DRB5*01:08N",
        )


if __name__ == "__main__":
    unittest.main()
