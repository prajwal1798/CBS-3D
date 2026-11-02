"""Regression tests for serial/MPI VTU material array naming.

Run:
    python3 -m unittest discover -s tests -p test_audit_cht_temperature_vtu.py -v
No production solver is involved.
"""
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import audit_cht_temperature_vtu as audit  # noqa: E402


class AuditVtuReaderTests(unittest.TestCase):
    def fixture(self, cell_data):
        return f"""<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
<UnstructuredGrid>
<Piece NumberOfPoints="4" NumberOfCells="2">
<Points>
<DataArray type="Float64" NumberOfComponents="3" format="ascii">
0 0 0  1 0 0  0 1 0  0 0 1
</DataArray>
</Points>
<Cells>
<DataArray Name="connectivity" type="Int32" format="ascii">0 1 2 3  0 1 2</DataArray>
<DataArray Name="offsets" type="Int32" format="ascii">4 7</DataArray>
<DataArray Name="types" type="UInt8" format="ascii">10 5</DataArray>
</Cells>
<PointData Scalars="temperature">
<DataArray Name="temperature" type="Float64" format="ascii">298.15 298.15 298.15 298.15</DataArray>
</PointData>
<CellData>
{cell_data}
<DataArray Name="bc_id" type="Int32" format="ascii">0 532</DataArray>
</CellData>
</Piece>
</UnstructuredGrid>
</VTKFile>
"""

    def parse(self, xml):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "tiny.vtu"
            path.write_text(xml, encoding="utf-8")
            return audit.read_piece(path)

    def test_serial_writer_mat_id(self):
        xml = self.fixture(
            '<DataArray Name="mat_id" type="Int32" format="ascii">0 -1</DataArray>'
        )
        piece = self.parse(xml)
        self.assertEqual(piece.mat_id, [0, -1])
        self.assertEqual(piece.bc_id, [0, 532])

    def test_mpi_writer_material_id(self):
        xml = self.fixture(
            '<DataArray Name="material_id" type="Int32" format="ascii">0 -1</DataArray>'
        )
        piece = self.parse(xml)
        self.assertEqual(piece.mat_id, [0, -1])
        self.assertEqual(piece.bc_id, [0, 532])

    def test_both_names_agree(self):
        xml = self.fixture(
            '<DataArray Name="mat_id">0 -1</DataArray>\n'
            '<DataArray Name="material_id">0 -1</DataArray>'
        )
        self.assertEqual(self.parse(xml).mat_id, [0, -1])

    def test_conflicting_names_rejected(self):
        xml = self.fixture(
            '<DataArray Name="mat_id">0 -1</DataArray>\n'
            '<DataArray Name="material_id">1 -1</DataArray>'
        )
        with self.assertRaisesRegex(RuntimeError, "conflicting material_id and mat_id"):
            self.parse(xml)

    def test_absent_material_array_rejected(self):
        with self.assertRaisesRegex(RuntimeError, "missing CellData material field"):
            self.parse(self.fixture(""))


if __name__ == "__main__":
    unittest.main()
