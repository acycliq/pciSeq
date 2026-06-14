import os
import sqlite3
import numpy as np
import pytest
from pciSeq.src.core.utils.io_utils import export_diagnostics

def test_export_diagnostics_includes_mrf(minimal_varbayes, tmp_path):
    """Test that diagnostics export includes the mrf column in cells table."""
    vb = minimal_varbayes
    vb.initialise_state()
    
    # Run one iteration to populate mrf
    vb.geneCount_upd()
    vb.gamma_upd()
    vb.cell_to_cellType()
    
    # Setup output directory
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    
    # Export diagnostics
    export_diagnostics(vb, str(output_dir))
    
    # Verify database exists
    db_path = output_dir / "diagnostics" / "diagnostics.db"
    assert db_path.exists()
    
    # Connect and check schema
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    
    # Check if mrf column exists in cells table
    cursor.execute("PRAGMA table_info(cells)")
    columns = [row[1] for row in cursor.fetchall()]
    assert "mrf" in columns, "mrf column missing from cells table"
    
    # Check if data is populated
    cursor.execute("SELECT mrf FROM cells LIMIT 1")
    row = cursor.fetchone()
    assert row is not None
    assert row[0] is not None
    
    # Verify the blob can be converted back to numpy and has correct size
    mrf_blob = row[0]
    mrf_array = np.frombuffer(mrf_blob, dtype=np.float32)
    assert mrf_array.shape == (vb.nK,), f"Expected shape ({vb.nK},), got {mrf_array.shape}"
    
    # Verify values match stashed mrf
    stashed_mrf = vb.cells.mrf[0].astype(np.float32)
    assert np.allclose(mrf_array, stashed_mrf)
    
    conn.close()
