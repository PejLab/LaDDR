from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from laddr.cli import CoverageConfig, get_sample_table
from laddr.coverage import base_covg_from_bigwigs, bin_covg_from_bigwigs, load_bins


class FakeBigWig:
    def __init__(self, path, data):
        self.path = path
        self.data = data[path]

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        pass

    def chroms(self):
        return {"chr1": 1000}

    def header(self):
        return {"sumData": self.data["sumData"]}

    def stats(self, chrom, start, end, type="mean", exact=True):
        return [self.data["stat"]]

    def values(self, chrom, start, end):
        return self.data["values"][start:end]


def fake_bigwig_module(data):
    return SimpleNamespace(open=lambda path: FakeBigWig(path, data))


def test_unstranded_manifest_parsing_is_unchanged(tmp_path):
    manifest = tmp_path / "coverage_manifest.tsv"
    manifest.write_text("dset1\tS1\tS1.bw\n")
    config = CoverageConfig(
        method="manifest",
        directory=Path("covg"),
        manifest=Path("coverage_manifest.tsv"),
        stranded=False,
    )

    table = get_sample_table(config, tmp_path)

    assert table.columns.tolist() == ["dataset", "sample", "path"]
    assert table.loc[0, "path"] == str((tmp_path / "covg" / "S1.bw").absolute())


def test_stranded_manifest_parsing_resolves_both_paths(tmp_path):
    manifest = tmp_path / "coverage_manifest.tsv"
    manifest.write_text("dset1\tS1\tS1.plus.bw\tS1.minus.bw\n")
    config = CoverageConfig(
        method="manifest",
        directory=Path("covg"),
        manifest=Path("coverage_manifest.tsv"),
        stranded=True,
    )

    table = get_sample_table(config, tmp_path)

    assert table.columns.tolist() == ["dataset", "sample", "plus_path", "minus_path"]
    assert table.loc[0, "plus_path"] == str((tmp_path / "covg" / "S1.plus.bw").absolute())
    assert table.loc[0, "minus_path"] == str((tmp_path / "covg" / "S1.minus.bw").absolute())


def test_stranded_directory_input_is_rejected(tmp_path):
    config = CoverageConfig(
        method="directory",
        directory=Path("covg"),
        manifest=Path("coverage_manifest.tsv"),
        stranded=True,
    )

    with pytest.raises(ValueError, match="manifest"):
        get_sample_table(config, tmp_path)


def test_load_bins_retains_bed_strand_column(tmp_path):
    bed = tmp_path / "batch_0.bed"
    bed.write_text("chr1\t10\t20\tGENE1_0_10\tadaptive\t+\n")

    bins = load_bins(bed)

    assert bins.loc[("GENE1", 5), "strand"] == "+"


def test_bin_covg_from_bigwigs_uses_gene_strand_and_combined_depth(monkeypatch):
    data = {
        "S1.plus.bw": {"sumData": 30.0, "stat": 30.0, "values": np.arange(10.0)},
        "S1.minus.bw": {"sumData": 10.0, "stat": 10.0, "values": np.arange(10.0)},
        "S2.plus.bw": {"sumData": 60.0, "stat": 60.0, "values": np.arange(10.0)},
        "S2.minus.bw": {"sumData": 20.0, "stat": 20.0, "values": np.arange(10.0)},
    }
    monkeypatch.setattr("laddr.coverage.pyBigWig", fake_bigwig_module(data))
    manifest = pd.DataFrame({
        "dataset": ["dset1", "dset1"],
        "sample": ["S1", "S2"],
        "plus_path": ["S1.plus.bw", "S2.plus.bw"],
        "minus_path": ["S1.minus.bw", "S2.minus.bw"],
    })
    bins = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "chrom_start": [0, 0],
            "chrom_end": [5, 5],
            "strand": ["+", "-"],
        },
        index=pd.MultiIndex.from_tuples([("gene_plus", 0), ("gene_minus", 0)], names=["gene_id", "pos"]),
    )

    covg = bin_covg_from_bigwigs(manifest, bins, median_coverage=60.0)

    assert covg.loc[("gene_plus", 0), "S1"] == 45.0
    assert covg.loc[("gene_minus", 0), "S1"] == 15.0
    assert covg.loc[("gene_plus", 0), "S2"] == 45.0
    assert covg.loc[("gene_minus", 0), "S2"] == 15.0


def test_base_covg_from_bigwigs_uses_requested_strand(monkeypatch):
    data = {
        "S1.plus.bw": {"sumData": 30.0, "stat": 0.0, "values": np.array([2.0, 4.0, 6.0])},
        "S1.minus.bw": {"sumData": 10.0, "stat": 0.0, "values": np.array([8.0, 10.0, 12.0])},
    }
    monkeypatch.setattr("laddr.coverage.pyBigWig", fake_bigwig_module(data))
    manifest = pd.DataFrame({
        "dataset": ["dset1"],
        "sample": ["S1"],
        "plus_path": ["S1.plus.bw"],
        "minus_path": ["S1.minus.bw"],
    })

    covg = base_covg_from_bigwigs(manifest, "chr1", 0, 3, median_coverage=20.0, strand="-")

    np.testing.assert_array_equal(covg[:, 0], np.array([4.0, 5.0, 6.0]))
