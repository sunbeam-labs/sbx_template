import csv
import os
import pytest
import shutil
import subprocess as sp
import tempfile


@pytest.fixture
def setup():
    temp_dir = tempfile.mkdtemp()

    reads_fp = os.path.abspath(".tests/data/reads/")

    project_dir = os.path.join(temp_dir, "project/")

    sp.check_output(["sunbeam", "init", "--data_fp", reads_fp, project_dir])

    config_fp = os.path.join(project_dir, "sunbeam_config.yml")

    config_str = f"sbx_template: {{example_rule_options: '--number'}}"

    sp.check_output(
        [
            "sunbeam",
            "config",
            "modify",
            "-i",
            "-s",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    yield temp_dir, project_dir

    shutil.rmtree(temp_dir)


@pytest.fixture
def run_sunbeam(setup):
    temp_dir, project_dir = setup

    # Run the test job.
    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            project_dir,
            "all_classify",
            "--directory",
            temp_dir,
        ]
    )

    output_fp = os.path.join(project_dir, "sunbeam_output")

    big_file_fp = os.path.join(output_fp, "qc/mush/big_file.txt")

    benchmarks_fp = os.path.join(project_dir, "stats/")

    yield big_file_fp, benchmarks_fp


def test_full_run(run_sunbeam):
    big_file_fp, benchmarks_fp = run_sunbeam

    # Check output
    assert os.path.exists(big_file_fp)


def test_benchmarks(run_sunbeam):
    big_file_fp, benchmarks_fp = run_sunbeam

    filename = os.listdir(benchmarks_fp)[0]
    with open(os.path.join(benchmarks_fp, filename)) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            assert r == 0
            # assert float(r["cpu_time"]) < 0.0
