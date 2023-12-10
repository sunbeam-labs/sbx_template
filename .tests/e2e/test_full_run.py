import csv
import os
import pytest
import shutil
import subprocess as sp
import tempfile
from pathlib import Path


@pytest.fixture
def setup():
    temp_dir = Path(tempfile.mkdtemp())

    reads_fp = Path(".tests/data/reads/").resolve()

    project_dir = temp_dir / "project/"

    sp.check_output(["sunbeam", "init", "--data_fp", reads_fp, project_dir])

    config_fp = project_dir / "sunbeam_config.yml"

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
            "all_template",
            "--directory",
            temp_dir,
        ]
    )

    output_fp = project_dir / "sunbeam_output"
    benchmarks_fp = project_dir / "stats/"

    yield output_fp, benchmarks_fp


def test_full_run(run_sunbeam):
    output_fp, benchmarks_fp = run_sunbeam

    big_file_fp = output_fp / "qc/mush/big_file.txt"

    # Check output
    assert big_file_fp.exists(), f"{big_file_fp} does not exist"


def test_benchmarks(run_sunbeam):
    output_fp, benchmarks_fp = run_sunbeam

    filename = os.listdir(benchmarks_fp)[0]
    with open(os.path.join(benchmarks_fp, filename)) as f:
        rd = csv.DictReader(f, delimiter="\t")
        print(rd.fieldnames)
        for r in rd:
            assert (
                float(r["cpu_time"]) < 1
            ), f"cpu_time for {r['rule']} is higher than 1: {r['cpu_time']}"
