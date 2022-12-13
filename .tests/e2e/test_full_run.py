import os
import shutil
import subprocess as sp
import tempfile
import unittest


class FullRunTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp_dir = tempfile.mkdtemp()

        self.reads_fp = ".tests/data/reads/"

        self.project_dir = os.path.join(self.temp_dir, "project/")

        sp.check_output(
            ["sunbeam", "init", "--data_fp", self.reads_fp, self.project_dir]
        )

        self.config_fp = os.path.join(self.project_dir, "sunbeam_config.yml")

        config_str = f"sbx_template: {{example_rule_options: '--number'}}"

        sp.check_output(
            [
                "sunbeam",
                "config",
                "modify",
                "-i",
                "-s",
                f"{config_str}",
                f"{self.config_fp}",
            ]
        )

        self.output_fp = os.path.join(self.project_dir, "sunbeam_output")

        self.big_file_fp = os.path.join(self.output_fp, "qc/mush/big_file.txt")

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_full_run(self):
        # Run the test job.
        sp.check_output(
            [
                "sunbeam",
                "run",
                "--profile",
                self.project_dir,
                "all_template",
                "--directory",
                self.temp_dir,
            ]
        )

        # Check output
        self.assertTrue(os.path.exists(self.big_file_fp))
