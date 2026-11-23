import copy
import tempfile
import unittest
from pathlib import Path

from profile_gpu import report_files, write_json, summary_tables


class ReportTests(unittest.TestCase):
    def record(self):
        def kernel(name, registers, waves):
            return {"shader": name + ".comp.spv", "executables": [{
                "name": "compute", "subgroup_size": 64, "statistics": [
                    {"name": "VGPRs", "value": registers, "description": "Native register count"},
                    {"name": "Subgroups per SIMD", "value": waves, "description": "Residency limit"},
                    {"name": "Hash", "value": 123, "description": "Not a pressure metric"},
                ]}]}
        return {"variant": "hw-wave", "scene": "A&B_1%", "integrator": "mispt",
                "width": 1024, "height": 1024, "spp": 256, "repeat": 0,
                "profile": {"supported": True, "device": {
                    "name": "Test GPU", "driver_name": "driver", "driver_info": "1.0"},
                    "kernels": [kernel("First_kernel", 128, 4), kernel("Second", 32, 12)]}}

    def test_independent_maxima_and_escaping(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            write_json(directory / "run.json", self.record())
            write_json(directory / "run.raw.json", {"not": "a run record"})
            report_files([directory], directory / "result.tex")
            result = (directory / "result.tex").read_text()
            self.assertIn(r"A\&B\_1\%", result)
            self.assertIn(r"First\_kernel", result)
            self.assertIn(r"\textbf{Maximum over kernels} & VGPRs & 128", result)
            self.assertIn(r"\textbf{Maximum over kernels} & Subgroups per SIMD & 12", result)
            self.assertNotIn("Hash &", result)
            self.assertIn("not measured achieved occupancy", result)

    def test_missing_metric_is_not_zero(self):
        record = copy.deepcopy(self.record())
        for kernel in record["profile"]["kernels"]:
            kernel["executables"][0]["statistics"] = []
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            write_json(directory / "run.json", record)
            report_files([directory / "run.json"], directory / "result.tex")
            result = (directory / "result.tex").read_text()
            self.assertIn("No recognized register or residency statistics", result)
            self.assertNotIn(r"\textbf{Maximum over kernels}", result)

    def test_raw_input_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            write_json(directory / "raw.json", {"supported": True})
            with self.assertRaisesRegex(ValueError, "not a run record"):
                report_files([directory / "raw.json"], directory / "result.tex")

    def test_gpu_columns_kernel_filter_repeats_and_absent_kernels(self):
        a = self.record()
        a["profile"]["device"]["name"] = "GPU A"
        a["profile"]["kernels"][0]["shader"] = "kernel_RayTrace2.comp.spv"
        b = copy.deepcopy(a)
        b["profile"]["device"]["name"] = "GPU B"
        b["profile"]["kernels"][0]["executables"][0]["statistics"][0]["value"] = 64
        # The other kernel affects the all-kernel maximum, never RayTrace2's count.
        b["profile"]["kernels"][1]["executables"][0]["statistics"][0]["value"] = 999
        repeat = copy.deepcopy(a)
        repeat["profile"]["kernels"][0]["executables"][0]["statistics"][0]["value"] = 160
        mega = copy.deepcopy(a)
        mega["variant"] = "hw-mega"
        mega["profile"]["kernels"][0]["shader"] = "PathTraceMega.comp.spv"
        output = "\n".join(summary_tables([a, b, repeat, mega]))
        self.assertIn("GPU 1 & GPU 2", output)
        self.assertIn("hw-wave & VGPRs & 160 & 999", output)
        self.assertIn("hw-wave & VGPRs & 160 & 64", output)
        self.assertIn("hw-mega & VGPRs & N/A & --", output)

    def test_spills_and_estimates_are_not_register_allocation(self):
        record = self.record()
        record["profile"]["kernels"][0]["executables"][0]["statistics"] += [
            {"name": "Pre-Sched VGPRs", "value": 9999},
            {"name": "Spilled VGPRs", "value": 8888},
        ]
        output = "\n".join(summary_tables([record]))
        self.assertNotIn("9999", output)
        self.assertNotIn("8888", output)
        self.assertIn("hw-wave & VGPRs & 128", output)


if __name__ == "__main__":
    unittest.main()
