import copy
import tempfile
import unittest
from pathlib import Path

from profile_gpu import ARC_OCCUPANCY, ARC_RESIDENCY, report_files, write_json, summary_tables


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
            self.assertIn(r"VK\_KHR\_pipeline\_executable\_properties", result)

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

    def test_gpu_columns_kernel_filter_and_repeats(self):
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
        mega["profile"]["kernels"][1]["executables"][0]["statistics"][0]["value"] = 999
        output = "\n".join(summary_tables([a, b, repeat, mega]))
        self.assertIn(r"\shortstack{(PT,\\GPU A)\\VGPRs}", output)
        self.assertIn(r"\shortstack{(PT,\\GPU B)\\VGPRs}", output)
        self.assertIn("hw-wave & 160 & 999", output)
        ray_table = output.split("Reported vector registers for PT kernels")[1]
        self.assertIn(r"hw-wave & kernel\_RayTrace2 & 160 & 64", ray_table)
        self.assertIn("hw-mega & PathTraceMega & 128 & --", ray_table)
        self.assertNotIn("Subgroup size (threads)", ray_table)
        self.assertNotIn("999", ray_table)
        self.assertNotIn("N/A", output)

    def test_pt_and_ao_share_maximum_table(self):
        pt = self.record()
        ao = copy.deepcopy(pt)
        ao["integrator"] = "rtao"
        ao["ao_radius"] = 2.0
        ao["profile"]["kernels"][0]["executables"][0]["statistics"][0]["value"] = 48
        output = "\n".join(summary_tables([pt, ao]))
        self.assertEqual(output.count("Maximum reported vector registers over rendering kernels"), 1)
        self.assertIn(r"\shortstack{(PT,\\Test GPU)\\VGPRs}", output)
        self.assertIn(r"\shortstack{(AO,\\Test GPU)\\VGPRs}", output)
        self.assertIn("hw-wave & 128 & 48", output)

    def test_spills_and_estimates_are_not_register_allocation(self):
        record = self.record()
        record["profile"]["kernels"][0]["executables"][0]["statistics"] += [
            {"name": "Pre-Sched VGPRs", "value": 9999},
            {"name": "Spilled VGPRs", "value": 8888},
        ]
        output = "\n".join(summary_tables([record]))
        self.assertNotIn("9999", output)
        self.assertNotIn("8888", output)
        self.assertIn("hw-wave & 128", output)

    def test_summary_merges_vendor_metrics_in_each_table(self):
        amd = self.record()
        amd["profile"]["device"]["name"] = "AMD GPU"
        amd["profile"]["kernels"][0]["shader"] = "kernel_RayTrace2.comp.spv"
        amd["profile"]["kernels"][0]["executables"][0]["statistics"].append(
            {"name": "SGPRs", "value": 96})
        intel = copy.deepcopy(amd)
        intel["profile"]["device"]["name"] = "Intel GPU"
        intel["profile"]["device"]["driver_name"] = "Intel open-source Mesa driver"
        intel_kernel = intel["profile"]["kernels"][0]
        intel_kernel["executables"][0]["statistics"] = [{"name": "Max live registers", "value": 42}]
        intel["profile"]["kernels"] = [intel_kernel]
        output = "\n".join(summary_tables([amd, intel]))
        tables = output.split(r"\begin{longtable}")[1:]
        self.assertEqual(len(tables), 2)
        for table, kernel in zip(tables, ("", r"kernel\_RayTrace2 & ")):
            self.assertIn(r"\shortstack{(PT,\\AMD GPU)\\VGPRs} & \shortstack{(PT,\\Intel GPU)\\Max live registers}", table)
            self.assertNotIn("SGPRs", table)
            self.assertNotIn("Subgroups per SIMD", table)
            self.assertIn(f"hw-wave & {kernel}128 & 42", table)

    def test_arc_residency_uses_workgroup_and_rounded_slm(self):
        record = self.record()
        record["profile"]["device"].update(name="Arc A380", vendor_id=0x8086, device_id=0x56A5,
                                            driver_name="Intel open-source Mesa driver")
        kernel = record["profile"]["kernels"][0]
        kernel["local_size"] = [256, 1, 1]
        executable = kernel["executables"][0]
        executable["subgroup_size"] = 16
        executable["statistics"] = [
            {"name": "Max live registers", "description": "Register pressure", "value": 317},
            {"name": "Workgroup Memory Size", "description": "Shared memory in bytes", "value": 16385},
        ]
        record["profile"]["kernels"] = [kernel]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            write_json(path / "arc.json", record)
            report_files([path / "arc.json"], path / "arc.tex")
            output = (path / "arc.tex").read_text()
        self.assertIn(f" & {ARC_RESIDENCY} & 64", output)
        self.assertIn(r" & " + ARC_OCCUPANCY.replace("%", r"\%") + " & 50", output)
        self.assertIn(r"\shortstack{(PT,\\Arc A380)\\Max live registers}", output)
        del kernel["local_size"]
        self.assertNotIn(ARC_RESIDENCY, "\n".join(summary_tables([record])))


if __name__ == "__main__":
    unittest.main()
