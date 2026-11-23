#!/usr/bin/env python3
"""Build independent Hydra variants, run them, and make a LaTeX supplement.

Standard library only. Run --help or see profiling/README.md.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import struct
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
VARIANTS = {
    "sw-wave": "kmake_wave_swrt_v1.json",
    "sw-mega": "kmake_mega_swrt_v1.json",
    "hw-wave": "kmake_wave.json",
    "hw-mega": "kmake_mega.json",
    "hw-mega-p": "kmake_mega_persistent.json",
}
METRIC = re.compile(r"register|\bregs?\b|\bgrfs?\b|gpr|spill|fill count|scratch|(?:shared|workgroup).*memory|lds|occupancy|"
                    r"subgroups per|\bwarps?\b|resident|threads per", re.I)
ARC_RESIDENCY = "Theoretical resident subgroups per Xe-core"
ARC_OCCUPANCY = "Theoretical Xe-core occupancy (%)"
PRESSURE = "Max live registers"
SUBGROUP_SIZE = "Subgroup size (threads)"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def revision(directory):
    return subprocess.check_output(["git", "-C", str(directory), "rev-parse", "HEAD"], text=True).strip()


def execute(command, cwd, log, env=None, timeout=None):
    command = [str(x) for x in command]
    print("Running:", subprocess.list2cmdline(command), flush=True)
    Path(log).parent.mkdir(parents=True, exist_ok=True)
    with Path(log).open("w", encoding="utf-8") as output:
        output.write(subprocess.list2cmdline(command) + "\n")
        output.flush()
        result = subprocess.run(command, cwd=cwd, env=env, stdout=output,
                                stderr=subprocess.STDOUT, timeout=timeout)
    if result.returncode:
        raise RuntimeError(f"exit {result.returncode}; see {log}")
    return Path(log).read_text(encoding="utf-8", errors="replace")


def build(args):
    slicer = args.slicer.resolve()
    executable = args.slicer_binary or slicer / "cmake-build-release/kslicer"
    if not executable.is_file():
        raise ValueError(f"kslicer not found: {executable}")
    for variant in args.variants:
        if variant.startswith("sw-") and not (ROOT / "external/HeavyRT/core/BVH2FatRT.cpp").is_file():
            raise ValueError("SW variants require HeavyRT in external/HeavyRT; see profiling/README.md")
    if any(v.startswith("sw-") for v in args.variants):
        heavy = ROOT / "external/HeavyRT"
        patch = ROOT / "profiling/heavyrt-compat.patch"
        installed = subprocess.run(["git", "-C", str(heavy), "apply", "--reverse", "--check", str(patch)],
                                   capture_output=True).returncode == 0
        if not installed:
            execute(["git", "-C", heavy, "apply", "--check", patch], ROOT,
                    args.build_root.resolve() / "heavyrt-check.log")
            execute(["git", "-C", heavy, "apply", patch], ROOT,
                    args.build_root.resolve() / "heavyrt-compat.log")
    for variant in args.variants:
        destination = args.build_root.resolve() / variant
        if (destination / "build.json").exists() and not args.rebuild:
            print(f"Already built: {variant} (use --rebuild to replace)")
            continue
        destination.mkdir(parents=True, exist_ok=True)
        # Invalidate the old manifest before changing any artifacts.
        (destination / "build.json").unlink(missing_ok=True)
        generated = ROOT / "shaders_generated"
        if generated.exists():
            shutil.rmtree(generated)
        config = ROOT / VARIANTS[variant]
        generation = execute([executable, config, "-selfdir", slicer], ROOT,
                             destination / "generate.log")
        if re.search(r"\berror:", generation, re.I):
            raise RuntimeError(f"translator reported errors; see {destination / 'generate.log'}")
        execute(["bash", "build_slang.sh"], generated, destination / "shaders.log",
                env=dict(os.environ, KSLICER_SHADER_JOBS=str(args.jobs)))
        source_shaders = list(generated.glob("*.slang"))
        if not source_shaders or any(not p.with_suffix(".comp.spv").is_file() for p in source_shaders):
            raise RuntimeError(f"missing SPIR-V after shader compilation for {variant}")
        cmake_dir = args.build_root.resolve() / "cmake"
        execute(["cmake", "-S", ROOT, "-B", cmake_dir, "-DCMAKE_BUILD_TYPE=Release",
                 "-DUSE_VULKAN=ON", "-DCAM_PLUGIN=OFF", "-DUSE_ENZYME=OFF",
                 "-DUSE_HEAVYRT=" + ("ON" if variant.startswith("sw-") else "OFF")],
                ROOT, destination / "configure.log")
        execute(["cmake", "--build", cmake_dir, "--config", "Release", "-j", args.jobs],
                ROOT, destination / "build.log")
        binaries = [ROOT / "bin-release/hydra", ROOT / "bin-release/Release/hydra.exe"]
        binary = next((p for p in binaries if p.is_file()), None)
        if binary is None:
            raise RuntimeError("CMake did not produce hydra")
        shutil.copy2(binary, destination / binary.name)
        shader_destination = destination / "shaders_generated"
        if shader_destination.exists():
            shutil.rmtree(shader_destination)
        shutil.copytree(generated, shader_destination)
        artifacts = [destination / binary.name, *sorted(shader_destination.glob("*.spv"))]
        manifest = {
            "variant": variant, "config": config.name, "config_sha256": digest(config),
            "source_revision": revision(ROOT), "slicer_revision": revision(slicer),
            "source_dirty": bool(subprocess.check_output(
                ["git", "-C", str(ROOT), "status", "--porcelain"], text=True).strip()),
            "heavyrt_patch_sha256": digest(ROOT / "profiling/heavyrt-compat.patch")
                if variant.startswith("sw-") else None,
            "built_at": datetime.now(timezone.utc).isoformat(),
            "binary": binary.name,
            "sha256": {str(p.relative_to(destination)): digest(p) for p in artifacts},
        }
        if variant.startswith("sw-"):
            manifest["heavyrt_revision"] = revision(ROOT / "external/HeavyRT")
        write_json(destination / "build.json", manifest)
        print(f"Built {variant}: {destination}", flush=True)


def parse_scene(value):
    name, separator, filename = value.partition("=")
    if not separator:
        filename = name
        name = Path(filename).parent.name
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", name):
        raise ValueError("scene labels must use letters, numbers, '.', '_' or '-'")
    path = Path(filename).resolve()
    if not path.is_file():
        raise ValueError(f"scene does not exist: {path}")
    return name, path


def arc_a380_anv(profile):
    device = profile["device"]
    return (device.get("vendor_id") == 0x8086 and device.get("device_id") == 0x56A5
            and device["driver_name"] == "Intel open-source Mesa driver")


def fixed_local_size(path):
    data = path.read_bytes()
    if len(data) < 20 or len(data) % 4:
        return None
    words = struct.unpack(f"<{len(data) // 4}I", data)
    if words[0] != 0x07230203:
        return None
    entries = []
    sizes = {}
    offset = 5
    while offset < len(words):
        count, opcode = words[offset] >> 16, words[offset] & 0xffff
        if not count or offset + count > len(words):
            return None
        if opcode == 15 and count >= 4:  # OpEntryPoint
            entries.append(words[offset + 1:offset + 3])
        if opcode == 16 and count == 6 and words[offset + 2] == 17:  # OpExecutionMode LocalSize
            sizes[words[offset + 1]] = words[offset + 3:offset + 6]
        offset += count
    # Hydra uses single-entry compute modules with a literal workgroup size.
    if len(entries) == 1 and entries[0][0] == 5:  # GLCompute
        size = sizes.get(entries[0][1])
        if size and all(size):
            return size
    return None


def run(args):
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    scenes = [parse_scene(s) for s in args.scene]
    if len({s[0] for s in scenes}) != len(scenes):
        raise ValueError("scene labels must be unique")
    if args.binary and (len(args.variants) != 1 or not args.resource_dir):
        raise ValueError("--binary requires exactly one --variants value and --resource-dir")
    failures = []
    for variant in args.variants:
        directory = args.build_root.resolve() / variant
        if args.binary:
            binary = args.binary.resolve()
            resources = args.resource_dir.resolve()
            manifest = {"variant": variant, "binary": str(binary), "sha256": digest(binary),
                        "external_build": True}
        else:
            manifest = read_json(directory / "build.json")
            if manifest["variant"] != variant:
                raise ValueError(f"wrong variant in {directory / 'build.json'}")
            for filename, expected in manifest["sha256"].items():
                if digest(directory / filename) != expected:
                    raise ValueError(f"artifact changed: {directory / filename}; rebuild it")
            binary = directory / manifest["binary"]
            resources = directory
        for scene, path in scenes:
            for repeat in range(args.repeats):
                stem = f"{variant}_{scene}_gpu{args.gpu_id}_{args.integrator}_{repeat}"
                raw = output / (stem + ".raw.json")
                record = output / (stem + ".json")
                # A failed rerun must not reuse an earlier successful record.
                raw.unlink(missing_ok=True)
                record.unlink(missing_ok=True)
                command = [binary, "--gpu", "-gpu_id", args.gpu_id, "-in", path,
                           "-integrator", args.integrator, "-width", args.width,
                           "-height", args.height, "-spp", args.spp, "-resource_dir", resources,
                           "-profile-json", raw, "-out", output / (stem + ".bmp")]
                if args.scene_dir:
                    command += ["-scn_dir", args.scene_dir.resolve()]
                if args.integrator == "rtao":
                    command += ["-ao-radius", args.ao_radius]
                try:
                    log = execute(command, ROOT, output / (stem + ".log"), timeout=args.timeout)
                    profile = read_json(raw)
                    if not profile["kernels"]:
                        raise RuntimeError("renderer recorded no compute kernels")
                    if not profile["supported"]:
                        raise RuntimeError("GPU/driver does not support pipeline executable statistics; "
                                           "raw metadata saved, use a vendor profiler on this machine")
                    if not any(e["statistics"] for k in profile["kernels"] for e in k["executables"]):
                        raise RuntimeError("driver returned no executable statistics")
                    if arc_a380_anv(profile):
                        for kernel in profile["kernels"]:
                            shader = Path(kernel["shader"])
                            if not shader.is_absolute():
                                shader = ROOT / shader
                            if shader.is_file() and (size := fixed_local_size(shader)):
                                kernel["local_size"] = size
                    if variant == "hw-mega-p":
                        widths = {e["subgroup_size"] for k in profile["kernels"] for e in k["executables"]}
                        if any(w == 0 or args.spp % w or (args.width * args.height) % w for w in widths):
                            raise RuntimeError("persistent rays require spp and pixel count divisible by subgroup size; "
                                               f"driver reported {sorted(widths)}; rerun with 1024x1024, 256 spp")
                    timings = re.findall(r"(\w+Block)\(exec\)\s*=\s*([\d.eE+-]+)\s*ms", log)
                    write_json(record, {
                        "schema_version": 1, "variant": variant, "scene": scene,
                        "scene_path": str(path), "scene_xml_sha256": digest(path),
                        "integrator": args.integrator, "width": args.width, "height": args.height,
                        "spp": args.spp, "repeat": repeat, "gpu_id": args.gpu_id,
                        "ao_radius": args.ao_radius if args.integrator == "rtao" else None,
                        "host": platform.node(), "platform": platform.platform(),
                        "gpu_label": args.gpu_label,
                        "captured_at": datetime.now(timezone.utc).isoformat(),
                        "command": [str(c) for c in command], "build": manifest,
                        "submit_wait_ms": {name: float(value) for name, value in timings},
                        "profile": profile,
                    })
                    print(f"Saved {record}", flush=True)
                except (OSError, RuntimeError, ValueError, subprocess.TimeoutExpired) as error:
                    failures.append(f"{stem}: {error}")
                    print(f"FAILED {failures[-1]}", file=sys.stderr, flush=True)
    records = sorted(p for p in output.glob("*.json") if not p.name.endswith(".raw.json"))
    if records:
        report_files(records, output / "profile.tex")
    if failures:
        raise RuntimeError("\n".join(failures))


def latex(value):
    escapes = {"\\": r"\textbackslash{}", "_": r"\_", "%": r"\%", "&": r"\&",
               "#": r"\#", "$": r"\$", "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}",
               "^": r"\textasciicircum{}"}
    return "".join(escapes.get(c, c) for c in str(value))


def number(value):
    if isinstance(value, bool):
        return str(value).lower()
    return f"{value:.6g}" if isinstance(value, float) else str(value)


def selected_statistics(executable):
    return [s for s in executable["statistics"] if METRIC.search(s["name"])]


def arc_residency_statistics(kernel, executable):
    size = kernel.get("local_size")
    width = executable["subgroup_size"]
    shared = next((s["value"] for s in executable["statistics"]
                   if s["name"] == "Workgroup Memory Size"), None)
    if (not size or not all(type(axis) is int and axis > 0 for axis in size)
            or width not in (8, 16, 32) or type(shared) is not int or not 0 <= shared <= 65536):
        return []
    # Intel's Xe-HPG occupancy model: 16 vector engines x 8 slots, 128 KiB SLM
    # per Xe-core, with power-of-two SLM allocation starting at 1 KiB.
    slots = 128
    subgroups = (math.prod(size) + width - 1) // width
    groups = slots // subgroups
    if shared:
        allocated = max(1024, 1 << (shared - 1).bit_length())
        groups = min(groups, 131072 // allocated)
    resident = groups * subgroups
    return [{"name": ARC_RESIDENCY, "value": resident,
             "description": "Modeled upper bound from SIMD width, SPIR-V workgroup size and "
                            "ANV shared memory, assuming enough workgroups to fill the Xe-core."},
            {"name": ARC_OCCUPANCY, "value": 100 * resident / slots,
             "description": "Modeled resident subgroups divided by 128 Xe-core hardware thread slots."}]


def gpu_key(record):
    d = record["profile"]["device"]
    return tuple(str(x) for x in (record.get("gpu_label", ""), d["name"], d.get("vendor_id", ""),
                 d.get("device_id", ""), d["driver_name"], d["driver_info"],
                 d.get("driver_version", ""), record.get("host", ""), record.get("gpu_id", 0)))


def summary_tables(records):
    lines = [r"\noindent\textbf{Register pressure summary}\par",
             "Repeated runs use the maximum for each cell. "
             "A dash indicates an unavailable count.\\\\"]
    devices = sorted({gpu_key(r) for r in records})
    labels = {}
    short_labels = {"RX 6800 XT": "6800XT", "Arc A380": "A380"}
    for key in devices:
        label = key[0] or re.sub(r"^(AMD Radeon |NVIDIA GeForce |NVIDIA )", "", key[1].split(" (")[0])
        labels[key] = short_labels.get(label, label)
    # Distinguish different driver/host experiments using the same GPU model.
    for key in devices:
        if sum(labels[other] == labels[key] for other in devices) > 1:
            label = labels[key]
            for i, other in enumerate(devices, 1):
                if labels[other] == label:
                    labels[other] += f" [GPU {i}]"
    for key in devices:
        lines.append(r"\noindent " + latex(labels[key]) + " identifies " + latex(
            f"{key[1]} ({key[4]} {key[5]}, driver version {key[6]}, host {key[7]}, device {key[8]})") + r"\par")
    groups = {}
    for r in records:
        key = (r["scene"], r.get("scene_xml_sha256", ""), r["width"], r["height"], r["spp"])
        groups.setdefault(key, []).append(r)
    # Compare each AO radius with the matching PT experiment without mixing AO settings.
    groups = {key + (radius,): [r for r in runs if r["integrator"] != "rtao" or r.get("ao_radius") == radius]
              for key, runs in groups.items()
              for radius in ({r.get("ao_radius") for r in runs if r["integrator"] == "rtao"} or {None})}

    def column_key(r):
        return r["integrator"], gpu_key(r), r.get("ao_radius")

    def column_order(column):
        workload = column[0]
        return (0 if workload == "mispt" else 1 if workload == "rtao" else 2, labels[column[1]], str(column))

    def header(column):
        workload = {"mispt": "PT", "rtao": "AO"}.get(column[0], column[0])
        return latex(f"{labels[column[1]]}/{workload}")

    variant_order = {variant: i for i, variant in enumerate(VARIANTS)}
    for group_id, (group, runs) in enumerate(sorted(groups.items(), key=lambda item: str(item[0]))):
        for target in ("maximum", "path-traversal"):
            cells = {}
            column_metrics = {}
            for r in runs:
                column = column_key(r)
                metric = PRESSURE if r["profile"]["device"]["driver_name"] == "Intel open-source Mesa driver" else "VGPRs"
                for k in r["profile"]["kernels"]:
                    name = Path(k["shader"]).name.removesuffix(".comp.spv")
                    if target == "path-traversal" and name not in (
                            "kernel_RayTrace2", "PathTraceMega", "NaivePathTraceMega"):
                        continue
                    row = (r["variant"], "" if target == "maximum" else name)
                    for e in k["executables"]:
                        for stat in e["statistics"]:
                            if stat["name"] == metric and type(stat["value"]) is int:
                                key = row, column
                                cells[key] = max(cells.get(key, stat["value"]), stat["value"])
                                column_metrics[column] = metric
            if not cells:
                continue
            rows = sorted({key[0] for key in cells}, key=lambda row: (variant_order[row[0]], row[1]))
            columns = sorted(column_metrics, key=column_order)
            caption = ("Peak register pressure across rendering kernels" if target == "maximum"
                       else "Register pressure for PT kernels that perform ray traversal")
            caption += f", {group[0]}, {group[2]} x {group[3]}, {group[4]} spp"
            radii = sorted({c[2] for c in columns if c[0] == "rtao" and c[2] is not None})
            if radii:
                caption += ", AO radius " + ", ".join(number(radius) for radius in radii)
            metrics = set(column_metrics.values())
            caption += "."
            if "VGPRs" in metrics:
                caption += " AMD RADV columns show allocated VGPRs per subgroup."
            if PRESSURE in metrics:
                caption += " Intel ANV columns show peak live GRFs before allocation and spilling (Max live registers)."
            if len(metrics) > 1:
                caption += " Each GPU uses its native counter and register units."
            target_id = f"{group_id}-{target.replace('_', '-')}"
            row_headers = ["Variant"] if target == "maximum" else ["Variant", "Kernel"]
            heading = r"\hline " + " & ".join(row_headers + [header(c) for c in columns]) + r" \\\hline"
            lines += [r"\begin{longtable}{" + "l" * len(row_headers) + "r" * len(columns) + "}",
                      r"\caption{" + latex(caption) + r"}\label{tab:hydra-summary-" + target_id + r"}\\",
                      heading, r"\endfirsthead", heading, r"\endhead"]
            for row in rows:
                values = [cells.get((row, c)) for c in columns]
                labels_for_row = [latex(value) for value in row if value]
                numbers = ["--" if value is None else latex(number(value)) for value in values]
                lines.append(" & ".join(labels_for_row + numbers) + r" \\")
            lines += [r"\hline", r"\end{longtable}"]
    lines.append(r"In the PT traversal summary, kernel\_RayTrace2 covers the standalone wavefront traversal stage. "
                 "Megakernel statistics cover the full fused shader, including traversal and shading. "
                 "RTAO appears in the rendering-kernel summary and the detailed tables.\\\\")
    return lines


def report_files(files, destination):
    records = []
    expanded = []
    for path in files:
        path = Path(path)
        expanded.extend(sorted(p for p in path.rglob("*.json") if not p.name.endswith(".raw.json"))
                        if path.is_dir() else [path])
    for path in expanded:
        value = read_json(path)
        if "profile" not in value or "variant" not in value:
            raise ValueError(f"not a run record: {path} (pass *.json excluding *.raw.json)")
        records.append(value)
    if not records:
        raise ValueError("no run records")
    lines = [r"% Generated by testing/profile_gpu.py. Requires \usepackage{longtable}.",
             r"\begingroup\scriptsize",
             r"Driver statistics come from the Vulkan extension \texttt{VK\_KHR\_pipeline\_executable\_properties}. "
             "The tables describe compute kernels recorded during rendering. "
             "Scene construction and pixel packing are excluded. Register counts and spills retain the "
             "driver's native units; they must not be compared across vendors as equivalent register files. "
             "Subgroups per SIMD is a static residency limit, not measured achieved occupancy. "
             "A maximum row is the independent maximum of one metric over the rendering kernels; "
             "it does not describe a single combined kernel.\\\\"]
    if any(arc_a380_anv(r["profile"]) for r in records):
        lines.append(r"On Arc A380, Max live registers is ANV's peak pre-allocation register pressure. "
                     r"The theoretical Xe-core occupancy ceiling uses Intel's workgroup packing model "
                     r"with 128 hardware thread slots and 128 KiB of shared local memory per Xe-core. "
                     r"The calculation assumes enough workgroups and uses the compiled SIMD width, "
                     r"SPIR-V workgroup size and ANV's workgroup shared memory size. "
                     r"Actual hardware occupancy requires runtime counters.\\\\")
    lines += summary_tables(records)
    for index, record in enumerate(records):
        p = record["profile"]
        device = p["device"]
        title = (f"{device['name']}; {device['driver_name']} {device['driver_info']}; "
                 f"{record['scene']}; {record['integrator']}; {record['variant']}; "
                 f"{record['width']} x {record['height']}; {record['spp']} spp; run {record['repeat'] + 1}")
        lines += [r"\begin{longtable}{p{0.30\linewidth}p{0.43\linewidth}r}",
                  r"\caption{" + latex(title) + r"}\label{tab:hydra-profile-" + str(index) + r"}\\",
                  r"\hline Kernel / executable & Metric (native units) & Value \\\hline",
                  r"\endfirsthead", r"\hline Kernel / executable & Metric (native units) & Value \\\hline",
                  r"\endhead"]
        maxima = {}
        definitions = {}
        for kernel in p["kernels"]:
            for executable in kernel["executables"]:
                kernel_name = Path(kernel["shader"]).name.removesuffix(".comp.spv")
                label = kernel_name + " / " + executable["name"]
                lines.append(latex(label) + " & " + SUBGROUP_SIZE + " & " + str(executable["subgroup_size"]) + r" \\")
                stats = selected_statistics(executable)
                if arc_a380_anv(p):
                    if kernel.get("local_size"):
                        lines.append(" & Workgroup size (invocations) & " +
                                     str(math.prod(kernel["local_size"])) + r" \\")
                    stats += arc_residency_statistics(kernel, executable)
                for stat in stats:
                    name, value = stat["name"], stat["value"]
                    lines.append(" & " + latex(name) + " & " + latex(number(value)) + r" \\")
                    definitions[name] = stat["description"]
                    if isinstance(value, (int, float)) and not isinstance(value, bool):
                        maxima[name] = max(maxima.get(name, value), value)
                lines.append(r"\hline")
        for name, value in sorted(maxima.items()):
            lines.append(r"\textbf{Maximum over kernels} & " + latex(name) + " & " + latex(number(value)) + r" \\")
        lines += [r"\hline", r"\end{longtable}"]
        for name, description in sorted(definitions.items()):
            lines.append(r"\noindent\textit{" + latex(name) + "}: " + latex(description) + r"\par")
        if not maxima:
            lines.append(r"\noindent No recognized register or residency statistics were returned.\par")
    lines.append(r"\endgroup")
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"LaTeX: {destination}")


def positive(value):
    value = int(value)
    if value <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    b = sub.add_parser("build", help="generate, compile and archive independent variants (Linux)")
    b.add_argument("--slicer", type=Path, required=True, help="kernel_slicer source directory")
    b.add_argument("--slicer-binary", type=Path)
    b.add_argument("--jobs", type=positive, default=min(os.cpu_count() or 4, 8))
    b.add_argument("--rebuild", action="store_true")
    r = sub.add_parser("run", help="run built variants and write JSON plus profile.tex")
    r.add_argument("--scene", action="append", required=True, help="Label=/path/to/scene.xml; repeatable")
    r.add_argument("--scene-dir", type=Path)
    r.add_argument("--gpu-id", type=int, default=0)
    r.add_argument("--gpu-label", default="", help="optional label for the GPU column in combined reports")
    r.add_argument("--integrator", choices=["mispt", "naivept", "shadowpt", "rtao"], default="mispt")
    r.add_argument("--ao-radius", type=float, default=2.0, help="RTAO radius in scene units")
    r.add_argument("--width", type=positive, default=1024)
    r.add_argument("--height", type=positive, default=1024)
    r.add_argument("--spp", type=positive, default=256)
    r.add_argument("--repeats", type=positive, default=1)
    r.add_argument("--timeout", type=positive, default=1800)
    r.add_argument("--output", type=Path, default=ROOT / "profiling_results/runs")
    r.add_argument("--binary", type=Path, help="profile an existing binary built with this profiler")
    r.add_argument("--resource-dir", type=Path, help="directory containing shaders_generated for --binary")
    for p in [b, r]:
        p.add_argument("--variants", nargs="+", choices=list(VARIANTS), default=list(VARIANTS))
        p.add_argument("--build-root", type=Path, default=ROOT / "profiling_results/builds")
    t = sub.add_parser("report", help="regenerate LaTeX from one or more run records")
    t.add_argument("records", type=Path, nargs="+", help="run JSON files or run directories (raw JSON is excluded)")
    t.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.action == "build":
            build(args)
        elif args.action == "run":
            if args.gpu_id < 0:
                raise ValueError("--gpu-id must be nonnegative")
            run(args)
        else:
            report_files(args.records, args.output)
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        parser.exit(1, f"ERROR: {error}\n")


if __name__ == "__main__":
    main()
