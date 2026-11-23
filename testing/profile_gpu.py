#!/usr/bin/env python3
"""Build independent Hydra variants, run them, and make a LaTeX supplement.

Standard library only. Run --help or see profiling/README.md.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
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
METRIC = re.compile(r"register|\bregs?\b|\bgrfs?\b|gpr|spill|scratch|shared.*memory|lds|occupancy|"
                    r"subgroups per|\bwarps?\b|resident|threads per", re.I)


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


def register_statistics(executable):
    # Allocation counts only: pre-scheduling estimates and spills have different meanings.
    return [s for s in executable["statistics"]
            if re.search(r"register|\bregs?\b|\bgrfs?\b|\b[sv]gprs?\b", s["name"], re.I)
            and not re.search(r"pre.?sched|spill", s["name"], re.I)
            and isinstance(s["value"], (int, float)) and not isinstance(s["value"], bool)]


def gpu_key(record):
    d = record["profile"]["device"]
    return tuple(str(x) for x in (record.get("gpu_label", ""), d["name"], d.get("vendor_id", ""),
                 d.get("device_id", ""), d["driver_name"], d["driver_info"],
                 d.get("driver_version", ""), record.get("host", ""), record.get("gpu_id", 0)))


def summary_tables(records):
    lines = [r"\noindent\textbf{Register allocation summary}\par",
             "Columns identify the workload and GPU. Counts retain the driver's native "
             "register-file names and units; SGPRs, VGPRs and GRFs are not added together. "
             "Repeated runs are reduced by taking the maximum for each cell. "
             "A dash means a measurement or metric was not reported. "
             "Subgroup widths and native metric definitions are listed in the detailed tables.\\\\"]
    devices = sorted({gpu_key(r) for r in records})
    labels = {}
    for i, key in enumerate(devices, 1):
        label = key[0] or re.sub(r"^(AMD Radeon |NVIDIA GeForce |NVIDIA )", "", key[1].split(" (")[0])
        labels[key] = label
        lines.append(r"\noindent " + latex(label) + ": " + latex(
            f"{key[1]}; {key[4]} {key[5]}; driver version {key[6]}; host {key[7]}; device {key[8]}") + r"\par")
    # Distinguish different driver/host experiments using the same GPU model.
    for key in devices:
        if sum(labels[other] == labels[key] for other in devices) > 1:
            label = labels[key]
            for i, other in enumerate(devices, 1):
                if labels[other] == label:
                    labels[other] += f" [GPU {i}]"
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
        return (0 if workload == "mispt" else 1 if workload == "rtao" else 2, str(column))

    def header(column):
        workload = {"mispt": "PT", "rtao": "AO"}.get(column[0], column[0])
        return r"\shortstack{(" + latex(workload) + r",\\" + latex(labels[column[1]]) + ")}"

    for group_id, (group, runs) in enumerate(sorted(groups.items(), key=lambda item: str(item[0]))):
        for target in ("maximum", "kernel_RayTrace2"):
            cells = {}
            for r in runs:
                column = column_key(r)
                for k in r["profile"]["kernels"]:
                    name = Path(k["shader"]).name.removesuffix(".comp.spv")
                    if target != "maximum" and name != target:
                        continue
                    for e in k["executables"]:
                        for stat in register_statistics(e):
                            key = r["variant"], stat["name"], column
                            cells[key] = max(cells.get(key, stat["value"]), stat["value"])
            if not cells:
                continue
            columns = sorted({key[2] for key in cells}, key=column_order)
            metrics = sorted({key[1] for key in cells})
            caption = ("Maximum register allocation over rendering kernels" if target == "maximum"
                       else "Register allocation of standalone kernel_RayTrace2")
            caption += f"; {group[0]}, {group[2]} x {group[3]}, {group[4]} spp"
            radii = sorted({c[2] for c in columns if c[0] == "rtao" and c[2] is not None})
            if radii:
                caption += "; AO radius " + ", ".join(number(radius) for radius in radii)
            for offset in range(0, len(columns), 4):
                panel = columns[offset:offset + 4]
                heading = r"\hline Variant & Register file / metric & " + " & ".join(header(c) for c in panel) + r" \\\hline"
                lines += [r"\begin{longtable}{ll" + "r" * len(panel) + "}",
                          r"\caption{" + latex(caption) + r"}\label{tab:hydra-register-" +
                          f"{group_id}-{target.replace('_', '-')}-{offset}" + r"}\\",
                          heading, r"\endfirsthead", heading, r"\endhead"]
                for variant in VARIANTS:
                    rows = []
                    for metric in metrics:
                        values = [cells.get((variant, metric, c)) for c in panel]
                        # In particular, omit fused megakernels from RayTrace2's table.
                        if all(value is None for value in values):
                            continue
                        rows.append(latex(variant) + " & " + latex(metric) + " & " +
                                    " & ".join("--" if value is None else latex(number(value)) for value in values) + r" \\")
                    if rows:
                        lines += rows + [r"\hline"]
                lines.append(r"\end{longtable}")
    lines.append(r"Standalone kernel\_RayTrace2 is available in PT wavefront variants. "
                 "Megakernels inline this stage; their register counts cannot be attributed to this function alone. "
                 "RTAO uses different kernel stages.\\\\")
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
             r"The numbers were obtained with the Vulkan extension \texttt{VK\_KHR\_pipeline\_executable\_properties}. "
             "These are compile-time driver statistics for compute kernels recorded during rendering. "
             "Scene construction and pixel packing are excluded. Register counts and spills retain the "
             "driver's native units; they must not be compared across vendors as equivalent register files. "
             "Subgroups per SIMD is a static residency limit, not measured achieved occupancy. "
             "A maximum row is the independent maximum of one metric over the rendering kernels; "
             "it does not describe a single combined kernel.\\\\"]
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
                lines.append(latex(label) + " & Subgroup size (threads) & " + str(executable["subgroup_size"]) + r" \\")
                for stat in selected_statistics(executable):
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
