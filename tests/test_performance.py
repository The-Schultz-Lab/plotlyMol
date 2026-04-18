"""
Performance Testing Suite for plotlyMol

This module provides quantitative performance testing to identify bottlenecks
and measure rendering/parsing performance under various conditions.

Usage:
    python tests/test_performance.py

    Or with pytest:
    pytest tests/test_performance.py -v --benchmark
"""

import sys
import time
import tracemalloc
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
import psutil

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from rdkit import Chem
from rdkit.Chem import AllChem

from plotlymol3d import create_vibration_animation, draw_3D_rep, parse_vibrations

# ============================================================================
# BENCHMARKING UTILITIES
# ============================================================================


def measure_performance(func: Callable, *args, **kwargs) -> tuple:
    """
    Measure execution time and memory usage of a function.

    Args:
        func: Function to benchmark
        *args, **kwargs: Arguments to pass to function

    Returns:
        (result, execution_time_ms, memory_increase_mb, peak_memory_mb)
    """
    # Start memory tracking
    tracemalloc.start()
    process = psutil.Process()
    mem_before = process.memory_info().rss / 1024 / 1024  # MB

    # Measure execution time
    start_time = time.perf_counter()
    result = func(*args, **kwargs)
    end_time = time.perf_counter()

    # Get memory usage
    mem_after = process.memory_info().rss / 1024 / 1024  # MB
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    execution_time_ms = (end_time - start_time) * 1000
    memory_increase_mb = mem_after - mem_before
    peak_memory_mb = peak / 1024 / 1024

    return result, execution_time_ms, memory_increase_mb, peak_memory_mb


def benchmark_multiple_runs(func: Callable, n_runs: int = 5, *args, **kwargs) -> dict:
    """
    Run a benchmark multiple times and compute statistics.

    Args:
        func: Function to benchmark
        n_runs: Number of runs
        *args, **kwargs: Arguments to pass to function

    Returns:
        Dictionary with timing and memory statistics
    """
    times = []
    memories = []
    peak_memories = []

    for _i in range(n_runs):
        _, exec_time, mem, peak = measure_performance(func, *args, **kwargs)
        times.append(exec_time)
        memories.append(mem)
        peak_memories.append(peak)

    return {
        "mean_time_ms": np.mean(times),
        "std_time_ms": np.std(times),
        "min_time_ms": np.min(times),
        "max_time_ms": np.max(times),
        "median_time_ms": np.median(times),
        "mean_memory_mb": np.mean(memories),
        "std_memory_mb": np.std(memories),
        "peak_memory_mb": np.mean(peak_memories),
    }


# ============================================================================
# RENDERING PERFORMANCE TESTS
# ============================================================================


def benchmark_rendering_by_size():
    """
    Benchmark rendering performance vs molecule size.

    Returns:
        DataFrame with results
    """
    print("\n" + "=" * 70)
    print("BENCHMARK: Rendering Performance vs Molecule Size")
    print("=" * 70)

    test_molecules = [
        ("Water", "O", 3),
        ("Ethanol", "CCO", 9),
        ("Benzene", "c1ccccc1", 12),
        ("Glucose", "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O", 24),
        ("Cholesterol", "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C", 74),
    ]

    results = []

    for name, smiles, approx_atoms in test_molecules:
        print(f"\n{name} ({approx_atoms} atoms):")

        for mode in ["ball+stick", "stick", "vdw"]:
            stats = benchmark_multiple_runs(
                draw_3D_rep, n_runs=3, smiles=smiles, mode=mode, resolution=32
            )

            results.append(
                {
                    "Molecule": name,
                    "Atoms": approx_atoms,
                    "Mode": mode,
                    "Mean Time (ms)": stats["mean_time_ms"],
                    "Std Time (ms)": stats["std_time_ms"],
                    "Peak Memory (MB)": stats["peak_memory_mb"],
                }
            )

            print(
                f"  {mode:12s}: {stats['mean_time_ms']:6.1f} ± {stats['std_time_ms']:4.1f} ms"
            )

    return pd.DataFrame(results)


def benchmark_resolution_impact():
    """
    Benchmark how sphere resolution affects performance.

    Returns:
        DataFrame with results
    """
    print("\n" + "=" * 70)
    print("BENCHMARK: Resolution Impact on Performance")
    print("=" * 70)

    test_smiles = "c1ccccc1"  # Benzene
    resolutions = [8, 16, 24, 32, 48, 64]

    results = []

    print("\nBenzene (12 atoms):")

    for resolution in resolutions:
        stats = benchmark_multiple_runs(
            draw_3D_rep,
            n_runs=3,
            smiles=test_smiles,
            mode="ball+stick",
            resolution=resolution,
        )

        results.append(
            {
                "Resolution": resolution,
                "Mean Time (ms)": stats["mean_time_ms"],
                "Std Time (ms)": stats["std_time_ms"],
                "Peak Memory (MB)": stats["peak_memory_mb"],
            }
        )

        print(
            f"  Resolution {resolution:2d}: {stats['mean_time_ms']:6.1f} ± {stats['std_time_ms']:4.1f} ms"
        )

    # Calculate performance ratios
    baseline_idx = [r["Resolution"] for r in results].index(32)
    baseline_time = results[baseline_idx]["Mean Time (ms)"]

    print("\nPerformance relative to default (resolution=32):")
    for res_data in results:
        ratio = res_data["Mean Time (ms)"] / baseline_time
        change = (ratio - 1) * 100
        print(
            f"  Resolution {res_data['Resolution']:2d}: {ratio:.2f}x ({change:+.0f}%)"
        )

    return pd.DataFrame(results)


# ============================================================================
# VIBRATION PERFORMANCE TESTS
# ============================================================================


def benchmark_vibration_parsing():
    """
    Benchmark vibration file parsing performance.

    Returns:
        DataFrame with results (or None if no fixtures available)
    """
    print("\n" + "=" * 70)
    print("BENCHMARK: Vibration File Parsing")
    print("=" * 70)

    # Look for test fixtures
    fixtures_dir = Path(__file__).parent / "fixtures"
    vib_files = {
        "Gaussian": fixtures_dir / "water_gaussian.log",
        "ORCA": fixtures_dir / "water_orca.out",
        "Molden": fixtures_dir / "water.molden",
    }

    results = []

    for program, filepath in vib_files.items():
        if not filepath.exists():
            print(f"\n{program}: File not found, skipping")
            continue

        print(f"\n{program} ({filepath.name}):")

        # Get file size
        file_size_kb = filepath.stat().st_size / 1024

        # Benchmark parsing
        stats = benchmark_multiple_runs(
            parse_vibrations, n_runs=5, filepath=str(filepath)
        )

        # Get mode count
        vib_data, _, _, _ = measure_performance(parse_vibrations, str(filepath))
        n_modes = len(vib_data.modes)
        n_atoms = len(vib_data.atomic_numbers)

        results.append(
            {
                "Program": program,
                "File Size (KB)": file_size_kb,
                "Atoms": n_atoms,
                "Modes": n_modes,
                "Parse Time (ms)": stats["mean_time_ms"],
                "Std (ms)": stats["std_time_ms"],
            }
        )

        print(
            f"  Parse time: {stats['mean_time_ms']:.1f} ± {stats['std_time_ms']:.1f} ms"
        )
        print(f"  {n_atoms} atoms, {n_modes} modes, {file_size_kb:.1f} KB")

    return pd.DataFrame(results) if results else None


def benchmark_vibration_visualization():
    """
    Benchmark vibration visualization modes.

    Returns:
        DataFrame with results (or None if no fixtures available)
    """
    print("\n" + "=" * 70)
    print("BENCHMARK: Vibration Visualization Modes")
    print("=" * 70)

    # Use water fixture
    fixtures_dir = Path(__file__).parent / "fixtures"
    vib_file = fixtures_dir / "water_gaussian.log"

    if not vib_file.exists():
        print("\nWater fixture not found, skipping")
        return None

    smiles = "O"
    vib_data = parse_vibrations(str(vib_file))

    results = []

    # Test static arrows
    print("\nStatic Arrows:")
    stats = benchmark_multiple_runs(
        draw_3D_rep,
        n_runs=3,
        smiles=smiles,
        vibration_file=str(vib_file),
        vibration_mode=1,
        vibration_display="arrows",
    )
    results.append(
        {
            "Mode": "Static Arrows",
            "Mean Time (ms)": stats["mean_time_ms"],
            "Std (ms)": stats["std_time_ms"],
            "Peak Memory (MB)": stats["peak_memory_mb"],
        }
    )
    print(f"  {stats['mean_time_ms']:.1f} ± {stats['std_time_ms']:.1f} ms")

    # Test heatmap
    print("\nHeatmap:")
    stats = benchmark_multiple_runs(
        draw_3D_rep,
        n_runs=3,
        smiles=smiles,
        vibration_file=str(vib_file),
        vibration_mode=1,
        vibration_display="heatmap",
    )
    results.append(
        {
            "Mode": "Heatmap",
            "Mean Time (ms)": stats["mean_time_ms"],
            "Std (ms)": stats["std_time_ms"],
            "Peak Memory (MB)": stats["peak_memory_mb"],
        }
    )
    print(f"  {stats['mean_time_ms']:.1f} ± {stats['std_time_ms']:.1f} ms")

    # Test animation
    print("\nAnimation (20 frames):")
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)

    stats = benchmark_multiple_runs(
        create_vibration_animation,
        n_runs=3,
        vib_data=vib_data,
        mode_number=1,
        mol=mol,
        amplitude=0.5,
        n_frames=20,
        mode="ball+stick",
    )
    results.append(
        {
            "Mode": "Animation (20 frames)",
            "Mean Time (ms)": stats["mean_time_ms"],
            "Std (ms)": stats["std_time_ms"],
            "Peak Memory (MB)": stats["peak_memory_mb"],
        }
    )
    print(f"  {stats['mean_time_ms']:.1f} ± {stats['std_time_ms']:.1f} ms")
    print(f"  ({stats['mean_time_ms']/20:.1f} ms per frame)")

    return pd.DataFrame(results)


def benchmark_animation_frame_count():
    """
    Benchmark animation performance vs frame count.

    Returns:
        DataFrame with results (or None if no fixtures available)
    """
    print("\n" + "=" * 70)
    print("BENCHMARK: Animation Frame Count Impact")
    print("=" * 70)

    # Use water fixture
    fixtures_dir = Path(__file__).parent / "fixtures"
    vib_file = fixtures_dir / "water_gaussian.log"

    if not vib_file.exists():
        print("\nWater fixture not found, skipping")
        return None

    smiles = "O"
    vib_data = parse_vibrations(str(vib_file))

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)

    frame_counts = [5, 10, 20, 30, 50]
    results = []

    for n_frames in frame_counts:
        stats = benchmark_multiple_runs(
            create_vibration_animation,
            n_runs=3,
            vib_data=vib_data,
            mode_number=1,
            mol=mol,
            amplitude=0.5,
            n_frames=n_frames,
            mode="ball+stick",
            resolution=32,
        )

        results.append(
            {
                "Frames": n_frames,
                "Total Time (ms)": stats["mean_time_ms"],
                "Std (ms)": stats["std_time_ms"],
                "Time per Frame (ms)": stats["mean_time_ms"] / n_frames,
                "Peak Memory (MB)": stats["peak_memory_mb"],
            }
        )

        print(
            f"  {n_frames:2d} frames: {stats['mean_time_ms']:7.1f} ms "
            f"({stats['mean_time_ms']/n_frames:5.1f} ms/frame)"
        )

    return pd.DataFrame(results)


# ============================================================================
# PERFORMANCE ANALYSIS
# ============================================================================


def analyze_performance_results(results: dict[str, pd.DataFrame]):
    """
    Analyze all benchmark results and provide recommendations.

    Args:
        results: Dictionary of benchmark name -> DataFrame
    """
    print("\n" + "=" * 70)
    print("PERFORMANCE ANALYSIS & RECOMMENDATIONS")
    print("=" * 70)

    # Rendering analysis
    if "rendering" in results:
        df = results["rendering"]
        print("\n1. RENDERING PERFORMANCE:")

        # By size
        small_mol = df[df["Atoms"] <= 20]["Mean Time (ms)"].mean()
        large_mol = df[df["Atoms"] > 50]["Mean Time (ms)"].mean()

        if not np.isnan(small_mol) and not np.isnan(large_mol):
            print(f"   • Small molecules (<20 atoms): {small_mol:.0f} ms average")
            print(f"   • Large molecules (>50 atoms): {large_mol:.0f} ms average")
            print(
                f"   • Scaling factor: {large_mol/small_mol:.1f}x slower for large molecules"
            )

        # By mode
        print("\n   Mode comparison:")
        for mode in df["Mode"].unique():
            mode_data = df[df["Mode"] == mode]
            avg = mode_data["Mean Time (ms)"].mean()
            print(f"     - {mode:12s}: {avg:.0f} ms average")

    # Resolution analysis
    if "resolution" in results:
        df = results["resolution"]
        print("\n2. RESOLUTION SETTINGS:")

        res_16 = df[df["Resolution"] == 16]["Mean Time (ms)"].values
        res_32 = df[df["Resolution"] == 32]["Mean Time (ms)"].values
        res_64 = df[df["Resolution"] == 64]["Mean Time (ms)"].values

        if len(res_16) > 0 and len(res_32) > 0 and len(res_64) > 0:
            print(f"   • Resolution 16 (Performance): {res_16[0]:.0f} ms")
            print(f"   • Resolution 32 (Balanced):    {res_32[0]:.0f} ms - DEFAULT")
            print(f"   • Resolution 64 (Quality):     {res_64[0]:.0f} ms")
            print(f"   • Speedup 64→32: {(1 - res_32[0]/res_64[0])*100:.0f}% faster")
            print(f"   • Speedup 32→16: {(1 - res_16[0]/res_32[0])*100:.0f}% faster")

    # Vibration parsing analysis
    if "parsing" in results and results["parsing"] is not None:
        df = results["parsing"]
        print("\n3. VIBRATION PARSING:")

        for _, row in df.iterrows():
            print(
                f"   • {row['Program']:8s}: {row['Parse Time (ms)']:5.1f} ms "
                f"({row['Modes']} modes, {row['File Size (KB)']:.0f} KB)"
            )

    # Vibration visualization analysis
    if "vib_viz" in results and results["vib_viz"] is not None:
        df = results["vib_viz"]
        print("\n4. VIBRATION VISUALIZATION:")

        for _, row in df.iterrows():
            print(f"   • {row['Mode']:25s}: {row['Mean Time (ms)']:6.1f} ms")

    # Animation analysis
    if "animation" in results and results["animation"] is not None:
        df = results["animation"]
        print("\n5. ANIMATION PERFORMANCE:")

        print("   Frame count recommendations:")
        for _, row in df.iterrows():
            print(
                f"     - {row['Frames']:2.0f} frames: {row['Total Time (ms)']:6.0f} ms "
                f"({row['Time per Frame (ms)']:.1f} ms/frame)"
            )

    # General recommendations
    print("\n6. GENERAL RECOMMENDATIONS:")
    print("   For GUI/Dash Performance:")
    print("   • Use resolution=16 for interactive preview (2-3x faster)")
    print("   • Switch to resolution=32 for final figures")
    print("   • Cache parsed vibration files using session state")
    print("   • Use 'stick' mode for molecules >50 atoms")
    print("   • Limit animations to 20-30 frames during development")
    print("   • Consider 'vdw' mode for very large systems (>100 atoms)")

    print("\n   Memory Management:")
    if "rendering" in results:
        avg_mem = results["rendering"]["Peak Memory (MB)"].mean()
        print(f"   • Average peak memory: {avg_mem:.1f} MB per render")
    print("   • Clear figure cache periodically in long sessions")
    print("   • Profile memory-intensive operations with tracemalloc")

    print("\n" + "=" * 70)


# ============================================================================
# MAIN EXECUTION
# ============================================================================


def run_all_benchmarks(save_results: bool = True):
    """
    Run all performance benchmarks and save results.

    Args:
        save_results: Whether to save results to CSV files

    Returns:
        Dictionary of DataFrames with results
    """
    results = {}

    # Run benchmarks
    results["rendering"] = benchmark_rendering_by_size()
    results["resolution"] = benchmark_resolution_impact()
    results["parsing"] = benchmark_vibration_parsing()
    results["vib_viz"] = benchmark_vibration_visualization()
    results["animation"] = benchmark_animation_frame_count()

    # Analyze results
    analyze_performance_results(results)

    # Save results
    if save_results:
        output_dir = Path(__file__).parent.parent / "benchmark_results"
        output_dir.mkdir(exist_ok=True)

        timestamp = time.strftime("%Y%m%d_%H%M%S")

        for name, df in results.items():
            if df is not None:
                output_file = output_dir / f"{name}_{timestamp}.csv"
                df.to_csv(output_file, index=False)
                print(f"\nSaved {name} results to: {output_file}")

    return results


if __name__ == "__main__":
    print("=" * 70)
    print("plotlyMol Performance Testing Suite")
    print("=" * 70)
    print("\nThis will run comprehensive performance benchmarks.")
    print("Estimated time: 2-5 minutes\n")

    results = run_all_benchmarks(save_results=True)

    print("\n" + "=" * 70)
    print("BENCHMARKING COMPLETE!")
    print("=" * 70)
    print("\nResults saved to: benchmark_results/")
    print("\nUse these results to optimize your GUI settings and identify bottlenecks.")
