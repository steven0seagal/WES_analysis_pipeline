#!/usr/bin/env python3

"""
Pipeline Performance Comparison Script
Compares GATK and DeepVariant pipelines on various metrics
"""

import os
import sys
import time
import psutil
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

class PipelineBenchmark:
    def __init__(self, sample_name, data_dir, results_dir):
        self.sample_name = sample_name
        self.data_dir = Path(data_dir)
        self.results_dir = Path(results_dir)
        self.metrics = {}

    def measure_resource_usage(self, func, *args, **kwargs):
        """Measure CPU, memory, and time usage of a function"""
        process = psutil.Process()
        start_time = time.time()
        start_cpu = process.cpu_times()
        start_memory = process.memory_info()

        result = func(*args, **kwargs)

        end_time = time.time()
        end_cpu = process.cpu_times()
        end_memory = process.memory_info()

        metrics = {
            'runtime_seconds': end_time - start_time,
            'cpu_user_seconds': end_cpu.user - start_cpu.user,
            'cpu_system_seconds': end_cpu.system - start_cpu.system,
            'memory_peak_mb': (end_memory.rss - start_memory.rss) / (1024 * 1024),
            'memory_percent': process.memory_percent()
        }

        return result, metrics

    def benchmark_gatk_pipeline(self):
        """Benchmark GATK pipeline"""
        print("Benchmarking GATK pipeline...")

        cmd = f"bash scripts/gatk_pipeline.sh -s {self.sample_name} -d {self.data_dir} -o {self.results_dir}/gatk"
        result, metrics = self.measure_resource_usage(os.system, cmd)

        self.metrics['gatk'] = metrics
        return metrics

    def benchmark_deepvariant_pipeline(self):
        """Benchmark DeepVariant pipeline"""
        print("Benchmarking DeepVariant pipeline...")

        cmd = f"bash scripts/deepvariant_pipeline.sh -s {self.sample_name} -d {self.data_dir} -o {self.results_dir}/deepvariant"
        result, metrics = self.measure_resource_usage(os.system, cmd)

        self.metrics['deepvariant'] = metrics
        return metrics

    def analyze_outputs(self):
        """Analyze pipeline outputs for quality metrics"""
        analysis = {}

        for pipeline in ['gatk', 'deepvariant']:
            pipeline_dir = self.results_dir / pipeline / self.sample_name

            if pipeline_dir.exists():
                # Count variants
                vcf_file = pipeline_dir / 'gvcf' / f'{self.sample_name}.g.vcf.gz'
                if vcf_file.exists():
                    # Count variants (simplified)
                    variant_count = 0
                    try:
                        import subprocess
                        result = subprocess.run(['bcftools', 'view', '-H', str(vcf_file)],
                                              capture_output=True, text=True)
                        variant_count = len(result.stdout.strip().split('\n')) if result.stdout.strip() else 0
                    except:
                        variant_count = 0

                    analysis[f'{pipeline}_variants'] = variant_count

                # Check file sizes
                bam_file = pipeline_dir / 'bam' / f'{self.sample_name}.analysis_ready.bam'
                if bam_file.exists():
                    analysis[f'{pipeline}_bam_size_mb'] = bam_file.stat().st_size / (1024 * 1024)

        return analysis

    def generate_report(self):
        """Generate performance comparison report"""
        print("\n" + "="*50)
        print("PIPELINE PERFORMANCE COMPARISON REPORT")
        print("="*50)

        if not self.metrics:
            print("No benchmark data available. Run benchmarks first.")
            return

        # Runtime comparison
        print("\nRUNTIME COMPARISON:")
        for pipeline, metrics in self.metrics.items():
            print(".2f")

        # Memory usage comparison
        print("\nMEMORY USAGE COMPARISON:")
        for pipeline, metrics in self.metrics.items():
            print(".2f")

        # CPU usage comparison
        print("\nCPU USAGE COMPARISON:")
        for pipeline, metrics in self.metrics.items():
            print(".2f")

        # Generate plots if matplotlib available
        try:
            self.generate_plots()
        except ImportError:
            print("\nMatplotlib not available for plotting")

    def generate_plots(self):
        """Generate comparison plots"""
        if not self.metrics:
            return

        # Prepare data
        pipelines = list(self.metrics.keys())
        runtimes = [self.metrics[p]['runtime_seconds'] for p in pipelines]
        memory_usage = [self.metrics[p]['memory_peak_mb'] for p in pipelines]

        # Create plots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        # Runtime plot
        ax1.bar(pipelines, runtimes, color=['blue', 'green'])
        ax1.set_title('Pipeline Runtime Comparison')
        ax1.set_ylabel('Time (seconds)')
        ax1.set_xlabel('Pipeline')

        # Memory plot
        ax2.bar(pipelines, memory_usage, color=['blue', 'green'])
        ax2.set_title('Pipeline Memory Usage Comparison')
        ax2.set_ylabel('Peak Memory (MB)')
        ax2.set_xlabel('Pipeline')

        plt.tight_layout()
        plt.savefig(self.results_dir / 'benchmark_comparison.png', dpi=300, bbox_inches='tight')
        print(f"\nComparison plot saved to: {self.results_dir / 'benchmark_comparison.png'}")

def main():
    parser = argparse.ArgumentParser(description='Compare WES pipeline performance')
    parser.add_argument('-s', '--sample', required=True, help='Sample name')
    parser.add_argument('-d', '--data-dir', default='data', help='Input data directory')
    parser.add_argument('-o', '--output-dir', default='benchmark_results', help='Output directory')
    parser.add_argument('--gatk-only', action='store_true', help='Benchmark only GATK')
    parser.add_argument('--deepvariant-only', action='store_true', help='Benchmark only DeepVariant')

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)

    # Initialize benchmark
    benchmark = PipelineBenchmark(args.sample, args.data_dir, output_dir)

    # Run benchmarks
    if not args.deepvariant_only:
        benchmark.benchmark_gatk_pipeline()

    if not args.gatk_only:
        benchmark.benchmark_deepvariant_pipeline()

    # Analyze outputs
    analysis = benchmark.analyze_outputs()

    # Generate report
    benchmark.generate_report()

    # Save detailed results
    results_file = output_dir / 'benchmark_results.json'
    import json
    with open(results_file, 'w') as f:
        json.dump({
            'metrics': benchmark.metrics,
            'analysis': analysis,
            'timestamp': time.time()
        }, f, indent=2)

    print(f"\nDetailed results saved to: {results_file}")

if __name__ == '__main__':
    main()