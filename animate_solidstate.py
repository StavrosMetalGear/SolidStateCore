#!/usr/bin/env python3
"""
animate_solidstate.py  —  Visualization helper for SolidStateCore CSV exports.

Usage:
    python animate_solidstate.py --file <csv_file> [--type band|dos|magnon|hall]

Requires: pandas, numpy, matplotlib
"""

import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def plot_band_structure(df):
    """Plot band structure from CSV with columns: k, band_0, band_1, ..."""
    fig, ax = plt.subplots(figsize=(8, 6))
    k = df.iloc[:, 0]
    for col in df.columns[1:]:
        ax.plot(k, df[col], label=col)
    ax.set_xlabel("k")
    ax.set_ylabel("E")
    ax.set_title("Band Structure")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


def plot_dos(df):
    """Plot density of states from CSV."""
    fig, ax = plt.subplots(figsize=(8, 6))
    if "E" in df.columns and "DOS" in df.columns:
        ax.plot(df["E"], df["DOS"])
    else:
        ax.plot(df.iloc[:, 0], df.iloc[:, 1])
    ax.set_xlabel("Energy")
    ax.set_ylabel("DOS")
    ax.set_title("Density of States")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


def plot_magnetization(df):
    """Plot magnetization vs temperature."""
    fig, ax = plt.subplots(figsize=(8, 6))
    if "T" in df.columns and "M" in df.columns:
        ax.plot(df["T"], df["M"])
    else:
        ax.plot(df.iloc[:, 0], df.iloc[:, 1])
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("Magnetization")
    ax.set_title("Magnetization vs Temperature")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


def plot_hall(df):
    """Plot quantum Hall effect data."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    B = df.iloc[:, 0]
    if "R_xy" in df.columns:
        ax1.plot(B, df["R_xy"])
    else:
        ax1.plot(B, df.iloc[:, 3])
    ax1.set_ylabel("R_xy (Ohm)")
    ax1.set_title("Quantum Hall Effect")
    ax1.grid(True, alpha=0.3)

    if "R_xx" in df.columns:
        ax2.plot(B, df["R_xx"])
    else:
        ax2.plot(B, df.iloc[:, 4])
    ax2.set_xlabel("B (T)")
    ax2.set_ylabel("R_xx (Ohm)")
    ax2.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


def plot_generic(df):
    """Generic CSV plotter."""
    fig, ax = plt.subplots(figsize=(8, 6))
    x = df.iloc[:, 0]
    for col in df.columns[1:]:
        ax.plot(x, df[col], label=col)
    ax.set_xlabel(df.columns[0])
    ax.set_ylabel("Value")
    ax.set_title("SolidStateCore Data")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="SolidStateCore CSV Visualizer")
    parser.add_argument("--file", required=True, help="Path to CSV file")
    parser.add_argument("--type", default="generic",
                        choices=["band", "dos", "mag", "hall", "generic"],
                        help="Plot type")
    args = parser.parse_args()

    df = pd.read_csv(args.file)
    print(f"Loaded {args.file}: {len(df)} rows, columns: {list(df.columns)}")

    plotters = {
        "band": plot_band_structure,
        "dos": plot_dos,
        "mag": plot_magnetization,
        "hall": plot_hall,
        "generic": plot_generic,
    }
    plotters[args.type](df)


if __name__ == "__main__":
    main()
