"""Command-line interface for taster.

Exposes the ``taster`` command, which wraps :func:`taster.run_partition_workflow`
and accepts all workflow parameters as command-line arguments.
"""
import argparse
import sys
from . import run_partition_workflow, available_solvents


def main():
    parser = argparse.ArgumentParser(
        prog='taster',
        description='Martini CG free energy partition coefficient calculations.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument('--itp', required=True,
                        help='Path to the molecule ITP file.')
    parser.add_argument('--structure', required=True,
                        help='Path to the input CG structure file.')
    parser.add_argument('--solvents', nargs='+', default=None,
                        help=f'Solvents to use. Available: {sorted(available_solvents())}. '
                             f'Defaults to all available.')
    parser.add_argument('--reference', default='water',
                        help='Reference solvent for partition coefficient.')
    parser.add_argument('--reps', type=int, default=1,
                        help='Number of replicates.')
    parser.add_argument('--ncores', type=int, default=None,
                        help='Number of parallel processes. Defaults to auto-detect.')
    parser.add_argument('--output-dir', default='./Partitions',
                        help='Root directory for output.')
    parser.add_argument('--gmx', default='gmx',
                        help='GROMACS executable name or path.')
    parser.add_argument('--temperature', type=float, default=298.0,
                        help='Temperature in Kelvin.')
    parser.add_argument('--cutoff', type=int, default=5000,
                        help='Number of initial frames to discard as equilibration.')
    parser.add_argument('--estimator', default='MBAR', choices=['TI', 'MBAR'],
                        help='Free energy estimator.')
    parser.add_argument('--nsteps', type=int, default=1250000,
                        help='Number of steps for the FEP production run (default 1250000 = 25 ns at dt=0.02 ps).')

    args = parser.parse_args()

    try:
        df = run_partition_workflow(
            itp=args.itp,
            structure=args.structure,
            solvents=args.solvents,
            reference=args.reference,
            reps=args.reps,
            ncores=args.ncores,
            output_dir=args.output_dir,
            gmx=args.gmx,
            T=args.temperature,
            cutoff=args.cutoff,
            estimator=args.estimator,
            nsteps=args.nsteps,
        )
        print(df.to_string(index=False))

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()