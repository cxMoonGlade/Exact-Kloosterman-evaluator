"""Direct cypari2 ordinary binary point counting; one fresh-process query."""
from time import perf_counter
ENTERED = perf_counter()
import argparse
import json
import resource
from direct_reference import MODULI


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("N", type=int)
    parser.add_argument("L", type=int)
    parser.add_argument("a", type=int)
    args = parser.parse_args()
    from cypari2 import Pari
    pari = Pari()
    start = perf_counter()
    coefficients = [(MODULI[args.N] >> i) & 1 for i in range(args.N + 1)]
    modulus = pari.Polrev(coefficients) * pari.Mod(1, 2)
    generator = pari.ffgen(modulus, "z")
    parameter = sum(generator ** i for i in range(args.N) if (args.a >> i) & 1)
    after_field = perf_counter()
    curve = pari.ellinit([1, 0, 0, parameter, 0], generator)
    after_curve = perf_counter()
    value = int(curve.ellcard()) - (1 << args.N) - 1
    if args.L == 3:
        value = value * value - (1 << args.N)
    finish = perf_counter()
    # Metadata assertions are outside the timed computation and preserve its output.
    actual_modulus = generator.mod()
    assert [int(actual_modulus.polcoef(i)) for i in range(args.N + 1)] == coefficients
    assert int(generator.p) == 2
    assert int(actual_modulus.poldegree()) == args.N
    print(json.dumps({"algorithm": "pari_direct", "N": args.N, "L": args.L,
                      "a": args.a, "modulus": MODULI[args.N], "value": value,
                      "compute_seconds": finish - start,
                      "entry_import_runtime_seconds": start - ENTERED,
                      "max_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                      "stages": {"field_parameter_seconds": after_field - start,
                                 "curve_seconds": after_curve - after_field,
                                 "point_count_recovery_seconds": finish - after_curve},
                      "pari_version": str(pari.version()), "status": "PASS"}))


if __name__ == "__main__":
    main()
