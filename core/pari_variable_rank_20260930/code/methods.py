"""PARI exact cyclic polynomial powers; independent of the ellcard shortcut."""
from pathlib import Path
from time import perf_counter
import struct
import sys

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT.parent
MEMORY_LIMIT = 91 * 1024**3


def memory_precheck(N):
    q = 1 << N
    word = struct.calcsize('l')
    minimum = word * (4*q + 2)
    return dict(required_bytes_lower_bound=minimum, memory_limit_bytes=MEMORY_LIMIT,
                feasible=minimum <= MEMORY_LIMIT, word_bytes=word,
                coefficient_slots=4*q-4, header_words=6,
                scope='PARI GEN polynomials F, first unreduced square F^2 and T; excludes integer payloads and temporary buffers',
                measured=False, certificate='PARI_SOURCE_AUDIT.md')


def setup_pari():
    from cypari2 import Pari
    pari = Pari()
    # Dynamic stack: reserve an upper address range; physical pages grow on use.
    # RLIMIT_AS also counts cypari2 heap copies and interpreter allocations.
    pari.allocatemem(8*1024**2, 90*1024**3, silent=True)
    pari.setrand(20260930)
    pari('qcbcharacter(z,g,m)={my(v=vector(m),u=g^0);'
         'for(j=1,m,v[j]=1-2*lift(trace(u));u*=g);Polrev(v,\'x)}')
    return pari


def _inputs(N, L, f, params):
    if type(N) is not int or N < 2 or type(L) is not int or L not in (4,6,8,12,16):
        raise ValueError('This registered adapter requires L in {4,6,8,12,16} and N>=2')
    if type(f) is not int or f.bit_length() != N+1:
        raise ValueError('Supplied polynomial must have degree N')
    if not params or any(type(a) is not int or not 0 < a < 1 << N for a in params):
        raise ValueError('Nonzero field labels required')


def field_parameter(pari, z, a):
    # Horner reads exactly the supplied polynomial-basis label; no isomorphism
    # or precomputed logarithm map is provided to the algorithm.
    out = z*0
    for i in range(a.bit_length()-1, -1, -1):
        out = out*z + ((a >> i) & 1)
    return out


def run_pari(pari, N, L, f, params, instrument=False):
    start = perf_counter()
    _inputs(N,L,f,params)
    check = memory_precheck(N)
    if not check['feasible']:
        return dict(status='memory_precheck_blocked', N=N,L=L,f=f,params=params,
                    outputs=None, compute_seconds=None, cold_seconds=None,
                    preprocess_seconds=None, timing_eligible=False, memory_precheck=check)
    coefficients = [(f >> i) & 1 for i in range(N+1)]
    modulus = pari.Polrev(coefficients, 'z') * pari.Mod(1, 2)
    if not modulus.polisirreducible():
        raise ValueError('Supplied polynomial is reducible')
    z = pari.ffgen(modulus, 'z')
    if [int(z.minpoly().polcoef(i)) for i in range(N+1)] != coefficients:
        raise ArithmeticError('PARI did not retain the supplied field polynomial')
    stages = dict(field_validation_construction=perf_counter()-start)
    t = perf_counter()
    g = z.ffprimroot()
    stages['primitive_search'] = perf_counter()-t
    q = 1 << N
    m = q-1
    t = perf_counter()
    F = pari('qcbcharacter')(z,g,m)
    T = pari('x')**m - 1
    stages['character_and_modulus'] = perf_counter()-t
    t = perf_counter()
    C = pari.lift(pari.Mod(F,T)**L)
    del F,T
    stages['cyclic_power'] = perf_counter()-t
    prep = perf_counter()-start
    outputs, queries = [], []
    for a in params:
        t = perf_counter()
        parameter = field_parameter(pari,z,a)
        tp = perf_counter()
        index = int(parameter.fflog(g,m))
        ti = perf_counter()
        value = int(C.polcoef(index))
        tf = perf_counter()
        outputs.append(value)
        queries.append(dict(a=a,seconds=tf-t,stages=dict(parameter=tp-t,
                            discrete_log=ti-tp, coefficient_extraction=tf-ti),
                            stage_gap_seconds=0.0))
    validation_table = None
    if instrument:
        validation_table = [int(C.polcoef(int(field_parameter(pari,z,a).fflog(g,m))))
                            for a in range(1,q)]
    t = perf_counter()
    del C, g, z, modulus, parameter
    stages['query_cleanup'] = perf_counter()-t
    total = perf_counter()-start
    return dict(status='ok',N=N,L=L,P=N*(L-1)+2,f=f,params=params,outputs=outputs,
                compute_seconds=total,preprocess_seconds=prep,
                cold_seconds=prep+queries[0]['seconds'],queries=queries,stages=stages,
                preprocess_stage_gap_seconds=prep-sum(v for k,v in stages.items() if k!='query_cleanup'),
                total_stage_gap_seconds=total-sum(stages.values())-sum(v['seconds'] for v in queries),
                timing_eligible=not instrument, full_table_entries=m, table_is_all_parameters=True,
                reusable_objects=['supplied PARI field', 'primitive element', 'full exact cyclic polynomial'],
                memory_precheck=check, validation_full_table=validation_table,
                output_scope='exact integers', pari_version=str(pari.version()),
                pari_stack_bytes=int(pari.stacksize()), pari_stack_max_bytes=int(pari.stacksizemax()),
                storage=dict(character_polynomial_coefficient_slots=m, output_polynomial_capacity=m,
                             modulus_coefficient_slots=q,
                             note='Derived slot counts, not measured payload or RSS'))


def run_shared(N,L,f,params):
    sys.path.insert(0, str(BASE/'query_optimization_20260930/code'))
    import candidate
    return candidate.run(N,L,f,params,'native_shared_powers',instrument=False)
