"""Exact paid inverse-permutation / additive Walsh recurrence for binary Kl_L."""
from pathlib import Path
import sys
from time import perf_counter
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'reference/code'))
import convolution as field


def hadamard(values):
    width=1
    while width<len(values):
        for start in range(0,len(values),2*width):
            for k in range(start,start+width):
                x,y=values[k],values[k+width]
                values[k],values[k+width]=x+y,x-y
        width*=2


def run(N,L,f,params,instrument=False):
    start=perf_counter(); stages={}
    if not field.irreducible(N,f):
        raise ValueError('reducible modulus')
    q=1<<N
    stages['field_validation']=perf_counter()-start
    t=perf_counter()
    # Linear trace form and its multiplication pairing, from the specified basis.
    trace_mask=sum(field.absolute_trace(1<<j,N,f)<<j for j in range(N))
    columns=[sum(((field.field_mul(1<<i,1<<j,N,f)&trace_mask).bit_count()%2)<<j for j in range(N)) for i in range(N)]
    pairing=[0]*q
    for a in range(1,q):
        bit=a&-a; pairing[a]=pairing[a^bit]^columns[bit.bit_length()-1]
    # Primitive search + cycle is paid; inverse uses index reflection, no free logs.
    primitive,exponents,logs=field._cyclic_index(N,f)
    inverses=[0]*q
    for i,a in enumerate(exponents):
        inverses[a]=exponents[(-i)%(q-1)]
    del exponents,logs
    stages['trace_pairing_and_inverse_index']=perf_counter()-t
    t=perf_counter()
    table=[1-2*((a&trace_mask).bit_count()%2) for a in range(q)]
    table[0]=0
    stages['character']=perf_counter()-t
    t=perf_counter()
    for rank in range(1,L):
        work=[table[inverses[x]] for x in range(q)]
        hadamard(work)
        table=[work[pairing[a]] for a in range(q)]
        table[0]=0
    if L>1:
        del work
    stages['inverse_fwht_recurrence']=perf_counter()-t
    prep=perf_counter()-start
    outputs=[]; query_seconds=[]
    for a in params:
        t=perf_counter(); value=table[a]; elapsed=perf_counter()-t
        outputs.append(value); query_seconds.append(elapsed)
    total=perf_counter()-start
    stages['lookup']=sum(query_seconds)
    result=dict(outputs=outputs,preprocess_seconds=prep,cold_seconds=prep+query_seconds[0],compute_seconds=total,query_seconds=query_seconds,stages_seconds=stages,stage_gap_seconds=total-sum(stages.values()),full_table_entries=q-1,primitive_element=primitive,reusable_objects=['inverse_index','trace_pairing_index','full_integer_table'],timing_eligible=not instrument)
    if instrument:
        result['storage']={'retained_slots':3*q,'retained_bytes':field._retained_python_bytes((inverses,pairing,table)), 'integer_bits_max':max(abs(x).bit_length() for x in table),'fwht_integer_additions':(L-1)*q*N}
    return result
