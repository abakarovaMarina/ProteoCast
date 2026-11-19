#!/usr/bin/env python3
"""MSA for a ProteoCast web job submitted with a UniProt accession (user 2026-11-19).

AlphaFold DB no longer serves its MSAs (403 on files/msa/*, "Anonymous caller does not have storage.objects.get access"
to bucket a3m-files-prod, seen 2026-11-19), so every UniProt job of proteocast.ijm.fr failed with an empty MSA. This script
runs on the server before the container and writes <job dir>/<ACC>.a3m (first header = the accession), which run.sh then
takes as an uploaded MSA (upload mode, structure still fetched from AlphaFold DB with the accession as 3rd argument).

Sources, in order:
  1. our own MSAs (ProteocastOthers runs: ColabFold search, UniRef30 2302 + ColabFold envdb 202108), in STORE as
     <ACC>.a3m.gz; used only if the query sequence is the current UniProt sequence (or UniProt can't be reached)
  2. the public ColabFold MMseqs2 server (api.colabfold.com, mode env: the same databases as our runs)
On success: <ACC>.a3m and msa_source.txt in the job dir, exit 0. On failure: msa_error.txt (message for the user), exit 1.
Standard library only (the host's python3).
Usage: get_msa.py <accession> <job dir>
"""
import gzip, io, json, os, re, sys, tarfile, time, urllib.error, urllib.parse, urllib.request

STORE = os.environ.get('MSA_STORE', '/data/msa_store')
CF = os.environ.get('COLABFOLD_API', 'https://api.colabfold.com')
UA = 'ProteoCast web server (https://proteocast.ijm.fr)'
CF_WAIT = int(os.environ.get('CF_WAIT', 5400))   # seconds before giving up on the ColabFold server
ACC_RE = re.compile(r'^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[0-9]+)?$')


def log(*a):
    print('[get_msa]', *a, flush=True)


def http(url, data=None, timeout=60):
    req = urllib.request.Request(url, data=urllib.parse.urlencode(data).encode() if data else None,
                                 headers={'User-Agent': UA})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def uniprot_seq(acc):
    try:
        txt = http(f'https://rest.uniprot.org/uniprotkb/{acc}.fasta').decode()
    except Exception as e:
        log('UniProt not reachable:', e)
        return None
    seq = ''.join(l.strip() for l in txt.splitlines() if l and not l.startswith('>'))
    return seq or ''


def a3m_query(lines):
    """sequence of the first record (the query), gaps/insertions removed"""
    seq, seen = [], False
    for l in lines:
        if l.startswith('>'):
            if seen:
                break
            seen = True
        elif seen:
            seq.append(l.strip())
    return re.sub(r'[^A-Z]', '', ''.join(seq))


def from_store(acc, useq):
    p = os.path.join(STORE, acc + '.a3m.gz')
    if not os.path.exists(p):
        return None
    lines = gzip.open(p, 'rt').read().splitlines()
    q = a3m_query(lines)
    if useq is not None and useq != q:
        log(f'stored MSA of {acc} is for another sequence version ({len(q)} aa, UniProt {len(useq)} aa): not used')
        return None
    return lines


def colabfold(acc, seq):
    """single query, mode env; returns a3m lines (uniref + environmental), query once"""
    q = f'>101\n{seq}\n'
    t0 = time.time()
    while True:   # submission; RATELIMIT / MAINTENANCE: wait and retry
        r = json.loads(http(f'{CF}/ticket/msa', {'q': q, 'mode': 'env'}))
        if r.get('status') not in ('RATELIMIT', 'MAINTENANCE', 'UNKNOWN'):
            break
        if time.time() - t0 > CF_WAIT:
            raise RuntimeError(f"ColabFold server busy ({r.get('status')})")
        log('ColabFold server:', r.get('status'), '- retrying in 60 s')
        time.sleep(60)
    tid = r.get('id')
    while r.get('status') in ('PENDING', 'RUNNING', 'UNKNOWN'):
        if time.time() - t0 > CF_WAIT:
            raise RuntimeError('ColabFold server did not answer in time')
        time.sleep(10)
        r = json.loads(http(f'{CF}/ticket/{tid}'))
    if r.get('status') != 'COMPLETE':
        raise RuntimeError(f"ColabFold server: {r.get('status')}")
    tf = tarfile.open(fileobj=io.BytesIO(http(f'{CF}/result/download/{tid}', timeout=600)))
    out = []
    for name in ('uniref.a3m', 'bfd.mgnify30.metaeuk30.smag30.a3m'):
        lines = tf.extractfile(name).read().decode().replace('\x00', '').splitlines()
        recs, cur = [], None   # records; drop the query (>101) of the second file, keep it once
        for l in lines:
            if l.startswith('>'):
                cur = [l]
                recs.append(cur)
            elif cur is not None and l.strip():
                cur.append(l)
        for rec in recs:
            if rec[0].split()[0] == '>101':
                if out:
                    continue
                rec = [f'>{acc}'] + rec[1:]
            out.extend(rec)
    return out


def main():
    acc, job = sys.argv[1].strip().upper(), sys.argv[2]
    if not ACC_RE.match(acc):
        open(os.path.join(job, 'msa_error.txt'), 'w').write(
            f'"{sys.argv[1]}" is not a UniProt accession (one accession per job, e.g. P38398 or P38398-2).')
        return 1
    useq = uniprot_seq(acc)
    if useq == '':
        open(os.path.join(job, 'msa_error.txt'), 'w').write(f'UniProt has no sequence for {acc}.')
        return 1
    src, lines = None, from_store(acc, useq)
    if lines:
        src = 'ProteoCast proteome runs (ColabFold search: UniRef30 2302 + ColabFold environmental database 202108)'
    elif useq:
        log(f'{acc}: no stored MSA, asking the ColabFold MMseqs2 server ({len(useq)} aa)')
        try:
            lines, src = colabfold(acc, useq), 'ColabFold MMseqs2 server (api.colabfold.com, UniRef30 + environmental database)'
        except Exception as e:
            log('ColabFold failed:', e)
            open(os.path.join(job, 'msa_error.txt'), 'w').write(
                f'No MSA could be built for {acc}: the AlphaFold Database no longer provides its MSAs and the ColabFold '
                f'MMseqs2 server failed ({e}). Please try again later or submit your own MSA.')
            return 1
    else:
        open(os.path.join(job, 'msa_error.txt'), 'w').write(
            f'No MSA could be built for {acc}: UniProt could not be reached. Please try again later or submit your own MSA.')
        return 1
    with open(os.path.join(job, acc + '.a3m'), 'w') as f:
        f.write('\n'.join(lines) + '\n')
    n = sum(l.startswith('>') for l in lines)
    open(os.path.join(job, 'msa_source.txt'), 'w').write(f'{src}\n{n} sequences\n')
    log(f'{acc}: MSA from {src}, {n} sequences')
    return 0


if __name__ == '__main__':
    sys.exit(main())
