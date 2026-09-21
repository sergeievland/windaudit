"""Additive difference constraints. An arc u -> v asserts g[v]-g[u] <= c.

Certificates prove incompatibility of the geometric model, not CT-verified
annotation errors. Radial monotonicity in the sampled neighbourhood is an
additional assumption beyond non-intersection alone.
"""
from dataclasses import dataclass
from collections import defaultdict
import hashlib
import json
import math
import time
from numbers import Integral
import numpy as np
from scipy.sparse import coo_matrix
from .frame import label_diff, node_tag, wrap_angle, TWO_PI
from .matching import build_links, point_slopes, projected_gap, _cell_of
from .repair import _LexProgram

@dataclass(frozen=True)
class Arc:
    u: tuple
    v: tuple
    c: int
    group: str
    witness: dict = None


def negative_cycle(arcs):
    """Super-source Bellman-Ford, with a checkable directed negative cycle."""
    nodes = sorted({x for a in arcs for x in (a.u, a.v)})
    for a in arcs:
        if not isinstance(a.c, Integral) or isinstance(a.c, bool):
            raise ValueError('integer arc bounds required')
    dist = dict.fromkeys(nodes, 0)
    pred = {}
    changed = None
    for _ in nodes:
        changed = None
        for i, a in enumerate(arcs):
            if dist[a.v] > dist[a.u] + a.c:
                dist[a.v] = dist[a.u] + a.c
                pred[a.v] = i
                changed = a.v
        if changed is None:
            return []
    if changed is None:
        return []
    for _ in nodes:
        changed = arcs[pred[changed]].u
    start, cycle = changed, []
    while True:
        i = pred[changed]
        cycle.append(i)
        changed = arcs[i].u
        if changed == start:
            break
    cycle.reverse()
    assert sum(arcs[i].c for i in cycle) < 0
    assert all(arcs[cycle[i]].v == arcs[cycle[(i+1) % len(cycle)]].u for i in range(len(cycle)))
    return cycle


def candidates(nodes, cfg):
    """Unique comparable pairs, using the legacy sector and projection rules.

    Order evidence has no same-wrap crowd veto: it requires a large radial
    separation instead. Every such modelling difference is null-controlled.
    """
    cells = defaultdict(list)
    nt = max(3, int(round(TWO_PI / cfg.sector_dtheta_rad)))
    for k in sorted(nodes):
        if nodes[k].seam_safe:
            for q, slope in zip(nodes[k].points, point_slopes(nodes[k], cfg.sense)):
                cells[_cell_of(q, cfg, nt)].append((k, q, slope))
    seen = set()
    out = []
    import itertools
    for zi, ti in sorted(cells):
        group = list(cells[(zi, ti)])
        for dz, dt in ((0, 1), (1, -1), (1, 0), (1, 1)):
            group += cells.get((zi+dz, (ti+dt) % nt), [])
        for (ka, a, sa), (kb, b, sb) in itertools.combinations(group, 2):
            if ka == kb or abs(a.z-b.z) > cfg.sector_dz or abs(wrap_angle(a.theta-b.theta)) > cfg.sector_dtheta_rad:
                continue
            if ka > kb:
                ka, a, sa, kb, b, sb = kb, b, sb, ka, a, sa
            key = (ka, a.pid, kb, b.pid)
            if key in seen:
                continue
            seen.add(key)
            gap = projected_gap(a, sa, b, sb)
            if gap is not None:
                out.append((ka, a, kb, b, gap))
    return out


def build_constraints(nodes, cfg, margin, match=None, pairs=None):
    if not math.isfinite(margin) or margin <= 0:
        raise ValueError('margin must be finite and positive')
    match = build_links(nodes, cfg) if match is None else match
    pairs = candidates(nodes, cfg) if pairs is None else pairs
    arcs = []
    evidence = []
    def witness(a, b, pa, pb, xa, xb, gap, kind):
        return dict(a=node_tag(a), b=node_tag(b), pid_a=pa, pid_b=pb,
                    xyz_a=xa, xyz_b=xb, projected_gap=gap, kind=kind)
    # Include each supported conflicting equality, not just consensus edges.
    equality_groups = defaultdict(list)
    for l in match.links:
        equality_groups[(l.a, l.b, l.offset)].append(l)
    for (a,b,c), links in sorted(equality_groups.items()):
        gid = 'equality:'+node_tag(a)+'/'+node_tag(b)+'/'+str(c)
        l = links[0]
        w = witness(a,b,l.pid_a,l.pid_b,l.xyz_a,l.xyz_b,l.gap,'equality')
        arcs.extend([Arc(a,b,c,gid,w), Arc(b,a,-c,gid,w)])
        evidence.extend((a,l.pid_a,b,l.pid_b,gid) for l in links)
    grouped = defaultdict(list)
    for a, qa, b, qb, gap in pairs:
        if abs(gap) <= margin * cfg.wrap_spacing_median:
            continue
        d = label_diff(qa, qb, cfg.sense)
        # a outside b => g[b]-g[a] <= label_a-label_b-1.
        u,v,c = (a,b,d-1) if gap > 0 else (b,a,-d-1)
        grouped[(u,v,c)].append((a,qa,b,qb,gap))
    for (u,v,c), ls in sorted(grouped.items()):
        if len(ls) < cfg.min_edge_support:
            continue
        gid = 'order:'+node_tag(u)+'/'+node_tag(v)+'/'+str(c)
        a,qa,b,qb,gap = ls[0]
        w = witness(a,b,qa.pid,qb.pid,qa.xyz,qb.xyz,gap,'order')
        arcs.append(Arc(u,v,c,gid,w))
        evidence.extend((a,qa.pid,b,qb.pid,gid) for a,qa,b,qb,gap in ls)
    return arcs, evidence


def configuration(cfg, margin):
    payload = dict(legacy=cfg.as_dict(), order_margin=margin,
                   evidence_cost='unit_per_relation_group', order_protocol=1)
    return dict(values=payload, sha256=hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest())


def cycle_certificate(arcs):
    ids = negative_cycle(arcs)
    return dict(consistent=not ids, total_bound=sum(arcs[i].c for i in ids),
                collections=sorted({node_tag(x) for i in ids for x in (arcs[i].u,arcs[i].v)}),
                arcs=[dict(u=node_tag(arcs[i].u),v=node_tag(arcs[i].v),bound=arcs[i].c,
                           group=arcs[i].group,witness=arcs[i].witness) for i in ids])


def repair(arcs, targets=None, edge_cut=False, method="big_m", time_budget=None):
    """Reuse the unchanged lexicographic solver and its all-optima ranges.

    Each relation group costs one, regardless of how many points support it.
    An equality's two arcs share one release variable. Quarantine minimises
    nodes, then lost groups. Cut minimises lost groups.
    """
    started = time.monotonic()
    nodes = sorted({x for a in arcs for x in (a.u,a.v)})
    wanted = nodes if targets is None else targets
    if not negative_cycle(arcs):
        return dict(removed=[],cut_groups=[],optimal=True,objective=0,groups_lost=0,
                    membership={node_tag(k):'consistent' for k in wanted})
    groups = sorted({a.group for a in arcs})
    ni,gi = {k:i for i,k in enumerate(nodes)}, {g:i for i,g in enumerate(groups)}
    n,m = len(nodes),len(groups)
    bound = max(1,n*max(abs(a.c) for a in arcs))
    if bound > 2**35:
        raise ValueError('graph exceeds supported numerical range')
    big_m = 2*bound+max(abs(a.c) for a in arcs)+1
    base = n if edge_cut else 2*n
    rows,cols,vals,upper = [],[],[],[]
    def row(cc,vv,ub):
        rows.extend([len(upper)]*len(cc)); cols.extend(cc); vals.extend(vv); upper.append(ub)
    ends = {}
    for a in arcs:
        endpoints = frozenset((a.u,a.v))
        if a.group in ends and ends[a.group] != endpoints:
            raise ValueError('a relation group must have fixed endpoints')
        ends[a.group] = endpoints
        row([ni[a.v],ni[a.u],base+gi[a.group]],[1,-1,-big_m],a.c)
    if not edge_cut:
        for g in groups:
            xs = [n+ni[k] for k in sorted(ends[g])]
            y = base+gi[g]
            for x in xs: row([x,y],[1,-1],0)
            row([y]+xs,[1]+[-1]*len(xs),0)
    # Valid short-cycle cuts strengthen the big-M relaxation without changing
    # feasible integer repairs or the all-optima membership definition.
    directed = defaultdict(list)
    for a in arcs: directed[(a.u,a.v)].append(a)
    cuts_seen = set()
    for a in arcs:
        for b in directed.get((a.v,a.u), []):
            if a.c+b.c >= 0: continue
            variables = tuple(sorted({base+gi[x.group] for x in (a,b)} if edge_cut
                                     else {n+ni[a.u],n+ni[a.v]}))
            if variables not in cuts_seen:
                row(list(variables),[-1]*len(variables),-1)
                cuts_seen.add(variables)
    # Add valid three-node cycle cuts from the tightest arc per direction.
    # These strengthen either formulation; they cannot exclude a feasible
    # repair. Remaining longer cycles are handled by gauges or separation.
    tight = {pair:min(aa,key=lambda x:x.c) for pair,aa in directed.items()}
    outgoing = defaultdict(list)
    for (u,v),a in tight.items(): outgoing[u].append(a)
    for (u,v),a in tight.items():
        for b in outgoing[v]:
            w=b.v
            if len({u,v,w}) != 3 or u != min(u,v,w): continue
            c=tight.get((w,u))
            if c is None or a.c+b.c+c.c >= 0: continue
            variables=tuple(sorted({base+gi[x.group] for x in (a,b,c)} if edge_cut
                                   else {n+ni[x] for x in (u,v,w)}))
            if variables not in cuts_seen and not any(set(pair)<=set(variables) for pair in cuts_seen if len(pair)<len(variables)):
                row(list(variables),[-1]*len(variables),-1)
                cuts_seen.add(variables)
    size = base+m
    A = coo_matrix((vals,(rows,cols)),shape=(len(upper),size)).tocsr()
    program_class = _LexProgram if method == "big_m" else _CycleProgram
    if method not in ("big_m", "cycles"):
        raise ValueError("unknown repair method")
    if method == "cycles":
        # Gauges are unnecessary in the exact cycle-cut formulation.
        A = A[len(arcs):]
        upper = upper[len(arcs):]
    prog = program_class(A,np.array(upper),np.r_[np.full(n,-bound),np.zeros(size-n)],
                       np.r_[np.full(n,bound),np.ones(size-n)],list(range(n,size)))
    if method == "cycles":
        prog.arc_data = (arcs, ni, gi, n, base, edge_cut)
        prog.deadline = None if time_budget is None else started + time_budget
    count = None
    if not edge_cut:
        stage = prog.minimise_and_fix(np.r_[np.ones(n),np.zeros(m)])
        if stage: count = stage[1]
    else: stage = True
    if stage:
        stage = prog.minimise_and_fix(np.ones(m) if edge_cut else np.r_[np.zeros(n),np.ones(m)])
    if not stage:
        return dict(removed=[],cut_groups=[],optimal=False,objective=None,groups_lost=None,
                    membership={node_tag(k):'unknown' for k in wanted},status='solver_did_not_prove_optimum')
    membership = {}
    if not edge_cut:
        for k in wanted:
            span = prog.value_range(n+ni[k]) if k in ni else (0,0)
            membership[node_tag(k)] = ('unknown' if span is None else 'certain' if span[0] else 'excluded' if not span[1] else 'possible')
    prog.select = list(range(base,size)) if edge_cut else list(range(n,2*n))
    final = prog.canonicalise(stage[0])
    removed = [] if edge_cut else [k for k in nodes if final[n+ni[k]]]
    cuts = [g for g in groups if final[base+gi[g]]]
    rest = [a for a in arcs if a.group not in cuts and a.u not in removed and a.v not in removed]
    return dict(removed=[node_tag(k) for k in removed],cut_groups=cuts,
                optimal=bool(prog.ok and not negative_cycle(rest)),objective=stage[1] if edge_cut else count,
                groups_lost=stage[1],membership=membership)


class _CycleProgram(_LexProgram):
    """Exact cycle separation; an independent scalable equivalent to big-M.

    Every feasible integer repair must hit each negative cycle. Whenever the
    relaxation proposes an inconsistent remainder, add its certificate as a
    valid cut and solve again. Finite binary selections ensure termination.
    The superclass still fixes objective optima and queries variable ranges.
    """
    def _solve(self, cost):
        from scipy.sparse import csr_matrix
        arcs, ni, gi, n, base, edge_cut = self.arc_data
        while True:
            if getattr(self, 'deadline', None) is None:
                res = super()._solve(cost)
            else:
                from scipy.optimize import Bounds, LinearConstraint, OptimizeResult, milp
                from scipy.sparse import vstack
                from .repair import _quiet_native_output, SOLVER_OPTIONS
                remaining = self.deadline - time.monotonic()
                if remaining <= 0:
                    return OptimizeResult(status=1, x=None, message='Order repair wall-time budget exhausted')
                with _quiet_native_output():
                    res = milp(cost, constraints=LinearConstraint(vstack(self.rows).tocsr(),
                               np.concatenate(self.row_lo),np.concatenate(self.row_hi)),
                               integrality=np.ones(self.n_vars),bounds=Bounds(self.lb,self.hb),
                               options=dict(SOLVER_OPTIONS,time_limit=remaining))
            if res.status != 0 or res.x is None: return res
            rest = [a for a in arcs if (res.x[base+gi[a.group]] < .5 if edge_cut
                    else res.x[n+ni[a.u]] < .5 and res.x[n+ni[a.v]] < .5)]
            cycle = negative_cycle(rest)
            if not cycle:
                # Quarantine group variables must exactly encode lost groups.
                # These OR rows are installed on the first solve below.
                return res
            vs = ({base+gi[rest[i].group] for i in cycle} if edge_cut else
                  {n+ni[k] for i in cycle for k in (rest[i].u,rest[i].v)})
            row = np.zeros((1,self.n_vars)); row[0,list(vs)] = -1
            self.rows.append(csr_matrix(row)); self.row_lo.append(np.array([-np.inf])); self.row_hi.append(np.array([-1.]))


def minimum_edge_cut(arcs, time_budget=None):
    """Exact SCC decomposition of the unit-cost group-cut objective.

    An arc between strongly connected components is on no directed cycle.
    Every group has fixed endpoints, so groups cannot couple cyclic SCCs.
    Objectives therefore add and component-wise canonical cuts give the same
    global canonical optimum. This decomposition is NOT used for quarantine,
    whose secondary cost can couple components through lost incident groups.
    """
    from scipy.sparse.csgraph import connected_components
    started=time.monotonic()
    if not negative_cycle(arcs): return repair(arcs,targets=[],edge_cut=True,method='cycles')
    ends={}
    for a in arcs:
        endpoints=frozenset((a.u,a.v))
        if a.group in ends and ends[a.group]!=endpoints:raise ValueError('a relation group must have fixed endpoints')
        ends[a.group]=endpoints
    nodes=sorted({k for a in arcs for k in (a.u,a.v)});idx={k:i for i,k in enumerate(nodes)}
    adj=coo_matrix((np.ones(len(arcs)),([idx[a.u] for a in arcs],[idx[a.v] for a in arcs])),shape=(len(nodes),len(nodes))).tocsr()
    _,labels=connected_components(adj,directed=True,connection='strong')
    components=defaultdict(list)
    for a in arcs:
        if labels[idx[a.u]]==labels[idx[a.v]]:components[int(labels[idx[a.u]])].append(a)
    cuts=[];objective=0.;optimal=True;details=[]
    for key,aa in sorted(components.items()):
        if not negative_cycle(aa):continue
        budget=None if time_budget is None else max(0.,time_budget-(time.monotonic()-started))
        result=repair(aa,targets=[],edge_cut=True,method='cycles',time_budget=budget)
        details.append(dict(nodes=sorted({node_tag(k) for a in aa for k in (a.u,a.v)}),optimal=result['optimal'],objective=result['objective']))
        cuts.extend(result['cut_groups']);optimal=optimal and result['optimal']
        if result['objective'] is not None:objective+=result['objective']
    cut_set=set(cuts)
    verified=not negative_cycle([a for a in arcs if a.group not in cut_set])
    return dict(removed=[],cut_groups=sorted(cut_set),optimal=bool(optimal and verified),
                objective=objective if optimal else None,groups_lost=len(cut_set),membership={},components=details)
