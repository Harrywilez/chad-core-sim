"""Structural proposals for the existing component-list design format.

The library contains component dictionaries, including optional placement/current.
No geometry algebra is invented here: ccsim.design builds and validates proposals.
"""
from __future__ import annotations
import copy
import math
import random


def propose(design, rng: random.Random, library, parameter_mutator, sigma=0.15,
            max_components=6, structural_rate=0.5, operation=None):
    """Return (independent candidate, move); geometry validity is checked by evaluator.

    A donor can be any component from any existing design. Replacing a family uses
    the donor's own parameters, preventing incompatible parameter carry-over.
    """
    if max_components < 1 or not 0 <= structural_rate <= 1:
        raise ValueError('max_components must be positive; structural_rate must be in [0,1]')
    if not design.get('components') or len(design['components']) > max_components:
        raise ValueError('seed component count must be in [1,max_components]')
    n = len(design['components'])
    allowed = ['parameter', 'current']
    if library:
        allowed.append('replace')
        if n < max_components:
            allowed.append('add')
    if n < max_components:
        allowed.append('duplicate')
    if n > 1:
        allowed.append('remove')
    if operation is None:
        structural = [x for x in allowed if x != 'parameter']
        operation = rng.choice(structural) if structural and rng.random() < structural_rate else 'parameter'
    if operation not in allowed:
        raise ValueError(f'operation {operation!r} unavailable at component count {n}')
    if operation == 'parameter':
        return parameter_mutator(design, sigma, rng), operation
    d = copy.deepcopy(design)
    cs = d['components']
    if operation == 'current':
        c = rng.choice(cs)
        old = float(c.get('current', 1.0))
        if not math.isfinite(old):
            raise ValueError('component current must be finite')
        # Multiplicative steps preserve existing senses, and span weak/strong donors.
        # A disabled coil can be reactivated at +/-0.1 relative current.
        new = old * 2 ** rng.uniform(-1., 1.) if old else rng.choice([-0.1, 0.1])
        if not math.isfinite(new):
            raise ValueError('current mutation overflowed')
        c['current'] = new
    elif operation == 'remove':
        del cs[rng.randrange(n)]
    elif operation == 'duplicate':
        c = copy.deepcopy(rng.choice(cs))
        t = list(c.get('translate', [0., 0., 0.]))
        j = rng.randrange(3)
        t[j] += rng.choice([-1, 1]) * rng.uniform(0.15, 0.4)
        c['translate'] = t
        c['current'] = c.get('current', 1.) * rng.choice([-1, 1])
        cs.append(c)
    else:
        idx = rng.randrange(n) if operation == 'replace' else None
        existing = {c['family'] for c in cs} if idx is None else {cs[idx]['family']}
        different = [c for c in library if c['family'] not in existing]
        c = copy.deepcopy(rng.choice(different or library))
        # Separate repeated donors and explore surrounding placements without changing senses.
        c['scale'] = max(0.1, c.get('scale', 1.) * rng.uniform(0.85, 1.15))
        c['translate'] = [float(t) + rng.uniform(-0.12, 0.12)
                          for t in c.get('translate', [0., 0., 0.])]
        if idx is None:
            cs.append(c)
        else:
            cs[idx] = c
    return d, operation


def component_library(designs):
    """Extract independent donors; preserve multi-circuit family internals/transforms."""
    return [copy.deepcopy(c) for d in designs for c in d['components']]


def choose_parent(best, population, rng, exploration=0.25):
    """Explore valid previous candidates, even below the incumbent; retain best separately."""
    if not 0 <= exploration <= 1:
        raise ValueError('exploration must be in [0,1]')
    return rng.choice(population) if population and rng.random() < exploration else best
