"""Compare every preview centerline sample with the Python design builder.

python3 ui/verify_torushelix.py --reference /path/to/frozen_lab
No field, particle, evaluator, browser, or model calls are made.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--reference', type=Path, default=Path(__file__).resolve().parent.parent)
parser.add_argument('--study-config', type=Path, help='Optional frozen current cap; never runs the evaluator')
parser.add_argument('--locked-dir', type=Path, help='Include the exact frozen direct and construction primary JSON files')
parser.add_argument('--out', type=Path, help='Save the geometry verification record')
args = parser.parse_args()
sys.path.insert(0, str(args.reference.resolve()))
from ccsim.design import build  # noqa: E402

ui = Path(__file__).resolve().parent
study_config = json.loads(args.study_config.read_text()) if args.study_config else None
cases = [
    {'name': 'defaults', 'components': [{'family': 'torushelix'}]},
    {'name': 'odd periods', 'components': [{'family': 'torushelix', 'params': {'periods': 3}}]},
    {'name': 'negative handedness', 'components': [{'family': 'torushelix', 'params': {'periods': -3, 'phase_deg': 37.5}}]},
    {'name': 'three toroidal turns', 'components': [{'family': 'torushelix', 'params': {'l': 3, 'periods': 2}}]},
    {'name': 'six distinct coils', 'components': [{'family': 'torushelix', 'params': {'l': 3, 'periods': 6}}]},
    {'name': 'zero periods', 'components': [{'family': 'torushelix', 'params': {'periods': 0}}]},
    {'name': 'single mode', 'components': [{'family': 'torushelix', 'params': {'l': 1, 'periods': 5, 'R0': 0.5, 'r': 0.2}}]},
    {'name': 'all transforms and signed current', 'components': [{'family': 'torushelix',
        'params': {'l': 2, 'periods': -4, 'R0': 0.53, 'r': 0.31, 'phase_deg': -28.5},
        'scale': 0.76, 'reflect': 'x', 'rotate': {'axis': [1, 2, 3], 'deg': 37},
        'translate': [0.1, -0.2, 0.03], 'current': -2.25}]},
    {'name': 'Python integer coercion', 'components': [{'family': 'torushelix', 'params': {'l': 2.9, 'periods': -3.7}}]},
    {'name': 'disabled coil omitted', 'components': [
        {'family': 'torushelix', 'current': 0},
        {'family': 'torushelix', 'params': {'periods': 5}, 'current': 1.75}]},
    {'name': 'study baseline woven mesh', 'wire_radius': 0.00448, 'components': [
        {'family': 'hopfmesh', 'params': {'eta': 0.6, 'circuits': 24, 'delta_eta': 0.015,
                                       'eps': 0.7, 'periods': 4, 'sense': -1}, 'current': 1}]},
    {'name': 'leading 48-circuit mesh', 'wire_radius': 0.00448, 'components': [
        {'family': 'hopfmesh', 'params': {'eta': 0.6, 'circuits': 48, 'delta_eta': 0.02,
                                       'eps': 0.7, 'periods': 4, 'sense': -1}, 'current': 1}]},
    {'name': 'leading two-revolution woven mirror', 'wire_radius': 0.00448, 'components': [
        {'family': 'hopfmirror', 'params': {'eta': 1.0, 'circuits': 24, 'delta_eta': 0.035,
                'mode': 'woven', 'revolutions': 2, 'sense': -1}, 'current': 1}]},
]
locked_inputs = []
if args.locked_dir:
    for arm in ('direct', 'construction'):
        source = args.locked_dir / f'{arm}-primary.json'
        raw = source.read_bytes()
        design = json.loads(raw)
        cases.append(design)
        locked_inputs.append({'arm': arm, 'path': str(source.resolve()),
                              'sha256': hashlib.sha256(raw).hexdigest(), 'name': design['name']})
script = r'''
const fs=require('fs'),vm=require('vm');
const request=JSON.parse(fs.readFileSync(0,'utf8'));
const context=vm.createContext({});
vm.runInContext(fs.readFileSync(request.engine,'utf8')+'\neval(ENGINE_SRC);',context);
const results=request.cases.map(design=>context.designParts(design).map(p=>({
  points:Array.from(p.pts),sign:p.sign,closed:p.closed
})));
// Verify dispatch separately: do not mistake generic approximate mirror geometry
// for the closed zone weave used by the leading two-revolution woven design.
const routing=[];
context.hopfMeshPair=(...args)=>{routing.push({kind:'zone_weave',args});return {A:new Float64Array(6),B:new Float64Array(6)};};
context.hopfMirrorPair=(...args)=>{routing.push({kind:'generic_mirror',args});return [new Float64Array(6),new Float64Array(6)];};
for(const params of [{mode:'woven',revolutions:2},{mode:'nested',revolutions:2},{mode:'woven',revolutions:1},{mode:'woven',revolutions:4}])context.designFamilyParts('hopfmirror',params);
process.stdout.write(JSON.stringify({results,routing}));
'''
output = subprocess.check_output(['node', '-e', script], input=json.dumps({
    'engine': str(ui / 'chad_core_lab_engine.js'), 'cases': cases}).encode())
javascript = json.loads(output)
actual = javascript['results']
assert [r['kind'] for r in javascript['routing']] == ['zone_weave', 'generic_mirror', 'generic_mirror', 'generic_mirror']
assert javascript['routing'][0]['args'][2:4] == [0, 1]
assert [r['args'][2:4] for r in javascript['routing'][1:]] == [[2, 'nested'], [1, 'woven'], [4, 'woven']]
max_error = 0.0
point_count = 0
coils = 0
per_case = []
for design, browser_parts in zip(cases, actual):
    reference_parts = build(design).windings
    case_error = 0.0
    tolerance = 1e-10 if any(c['family'] in ('hopfmesh', 'hopfmirror') for c in design['components']) else 1e-12
    assert len(browser_parts) == len(reference_parts), design['name']
    for got, reference in zip(browser_parts, reference_parts):
        points = np.asarray(got['points']).reshape(-1, 3)
        expected = reference.points
        assert points.shape == expected.shape, design['name']
        assert got['sign'] == reference.current_A, design['name']
        assert got['closed'] is True, design['name']
        assert np.array_equal(points[0], points[-1]), design['name']
        error = float(np.max(np.abs(points - expected)))
        assert error < tolerance, (design['name'], error, tolerance)
        max_error = max(error, max_error)
        case_error = max(error, case_error)
        point_count += len(points)
        coils += 1
    per_case.append({'name': design['name'], 'coils': len(browser_parts),
                     'points_per_coil': [len(p['points']) // 3 for p in browser_parts],
                     'relative_signed_currents': [p['sign'] for p in browser_parts],
                     'max_absolute_error': case_error, 'tolerance': tolerance})
    if study_config:
        cap = study_config['I_max_A']
        browser_relative = np.asarray([p['sign'] for p in browser_parts])
        reference_relative = np.asarray([p.current_A for p in reference_parts])
        browser_physical = browser_relative * cap / np.max(np.abs(browser_relative))
        reference_physical = reference_relative * cap / np.max(np.abs(reference_relative))
        assert np.array_equal(browser_physical, reference_physical), design['name']
        per_case[-1]['physical_signed_currents_at_frozen_cap_A'] = browser_physical.tolist()
for name in ('chad_core_lab.html', 'chad_core_lab_standalone.html'):
    assert (ui / 'chad_core_lab_engine.js').read_text() in (ui / name).read_text(), name
report = {'cases': len(cases), 'coils': coils, 'centerline_points': point_count,
                  'max_absolute_error': max_error,
                  'tolerances': {'torushelix': 1e-12, 'iterated_weave': 1e-10},
                  'reference_geometry_sha256': hashlib.sha256((args.reference / 'ccsim/geometry.py').read_bytes()).hexdigest(),
                  'reference_design_sha256': hashlib.sha256((args.reference / 'ccsim/design.py').read_bytes()).hexdigest(),
                  'per_case': per_case,
                  'locked_primary_inputs': locked_inputs,
                  'mirror_mode_revolution_dispatch_checks': javascript['routing'],
                  'scope': 'Point geometry qualified for the listed designs. Other generic mirror modes/revolutions remain explicitly approximate; dispatch tests only establish that these parameters select distinct branches.',
                  'physical_current_scope': 'Derived from matched relative signed coil order and the frozen evaluator cap; the browser preview itself uses relative currents.',
                  'assembled_pages_match_engine': True}
if args.study_config:
    report['study_config_sha256'] = hashlib.sha256(args.study_config.read_bytes()).hexdigest()
if args.out:
    args.out.write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
