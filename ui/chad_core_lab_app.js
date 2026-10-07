// ===================== main thread: geometry, rendering, UI =====================
eval(ENGINE_SRC);   // geometry + numerics available here; the worker gets the same source

const $ = id => document.getElementById(id);
const fmt = (v, d=3) => (v===undefined||v===null||!isFinite(v)) ? '–' : (Math.abs(v)>=1e4||Math.abs(v)<1e-3&&v!==0 ? v.toExponential(d-1) : v.toFixed(d));
const state = { circuits: [], active: [], returns: [], clearance: NaN, length: 0, activeLength: 0, grid: null, lines: null, particles: null, name: '', busy: false };

// ---------- geometry parameter forms ----------
const FAMILY_PARAMS = {
  hopfdrift: [['variant','select',['s64_150','s64_180','s64_210','d48_150','d48_180','d48_210'],'s64_180','Codex variant (η₁, sweep)'],['circuits','number',[4,40,1],10,'circuits']],
  hopfcont: [['eta0','number',[0.4,0.88,0.01],0.82,'η start'],['eta1','number',[0.2,0.88,0.01],0.64,'η end'],['sweep','number',[30,3600,10],720,'φ sweep (deg)'],['circuits','number',[4,200,1],48,'circuits']],
  hopftorus: [['eta','number',[0.2,0.88,0.01],0.70,'η (torus)'],['circuits','number',[4,200,1],24,'circuits (fibres)'],['revs','number',[0.25,8,0.25],1,'revolutions of φ']],
  precess: [['inward','number',[0.05,0.4,0.01],0.22,'inward step'],['rotation','number',[0,60,1],24,'rotation per circuit (deg)'],['slip','number',[0,20,0.5],3,'minor slip (deg)']],
  phase: [['inward','number',[0.05,0.4,0.01],0.32,'inward step'],['slip','number',[0,20,0.5],3,'minor slip (deg)']],
  codexhopf: [['inward','number',[0,0.5,0.01],0.18,'inward'],['rotation','number',[0,60,1],20,'rotation (deg)']],
  recursive: [['baseR','number',[0.3,0.9,0.01],0.72,'base circle radius'],['r1','number',[0.05,0.6,0.01],0.24,'level-1 radius ratio'],['n1','number',[3,40,1],12,'level-1 turns'],['r2','number',[0,0.8,0.01],0.30,'level-2 radius ratio (0 = off)'],['n2','number',[1,20,1],6,'level-2 turns per turn'],['alt','select',['same','opposite'],'same','level-2 handedness']],
  solenoid: [['R','number',[0.1,0.8,0.01],0.45,'radius'],['L','number',[0.2,1.6,0.05],1.1,'length'],['turns','number',[2,60,1],12,'turns']],
  hopfmirror: [['eta','number',[0.3,0.88,0.01],0.70,'η (torus)'],['circuits','number',[4,80,1],24,'circuits per strand'],['revs','number',[0.25,8,0.25],1,'revolutions of φ'],['mode','select',['woven','nested'],'woven','fit: woven (interlocked) / nested'],['deta','number',[0.005,0.08,0.005],0.025,'δη over/under offset'],['sense','select',['same','opposed'],'opposed','mirror current']],
  hopfmesh: [['eta','number',[0.3,0.88,0.01],0.60,'η (torus; smaller = thinner tube)'],['circuits','number',[6,60,1],24,'circuits per strand (closed knots)'],['eps','number',[0,0.9,0.05],0.7,'helical modulation ε (0 = pure mirror pair)'],['periods','number',[1,9,1],4,'field periods n'],['deta','number',[0.005,0.08,0.005],0.03,'δη over/under offset (≤ 0.015 weaves, ≥ 0.03 nests)'],['hmode','select',['toroidal','poloidal'],'toroidal','modulated current: toroidal (stellarator) / poloidal (ripple)'],['sense','select',['opposed','same'],'opposed','mirror current']],
  baseball: [['radius','number',[0.3,0.9,0.01],0.72,'sphere radius'],['amp','number',[5,80,1],40,'seam amplitude (deg)'],['turns','number',[1,16,1],6,'turns'],['pitch','number',[0.01,0.08,0.005],0.035,'radial pitch per turn']],
  yinyang: [['router','number',[0.4,0.9,0.01],0.76,'outer radius'],['rinner','number',[0.3,0.85,0.01],0.60,'inner radius'],['amp','number',[5,80,1],40,'seam amplitude (deg)'],['turns','number',[1,12,1],4,'turns each'],['sense','select',['same','opposed'],'same','inner current']],
  picket: [['rings','number',[2,15,1],5,'rings'],['radius','number',[0.3,0.9,0.01],0.72,'sphere radius'],['alt','select',['alternating','aligned'],'alternating','ring currents']],
  sphere: [['radius','number',[0.3,0.9,0.01],0.76,'sphere radius'],['turns','number',[4,80,1],24,'turns (uniform in z)']],
  link: [['radius','number',[0.3,0.8,0.01],0.62,'ring radius'],['sense','select',['same','opposed'],'same','second ring current']],
  cone: [['rbase','number',[0.05,0.9,0.01],0.32,'base radius'],['rtip','number',[0.01,0.5,0.01],0.06,'tip radius'],['height','number',[0.1,1.6,0.05],0.5,'height'],['turns','number',[2,40,1],8,'turns'],['tip','select',['up','down'],'up','tip points']],
  design: [],
};
function renderParams(){
  const fam=$('family').value; const host=$('geo-params'); host.innerHTML='';
  if(fam==='design'){ const ta=document.createElement('textarea'); ta.id='p_designjson'; ta.rows=10; ta.spellcheck=false; ta.style.cssText='width:100%;box-sizing:border-box;font-family:inherit;font-size:11px;background:#0b1015;color:var(--text);border:1px solid var(--line-2);border-radius:3px;padding:6px;resize:vertical';
    ta.value=state.designSpec?JSON.stringify(state.designSpec,null,1):''; ta.placeholder='paste a ccsim design JSON, or pick one in the Designs tab'; host.appendChild(ta);
    const n=document.createElement('p'); n.className='note'; n.textContent='components with ccsim family names and parameters; edit and press Build geometry'; host.appendChild(n); return; }
  for(const [key,type,range,def,label] of FAMILY_PARAMS[fam]){
    const row=document.createElement('div'); row.className='row'; const lab=document.createElement('label'); lab.textContent=label; row.appendChild(lab);
    let inp; if(type==='select'){inp=document.createElement('select'); range.forEach(o=>{const op=document.createElement('option');op.value=o;op.textContent=o;inp.appendChild(op);}); inp.value=def;}
    else {inp=document.createElement('input'); inp.type='number'; inp.min=range[0]; inp.max=range[1]; inp.step=range[2]; inp.value=def;}
    inp.id='p_'+key; row.appendChild(inp); host.appendChild(row);
  }
}
const P = k => { const el=$('p_'+k); return el ? (el.tagName==='SELECT'?el.value:parseFloat(el.value)) : undefined; };
const VARIANTS={s64_150:[0.64,150],s64_180:[0.64,180],s64_210:[0.64,210],d48_150:[0.48,150],d48_180:[0.48,180],d48_210:[0.48,210]};
function makeCore(){
  const fam=$('family').value;
  switch(fam){
    case 'hopfdrift': {const [e1,sw]=VARIANTS[P('variant')]; return {pts:hopfDriftPts(0.82,e1,sw,P('circuits'),Math.PI/2,true), name:`hopf-drift ${P('variant')} × ${P('circuits')}`};}
    case 'hopfcont': return {pts:hopfDriftPts(P('eta0'),P('eta1'),P('sweep'),P('circuits'),Math.PI/2,false), name:`hopf continued ${P('sweep')}° / ${P('circuits')} (η ${P('eta0')}→${P('eta1')})`};
    case 'hopftorus': return {pts:hopfTorusPts(P('eta'),P('circuits'),P('revs')), name:`hopf torus η=${P('eta')}, ${P('circuits')}×${P('revs')}`};
    case 'precess': return {pts:codexToroidal('precess',P('inward'),P('rotation'),P('slip')), name:`precessing loops (in ${P('inward')}, rot ${P('rotation')}°)`};
    case 'phase': return {pts:codexToroidal('phase',P('inward'),3,P('slip')), name:`phase-slip spiral (in ${P('inward')})`};
    case 'codexhopf': return {pts:codexHopf(P('inward'),P('rotation')), name:`Hopf-coordinate drift (in ${P('inward')}, rot ${P('rotation')}°)`};
    case 'recursive': {const lv=[[P('r1'),P('n1')]]; if(P('r2')>0)lv.push([P('r2'),P('n2')]); return {pts:recursiveWinding(P('baseR'),lv,lv.length>1?24:48,P('alt')==='opposite'), name:`recursive ${lv.length+1} levels (${lv.map(l=>l[1]).join('×')})`};}
    case 'solenoid': return {pts:solenoidPts(P('R'),P('L'),P('turns'),60), name:`solenoid ${P('turns')} turns`};
    case 'hopfmirror': {const s=P('sense')==='opposed'?-1:1; if(Math.abs(P('revs')-2)<1e-9&&P('mode')==='woven'){const m=hopfMeshPair(P('eta'),P('circuits'),0,1,P('deta'),'toroidal'); return {parts:[{pts:m.A,sign:1,closed:true},{pts:m.B,sign:s,closed:true}], name:`hopf mirror pair η=${P('eta')} ${P('circuits')}×2 closed knots ${P('sense')} (${m.crossings} crossings, ${m.weave})`};}
      const [A,B]=hopfMirrorPair(P('eta'),P('circuits'),P('revs'),P('mode'),P('deta')); return {parts:[{pts:A,sign:1},{pts:B,sign:s}], name:`hopf mirror ${P('mode')} η=${P('eta')} ${P('circuits')}×${P('revs')} ${P('sense')} (session-5 scheme)`};}
    case 'hopfmesh': {const m=hopfMeshPair(P('eta'),P('circuits'),P('eps'),P('periods'),P('deta'),P('hmode')); const s=P('sense')==='opposed'?-1:1; return {parts:[{pts:m.A,sign:1,closed:true},{pts:m.B,sign:s,closed:true}], name:`meshed torus II η=${P('eta')} ${P('circuits')}×2 ε=${P('eps')} n=${P('periods')} ${P('hmode')} ${P('sense')} (${m.crossings} crossings, ${m.weave}${m.weave!=='nested'?', '+Math.round(100*(1-m.zoneFraction))+'% free to weave':''})`};}
    case 'baseball': return {pts:baseballSeam(P('radius'),P('amp'),P('turns'),P('pitch')), name:`baseball seam (${P('turns')} turns, ${P('amp')}°)`};
    case 'yinyang': {const o=baseballSeam(P('router'),P('amp'),P('turns'),0.03); const i=baseballSeam(P('rinner'),P('amp'),P('turns'),0.03); const R=rotMat([0,0,1],Math.PI/2); const ir=new Float64Array(i.length); for(let k=0;k<i.length;k+=3){ir[k]=R[0][0]*i[k]+R[0][1]*i[k+1];ir[k+1]=R[1][0]*i[k]+R[1][1]*i[k+1];ir[k+2]=i[k+2];}
      return {parts:[{pts:o,sign:1},{pts:ir,sign:P('sense')==='opposed'?-1:1}], name:`yin-yang pair (${P('turns')}+${P('turns')} turns)`};}
    case 'picket': {const parts=picketFence(P('rings'),P('radius')); if(P('alt')==='aligned')parts.forEach(p=>p.sign=1); return {parts, name:`picket fence ${P('rings')} rings${P('alt')==='aligned'?' aligned':''}`};}
    case 'sphere': return {pts:sphereWinding(P('radius'),P('turns')), name:`sphere winding ${P('turns')} turns (uniform B)`};
    case 'link': {const parts=hopfLink(P('radius')); if(P('sense')==='opposed')parts[1].sign=-1; return {parts, name:`Hopf link${P('sense')==='opposed'?', opposed':''}`};}
    case 'cone': {let pts=conePts(P('rbase'),P('rtip'),P('height'),P('turns'),90); if(P('tip')==='down'){for(let i=2;i<pts.length;i+=3)pts[i]=-pts[i];} return {pts, name:`tornado coil ${P('turns')} turns (${P('rbase')}→${P('rtip')}, tip ${P('tip')})`};}
    case 'design': { const txt=($('p_designjson')||{}).value; let spec=state.designSpec; if(txt&&txt.trim()){spec=JSON.parse(txt);} if(!spec||!spec.components)throw new Error('no design loaded — pick one in the Designs tab or paste JSON');
      state.designSpec=spec; const dp=designParts(spec); const approx=dp.some(p=>p.approx); const parts=dp.map(p=>({pts:p.pts,sign:p.sign,closed:p.closed}));
      return {parts, name:(spec.name||'design')+' ['+designSummary(spec)+']'+(approx?' (approx. weave)':'')}; }
  }
}
function buildGeometry(){
  try{
    const core=makeCore(); const layout=$('layout').value; const cs=parseFloat($('corescale').value);
    let cores=[]; if(core.parts){ for(const part of core.parts){ for(const c of assemble(part.pts,layout,cs)) cores.push({pts:c.pts,sign:c.sign*part.sign,closed:!!part.closed}); } } else cores=assemble(core.pts,layout,cs); state.circuits=[]; state.active=[]; state.returns=[]; let L=0,La=0; let clr=Infinity;
    for(const c of cores){let closed, ret; if(c.closed){closed=Float64Array.from(c.pts); ret=new Float64Array(0);} else {ret=remoteReturn(c.pts,3.8,300); closed=new Float64Array(c.pts.length+ret.length-3); closed.set(c.pts); closed.set(ret.subarray(3),c.pts.length);}
      state.circuits.push({pts:closed,sign:c.sign}); state.active.push({pts:c.pts,sign:c.sign}); state.returns.push(ret); L+=polyLength(closed); La+=polyLength(c.pts);
      clr=Math.min(clr,minNonlocal(c.pts,12,cores.length>2||c.pts.length>6000?3:1));}
    for(let i=0;i<cores.length;i++)for(let j=i+1;j<cores.length;j++){const a=cores[i].pts,b=cores[j].pts; let best=Infinity; for(let s=0;s<a.length;s+=9)for(let t=0;t<b.length;t+=9){const dx=a[s]-b[t],dy=a[s+1]-b[t+1],dz=a[s+2]-b[t+2];const d=dx*dx+dy*dy+dz*dz;if(d<best)best=d;} clr=Math.min(clr,Math.sqrt(best));}
    state.clearance=clr; state.length=L; state.activeLength=La; state.name=core.name+(layout!=='single'?` × ${layout}`:''); state.grid=null; state.lines=null; state.particles=null;
    drawConductors(); clearLines(); clearParticles(); drawSlice(null); updateGeoReadout(); updateScaleReadout(); state.poincare=null; $('poincare-canvas').hidden=true; $('ro-poinc').innerHTML='<span class="note">run a Poincaré section (Field block) — toroidal families only</span>'; setStatus('geometry built','ok');
    if($('family').value!=='design'){state.designRow=null;} updateDesignReadout();
    $('geo-note').textContent=`${cores.length} circuit${cores.length>1?'s':''}, ${state.active.reduce((a,c)=>a+c.pts.length/3,0)} conductor points. Clearance ${fmt(clr,4)} ball units.`;
  }catch(e){setStatus(e.message,'bad'); $('geo-note').textContent=e.message;}
}

// ---------- worker ----------
const workerSrc = ENGINE_SRC + `
self.onmessage = function(ev){ const m=ev.data;
  if(m.type==='field'){ const B=biotSavartGrid(m.circuits,m.n,m.half,m.a); const D=wireDistanceGrid(m.circuits,m.n,m.half,2); const stats=fieldStats(B,D,m.n,m.half); let seed=99; const rng=()=>{seed=(seed*1664525+1013904223)>>>0;return seed/4294967296;}; stats.exact=maxwellChecks(m.circuits,m.a,rng); self.postMessage({type:'field',B,D,stats},[B.buffer,D.buffer]); }
  else if(m.type==='lines'){ const res=traceLines(m.B,m.D,m.n,m.half,m.seeds,m.wall,m.step,m.maxLen,m.wireClear); self.postMessage({type:'lines',res}); }
  else if(m.type==='particles'){ const res=runParticles(m.B,m.D,m.n,m.half,m.params); self.postMessage({type:'particles',res}); }
  else if(m.type==='poincare'){ const res=poincareAnalyse(m.circuits,m.a,m.R0,m.rTube,m.nSeeds,m.turns,m.step,m.wireClear); self.postMessage({type:'poincare',res}); }
};`;
const worker = new Worker(URL.createObjectURL(new Blob([workerSrc],{type:'text/javascript'})));
let pending=null;
worker.onmessage = ev => { const m=ev.data; state.busy=false;
  if(m.type==='field'){ state.grid={B:m.B,D:m.D,n:pending.n,half:pending.half,stats:m.stats,a:pending.a}; drawSlice(state.grid); updateFieldReadout(); updateScaleReadout(); updateLiveChecks(); setStatus(`field ready · B_rms ${m.stats.brms.toExponential(2)} T/A · Ampère ${m.stats.exact.ampere.toFixed(4)} μ₀I`,'ok'); }
  else if(m.type==='lines'){ state.lines=m.res; drawLines(m.res); updateLinesReadout(); setStatus('field lines traced','ok'); }
  else if(m.type==='particles'){ state.particles=m.res; setupPlayback(m.res); updatePartReadout(); setStatus(`particles: ${m.res.retained}/${m.res.N} retained`,'ok'); }
  else if(m.type==='poincare'){ state.poincare=m.res; drawPoincare(m.res); updatePoincareReadout(m.res); const ns=m.res.lines.filter(l=>l.kind==='surface').length; setStatus(`Poincaré: ${ns}/${m.res.lines.length} seeds on closed surfaces`,'ok'); }
  $('compute').disabled=$('lines').disabled=$('run').disabled=$('poincare').disabled=false;
};
function setStatus(msg,kind){ $('status').textContent=msg; $('dot').className='dot'+(kind==='busy'?' busy':kind==='ok'?' ok':''); }
function post(msg,transfer){ state.busy=true; $('compute').disabled=$('lines').disabled=$('run').disabled=$('poincare').disabled=true; worker.postMessage(msg,transfer||[]); }
function runPoincare(){ if(!state.circuits.length)buildGeometry(); const a=parseFloat($('wirer').value); const act=state.active[0].pts; let rmin=Infinity,rmax=0; for(let i=0;i<act.length;i+=3){const r=Math.hypot(act[i],act[i+1]); if(r<rmin)rmin=r; if(r>rmax)rmax=r;}
  const R0=(rmin+rmax)/2, rTube=(rmax-rmin)/2; setStatus('Poincaré section on the exact field (≈ ½–2 min)…','busy');
  post({type:'poincare',circuits:state.circuits.map(c=>({pts:Float64Array.from(c.pts),sign:c.sign})),a,R0,rTube,nSeeds:parseInt($('npoinc').value),turns:parseInt($('nturns').value),step:0.012,wireClear:Math.max(0.02,2*a)}); }
function drawPoincare(res){ const cv=$('poincare-canvas'); cv.hidden=false; const ctx=cv.getContext('2d'); const W=cv.width,H=cv.height; ctx.fillStyle='rgba(11,16,21,0.92)'; ctx.fillRect(0,0,W,H);
  // frame: window around the axis
  const act=state.active[0].pts; const pts=[]; for(let i=0;i<act.length;i+=3){const phi=Math.atan2(act[i+1],act[i]); if(Math.abs(phi)<0.06)pts.push([Math.hypot(act[i],act[i+1]),act[i+2]]);}
  let Rs=pts.map(p=>p[0]),zs=pts.map(p=>p[1]); const Rc=(Math.min(...Rs)+Math.max(...Rs))/2, half=Math.max(Math.max(...Rs)-Math.min(...Rs),Math.max(...zs)-Math.min(...zs))/2*1.1||0.5;
  const X=R=>W/2+(R-Rc)/half*(W/2-8), Y=z=>H/2-z/half*(H/2-8);
  ctx.fillStyle='#e6b34f'; for(const [R,z] of pts){ctx.fillRect(X(R)-1,Y(z)-1,2,2);}
  const col={surface:'#45b8a6','island/chaotic':'#8a7bd6',open:'#e05d57'};
  for(const l of res.lines){ctx.fillStyle=col[l.kind]||'#888'; for(let j=0;j<l.punct.length;j+=2){ctx.fillRect(X(l.punct[j])-1.2,Y(l.punct[j+1])-1.2,2.4,2.4);}}
  ctx.strokeStyle='#dde5ea'; ctx.beginPath(); ctx.moveTo(X(res.axis[0])-5,Y(res.axis[1])); ctx.lineTo(X(res.axis[0])+5,Y(res.axis[1])); ctx.moveTo(X(res.axis[0]),Y(res.axis[1])-5); ctx.lineTo(X(res.axis[0]),Y(res.axis[1])+5); ctx.stroke();
  ctx.fillStyle='#8695a0'; ctx.font='10px IBM Plex Mono, monospace'; ctx.fillText('Poincaré section φ = 0  (R, z)',8,12); }
function updatePoincareReadout(res){ const surf=res.lines.filter(l=>l.kind==='surface').sort((a,b)=>a.r-b.r); const rows=[['magnetic axis (R, z)',`${res.axis[0].toFixed(3)}, ${res.axis[1].toFixed(3)}`],['seeds on closed surfaces',`${surf.length} / ${res.lines.length}`],['islands / chaotic',`${res.lines.filter(l=>l.kind==='island/chaotic').length}`],['last closed surface r',surf.length?surf[surf.length-1].r.toFixed(3):'—'],['ι axis → edge',surf.length?`${surf[0].iota.toFixed(3)} → ${surf[surf.length-1].iota.toFixed(3)}`:'—']];
  $('ro-poinc').innerHTML='<table>'+rows.map(r=>`<tr><td>${r[0]}</td><td>${r[1]}</td></tr>`).join('')+(surf.length?`<tr><td colspan=2 style="color:var(--muted)">ι profile: ${surf.map(l=>`${l.r.toFixed(2)}:${l.iota.toFixed(2)}`).join(' · ')}</td></tr>`:'')+'</table>'; }
function computeField(){ if(!state.circuits.length)buildGeometry(); const n=parseInt($('gridn').value); const a=parseFloat($('wirer').value); pending={n,half:0.87,a};
  setStatus(`computing ${n}³ Biot–Savart grid…`,'busy'); post({type:'field',circuits:state.circuits.map(c=>({pts:Float64Array.from(c.pts),sign:c.sign})),n,half:0.87,a}); }
function seedsInBall(count,rng,D,n,half){ const dist=D?makeInterpScalar(D,n,half):null; const out=[]; let guard=0;
  while(out.length<count*3&&guard<100000){guard++; const x=rng()*1.16-0.58,y=rng()*1.16-0.58,z=rng()*1.16-0.58; const r=vlen(x,y,z); if(r<0.18||r>0.58)continue; if(dist&&dist(x,y,z)<0.035)continue; out.push(x,y,z);} return Float64Array.from(out); }
function mulberry(seed){ let a=seed>>>0; return function(){ a|=0; a=a+0x6D2B79F5|0; let t=Math.imul(a^a>>>15,1|a); t=t+Math.imul(t^t>>>7,61|t)^t; return ((t^t>>>14)>>>0)/4294967296; }; }
function traceFieldLines(){ if(!state.grid){setStatus('compute the field first','bad');return;} const g=state.grid; const rng=mulberry(parseInt($('seed').value)+7); const seeds=seedsInBall(parseInt($('nlines').value),rng,g.D,g.n,g.half);
  setStatus('tracing field lines…','busy'); post({type:'lines',B:g.B,D:g.D,n:g.n,half:g.half,seeds,wall:0.82,step:0.006,maxLen:40,wireClear:Math.max(0.008,g.a)}); }
function runParticlesUI(){ if(!state.grid){setStatus('compute the field first','bad');return;} const g=state.grid; const N=parseInt($('npart').value); const rng=mulberry(parseInt($('seed').value));
  const x0=seedsInBall(N,rng,g.D,g.n,g.half); const speed=parseFloat($('speed').value); const v0=new Float64Array(N*3); const q=new Float64Array(N); const qs=$('qsign').value;
  for(let i=0;i<N;i++){let dx,dy,dz,m; do{dx=rng()*2-1;dy=rng()*2-1;dz=rng()*2-1;m=vlen(dx,dy,dz);}while(m>1||m<1e-3); v0[3*i]=speed*dx/m;v0[3*i+1]=speed*dy/m;v0[3*i+2]=speed*dz/m; q[i]=qs==='0'?(i%2?-1:1):parseFloat(qs);}
  // normalise the field to B_rms = 1 in the ROI so Ω = 1 there (code units): scale B by 1/brms
  const gain=8/g.stats.brms; /* Codex normalisation: B_rms = 8 in the ROI */ const Bn=new Float32Array(g.B.length); for(let i=0;i<Bn.length;i++)Bn[i]=g.B[i]*gain;
  const transits=parseFloat($('transits').value); const tMax=transits*0.82/speed; const omdt=parseFloat($('omdt').value); const dt=Math.min(omdt/8.0, tMax/400);
  const sampleEvery=Math.max(1,Math.floor(tMax/dt/240));
  setStatus(`pushing ${N} particles for ${transits} transits…`,'busy'); state.partParams={speed,tMax,dt,transits};
  post({type:'particles',B:Bn,D:g.D,n:g.n,half:g.half,params:{x0,v0,q,tMax,dt,wall:0.82,wireR:g.a,omdt,sampleEvery}},[Bn.buffer]); }

// ---------- three.js scene ----------
const view=$('view'); const renderer=new THREE.WebGLRenderer({antialias:true,alpha:true}); renderer.setPixelRatio(Math.min(window.devicePixelRatio,2)); view.appendChild(renderer.domElement);
const scene=new THREE.Scene(); const camera=new THREE.PerspectiveCamera(38,1,0.01,100); let cam={theta:0.7,phi:1.1,r:3.2}; let autoRot=true;
function placeCamera(){ camera.position.set(cam.r*Math.sin(cam.phi)*Math.cos(cam.theta),cam.r*Math.sin(cam.phi)*Math.sin(cam.theta),cam.r*Math.cos(cam.phi)); camera.up.set(0,0,1); camera.lookAt(0,0,0); }
function resize(){ const w=view.clientWidth,h=view.clientHeight; renderer.setSize(w,h,false); camera.aspect=w/h; camera.updateProjectionMatrix(); } window.addEventListener('resize',resize);
let drag=null; view.addEventListener('mousedown',e=>{drag={x:e.clientX,y:e.clientY};autoRot=false;$('turntable').checked=false;}); window.addEventListener('mouseup',()=>drag=null);
window.addEventListener('mousemove',e=>{if(!drag)return; cam.theta-=(e.clientX-drag.x)*0.008; cam.phi=Math.min(Math.max(cam.phi-(e.clientY-drag.y)*0.008,0.05),3.09); drag={x:e.clientX,y:e.clientY};});
view.addEventListener('wheel',e=>{e.preventDefault(); cam.r=Math.min(Math.max(cam.r*Math.exp(e.deltaY*0.001),0.6),12);},{passive:false});
view.addEventListener('dblclick',()=>{cam={theta:0.7,phi:1.1,r:3.2};});
$('turntable').addEventListener('change',e=>autoRot=e.target.checked);
const wallGeo=new THREE.SphereGeometry(0.82,48,32); const wallMesh=new THREE.Mesh(wallGeo,new THREE.MeshBasicMaterial({color:0x2a3a48,wireframe:true,transparent:true,opacity:0.18})); scene.add(wallMesh);
$('showwall').addEventListener('change',e=>wallMesh.visible=e.target.checked);
const groups={cond:new THREE.Group(),ret:new THREE.Group(),lines:new THREE.Group(),part:new THREE.Group(),trails:new THREE.Group(),slice:new THREE.Group()}; Object.values(groups).forEach(g=>scene.add(g));
function clearGroup(g){ while(g.children.length){const c=g.children.pop(); if(c.geometry)c.geometry.dispose(); if(c.material)c.material.dispose();} }
function lineFrom(pts,color,opacity){ const geo=new THREE.BufferGeometry(); geo.setAttribute('position',new THREE.BufferAttribute(Float32Array.from(pts),3)); return new THREE.Line(geo,new THREE.LineBasicMaterial({color,transparent:opacity<1,opacity})); }
function drawConductors(){ clearGroup(groups.cond); clearGroup(groups.ret);
  state.active.forEach(c=>groups.cond.add(lineFrom(c.pts,c.sign>0?0xd4884a:0x6fa8dc,1)));
  state.returns.forEach((r,i)=>groups.ret.add(lineFrom(r,state.active[i].sign>0?0x7a5a3a:0x3a5a7a,0.5))); groups.ret.visible=$('showreturn').checked; }
$('showreturn').addEventListener('change',e=>groups.ret.visible=e.target.checked);
function clearLines(){ clearGroup(groups.lines); }
function drawLines(res){ clearLines(); for(const l of res){ for(const d of l.dirs){ const color=d.end==='wall'?0xe05d57:d.end==='wire'?0xe6b34f:d.end==='closed'?0x6fc47c:0x45b8a6; groups.lines.add(lineFrom(d.path,color,0.55)); } } }
let sliceMesh=null;
function drawSlice(g){ clearGroup(groups.slice); sliceMesh=null; if(!g||!$('showslice').checked)return; const n=g.n,half=g.half; const cv=document.createElement('canvas'); cv.width=cv.height=n; const ctx=cv.getContext('2d'); const img=ctx.createImageData(n,n);
  const k=(n-1)>>1; let lo=Infinity,hi=-Infinity; const vals=new Float32Array(n*n); for(let i=0;i<n;i++)for(let j=0;j<n;j++){const id=3*((i*n+j)*n+k); const b=Math.log10(vlen(g.B[id],g.B[id+1],g.B[id+2])+1e-12); vals[j*n+i]=b; if(b<lo)lo=b; if(b>hi)hi=b;}
  hi=Math.min(hi,lo+3.5); for(let p=0;p<n*n;p++){const t=Math.min(Math.max((vals[p]-lo)/(hi-lo),0),1); // teal->copper ramp
    const r=Math.round(20+t*t*215),gg=Math.round(60+t*110-t*t*40),b=Math.round(80+(1-t)*90); img.data[4*p]=r;img.data[4*p+1]=gg;img.data[4*p+2]=b;img.data[4*p+3]=150;}
  ctx.putImageData(img,0,0); const tex=new THREE.CanvasTexture(cv); tex.magFilter=THREE.LinearFilter; const mesh=new THREE.Mesh(new THREE.PlaneGeometry(2*half,2*half),new THREE.MeshBasicMaterial({map:tex,transparent:true,side:THREE.DoubleSide,depthWrite:false})); groups.slice.add(mesh); sliceMesh=mesh; }
$('showslice').addEventListener('change',()=>drawSlice(state.grid));
// particles
let pb=null; // playback state
function clearParticles(){ clearGroup(groups.part); clearGroup(groups.trails); pb=null; $('scrub').max=0; $('tnow').textContent='0.000'; $('tmax').textContent='0'; $('alive').textContent='–'; drawSurvival(null); }
function setupPlayback(res){ clearGroup(groups.part); clearGroup(groups.trails); const N=res.N; const geo=new THREE.BufferGeometry(); geo.setAttribute('position',new THREE.BufferAttribute(new Float32Array(N*3),3)); geo.setAttribute('color',new THREE.BufferAttribute(new Float32Array(N*3),3));
  const pts=new THREE.Points(geo,new THREE.PointsMaterial({size:0.03,vertexColors:true,sizeAttenuation:true})); groups.part.add(pts);
  // trails: one line per particle for the first 40 particles
  const nt=Math.min(N,60); const trails=[]; for(let i=0;i<nt;i++){const g=new THREE.BufferGeometry(); g.setAttribute('position',new THREE.BufferAttribute(new Float32Array(res.frames.length*3),3)); g.setDrawRange(0,0); const l=new THREE.Line(g,new THREE.LineBasicMaterial({color:i%2?0x7ad9ca:0xf0a868,transparent:true,opacity:0.5})); groups.trails.add(l); trails.push(l);}
  pb={res,frame:0,playing:true,pts,trails,nt,lastT:performance.now()}; $('scrub').max=res.frames.length-1; $('scrub').value=0; $('tmax').textContent=res.tMax.toFixed(2); $('play').textContent='⏸ Pause'; drawSurvival(res); showFrame(0); }
function showFrame(f){ if(!pb)return; const {res,pts,trails,nt}=pb; pb.frame=f; const pos=pts.geometry.attributes.position.array, col=pts.geometry.attributes.color.array; const fr=res.frames[f], al=res.aliveFrames[f]; let alive=0;
  for(let i=0;i<res.N;i++){ if(al[i]){pos[3*i]=fr[3*i];pos[3*i+1]=fr[3*i+1];pos[3*i+2]=fr[3*i+2]; col[3*i]=0.94;col[3*i+1]=0.66;col[3*i+2]=0.41; alive++;} else { const k=res.lossK[i]; pos[3*i]=fr[3*i];pos[3*i+1]=fr[3*i+1];pos[3*i+2]=fr[3*i+2]; if(k===1){col[3*i]=0.88;col[3*i+1]=0.36;col[3*i+2]=0.34;} else if(k===2){col[3*i]=0.9;col[3*i+1]=0.7;col[3*i+2]=0.31;} else {col[3*i]=0.54;col[3*i+1]=0.48;col[3*i+2]=0.84;} } }
  pts.geometry.attributes.position.needsUpdate=true; pts.geometry.attributes.color.needsUpdate=true;
  if($('trails').checked){ for(let i=0;i<nt;i++){ const arr=trails[i].geometry.attributes.position.array; let n=0; for(let k=0;k<=f;k++){ if(!res.aliveFrames[k][i]&&k>0&&!res.aliveFrames[k-1][i])break; arr[3*n]=res.frames[k][3*i];arr[3*n+1]=res.frames[k][3*i+1];arr[3*n+2]=res.frames[k][3*i+2]; n++; } trails[i].geometry.setDrawRange(0,n); trails[i].geometry.attributes.position.needsUpdate=true; } groups.trails.visible=true; } else groups.trails.visible=false;
  $('scrub').value=f; $('tnow').textContent=res.times[f].toFixed(3); $('alive').textContent=`${alive}/${res.N}`; drawSurvival(res,f); }
$('scrub').addEventListener('input',e=>{ if(pb){pb.playing=false;$('play').textContent='▶ Play';showFrame(parseInt(e.target.value));} });
$('play').addEventListener('click',()=>{ if(!pb)return; pb.playing=!pb.playing; $('play').textContent=pb.playing?'⏸ Pause':'▶ Play'; if(pb.playing&&pb.frame>=pb.res.frames.length-1)pb.frame=0; });
$('clearp').addEventListener('click',clearParticles); $('trails').addEventListener('change',()=>{ if(pb)showFrame(pb.frame); });
function drawSurvival(res,f){ const cv=$('survival'); const ctx=cv.getContext('2d'); const W=cv.width,H=cv.height; ctx.clearRect(0,0,W,H); ctx.fillStyle='#0b1015'; ctx.fillRect(0,0,W,H); ctx.strokeStyle='#25303a'; ctx.lineWidth=1; for(let k=0;k<=4;k++){const y=8+(H-20)*k/4; ctx.beginPath();ctx.moveTo(34,y);ctx.lineTo(W-6,y);ctx.stroke();}
  ctx.fillStyle='#8695a0'; ctx.font='16px IBM Plex Mono, monospace'; ctx.fillText('1',10,16); ctx.fillText('0',10,H-8); ctx.fillText('survival',40,H-4);
  if(!res)return; const T=res.tMax; const pts=[]; for(let k=0;k<res.times.length;k++){let a=0; const al=res.aliveFrames[k]; for(let i=0;i<res.N;i++)a+=al[i]; pts.push([34+(W-40)*res.times[k]/T, 8+(H-20)*(1-a/res.N)]);}
  ctx.strokeStyle='#45b8a6'; ctx.lineWidth=2.5; ctx.beginPath(); pts.forEach((p,i)=>i?ctx.lineTo(p[0],p[1]):ctx.moveTo(p[0],p[1])); ctx.stroke();
  if(f!==undefined){ctx.strokeStyle='#d4884a'; ctx.lineWidth=2; ctx.beginPath(); const x=34+(W-40)*res.times[f]/T; ctx.moveTo(x,6);ctx.lineTo(x,H-12);ctx.stroke();} }
let fpsT=performance.now(),fpsN=0;
function animate(now){ requestAnimationFrame(animate); if(autoRot)cam.theta+=0.0035; placeCamera(); if(pb&&pb.playing){ const dtms=now-pb.lastT; if(dtms>33){ pb.lastT=now; let f=pb.frame+1; if(f>=pb.res.frames.length){f=0;} showFrame(f);} } renderer.render(scene,camera); fpsN++; if(now-fpsT>1000){$('fps').textContent=`${fpsN} fps`;fpsN=0;fpsT=now;} }

// ---------- readouts ----------
const SPECIES={code:{m:1,q:1,label:'code'},D:{m:2.01355*1.66053906660e-27,q:1.602176634e-19},e:{m:9.1093837015e-31,q:1.602176634e-19},alpha:{m:4.0015*1.66053906660e-27,q:2*1.602176634e-19},Ar:{m:39.948*1.66053906660e-27,q:1.602176634e-19}};
const SCALES={dimensionless:{L:1,I:1,a:0.005,label:'unit ball, B_rms = 1'},benchtop:{L:0.18,I:1000,a:0.001,label:'18 cm ball, 1 kA, 1 mm Cu'},reactor:{L:1.0,I:5.0*1.1/(1.25663706212e-6*12),a:0.01,label:'1 m ball, 364 kA, 1 cm conductor'}};
function updateGeoReadout(){ const nC=state.circuits.length; $('ro-geo').innerHTML=`<table><tr><td>configuration</td><td>${state.name}</td></tr><tr><td>circuits</td><td>${nC}</td></tr><tr><td>active length</td><td>${fmt(state.activeLength,2)} ball units</td></tr><tr><td>total incl. returns</td><td>${fmt(state.length,2)}</td></tr><tr><td>conductor clearance</td><td>${fmt(state.clearance,4)}</td></tr><tr><td>max conductor radius</td><td>${fmt(state.clearance*0.45,4)}</td></tr></table>`; }
function updateFieldReadout(){ const g=state.grid; if(!g){$('ro-field').innerHTML='<span class="pill">no field yet</span>';return;} const s=g.stats;
  $('ro-field').innerHTML=`<table><tr><td>grid</td><td>${g.n}³ · ±${g.half}</td></tr><tr><td>B_rms in ROI (per A, unit ball)</td><td>${s.brms.toExponential(3)} T</td></tr><tr><td>B_max on grid</td><td>${s.bmax.toExponential(3)} T</td></tr><tr><td>B_max / B_rms</td><td>${fmt(s.bmax/s.brms,1)}</td></tr><tr><td>core / shell |B| (median)</td><td>${fmt(s.bcore/s.bshell,2)}</td></tr><tr><td>Ampère ∮B·dl / μ₀I (exact field)</td><td>${s.exact.ampere.toFixed(5)} <span class="pill ${Math.abs(s.exact.ampere-1)<1e-3?'ok':'bad'}">${Math.abs(s.exact.ampere-1)<1e-3?'pass':'fail'}</span></td></tr><tr><td>∇·B, ∇×B (exact, rel)</td><td>${s.exact.divRelMax.toExponential(1)} · ${s.exact.curlRelMax.toExponential(1)}</td></tr></table>`; }
function updateLinesReadout(){ const L=state.lines; if(!L){$('ro-lines').innerHTML='<span class="pill">not traced</span>';return;} const n=L.length; const ends=[]; L.forEach(l=>l.dirs.forEach(d=>ends.push(d.end)));
  const wallBoth=L.filter(l=>l.dirs[0].end==='wall'&&l.dirs[1].end==='wall').length/n, closed=L.filter(l=>l.closed).length/n, wire=L.filter(l=>l.dirs.some(d=>d.end==='wire')).length/n;
  const R=L.map(l=>l.mirror).filter(isFinite).sort((a,b)=>a-b); const med=R[R.length>>1], p90=R[Math.floor(R.length*0.9)]; const pred=L.reduce((a,l)=>a+l.pred,0)/n;
  $('ro-lines').innerHTML=`<table><tr><td>seeds</td><td>${n}</td></tr><tr><td>reach wall both ways</td><td>${(100*wallBoth).toFixed(0)} %</td></tr><tr><td>hit a conductor</td><td>${(100*wire).toFixed(0)} %</td></tr><tr><td>closed / very long</td><td>${(100*closed).toFixed(0)} %</td></tr><tr><td>mirror ratio median · p90</td><td>${fmt(med,2)} · ${fmt(p90,2)}</td></tr><tr><td><b>predicted adiabatic retention</b></td><td><b>${(100*pred).toFixed(1)} %</b></td></tr></table><p class="note">√(1 − B_seed/B_mirror) per line, isotropic pitch. This is the ceiling for particles with ρ_L ≪ L_B; drifts then empty the trap on ~(L_B/ρ_L) transits.</p>`; }
function updatePartReadout(){ const r=state.particles; if(!r){$('ro-part').innerHTML='<span class="pill">no run</span>';return;} const p=state.partParams; const N=r.N; const z=1.96; const pr=r.retained/N; const den=1+z*z/N; const c=(pr+z*z/(2*N))/den; const h=z*Math.sqrt(pr*(1-pr)/N+z*z/(4*N*N))/den;
  const lostT=r.lossT.filter(isFinite).sort((a,b)=>a-b); const med=r.retained>N/2?`> ${r.tMax.toFixed(2)}`:fmt(lostT[Math.floor(N/2)-1]||r.tMax,2);
  $('ro-part').innerHTML=`<table><tr><td>particles</td><td>${N}</td></tr><tr><td>code speed · ρ_L at B_rms</td><td>${p.speed} · ${(p.speed/8).toFixed(4)}</td></tr><tr><td>window</td><td>${p.transits} transits = ${r.tMax.toFixed(2)}</td></tr><tr><td><b>retained</b></td><td><b>${(100*pr).toFixed(1)} %</b> [${(100*Math.max(0,c-h)).toFixed(0)}, ${(100*Math.min(1,c+h)).toFixed(0)}]</td></tr><tr><td>wall · conductor · grid</td><td>${r.nwall} · ${r.nwire} · ${r.nexit}</td></tr><tr><td>median loss time</td><td>${med}</td></tr><tr><td>μ variation (median)</td><td>${fmt(r.muVarMedian,3)}</td></tr><tr><td>max sub-steps</td><td>${r.subMax}</td></tr></table>`; }
function updateScaleReadout(){ const sc=SCALES[$('scale').value]; const sp=SPECIES[$('species').value]; const g=state.grid; let html=`<table><tr><td>preset</td><td>${sc.label}</td></tr>`;
  if(state.circuits.length){ const a=Math.min(sc.a,0.45*state.clearance*sc.L); const R=1.68e-8*state.length*sc.L/(Math.PI*a*a); html+=`<tr><td>conductor radius used</td><td>${(a*1e3).toFixed(2)} mm${a<sc.a?' <span class="pill warn">thinned to fit clearance</span>':''}</td></tr><tr><td>resistance (Cu)</td><td>${fmt(R,4)} Ω</td></tr><tr><td>Joule power</td><td>${(sc.I*sc.I*R).toExponential(2)} W</td></tr><tr><td>current density</td><td>${(sc.I/(Math.PI*a*a)).toExponential(2)} A/m²</td></tr>`; }
  if(g){ const brms=g.stats.brms*sc.I/sc.L, bmax=g.stats.bmax*sc.I/sc.L; html+=`<tr><td>B_rms (ROI)</td><td>${brms.toExponential(3)} T</td></tr><tr><td>B_max</td><td>${bmax.toExponential(3)} T</td></tr>`;
    if($('species').value!=='code'){ const E=parseFloat($('energy').value)*1.602176634e-19; const v=Math.sqrt(2*E/sp.m); const rho=sp.m*v/(sp.q*brms); const code=8*rho/sc.L; /* B_rms = 8 in code units */ html+=`<tr><td>${$('species').value} speed</td><td>${v.toExponential(2)} m/s</td></tr><tr><td>gyroradius at B_rms</td><td>${(rho*1e3).toFixed(2)} mm → code speed ${code.toExponential(2)}</td></tr><tr><td>gyroperiod at B_rms</td><td>${(2*Math.PI*sp.m/(sp.q*brms)).toExponential(2)} s</td></tr><tr><td>wall transit</td><td>${(0.82*sc.L/v).toExponential(2)} s</td></tr>`; $('speed').value=code.toPrecision(3); $('part-note').textContent=`Code speed set from ${$('species').value} at ${$('energy').value} eV in this field: ρ_L/R_ball = ${code.toExponential(2)}.`; } }
  html+='</table>'; $('ro-scale').innerHTML=html; }
['scale','species','energy'].forEach(id=>$(id).addEventListener('change',updateScaleReadout));
$('species').addEventListener('change',()=>{ const s=$('species').value; const sc=$('scale').value; if(s==='code')$('part-note').textContent=''; });
function updateLiveChecks(){ const g=state.grid; if(!g)return; const e=g.stats.exact; $('livechecks').innerHTML=`<span>∇·B · L / |B| (exact field, FD)</span><span>${e.divRelMax.toExponential(2)}</span><span>∇×B · L / |B| outside conductor</span><span>${e.curlRelMax.toExponential(2)}</span><span>Ampère loop / μ₀I</span><span>${e.ampere.toFixed(6)}</span><span>check points</span><span>${e.nPoints}</span><span>trilinear grid ∇·B (resolution indicator)</span><span>${g.stats.divRelMax.toExponential(2)}</span><span>grid</span><span>${g.n}³</span><span>conductor model</span><span>uniform current, a = ${g.a}</span>`; }

// ---------- ledger from precomputed Python results ----------
function buildLedger(){ const rows=DATA.configs; const el=$('ledger'); const bar=v=>`<span class="bar" style="width:${Math.round(60*(v||0))}px"></span>`;
  el.innerHTML=`<table><thead><tr><th>configuration</th><th>clear.</th><th>R_mir p90</th><th>predicted</th><th>e⁻ S4</th><th>D⁺ S4</th><th>D⁺ S20</th><th>dimless</th><th>B_rms bench</th><th>P bench</th></tr></thead><tbody>`+rows.map((r,i)=>`<tr data-i="${i}"><td>${r.name.replace(' (10 circuits)','')}</td><td>${fmt(r.clearance,3)}</td><td>${fmt(r.topology.mirror_ratio_p90,1)}</td><td>${bar(r.topology.predicted_adiabatic_retention)}${(100*r.topology.predicted_adiabatic_retention).toFixed(0)}%</td><td>${bar(r.S4.bench_e)}${(100*r.S4.bench_e).toFixed(0)}%</td><td>${bar(r.S4.reactor_D)}${(100*r.S4.reactor_D).toFixed(0)}%</td><td>${(100*r.S20.reactor_D).toFixed(0)}%</td><td>${(100*r.retention_end.dimless).toFixed(0)}%</td><td>${(1e3*r.scales.benchtop.B_rms).toFixed(1)} mT</td><td>${(r.scales.benchtop.joule_W/1e3).toFixed(0)} kW</td></tr>`).join('')+'</tbody></table>';
  el.querySelectorAll('tr[data-i]').forEach(tr=>tr.addEventListener('click',()=>{ loadFromLedger(rows[parseInt(tr.dataset.i)]); el.querySelectorAll('tr').forEach(t=>t.classList.remove('sel')); tr.classList.add('sel'); }));
  const A=DATA.assemblies||[]; $('assemblies').innerHTML=`<table><thead><tr><th>configuration</th><th>predicted</th><th>dimless</th><th>S4</th><th>S8</th><th>S20</th></tr></thead><tbody>`+A.map((r,i)=>`<tr data-a="${i}"><td>${r.name.replace(' (10 circuits)','')}</td><td>${bar(r.pred)}${(100*r.pred).toFixed(0)}%</td><td>${(100*r.dimless).toFixed(0)}%</td><td>${bar(r.S4)}${(100*r.S4).toFixed(0)}%</td><td>${(100*r.S8).toFixed(0)}%</td><td>${(100*r.S20).toFixed(0)}%</td></tr>`).join('')+'</tbody></table>';
  $('assemblies').querySelectorAll('tr[data-a]').forEach(tr=>tr.addEventListener('click',()=>{ const r=A[parseInt(tr.dataset.a)]; const n=r.name; const fam=$('family');
    if(n.includes('precess')){fam.value='precess';renderParams();} else if(n.startsWith('hopf continued')){fam.value='hopfcont';renderParams();const m=n.match(/(\d+)° \/ (\d+)/); if(m){$('p_sweep').value=m[1];$('p_circuits').value=m[2];} if(n.includes('0.50'))$('p_eta1').value=0.5;} else if(n.startsWith('hopf torus')){fam.value='hopftorus';renderParams();} else if(n.includes('hopf-drift')||n.includes('Shallow 180')){fam.value='hopfdrift';renderParams();$('p_variant').value='s64_180';} else {fam.value='solenoid';renderParams();}
    const lay=n.includes(' × ')?n.split(' × ')[1]:'single'; $('layout').value=LAYOUTS[lay]?lay:'single'; buildGeometry(); document.querySelector('.tabs button[data-tab=readout]').click(); }));
  const PB=(DATA.playbook||[]).slice().sort((a,b)=>b.S20-a.S20); $('playbook').innerHTML=`<table><thead><tr><th>configuration</th><th>closed</th><th>S4</th><th>S8</th><th>S20</th><th>core</th><th>B_rms</th></tr></thead><tbody>`+PB.map((r,i)=>`<tr data-p="${i}"><td>${r.name}</td><td>${(100*r.closed).toFixed(0)}%</td><td>${bar(r.S4)}${(100*r.S4).toFixed(0)}%</td><td>${(100*r.S8).toFixed(0)}%</td><td>${bar(r.S20)}${(100*r.S20).toFixed(0)}%</td><td>${(100*r.core).toFixed(0)}%</td><td>${r.brms.toFixed(1)} T</td></tr>`).join('')+'</tbody></table>';
  $('playbook').querySelectorAll('tr[data-p]').forEach(tr=>tr.addEventListener('click',()=>{ const r=PB[parseInt(tr.dataset.p)]; if(!r.ui)return; const fam=$('family'); fam.value=r.ui.family; renderParams(); for(const [k,v] of Object.entries(r.ui.params||{})){const el=$('p_'+k); if(el)el.value=v;} $('layout').value=LAYOUTS[r.ui.layout]?r.ui.layout:'single'; if(r.ui.corescale)$('corescale').value=r.ui.corescale; buildGeometry(); $('playbook').querySelectorAll('tr').forEach(t=>t.classList.remove('sel')); tr.classList.add('sel'); document.querySelector('.tabs button[data-tab=readout]').click(); }));
  const SF=DATA.surfaces||[]; $('surfaces').innerHTML=`<table><thead><tr><th>configuration</th><th>surfaces</th><th>r_out</th><th>ι axis→edge</th><th>well</th></tr></thead><tbody>`+SF.map((r,i)=>`<tr data-s="${i}"><td>${r.name}</td><td>${bar(r.surf)}${(100*r.surf).toFixed(0)}%</td><td>${r.r_out.toFixed(3)}</td><td>${isFinite(r.iota_axis)?r.iota_axis.toFixed(3)+' → '+r.iota_edge.toFixed(3):'—'}</td><td>${isFinite(r.well)?(r.well>0?'+':'')+r.well.toFixed(2):'—'}</td></tr>`).join('')+'</tbody></table>';
  $('surfaces').querySelectorAll('tr[data-s]').forEach(tr=>tr.addEventListener('click',()=>{ const r=SF[parseInt(tr.dataset.s)]; if(!r.ui||!r.ui.family)return; $('family').value=r.ui.family; renderParams(); for(const [k,v] of Object.entries(r.ui.params||{})){const el=$('p_'+k); if(el)el.value=v;} $('layout').value='single'; buildGeometry(); $('surfaces').querySelectorAll('tr').forEach(t=>t.classList.remove('sel')); tr.classList.add('sel'); document.querySelector('.tabs button[data-tab=readout]').click(); }));
  const TR=DATA.transport||[]; if(TR.length) $('transport').innerHTML=`<table><thead><tr><th>run (guiding centre, exact field)</th><th>ν/transit</th><th>S20</th><th>S40</th><th>S60</th></tr></thead><tbody>`+TR.map(r=>`<tr><td>${r.name}</td><td>${r.nu}</td><td>${(100*r.S20).toFixed(0)}%</td><td>${(100*r.S40).toFixed(0)}%</td><td>${bar(r.S60)}${(100*r.S60).toFixed(0)}%</td></tr>`).join('')+'</tbody></table>';
  const sw=DATA.sweep.sweep; $('sweep').innerHTML=`<table><thead><tr><th>speed</th><th>ρ_L/L_B</th><th>retention</th><th>μ var</th></tr></thead><tbody>`+sw.map(s=>`<tr><td>${s.speed}</td><td>${s.adiabaticity.toExponential(2)}</td><td>${bar(s.retention)}${(100*s.retention).toFixed(1)}%</td><td>${s.mu_variation.toFixed(2)}</td></tr>`).join('')+`</tbody></table><p class="note">adiabatic prediction ${(100*DATA.sweep.predicted_adiabatic).toFixed(1)} % for ${DATA.sweep.winding}</p>`;
  const cx=rows.filter(r=>r.codex); $('codex').innerHTML=`<table><thead><tr><th>variant</th><th>Codex</th><th>ccsim</th></tr></thead><tbody>`+cx.map(r=>`<tr><td>${r.name.replace('hopf-drift ','').replace(' (10 circuits)','')}</td><td>${(100*r.codex.retention).toFixed(1)}% [${(100*r.codex.ci[0]).toFixed(0)},${(100*r.codex.ci[1]).toFixed(0)}]</td><td>${(100*r.codex.ccsim).toFixed(1)}% [${(100*r.codex.ccsim_wilson[0]).toFixed(0)},${(100*r.codex.ccsim_wilson[1]).toFixed(0)}]</td></tr>`).join('')+'</tbody></table>';
  const v=DATA.verify; const kv=$('verify'); const items=[['loop off-axis (Biot–Savart vs elliptic)',v.fields.loop_offaxis_rel_max],['solenoid on-axis',v.fields.solenoid_axis_rel_max],['straight wire',v.fields.straight_wire_rel_max],['conductor interior',v.fields.conductor_interior_rel_max],['∇·B on Hopf winding',v.maxwell.div_B_rel_max],['∇×B outside conductor',v.maxwell.curl_B_rel_max],['∇×A − B',v.maxwell.curlA_minus_B_rel_max],['Ampère loop',v.maxwell.ampere_rel_err_max],['Faraday EMF',v.faraday],['gyroradius',v.particles.gyroradius_rel_err],['gyroperiod',v.particles.gyroperiod_rel_err],['μ conservation (uniform B)',v.particles.mu_variation],['mirror loss-cone agreement',1-v.mirror.per_particle_agreement],['∇B drift',v.gradB.rel_err],['D–T ⟨σv⟩ 10 keV vs table',v.fusion.DT_10keV_rel_err],['beam–target MC vs nσv',Math.abs(v.fusion.beam_target_mc_vs_n_sigma_v_rel)],['EEDF integral',v.chemistry.eedf_constant_sigma_rel_err]]; if(v.sphere)items.splice(4,0,['sphere winding centre B (closed form)',v.sphere.B_centre_rel_err]);
  kv.innerHTML=items.map(([k,val])=>`<span>${k}</span><span>${val.toExponential(1)}</span>`).join(''); }
function loadFromLedger(r){ const n=r.name; const fam=$('family');
  if(n.startsWith('hopf-drift')){fam.value='hopfdrift';renderParams();$('p_variant').value={'Shallow 150°':'s64_150','Shallow 180°':'s64_180','Shallow 210°':'s64_210','Deep 150°':'d48_150','Deep 180°':'d48_180','Deep 210°':'d48_210'}[n.split('hopf-drift ')[1].split(' (')[0]];}
  else if(n.startsWith('codex phase')){fam.value='phase';renderParams();}
  else if(n.startsWith('codex precess')){fam.value='precess';renderParams();}
  else if(n.startsWith('codex Hopf')){fam.value='codexhopf';renderParams();}
  else if(n.startsWith('recursive-L2')){fam.value='recursive';renderParams();$('p_r2').value=0;}
  else if(n.startsWith('recursive-L3')){fam.value='recursive';renderParams();if(n.includes('alt'))$('p_alt').value='opposite';}
  else {fam.value='solenoid';renderParams();}
  $('layout').value='single'; buildGeometry(); document.querySelector('.tabs button[data-tab=readout]').click(); }

// ---------- Designs: the design loop's ledger, live feeds, score card ----------
const DESIGNS=new Map();           // id → compact row (see ui/ledger_rows.py)
let feed={mode:'snapshot',label:'',lines:0,newIds:new Set()};
function addRows(rows,isNew){ let added=0; for(const r of rows||[]){ if(!r||!r.id)continue; if(!DESIGNS.has(r.id)){added++; if(isNew)feed.newIds.add(r.id);} DESIGNS.set(r.id,r);} return added; }
function setFeed(mode,label,live){ feed.mode=mode; feed.label=label; const el=$('feed-status'); el.textContent=label; el.className='feed'+(live?' live':''); $('livedot').className='livedot'+(live?' on':''); }
function rowsSorted(){ const key=$('dsort').value; const q=$('dfilter').value.trim().toLowerCase(); let rows=[...DESIGNS.values()];
  if($('dlatest').checked){ const byHash=new Map(); for(const r of rows){ const o=byHash.get(r.hash); if(!o||(r.t||'')>(o.t||''))byHash.set(r.hash,r);} const counts={}; for(const r of rows)counts[r.hash]=(counts[r.hash]||0)+1; rows=[...byHash.values()].map(r=>Object.assign({},r,{evals:counts[r.hash]})); }
  if(q) rows=rows.filter(r=>((r.name||'')+' '+(r.fam||'')+' '+(r.hash||'')).toLowerCase().includes(q));
  const num=v=>(v===null||v===undefined||isNaN(v))?-Infinity:v;
  rows.sort(key==='t'?((a,b)=>(b.t||'').localeCompare(a.t||'')):((a,b)=>num(b[key])-num(a[key])||(b.t||'').localeCompare(a.t||'')));
  return rows; }
function renderDesigns(){ const rows=rowsSorted(); const el=$('designs'); const bar=v=>`<span class="bar" style="width:${Math.round(60*Math.max(0,Math.min(1,v||0)))}px"></span>`; const best=Math.max(1e-9,...[...DESIGNS.values()].map(r=>r.score||0));
  const pct=v=>(v===null||v===undefined)?'–':(100*v).toFixed(0)+'%'; const f1=v=>(v===null||v===undefined)?'–':(+v).toFixed(1);
  el.innerHTML=`<table><thead><tr><th>#</th><th>design</th><th>score</th><th>τ_all</th><th>τ_pass</th><th>build</th><th>S20</th></tr></thead><tbody>`+rows.map((r,i)=>`<tr data-id="${r.id}" class="${feed.newIds.has(r.id)?'new':''}${state.designRow&&state.designRow.id===r.id?' sel':''}"><td>${i+1}</td><td title="${(r.name||'')+' — '+(r.fam||'')+(r.closed!==undefined&&r.closed!==null?' · closed lines '+(100*r.closed).toFixed(0)+'%':'')+(r.I_kA?' · '+Math.round(r.I_kA)+' kA for 2 T':'')+(r.S_end!==undefined&&r.S_end!==null?' · S_end '+(100*r.S_end).toFixed(0)+'% at '+r.T+' transits':'')+' · '+(r.t||'').replace('T',' ').replace('Z','')}">${r.name||'?'}<br><span style="color:var(--muted)">${r.fam||''}${r.evals>1?`<span class="badge">×${r.evals}</span>`:''}${r.valid===false?'<span class="badge" style="background:var(--wall)">invalid</span>':''}${r.score_plasma!==undefined&&r.score_plasma!==null?`<span class="badge" title="score_plasma (with --surfaces)">plasma ${f1(r.score_plasma)}</span>`:''}</span></td><td>${bar((r.score||0)/best)}${f1(r.score)}</td><td>${f1(r.tau_all)}</td><td>${f1(r.tau_pass)}</td><td>${r.build===undefined||r.build===null?'–':(+r.build).toFixed(2)}</td><td>${pct(r.S20)}</td></tr>`).join('')+'</tbody></table>';
  el.querySelectorAll('tr[data-id]').forEach(tr=>tr.addEventListener('click',()=>loadDesignRow(DESIGNS.get(tr.dataset.id))));
  feed.newIds.clear(); drawScoreHistory(); }
function loadDesignRow(r){ if(!r||!r.design)return; state.designRow=r; state.designSpec=r.design; $('family').value='design'; renderParams(); $('layout').value='single'; buildGeometry(); $('designs').querySelectorAll('tr').forEach(t=>t.classList.toggle('sel',t.dataset.id===r.id)); document.querySelector('.tabs button[data-tab=readout]').click(); }
function drawScoreHistory(){ const cv=$('score-history'); const ctx=cv.getContext('2d'); const W=cv.width,H=cv.height; ctx.fillStyle='#0b1015'; ctx.fillRect(0,0,W,H); const rows=[...DESIGNS.values()].filter(r=>r.t).sort((a,b)=>a.t.localeCompare(b.t)); if(!rows.length)return;
  const n=rows.length; const cur=rows.filter(r=>(r.T||0)>=40).map(r=>r.score||0); const smax=Math.max(1,1.08*(cur.length?Math.max(...cur):Math.max(...rows.map(r=>r.score||0)))); const X=i=>28+(W-40)*(n>1?i/(n-1):0.5), Y=s=>H-16-(H-28)*Math.min(1,(s||0)/smax);
  ctx.strokeStyle='#25303a'; ctx.lineWidth=1; for(let k=0;k<=4;k++){const y=H-16-(H-28)*k/4; ctx.beginPath();ctx.moveTo(28,y);ctx.lineTo(W-10,y);ctx.stroke();}
  ctx.fillStyle='#8695a0'; ctx.font='11px IBM Plex Mono, monospace'; ctx.fillText(smax.toFixed(0),4,16); ctx.fillText('0',4,H-14); ctx.fillText(`score vs evaluation order (${n}; clipped to the current score's best — the early 300s were score v1)`,34,H-4);
  let best=0; ctx.strokeStyle='#e6b34f'; ctx.lineWidth=1.5; ctx.beginPath(); rows.forEach((r,i)=>{best=Math.max(best,Math.min(r.score||0,smax)); const x=X(i),y=Y(best); i?ctx.lineTo(x,y):ctx.moveTo(x,y);}); ctx.stroke();
  const famColor=f=>({hopfmesh:'#45b8a6',hopfmirror:'#7ad9ca',hopftorus:'#8a7bd6',precess:'#f0a868',yinyang:'#e05d57',solenoid:'#8695a0'})[(f||'').split(/[ ×+]/)[0]]||'#dde5ea';
  rows.forEach((r,i)=>{ctx.fillStyle=famColor(r.fam); const x=X(i),y=Y(r.score); ctx.fillRect(x-1.5,y-1.5,3,3); if(state.designRow&&state.designRow.id===r.id){ctx.strokeStyle='#dde5ea'; ctx.beginPath(); ctx.arc(x,y,4,0,2*Math.PI); ctx.stroke();}}); }
function updateDesignReadout(){ const r=state.designRow; const sec=$('sec-design'); if(!r){sec.hidden=true;return;} sec.hidden=false; const f1=v=>(v===null||v===undefined)?'–':(+v).toFixed(1); const pct=v=>(v===null||v===undefined)?'–':(100*v).toFixed(0)+' %';
  let rows=[['design',r.name],['components',r.fam],['<b>score</b>',`<b>${f1(r.score)}</b> = τ_c ${f1(r.tau_c)} × build ${r.build===undefined?'–':(+r.build).toFixed(2)}`],['τ_all · τ_pass (transits)',`${f1(r.tau_all)} · ${f1(r.tau_pass)}`],['S4 · S20 · S_end (T)',`${pct(r.S4)} · ${pct(r.S20)} · ${pct(r.S_end)} (${r.T})`],['closed lines · adiabatic pred.',`${pct(r.closed)} · ${pct(r.pred)}`],['current for 2 T (reactor)',r.I_kA?`${Math.round(r.I_kA)} kA (${(r.brms364||0).toFixed(2)} T at 364 kA)`:'–'],['clearance · conductor radius',`${r.clr===undefined?'–':(+r.clr).toFixed(4)} · ${r.wire===undefined?'–':(+r.wire).toFixed(4)}`],['markers · evaluated',`${r.n||'–'} · ${(r.t||'').replace('T',' ').replace('Z','')} (${r.secs?Math.round(r.secs)+' s':'–'})`]];
  if(r.surf) rows.push(['surfaces (Poincaré)',`${pct(r.surf.surf)} of seeds · r_out ${(r.surf.r_out||0).toFixed(3)} · ι ${(r.surf.iota_axis||0).toFixed(3)} → ${(r.surf.iota_edge||0).toFixed(3)} · well ${(r.surf.well||0).toFixed(2)}`],['<b>score_plasma</b>',`<b>${f1(r.score_plasma)}</b>${r.plasma?` (τ_pol ${f1(r.plasma.tau_pol)} transits; ι needed ${(r.plasma.iota_needed||0).toFixed(3)}; ${r.plasma.ok?'surfaces shorted — exempt':'capped'})`:''}`]);
  if(r.valid===false) rows.push(['<span style="color:var(--wall)">invalid</span>',(r.problems||[]).join('; ')]);
  $('ro-design').innerHTML='<table>'+rows.map(x=>`<tr><td>${x[0]}</td><td>${x[1]}</td></tr>`).join('')+'</table><p class="note">Numbers from the Python evaluation (ccsim.evaluate, numba Biot–Savart + Boris, 160 markers unless stated); the geometry here is rebuilt from the same design JSON.</p>'; }
function followLatest(){ if(!$('follow').checked||state.busy)return; const rows=[...DESIGNS.values()].filter(r=>r.t&&r.design).sort((a,b)=>b.t.localeCompare(a.t)); if(rows.length&&(!state.designRow||state.designRow.id!==rows[0].id))loadDesignRow(rows[0]); }
async function initFeeds(){
  addRows(DATA.designs||[],false); setFeed('snapshot',`built-in snapshot: ${DESIGNS.size} evaluations (${DATA.designs_built||'?'})`,false); renderDesigns();
  // 1. the published page: shared database written by the design loop's runs (artifact `db` capability)
  if(window.claude&&typeof window.claude.use==='function'){ try{ const db=await window.claude.use('db'); if(db){
      db.collection('designs').orderBy('t','desc').limit(600).onSnapshot(snap=>{ const rows=snap.docs.filter(d=>d.exists).map(d=>d.data()); const added=addRows(rows,!snap.metadata.fromCache); setFeed('db',`live feed: ${DESIGNS.size} evaluations · ${rows.length} in the shared feed${snap.metadata.fromCache?' (cached)':''}`,true); renderDesigns(); if(added)followLatest(); },
        e=>setFeed('db','live feed unavailable ('+(e&&e.code||'error')+') — showing the built-in snapshot',false));
      db.doc('status/run').onSnapshot(snap=>{ const st=snap.exists?snap.data():null; const el=$('run-status'); if(st&&st.label){el.hidden=false; el.textContent=`${st.running?'● running':'○ finished'}: ${st.label}${st.best?' · best '+st.best:''}${st.note?' — '+st.note:''} (${(st.updated||'').replace('T',' ').replace('Z','')})`;} else el.hidden=true; },()=>{});
      return; } }catch(e){} }
  // 2. served locally by ui/serve.py: poll the ledger file
  if(location.protocol.startsWith('http')){ try{ const r=await fetch('/api/ledger?after=0',{cache:'no-store'}); if(r.ok){ const j=await r.json(); addRows(j.rows,false); feed.lines=j.count; setFeed('local',`local ledger: ${DESIGNS.size} evaluations · polling results/design_ledger.jsonl (${feed.lines} lines)`,true); renderDesigns();
        const poll=async()=>{ try{ const r2=await fetch(`/api/ledger?after=${feed.lines}`,{cache:'no-store'}); if(r2.ok){ const j2=await r2.json(); if(j2.count!==feed.lines||j2.rows.length){ feed.lines=j2.count; const added=addRows(j2.rows,true); setFeed('local',`local ledger: ${DESIGNS.size} evaluations · live (${feed.lines} lines)`,true); if(added){renderDesigns(); followLatest();} } }
            const s=await fetch('/api/status',{cache:'no-store'}); if(s.ok){const st=await s.json(); const el=$('run-status'); if(st&&st.label){el.hidden=false; el.textContent=`${st.running?'● running':'○ finished'}: ${st.label}${st.best?' · best '+st.best:''}${st.note?' — '+st.note:''}`;} else el.hidden=true;} }catch(e){} setTimeout(poll,3000); }; setTimeout(poll,3000); return; } }catch(e){} }
}
['dsort','dfilter','dlatest'].forEach(id=>$(id).addEventListener(id==='dfilter'?'input':'change',renderDesigns));
$('dload').addEventListener('click',()=>{ try{ const spec=JSON.parse($('djson').value); state.designRow=null; state.designSpec=spec; $('family').value='design'; renderParams(); $('layout').value='single'; buildGeometry(); document.querySelector('.tabs button[data-tab=readout]').click(); }catch(e){setStatus('design JSON: '+e.message,'bad');} });
$('dcopy').addEventListener('click',()=>{ const spec=state.designSpec; if(!spec){setStatus('no design loaded','bad');return;} $('djson').value=JSON.stringify(spec,null,1); try{navigator.clipboard.writeText($('djson').value);}catch(e){} setStatus('design JSON copied','ok'); });
$('follow').addEventListener('change',followLatest);

// ---------- wiring ----------
document.querySelectorAll('.tabs button').forEach(b=>b.addEventListener('click',()=>{document.querySelectorAll('.tabs button').forEach(x=>x.classList.remove('on'));b.classList.add('on');document.querySelectorAll('.tabpane').forEach(p=>p.classList.remove('on'));$('tab-'+b.dataset.tab).classList.add('on');}));
Object.keys(LAYOUTS).forEach(k=>{const o=document.createElement('option');o.value=k;o.textContent=k;$('layout').appendChild(o);});
$('family').addEventListener('change',renderParams); $('build').addEventListener('click',buildGeometry); $('compute').addEventListener('click',computeField); $('lines').addEventListener('click',traceFieldLines); $('run').addEventListener('click',runParticlesUI); $('poincare').addEventListener('click',runPoincare);
renderParams(); buildLedger(); initFeeds(); resize(); placeCamera(); requestAnimationFrame(animate);
// Study links select the recorded JSON through the existing design loader.
// Geometry preview does not recompute or replace the study's external metrics.
const studySelection = new URLSearchParams(location.search).get('study');
if(studySelection){
  $('follow').checked=false;
  setStatus('loading study geometry','busy');
  fetch('/results/confinement_study.json',{cache:'no-store'}).then(r=>{if(!r.ok)throw new Error('Study data unavailable');return r.json();}).then(study=>{
    const row=(study.evaluations||[]).find(r=>String(r.id)===studySelection);
    if(!row||!row.design||!Array.isArray(row.design.components))throw new Error('Study candidate not found');
    $('djson').value=JSON.stringify(row.design,null,2);
    $('dload').click();
  }).catch(error=>setStatus(error.message,'bad'));
}else{
  // Ordinary opening state retains the existing default simulation.
  buildGeometry(); computeField();
  const autoChain=ev=>{ if(ev.data.type==='field'){ setTimeout(()=>{traceFieldLines();},50); } else if(ev.data.type==='lines'&&!state.particles){ setTimeout(()=>{ runParticlesUI(); },50); worker.removeEventListener('message',autoChain);} };
  worker.addEventListener('message',autoChain);
}
