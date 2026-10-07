// ===================== numerics shared by main thread and worker =====================
const ENGINE_SRC = String.raw`
function vlen(x,y,z){return Math.sqrt(x*x+y*y+z*z);}
// --- polyline utilities ---
function resample(pts, count){
  const n=pts.length/3; const s=new Float64Array(n); s[0]=0;
  for(let i=1;i<n;i++){const dx=pts[3*i]-pts[3*i-3],dy=pts[3*i+1]-pts[3*i-2],dz=pts[3*i+2]-pts[3*i-1]; s[i]=s[i-1]+vlen(dx,dy,dz);}
  const L=s[n-1]; const out=new Float64Array(count*3); let j=0;
  for(let k=0;k<count;k++){const q=L*k/(count-1); while(j<n-2 && s[j+1]<q) j++; const seg=s[j+1]-s[j]; const w=seg>0?(q-s[j])/seg:0;
    for(let c=0;c<3;c++) out[3*k+c]=pts[3*j+c]*(1-w)+pts[3*(j+1)+c]*w;}
  return out;
}
function polyLength(p){let L=0;for(let i=3;i<p.length;i+=3)L+=vlen(p[i]-p[i-3],p[i+1]-p[i-2],p[i+2]-p[i-1]);return L;}
function smoothstep5(u){return u*u*u*(10+u*(-15+6*u));}
function unit3(v){const n=vlen(v[0],v[1],v[2])||1e-300;return [v[0]/n,v[1]/n,v[2]/n];}
function cross(a,b){return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];}
function dot(a,b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
function bishopFrame(p){
  const n=p.length/3; const t=new Array(n),n1=new Array(n),n2=new Array(n);
  for(let i=0;i<n-1;i++) t[i]=unit3([p[3*i+3]-p[3*i],p[3*i+4]-p[3*i+1],p[3*i+5]-p[3*i+2]]); t[n-1]=t[n-2];
  for(let i=n-2;i>=1;i--) t[i]=unit3([t[i-1][0]+t[i][0],t[i-1][1]+t[i][1],t[i-1][2]+t[i][2]]);
  let trial=Math.abs(t[0][2])<0.9?[0,0,1]:[1,0,0]; const d0=dot(trial,t[0]);
  n1[0]=unit3([trial[0]-d0*t[0][0],trial[1]-d0*t[0][1],trial[2]-d0*t[0][2]]);
  for(let i=0;i<n-1;i++){
    const v1=[p[3*i+3]-p[3*i],p[3*i+4]-p[3*i+1],p[3*i+5]-p[3*i+2]]; const c1=dot(v1,v1);
    if(c1<1e-30){n1[i+1]=n1[i];continue;}
    const k=2/c1; const dn=dot(v1,n1[i]),dt=dot(v1,t[i]);
    const rL=[n1[i][0]-k*dn*v1[0],n1[i][1]-k*dn*v1[1],n1[i][2]-k*dn*v1[2]];
    const tL=[t[i][0]-k*dt*v1[0],t[i][1]-k*dt*v1[1],t[i][2]-k*dt*v1[2]];
    const v2=[t[i+1][0]-tL[0],t[i+1][1]-tL[1],t[i+1][2]-tL[2]]; const c2=dot(v2,v2);
    let r=rL; if(c2>=1e-30){const k2=2/c2,dr=dot(v2,rL); r=[rL[0]-k2*dr*v2[0],rL[1]-k2*dr*v2[1],rL[2]-k2*dr*v2[2]];}
    const dd=dot(r,t[i+1]); n1[i+1]=unit3([r[0]-dd*t[i+1][0],r[1]-dd*t[i+1][1],r[2]-dd*t[i+1][2]]);
  }
  for(let i=0;i<n;i++) n2[i]=cross(t[i],n1[i]);
  return {t,n1,n2};
}
function helixAround(parent, radius, turns, nPerTurn, phase, hand){
  const n=Math.max(Math.floor(turns*nPerTurn)+1,8); const base=resample(parent,n); const f=bishopFrame(base);
  const out=new Float64Array(n*3);
  for(let i=0;i<n;i++){const a=phase+hand*2*Math.PI*turns*i/(n-1); const c=Math.cos(a)*radius,s=Math.sin(a)*radius;
    for(let k=0;k<3;k++) out[3*i+k]=base[3*i+k]+c*f.n1[i][k]+s*f.n2[i][k];}
  return out;
}
function circleOpen(R, n, turns0){const th=2*Math.PI*(1-0.5/turns0); const out=new Float64Array(n*3);
  for(let i=0;i<n;i++){const a=th*i/(n-1); out[3*i]=R*Math.cos(a);out[3*i+1]=R*Math.sin(a);out[3*i+2]=0;} return out;}
function recursiveWinding(baseR, levels, nPerTurn, alt){
  let curve=circleOpen(baseR, Math.max(360,24*Math.floor(levels[0][1])), levels[0][1]); let radius=baseR, turns=1;
  levels.forEach((lv,k)=>{radius*=lv[0]; turns*=lv[1]; curve=helixAround(curve,radius,turns,nPerTurn,0,(alt&&k%2===1)?-1:1);});
  return curve;
}
function solenoidPts(R,L,turns,nPerTurn){const n=turns*nPerTurn+1; const out=new Float64Array(n*3);
  for(let i=0;i<n;i++){const th=2*Math.PI*turns*i/(n-1); out[3*i]=R*Math.cos(th);out[3*i+1]=R*Math.sin(th);out[3*i+2]=-L/2+L*i/(n-1);} return out;}
function rotZ(p,ang){const c=Math.cos(ang),s=Math.sin(ang);return [c*p[0]-s*p[1],s*p[0]+c*p[1],p[2]];}
function rotY(p,ang){const c=Math.cos(ang),s=Math.sin(ang);return [c*p[0]+s*p[2],p[1],-s*p[0]+c*p[2]];}
function normalizeBall(p){let m=0;for(let i=0;i<p.length;i+=3)m=Math.max(m,vlen(p[i],p[i+1],p[i+2]));for(let i=0;i<p.length;i++)p[i]/=m;return p;}
function codexToroidal(kind, inward, rotationDeg, slipDeg){
  const N=10,K=11,nhi=8193; const raw=new Float64Array(nhi*3);
  for(let i=0;i<nhi;i++){const u=i/(nhi-1),h=smoothstep5(u),loops=N*u,th=2*Math.PI*loops;
    const r=4.45+((4.45-inward*N)-4.45)*h, env=0.52+(0.52*0.88-0.52)*h, psi=K*th+slipDeg*Math.PI/180*loops;
    let v=[(r+env*Math.cos(psi))*Math.cos(th),(r+env*Math.cos(psi))*Math.sin(th),env*Math.sin(psi)];
    if(kind==='precess'){const chi=rotationDeg*Math.PI/180*loops; v=rotZ(rotY(rotZ(v,-chi),28*Math.PI/180),chi);}
    raw[3*i]=v[0];raw[3*i+1]=v[1];raw[3*i+2]=v[2];}
  return normalizeBall(resample(raw,1801));
}
function hopfPt(eta,phi,tau){const a1=tau+phi/2,a2=tau-phi/2; const q1=Math.cos(eta)*Math.cos(a1),q2=Math.cos(eta)*Math.sin(a1),q3=Math.sin(eta)*Math.cos(a2),q4=Math.sin(eta)*Math.sin(a2);
  const den=1-q4; if(den<=1e-3) return null; return [q1/den,q2/den,q3/den];}
function codexHopf(inward, rotationDeg){
  const N=10,nhi=8193; const raw=new Float64Array(nhi*3); const eta0=0.82, eta1=Math.max(0.34,eta0-(0.16+0.17*inward));
  for(let i=0;i<nhi;i++){const u=i/(nhi-1),h=smoothstep5(u),loops=N*u; const p=hopfPt(eta0+(eta1-eta0)*h,0.18+rotationDeg*Math.PI/180*loops,2*Math.PI*loops); if(!p) throw new Error('pole');
    raw[3*i]=p[0];raw[3*i+1]=p[1];raw[3*i+2]=p[2];}
  return normalizeBall(resample(raw,1801));
}
function hopfArclengthEta(eta0,eta1,sweepRad,nhi){ // v(u) that makes the (eta,phi) drift arclength-uniform on S3
  const M=4097; const vt=new Float64Array(M),cum=new Float64Array(M); const de=eta1-eta0;
  for(let i=0;i<M;i++){vt[i]=i/(M-1);}
  let acc=0; cum[0]=0;
  for(let i=1;i<M;i++){const e0=eta0+de*vt[i-1],e1=eta0+de*vt[i]; const s0=Math.sqrt(4*de*de+Math.pow(Math.sin(2*e0)*sweepRad,2)),s1=Math.sqrt(4*de*de+Math.pow(Math.sin(2*e1)*sweepRad,2)); acc+=0.5*(s0+s1)*(vt[i]-vt[i-1]); cum[i]=acc;}
  const out=new Float64Array(nhi); let j=0;
  for(let i=0;i<nhi;i++){const u=i/(nhi-1)*acc; while(j<M-2&&cum[j+1]<u)j++; const w=(cum[j+1]-cum[j])>0?(u-cum[j])/(cum[j+1]-cum[j]):0; out[i]=vt[j]*(1-w)+vt[j+1]*w;}
  return out;
}
function hopfDriftPts(eta0,eta1,sweepDeg,circuits,tau0,scaleByEta0){
  const nhi=Math.max(4097,circuits*400+1); const v=hopfArclengthEta(eta0,eta1,sweepDeg*Math.PI/180,nhi); const raw=new Float64Array(nhi*3);
  for(let i=0;i<nhi;i++){const u=i/(nhi-1); const p=hopfPt(eta0+(eta1-eta0)*v[i],sweepDeg*Math.PI/180*v[i],tau0+2*Math.PI*circuits*u); if(!p) throw new Error('stereographic pole reached — keep eta below ~0.9');
    raw[3*i]=p[0];raw[3*i+1]=p[1];raw[3*i+2]=p[2];}
  const sc=scaleByEta0?(1-Math.sin(eta0))/Math.cos(eta0):1; for(let i=0;i<raw.length;i++)raw[i]*=sc;
  const out=resample(raw,circuits*150+1); return scaleByEta0?out:normalizeBall(out);
}
function hopfTorusPts(eta,circuits,revs){
  const nhi=Math.max(4097,circuits*400+1); const raw=new Float64Array(nhi*3);
  for(let i=0;i<nhi;i++){const u=i/(nhi-1); const p=hopfPt(eta,2*Math.PI*revs*u,Math.PI/2+2*Math.PI*circuits*u); if(!p) throw new Error('pole'); raw[3*i]=p[0];raw[3*i+1]=p[1];raw[3*i+2]=p[2];}
  return normalizeBall(resample(raw,circuits*150+1));
}
function weaveProfile(uc, sgn, u){ // smooth ±1 profile through the crossing values (cosine ramps)
  const idx=Array.from(uc.keys()).sort((a,b)=>uc[a]-uc[b]); const us=idx.map(i=>uc[i]), ss=idx.map(i=>sgn[i]);
  const uu=[us[0]-(us[1]-us[0]),...us,us[us.length-1]+(us[us.length-1]-us[us.length-2])]; const s2=[-ss[0],...ss,-ss[ss.length-1]];
  const out=new Float64Array(u.length); let k=0;
  for(let i=0;i<u.length;i++){while(k<uu.length-2&&uu[k+1]<=u[i])k++; const t=Math.min(1,Math.max(0,(u[i]-uu[k])/(uu[k+1]-uu[k]))); out[i]=s2[k]+(s2[k+1]-s2[k])*(1-Math.cos(Math.PI*t))/2;}
  return out;
}
function hopfMirrorPair(eta, circuits, revs, mode, dEta){ // Hopf torus + its z-mirror (anti-fibres); woven = interlocked at every crossing
  const N=circuits, r=revs, nhi=Math.max(4097,N*400+1); const u=new Float64Array(nhi); for(let i=0;i<nhi;i++)u[i]=i/(nhi-1);
  const curve=(etaArr,mir)=>{const raw=new Float64Array(nhi*3); for(let i=0;i<nhi;i++){const p=hopfPt(etaArr[i],2*Math.PI*r*u[i],Math.PI/2+2*Math.PI*N*u[i]); if(!p) throw new Error('pole: keep η below ~0.9'); raw[3*i]=p[0];raw[3*i+1]=p[1];raw[3*i+2]=mir?-p[2]:p[2];} return raw;};
  let A,B;
  if(mode==='nested'){A=curve(new Float64Array(nhi).fill(eta+dEta),false); B=curve(new Float64Array(nhi).fill(eta-dEta),true);}
  else { const p=N+r/2,q=N-r/2; const uc=[],vc=[],sg=[];
    for(let j=-Math.floor(p)-2;j<=Math.floor(p)+2;j++)for(let k=-1;k<=Math.floor(2*q)+2;k++){const uu=(j/p+k/q)/2,vv=(k/q-j/p)/2; if(uu>=0&&uu<=1&&vv>=0&&vv<=1){uc.push(uu);vc.push(vv);sg.push(k%2?-1:1);}}
    const sA=weaveProfile(uc,sg,u), sB=weaveProfile(vc,sg.map(s=>-s),u);
    const eA=new Float64Array(nhi), eB=new Float64Array(nhi); for(let i=0;i<nhi;i++){eA[i]=eta+dEta*sA[i]; eB[i]=eta+dEta*sB[i];}
    A=curve(eA,false); B=curve(eB,true); }
  let m=0; for(let i=0;i<A.length;i+=3){m=Math.max(m,vlen(A[i],A[i+1],A[i+2]),vlen(B[i],B[i+1],B[i+2]));}
  for(let i=0;i<A.length;i++){A[i]/=m;B[i]/=m;}
  return [resample(A,N*150+1),resample(B,N*150+1)];
}
function weaveProfilePeriodic(uc, sgn, u){ // periodic (closed strand) version: crossings repeated at u±1
  const idx=Array.from(uc.keys()).sort((a,b)=>uc[a]-uc[b]); const us=idx.map(i=>uc[i]), ss=idx.map(i=>sgn[i]);
  const uu=[...us.map(x=>x-1),...us,...us.map(x=>x+1)], s2=[...ss,...ss,...ss];
  const out=new Float64Array(u.length); let k=0;
  for(let i=0;i<u.length;i++){while(k<uu.length-2&&uu[k+1]<=u[i])k++; const t=Math.min(1,Math.max(0,(u[i]-uu[k])/(uu[k+1]-uu[k]))); out[i]=s2[k]+(s2[k+1]-s2[k])*(1-Math.cos(Math.PI*t))/2;}
  return out;
}
function modulatedStrand(p, q, eps, n, phase, handed, nhi){ // turns = level sets of a2 + (eps/2) sin(2a2 - handed n a1 + phase) - (q/p) a1 (fixed point iterated to convergence, as in ccsim)
  const a1=new Float64Array(nhi), a2=new Float64Array(nhi); const L=2*Math.PI*p;
  for(let i=0;i<nhi;i++){a1[i]=L*i/(nhi-1); let y=(q/p)*a1[i]; for(let it=0;it<2000;it++){const yn=(q/p)*a1[i]-0.5*eps*Math.sin(2*y-handed*n*a1[i]+phase); const d=Math.abs(yn-y); y=yn; if(d<1e-13)break;} a2[i]=y;}
  return [a1,a2];
}
function hopfMeshPair(eta, N, eps, periods, dEta, helicalMode){ // meshed torus II: closed (N+1,N-1) torus knots, exact rotated mirror, optional chiral density modulation (same algorithm and sampling as ccsim.geometry.hopf_helical_pair)
  const p=N+1, q=N-1, alpha=Math.PI, nhi=20000, L=2*Math.PI*p;
  const [a1A,a2A]=modulatedStrand(p,q,eps,periods,0,+1,nhi);
  const epsB=helicalMode==='poloidal'?eps:-eps; const [a1B,a2B]=modulatedStrand(p,q,epsB,periods,-periods*alpha,-1,nhi);
  const interpB=t=>{ // periodic extension a2B(t+L)=a2B(t)+2πq
    let m=Math.floor(t/L); let tt=t-m*L; const x=tt/L*(nhi-1); const i=Math.min(nhi-2,Math.max(0,Math.floor(x))); const f=x-i; return a2B[i]+(a2B[i+1]-a2B[i])*f+m*2*Math.PI*q; };
  // scan A one sample past each end (periodic extension) so a crossing exactly at u = 0 is never missed; branches m and m+p coincide, so scan the p distinct ones
  const a1s=new Float64Array(nhi+2), a2s=new Float64Array(nhi+2); a1s[0]=a1A[nhi-2]-L; a2s[0]=a2A[nhi-2]-2*Math.PI*q; a1s[nhi+1]=a1A[1]+L; a2s[nhi+1]=a2A[1]+2*Math.PI*q; for(let i=0;i<nhi;i++){a1s[i+1]=a1A[i];a2s[i+1]=a2A[i];}
  const uc=[], vc=[];
  for(let m=0;m<Math.round(p);m++){
    let Dprev=null;
    for(let i=0;i<nhi+2;i++){const t=a1s[i]-alpha+2*Math.PI*m; let D=(a2s[i]+interpB(t))%(2*Math.PI); if(D<0)D+=2*Math.PI; D-=Math.PI;
      if(i>0&&Math.sign(D)!==Math.sign(Dprev)&&Math.abs(D-Dprev)<Math.PI){const f=Dprev/(Dprev-D); const s=a1s[i-1]+f*(a1s[i]-a1s[i-1]); let u=(s/L)%1; if(u<0)u+=1; if(u>1-1e-9)u=0; let v=((s-alpha+2*Math.PI*m)/L)%1; if(v<0)v+=1; uc.push(u);vc.push(v);}
      Dprev=D;}
  }
  // dedupe on u: merge anything closer than 1e-6 (crossings are ~1/(2N²) apart), including across the wrap
  const order=Array.from(uc.keys()).sort((a,b)=>uc[a]-uc[b]); const U=[],V=[]; let last=-1;
  for(const i of order){if(uc[i]-last>1e-6){U.push(uc[i]);V.push(vc[i]);last=uc[i];}}
  if(U.length>1&&U[0]+1-U[U.length-1]<1e-6){U.pop();V.pop();}
  if(U.length%2) throw new Error('odd crossing count '+U.length+': weave parity broken');
  const sign=U.map((_,i)=>(i%2?-1:1)); // plain-weave checkerboard: alternate by rank along A
  const ca=Math.cos(alpha), sa=Math.sin(alpha);
  const curveA=etaOf=>{const P=new Float64Array(nhi*3); for(let i=0;i<nhi;i++){const pa=hopfPt(etaOf(i), a1A[i]-a2A[i], 0.5*(a1A[i]+a2A[i])); if(!pa) throw new Error('pole: keep η below ~0.9'); P[3*i]=pa[0];P[3*i+1]=pa[1];P[3*i+2]=pa[2];} return P;};
  const curveB=etaOf=>{const P=new Float64Array(nhi*3); for(let i=0;i<nhi;i++){const pb=hopfPt(etaOf(i), a1B[i]-a2B[i], 0.5*(a1B[i]+a2B[i])); if(!pb) throw new Error('pole'); P[3*i]=ca*pb[0]-sa*pb[1]; P[3*i+1]=sa*pb[0]+ca*pb[1]; P[3*i+2]=-pb[2];} return P;};
  const A0=curveA(()=>eta), B0=curveB(()=>eta); const h=1e-4; const Ap=curveA(()=>eta+h), Am=curveA(()=>eta-h), Bp=curveB(()=>eta+h), Bm=curveB(()=>eta-h);
  const rateA=new Float64Array(nhi), rateB=new Float64Array(nhi); for(let i=0;i<nhi;i++){rateA[i]=vlen(Ap[3*i]-Am[3*i],Ap[3*i+1]-Am[3*i+1],Ap[3*i+2]-Am[3*i+2])/(2*h); rateB[i]=vlen(Bp[3*i]-Bm[3*i],Bp[3*i+1]-Bm[3*i+1],Bp[3*i+2]-Bm[3*i+2])/(2*h);}
  const zw=zoneWeave(A0,B0,rateA,rateB,dEta,U,V,sign,nhi);
  const A=curveA(i=>eta+dEta*zw.sA[i]), B=curveB(i=>eta+dEta*zw.sB[i]);
  let mx=0; for(let i=0;i<A.length;i+=3){mx=Math.max(mx,vlen(A[i],A[i+1],A[i+2]),vlen(B[i],B[i+1],B[i+2]));}
  for(let i=0;i<A.length;i++){A[i]/=mx;B[i]/=mx;}
  const rA=resample(A,N*150+1), rB=resample(B,N*150+1);
  for(const P of [rA,rB]){P[P.length-3]=P[0];P[P.length-2]=P[1];P[P.length-1]=P[2];}
  return {A:rA,B:rB,crossings:U.length,defects:zw.defects,zoneFraction:zw.zoneFraction,components:zw.components,weave:zw.zoneFraction>0.99?'nested':(zw.zoneFraction<0.5?'woven':'mixed')};
}
// ---- zone weave (ccsim.geometry._zone_weave): profiles that keep two strands on one surface apart everywhere ----
function gridIndex(P,cell){ const m=new Map(); const key=(x,y,z)=>x+','+y+','+z; for(let i=0;i<P.length/3;i++){const k=key(Math.floor(P[3*i]/cell),Math.floor(P[3*i+1]/cell),Math.floor(P[3*i+2]/cell)); let a=m.get(k); if(!a){a=[];m.set(k,a);} a.push(i);} return {m,cell,key}; }
function ballQuery(g,P,x,y,z,r,out){ // indices of P within r of (x,y,z); r must be ≤ g.cell
  out.length=0; const cx=Math.floor(x/g.cell),cy=Math.floor(y/g.cell),cz=Math.floor(z/g.cell); const r2=r*r;
  for(let dx=-1;dx<=1;dx++)for(let dy=-1;dy<=1;dy++)for(let dz=-1;dz<=1;dz++){const a=g.m.get(g.key(cx+dx,cy+dy,cz+dz)); if(!a)continue; for(const j of a){const ex=P[3*j]-x,ey=P[3*j+1]-y,ez=P[3*j+2]-z; if(ex*ex+ey*ey+ez*ez<r2)out.push(j);}} return out; }
function runsPeriodic(mask){ const n=mask.length; let all=true,any=false; for(let i=0;i<n;i++){if(mask[i])any=true; else all=false;} if(all)return [[0,n]]; if(!any)return [];
  let k=0; while(mask[k])k++; const runs=[]; let start=-1; for(let t=0;t<n;t++){const v=mask[(t+k)%n]; if(v&&start<0)start=t; if(!v&&start>=0){runs.push([start+k,t+k]);start=-1;}} if(start>=0)runs.push([start+k,n+k]); return runs; }
function zoneWeave(A0,B0,rateA,rateB,dEta,uc,vc,sign,n){
  const margin=1.15; const hA=new Float64Array(n), hB=new Float64Array(n); let hmax=0; for(let i=0;i<n;i++){hA[i]=2*Math.abs(dEta)*rateA[i]; hB[i]=2*Math.abs(dEta)*rateB[i]; hmax=Math.max(hmax,hA[i],hB[i]);}
  const cell=Math.max(1e-6,margin*hmax); const gA=gridIndex(A0,cell), gB=gridIndex(B0,cell); const buf=[];
  const maskA=new Uint8Array(n), maskB=new Uint8Array(n);
  for(let i=0;i<n;i++){ maskA[i]=ballQuery(gB,B0,A0[3*i],A0[3*i+1],A0[3*i+2],margin*hA[i],buf).length?1:0; }
  for(let j=0;j<n;j++){ maskB[j]=ballQuery(gA,A0,B0[3*j],B0[3*j+1],B0[3*j+2],margin*hB[j],buf).length?1:0; }
  const zonesA=runsPeriodic(maskA), zonesB=runsPeriodic(maskB); const zidA=new Int32Array(n).fill(-1), zidB=new Int32Array(n).fill(-1);
  zonesA.forEach(([a,b],k)=>{for(let t=a;t<b;t++)zidA[t%n]=k;}); zonesB.forEach(([a,b],k)=>{for(let t=a;t<b;t++)zidB[t%n]=k;});
  const nA=zonesA.length, nB=zonesB.length; const parent=new Int32Array(nA+nB); for(let i=0;i<nA+nB;i++)parent[i]=i;
  const find=x=>{while(parent[x]!==x){parent[x]=parent[parent[x]]; x=parent[x];} return x;}; const union=(x,y)=>{const rx=find(x),ry=find(y); if(rx!==ry)parent[Math.max(rx,ry)]=Math.min(rx,ry);};
  for(let i=0;i<n;i++){ if(zidA[i]<0)continue; ballQuery(gB,B0,A0[3*i],A0[3*i+1],A0[3*i+2],hA[i],buf); for(const j of buf){if(zidB[j]>=0)union(zidA[i],nA+zidB[j]);} }
  for(let j=0;j<n;j++){ if(zidB[j]<0)continue; ballQuery(gA,A0,B0[3*j],B0[3*j+1],B0[3*j+2],hB[j],buf); for(const i of buf){if(zidA[i]>=0)union(zidA[i],nA+zidB[j]);} }
  const compSign=new Map(); const order=Array.from(uc.keys()).sort((a,b)=>uc[a]-uc[b]); const idxOf=u=>Math.min(n-1,Math.round(u*(n-1)));
  for(const k of order){ const za=zidA[idxOf(uc[k])]; if(za<0)continue; const c=find(za); if(!compSign.has(c))compSign.set(c,sign[k]); }
  let defects=0; const crossingSign=uc.map((u,k)=>{const za=zidA[idxOf(u)]; const s=za>=0&&compSign.has(find(za))?compSign.get(find(za)):sign[k]; if(s!==sign[k])defects++; return s;});
  function profile(zones,compOf,factor){ const zs=new Float64Array(zones.length); let prev=1; for(let k=0;k<zones.length;k++){const c=compOf(k); if(compSign.has(c))zs[k]=compSign.get(c); else {zs[k]=prev; compSign.set(c,prev);} prev=zs[k];}
    const out=new Float64Array(n).fill(factor); if(!zones.length)return out; if(zones.length===1&&zones[0][1]-zones[0][0]>=n){out.fill(zs[0]*factor);return out;}
    for(let k=0;k<zones.length;k++){const [a,b]=zones[k]; for(let t=a;t<b;t++)out[t%n]=zs[k]*factor; let a2=zones[(k+1)%zones.length][0]; if(a2<=b)a2+=n; const s0=zs[k],s1=zs[(k+1)%zones.length]; const len=a2-b;
      for(let t=b;t<a2;t++){const tt=(t-b+1)/(len+1); out[t%n]=(s0!==s1?s0+(s1-s0)*(1-Math.cos(Math.PI*tt))/2:s0)*factor;}}
    return out; }
  const sA=profile(zonesA,k=>find(k),1), sB=profile(zonesB,k=>find(nA+k),-1);
  let zoneCount=0; for(let i=0;i<n;i++)if(zidA[i]>=0)zoneCount++; const comps=new Set(); for(let x=0;x<nA+nB;x++)comps.add(find(x));
  return {sA,sB,defects,zoneFraction:zoneCount/n,components:comps.size,crossingSign};
}
function baseballSeam(radius, ampDeg, turns, pitch){ // latitude A·sin(2t), longitude t, radius stepping inward each turn
  const A=ampDeg*Math.PI/180, n=turns*720+1, out=new Float64Array(n*3);
  for(let i=0;i<n;i++){const t=2*Math.PI*turns*i/(n-1); const lam=A*Math.sin(2*t); const r=radius-pitch*t/(2*Math.PI); out[3*i]=r*Math.cos(lam)*Math.cos(t); out[3*i+1]=r*Math.cos(lam)*Math.sin(t); out[3*i+2]=r*Math.sin(lam);}
  return out;
}
function circlePts(R,n,cx,cy,cz,axis){const out=new Float64Array(n*3); for(let i=0;i<n;i++){const t=2*Math.PI*i/(n-1); const c=R*Math.cos(t),s=R*Math.sin(t); // same orientation conventions as ccsim.geometry.circular_loop
  if(axis==='z'){out[3*i]=cx+c;out[3*i+1]=cy+s;out[3*i+2]=cz;} else if(axis==='y'){out[3*i]=cx+s;out[3*i+1]=cy;out[3*i+2]=cz+c;} else {out[3*i]=cx;out[3*i+1]=cy+c;out[3*i+2]=cz+s;}} return out;}
function conePts(rb,rt,h,turns,nPerTurn){ // 'tornado' coil: conical helix, base at -h/2, tip at +h/2 (ccsim.geometry.cone_helix)
  const n=turns*nPerTurn+1, out=new Float64Array(n*3); for(let i=0;i<n;i++){const t=i/(n-1); const r=rb+(rt-rb)*t, a=2*Math.PI*turns*t; out[3*i]=r*Math.cos(a);out[3*i+1]=r*Math.sin(a);out[3*i+2]=-h/2+h*t;} return out; }
function torusHelicalWindings(R0,r,l,periods,phaseDeg){
  // Match ccsim.geometry.torus_helical_windings: 2*gcd(|n|,l) distinct
  // closed coils, each spanning l/gcd(|n|,l) toroidal turns at 720 samples/turn.
  l=Math.trunc(l); periods=Math.trunc(periods);
  if(!Number.isFinite(l)||l<1||!Number.isFinite(periods)||![R0,r,phaseDeg].every(Number.isFinite))throw new Error('invalid torushelix parameters');
  let a=Math.abs(periods),b=l; while(b){const remainder=a%b;a=b;b=remainder;}
  const g=a,m=l/g,n=m*720+1,phase=phaseDeg*Math.PI/180,parts=[];
  for(let k=0;k<2*g;k++){
    const pts=new Float64Array(n*3);
    for(let i=0;i<n;i++){
      const phi=(2*Math.PI*m/(n-1))*i,theta=(periods/l)*phi+k*Math.PI/l+phase;
      const rho=R0+r*Math.cos(theta);
      pts[3*i]=rho*Math.cos(phi);pts[3*i+1]=rho*Math.sin(phi);pts[3*i+2]=r*Math.sin(theta);
    }
    pts[n*3-3]=pts[0];pts[n*3-2]=pts[1];pts[n*3-1]=pts[2];
    parts.push({pts,sign:k%2?-1:1,closed:true});
  }
  return parts;
}
// ---- declarative designs (ccsim.design): components with Python parameter names, transform scale → reflect → rotate → translate ----
var DESIGN_DEFAULTS={circle:{radius:0.7,axis:'z'},solenoid:{radius:0.45,length:1.1,turns:12},precess:{inward:0.22,rotation_deg:24,slip_deg:3},phase:{inward:0.32,slip_deg:3},codexhopf:{inward:0.18,rotation_deg:20},
  hopfdrift:{variant:'s64_180',circuits:10},hopfcont:{eta0:0.82,eta1:0.64,sweep_deg:720,circuits:48},hopftorus:{eta:0.70,circuits:24,revolutions:1.0},recursive:{base_radius:0.72,levels:[[0.24,12]],alternate:false},
  baseball:{radius:0.72,amplitude_deg:40,turns:6,radial_pitch:0.035},sphere:{radius:0.76,turns:24},hopfmirror:{eta:0.70,circuits:24,revolutions:2.0,mode:'woven',delta_eta:0.03,sense:-1},
  hopfmesh:{eta:0.60,circuits:24,eps:0.70,periods:4,delta_eta:0.03,helical_mode:'toroidal',sense:-1},yinyang:{radius_outer:0.76,radius_inner:0.60,amplitude_deg:40,turns:4,sense:1},picket:{rings:5,radius:0.72,alternating:true},
  link:{radius:0.62,sense:1},cone:{radius_base:0.32,radius_tip:0.06,height:0.5,turns:8},torushelix:{R0:0.60,r:0.40,l:2,periods:4,phase_deg:0.0}};
var DESIGN_VARIANTS={s64_150:[0.64,150],s64_180:[0.64,180],s64_210:[0.64,210],d48_150:[0.48,150],d48_180:[0.48,180],d48_210:[0.48,210]};
function designFamilyParts(family, params){ // → [{pts,sign,closed}] in the unit ball, mirroring ccsim.design._parts
  const d=Object.assign({},DESIGN_DEFAULTS[family]||{},params||{}); const one=(pts,closed)=>[{pts,sign:1,closed:!!closed}];
  switch(family){
    case 'circle': return one(circlePts(d.radius,721,0,0,0,d.axis||'z'),true);
    case 'solenoid': return one(solenoidPts(d.radius,d.length,Math.round(d.turns),60),false);
    case 'precess': return one(codexToroidal('precess',d.inward,d.rotation_deg,d.slip_deg),false);
    case 'phase': return one(codexToroidal('phase',d.inward,3,d.slip_deg),false);
    case 'codexhopf': return one(codexHopf(d.inward,d.rotation_deg),false);
    case 'hopfdrift': {const [e1,sw]=DESIGN_VARIANTS[d.variant]||DESIGN_VARIANTS.s64_180; return one(hopfDriftPts(0.82,e1,sw,Math.round(d.circuits),Math.PI/2,true),false);}
    case 'hopfcont': return one(hopfDriftPts(d.eta0,d.eta1,d.sweep_deg,Math.round(d.circuits),Math.PI/2,false),false);
    case 'hopftorus': {const rv=d.revolutions; const closed=Math.abs(rv-Math.round(rv))<1e-9&&Math.round(rv)%2===0; const pts=hopfTorusPts(d.eta,Math.round(d.circuits),rv); if(closed){pts[pts.length-3]=pts[0];pts[pts.length-2]=pts[1];pts[pts.length-1]=pts[2];} return one(pts,closed);}
    case 'recursive': {const lv=d.levels.map(l=>[+l[0],Math.round(l[1])]); return one(recursiveWinding(d.base_radius,lv,lv.length>1?24:48,!!d.alternate),false);}
    case 'baseball': return one(baseballSeam(d.radius,d.amplitude_deg,Math.round(d.turns),d.radial_pitch),false);
    case 'sphere': return one(sphereWinding(d.radius,Math.round(d.turns)),false);
    case 'hopfmirror': { // ccsim default scheme (τ₀ = 0, mirror rotated by π): closed (N+1, N−1) knots at 2 revolutions = the unmodulated meshed torus
      if(Math.abs(d.revolutions-2)<1e-9&&d.mode!=='nested'){const m=hopfMeshPair(d.eta,Math.round(d.circuits),0,1,d.delta_eta,'toroidal'); return [{pts:m.A,sign:1,closed:true},{pts:m.B,sign:Math.sign(d.sense)||-1,closed:true}];}
      const [A,B]=hopfMirrorPair(d.eta,Math.round(d.circuits),d.revolutions,d.mode,d.delta_eta); return [{pts:A,sign:1,closed:false,approx:true},{pts:B,sign:Math.sign(d.sense)||-1,closed:false,approx:true}]; }
    case 'hopfmesh': {const m=hopfMeshPair(d.eta,Math.round(d.circuits),d.eps,Math.round(d.periods),d.delta_eta,d.helical_mode||'toroidal'); return [{pts:m.A,sign:1,closed:true},{pts:m.B,sign:Math.sign(d.sense)||-1,closed:true}];}
    case 'yinyang': {const o=baseballSeam(d.radius_outer,d.amplitude_deg,Math.round(d.turns),0.03); const i=baseballSeam(d.radius_inner,d.amplitude_deg,Math.round(d.turns),0.03); const R=rotMat([0,0,1],Math.PI/2); const ir=new Float64Array(i.length); for(let k=0;k<i.length;k+=3){ir[k]=R[0][0]*i[k]+R[0][1]*i[k+1];ir[k+1]=R[1][0]*i[k]+R[1][1]*i[k+1];ir[k+2]=i[k+2];} return [{pts:o,sign:1,closed:false},{pts:ir,sign:Math.sign(d.sense)||1,closed:false}];}
    case 'picket': {const parts=picketFence(Math.round(d.rings),d.radius); return parts.map(p=>({pts:p.pts,sign:d.alternating===false?1:p.sign,closed:true}));}
    case 'link': {const parts=hopfLink(d.radius); return [{pts:parts[0].pts,sign:1,closed:true},{pts:parts[1].pts,sign:Math.sign(d.sense)||1,closed:true}];}
    case 'cone': return one(conePts(d.radius_base,d.radius_tip,d.height,Math.round(d.turns),90),false);
    case 'torushelix': return torusHelicalWindings(d.R0,d.r,d.l,d.periods,d.phase_deg);
  }
  throw new Error('unknown family '+family);
}
function designTransform(pts, comp){ const sc=comp.scale===undefined?1:+comp.scale; const m=[1,1,1]; if(comp.reflect){m['xyz'.indexOf(comp.reflect)]=-1;}
  const R=comp.rotate?rotMat(comp.rotate.axis||[0,0,1],(comp.rotate.deg||0)*Math.PI/180):null; const tr=comp.translate||[0,0,0]; const out=new Float64Array(pts.length);
  for(let i=0;i<pts.length;i+=3){let x=pts[i]*sc*m[0],y=pts[i+1]*sc*m[1],z=pts[i+2]*sc*m[2]; if(R){const X=R[0][0]*x+R[0][1]*y+R[0][2]*z,Y=R[1][0]*x+R[1][1]*y+R[1][2]*z,Z=R[2][0]*x+R[2][1]*y+R[2][2]*z; x=X;y=Y;z=Z;} out[i]=x+tr[0];out[i+1]=y+tr[1];out[i+2]=z+tr[2];}
  return out; }
function designParts(design){ // → [{pts,sign,closed,component,family}] world coordinates; closed knots stay closed, open strands get a remote return (in buildGeometry)
  const parts=[]; (design.components||[]).forEach((comp,ci)=>{ const cur=comp.current===undefined?1:+comp.current; if(cur===0)return; for(const p of designFamilyParts(comp.family,comp.params)){ parts.push({pts:designTransform(p.pts,comp),sign:cur*p.sign,closed:p.closed,approx:!!p.approx,component:ci,family:comp.family}); } }); return parts; }
function designSummary(design){ const c=design.components||[]; const fams={}; c.forEach(x=>{fams[x.family]=(fams[x.family]||0)+1;}); return Object.entries(fams).map(([f,n])=>n>1?f+'×'+n:f).join(' + '); }
function picketFence(rings, radius){ const parts=[]; for(let k=0;k<rings;k++){const z=(rings>1?(-0.62+1.24*k/(rings-1)):0)*radius; const rk=Math.sqrt(Math.max(radius*radius-z*z,1e-6)); parts.push({pts:circlePts(rk,720,0,0,z,'z'),sign:k%2?-1:1});} return parts; }
function sphereWinding(radius, turns){ const n=turns*90+1, out=new Float64Array(n*3); for(let i=0;i<n;i++){const t=i/(n-1); const z=radius*0.94*(2*t-1); const rho=Math.sqrt(Math.max(radius*radius-z*z,0)); const a=2*Math.PI*turns*t; out[3*i]=rho*Math.cos(a);out[3*i+1]=rho*Math.sin(a);out[3*i+2]=z;} return out; }
function hopfLink(radius){ return [{pts:circlePts(radius,720,-radius/2,0,0,'z'),sign:1},{pts:circlePts(radius,720,radius/2,0,0,'y'),sign:1}]; }
function bezier(p0,p1,p2,p3,n,out,off){for(let i=0;i<n;i++){const t=i/(n-1),a=(1-t)**3,b=3*(1-t)**2*t,c=3*(1-t)*t*t,d=t**3;for(let k=0;k<3;k++)out[off+3*i+k]=a*p0[k]+b*p1[k]+c*p2[k]+d*p3[k];}}
function remoteReturn(active, radius, nPoints){
  const n=active.length/3; const start=[active[0],active[1],active[2]], end=[active[3*n-3],active[3*n-2],active[3*n-1]];
  const t0=unit3([active[3]-active[0],active[4]-active[1],active[5]-active[2]]), t1=unit3([end[0]-active[3*n-6],end[1]-active[3*n-5],end[2]-active[3*n-4]]);
  const n0=unit3(start),n1=unit3(end); const q0=n0.map(x=>x*radius),q1=n1.map(x=>x*radius);
  let axis=cross(n1,n0); if(vlen(axis[0],axis[1],axis[2])<0.12){let trial=[0,0,1]; if(Math.abs(dot(trial,n0))>0.85) trial=[0,1,0]; axis=cross(n0,trial);} axis=unit3(axis); if(axis[2]<0) axis=axis.map(x=>-x);
  const apex=axis.map(x=>x*radius*1.24); const tang=unit3([q0[0]-q1[0],q0[1]-q1[1],q0[2]-q1[2]]); const leg=0.09*radius,outer=0.22*radius;
  const raw=new Float64Array((80+95+95+80)*3);
  bezier(end,[end[0]+leg*t1[0],end[1]+leg*t1[1],end[2]+leg*t1[2]],[q1[0]-outer*n1[0],q1[1]-outer*n1[1],q1[2]-outer*n1[2]],q1,80,raw,0);
  bezier(q1,[q1[0]+outer*n1[0],q1[1]+outer*n1[1],q1[2]+outer*n1[2]],[apex[0]-outer*tang[0],apex[1]-outer*tang[1],apex[2]-outer*tang[2]],apex,95,raw,240);
  bezier(apex,[apex[0]+outer*tang[0],apex[1]+outer*tang[1],apex[2]+outer*tang[2]],[q0[0]+outer*n0[0],q0[1]+outer*n0[1],q0[2]+outer*n0[2]],q0,95,raw,525);
  bezier(q0,[q0[0]-outer*n0[0],q0[1]-outer*n0[1],q0[2]-outer*n0[2]],[start[0]-leg*t0[0],start[1]-leg*t0[1],start[2]-leg*t0[2]],start,80,raw,810);
  return resample(raw,nPoints+1);
}
function rotMat(axis,ang){const a=unit3(axis),c=Math.cos(ang),s=Math.sin(ang),[x,y,z]=a;
  return [[c+x*x*(1-c),x*y*(1-c)-z*s,x*z*(1-c)+y*s],[y*x*(1-c)+z*s,c+y*y*(1-c),y*z*(1-c)-x*s],[z*x*(1-c)-y*s,z*y*(1-c)+x*s,c+z*z*(1-c)]];}
var LAYOUTS={
  'single':[[[0,0,0],0,[0,0,1],0,1]],
  'pair-axial-same':[[[0,0,1],0.5,[0,0,1],0,1],[[0,0,-1],0.5,[0,0,1],0,1]],
  'pair-axial-opposed':[[[0,0,1],0.5,[0,0,1],0,1],[[0,0,-1],0.5,[0,0,1],0,-1]],
  'pair-facing-flipped':[[[0,0,1],0.5,[1,0,0],Math.PI,1],[[0,0,-1],0.5,[0,0,1],0,1]],
  'pair-side-by-side':[[[1,0,0],0.5,[0,0,1],0,1],[[-1,0,0],0.5,[0,0,1],0,1]],
  'pair-side-crossed':[[[1,0,0],0.5,[0,0,1],0,1],[[-1,0,0],0.5,[1,0,0],Math.PI/2,1]],
  'triad-120':[0,1,2].map(k=>[[Math.cos(k*2*Math.PI/3),Math.sin(k*2*Math.PI/3),0],0.55,[0,0,1],k*2*Math.PI/3,1]),
  'triad-120-alternating':[0,1,2].map(k=>[[Math.cos(k*2*Math.PI/3),Math.sin(k*2*Math.PI/3),0],0.55,[0,0,1],k*2*Math.PI/3,k%2?-1:1]),
  'tetra-4':[[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]].map(d=>[d,0.5,[0,0,1],0,1]),
  'octa-6-cusp':[[[1,0,0],0.55,[0,1,0],Math.PI/2,1],[[-1,0,0],0.55,[0,1,0],-Math.PI/2,1],[[0,1,0],0.55,[1,0,0],-Math.PI/2,-1],[[0,-1,0],0.55,[1,0,0],Math.PI/2,-1],[[0,0,1],0.55,[0,0,1],0,1],[[0,0,-1],0.55,[1,0,0],Math.PI,1]],
  'octa-6-aligned':[[[1,0,0],0.55,[0,1,0],Math.PI/2,1],[[-1,0,0],0.55,[0,1,0],-Math.PI/2,1],[[0,1,0],0.55,[1,0,0],-Math.PI/2,1],[[0,-1,0],0.55,[1,0,0],Math.PI/2,1],[[0,0,1],0.55,[0,0,1],0,1],[[0,0,-1],0.55,[1,0,0],Math.PI,1]],
  'cube-8':[[1,1,1],[1,1,-1],[1,-1,1],[1,-1,-1],[-1,1,1],[-1,1,-1],[-1,-1,1],[-1,-1,-1]].map(d=>[d,0.55,[0,0,1],0,1]),
  'cube-8-checker':[[1,1,1],[1,1,-1],[1,-1,1],[1,-1,-1],[-1,1,1],[-1,1,-1],[-1,-1,1],[-1,-1,-1]].map(d=>[d,0.55,[0,0,1],0,d[0]*d[1]*d[2]>0?1:-1]),
};
function assemble(core, layout, coreScale){
  const spec=LAYOUTS[layout]; const cores=[];
  for(const [dir,dist,axis,ang,sign] of spec){const d=vlen(dir[0],dir[1],dir[2])>0?unit3(dir):[0,0,0]; const R=rotMat(axis,ang); const sc=spec.length>1?coreScale:1;
    const out=new Float64Array(core.length);
    for(let i=0;i<core.length;i+=3){const x=core[i]*sc,y=core[i+1]*sc,z=core[i+2]*sc;
      out[i]=R[0][0]*x+R[0][1]*y+R[0][2]*z+dist*d[0]; out[i+1]=R[1][0]*x+R[1][1]*y+R[1][2]*z+dist*d[1]; out[i+2]=R[2][0]*x+R[2][1]*y+R[2][2]*z+dist*d[2];}
    cores.push({pts:out,sign});}
  return cores;
}
function minNonlocal(p, skip, stride){ // conductor clearance on a subsampled polyline
  const n=Math.floor(p.length/3/stride); let best=Infinity;
  for(let i=0;i<n;i++){const xi=p[3*i*stride],yi=p[3*i*stride+1],zi=p[3*i*stride+2];
    for(let j=i+skip;j<n-skip;j++){if(n-j<=skip&&i<skip)continue; const dx=p[3*j*stride]-xi,dy=p[3*j*stride+1]-yi,dz=p[3*j*stride+2]-zi; const d2=dx*dx+dy*dy+dz*dz; if(d2<best)best=d2;}}
  return Math.sqrt(best);
}
// --- exact finite-segment Biot-Savart on a grid (unit current per circuit, sign per circuit) ---
function biotSavartGrid(circuits, n, half, a){
  const N=n*n*n; const B=new Float32Array(N*3); const axis=new Float64Array(n); for(let i=0;i<n;i++)axis[i]=-half+2*half*i/(n-1);
  const pref=1e-7; const a2=a*a;
  for(const {pts,sign} of circuits){
    const m=pts.length/3-1;
    const p1x=new Float64Array(m),p1y=new Float64Array(m),p1z=new Float64Array(m),p2x=new Float64Array(m),p2y=new Float64Array(m),p2z=new Float64Array(m),tx=new Float64Array(m),ty=new Float64Array(m),tz=new Float64Array(m),Ls=new Float64Array(m);
    for(let s=0;s<m;s++){p1x[s]=pts[3*s];p1y[s]=pts[3*s+1];p1z[s]=pts[3*s+2];p2x[s]=pts[3*s+3];p2y[s]=pts[3*s+4];p2z[s]=pts[3*s+5];const dx=p2x[s]-p1x[s],dy=p2y[s]-p1y[s],dz=p2z[s]-p1z[s];const L=vlen(dx,dy,dz)||1e-300;Ls[s]=L;tx[s]=dx/L;ty[s]=dy/L;tz[s]=dz/L;}
    let idx=0;
    for(let i=0;i<n;i++)for(let j=0;j<n;j++)for(let k=0;k<n;k++){
      const x=axis[i],y=axis[j],z=axis[k]; let bx=0,by=0,bz=0;
      for(let s=0;s<m;s++){
        const Rix=x-p1x[s],Riy=y-p1y[s],Riz=z-p1z[s],Rfx=x-p2x[s],Rfy=y-p2y[s],Rfz=z-p2z[s];
        const ri=Math.sqrt(Rix*Rix+Riy*Riy+Riz*Riz),rf=Math.sqrt(Rfx*Rfx+Rfy*Rfy+Rfz*Rfz);
        const rirf=ri*rf; const denom=rirf*(rirf+Rix*Rfx+Riy*Rfy+Riz*Rfz); if(denom<=1e-300)continue;
        const s0=Rix*tx[s]+Riy*ty[s]+Riz*tz[s]; const rho2=Math.max(ri*ri-s0*s0,0);
        let coef=(ri+rf)/denom; if(rho2<a2&&s0>=-a&&s0<=Ls[s]+a) coef*=rho2/a2;
        bx+=coef*(Riy*Rfz-Riz*Rfy); by+=coef*(Riz*Rfx-Rix*Rfz); bz+=coef*(Rix*Rfy-Riy*Rfx);
      }
      B[3*idx]+=pref*sign*bx; B[3*idx+1]+=pref*sign*by; B[3*idx+2]+=pref*sign*bz; idx++;
    }
  }
  return B;
}
var _segCache=null;
function segmentTable(circuits){ // flattened segment arrays (cached per circuit set): p1, p2, unit tangent, length, signed current
  if(_segCache&&_segCache.key===circuits) return _segCache;
  let S=0; for(const c of circuits)S+=c.pts.length/3-1;
  const p1=new Float64Array(S*3),p2=new Float64Array(S*3),th=new Float64Array(S*3),L=new Float64Array(S),I=new Float64Array(S); let k=0;
  for(const c of circuits){const pts=c.pts,m=pts.length/3-1; for(let s=0;s<m;s++){const dx=pts[3*s+3]-pts[3*s],dy=pts[3*s+4]-pts[3*s+1],dz=pts[3*s+5]-pts[3*s+2]; const l=Math.sqrt(dx*dx+dy*dy+dz*dz); if(l<1e-15)continue;
    p1[3*k]=pts[3*s];p1[3*k+1]=pts[3*s+1];p1[3*k+2]=pts[3*s+2];p2[3*k]=pts[3*s+3];p2[3*k+1]=pts[3*s+4];p2[3*k+2]=pts[3*s+5];th[3*k]=dx/l;th[3*k+1]=dy/l;th[3*k+2]=dz/l;L[k]=l;I[k]=c.sign;k++;}}
  _segCache={key:circuits,S:k,p1,p2,th,L,I}; return _segCache;
}
function biotSavartPoints(circuits, P, a){
  const T=segmentTable(circuits); const S=T.S,p1=T.p1,p2=T.p2,th=T.th,L=T.L,I=T.I; const M=P.length/3; const out=new Float64Array(M*3); const pref=1e-7,a2=a*a;
  for(let q=0;q<M;q++){const x=P[3*q],y=P[3*q+1],z=P[3*q+2]; let bx=0,by=0,bz=0;
    for(let s=0;s<S;s++){const s3=3*s; const Rix=x-p1[s3],Riy=y-p1[s3+1],Riz=z-p1[s3+2],Rfx=x-p2[s3],Rfy=y-p2[s3+1],Rfz=z-p2[s3+2];
      const ri=Math.sqrt(Rix*Rix+Riy*Riy+Riz*Riz),rf=Math.sqrt(Rfx*Rfx+Rfy*Rfy+Rfz*Rfz); const rirf=ri*rf; const denom=rirf*(rirf+Rix*Rfx+Riy*Rfy+Riz*Rfz); if(denom<=1e-300)continue;
      let coef=(ri+rf)/denom*I[s]; const s0=Rix*th[s3]+Riy*th[s3+1]+Riz*th[s3+2]; const rho2=ri*ri-s0*s0; if(rho2<a2&&s0>=-a&&s0<=L[s]+a)coef*=(rho2>0?rho2:0)/a2;
      bx+=coef*(Riy*Rfz-Riz*Rfy);by+=coef*(Riz*Rfx-Rix*Rfz);bz+=coef*(Rix*Rfy-Riy*Rfx);}
    out[3*q]=pref*bx;out[3*q+1]=pref*by;out[3*q+2]=pref*bz;}
  return out;
}
function maxwellChecks(circuits, a, rng){
  // 24 random points away from conductors: div B and curl B by central differences on the EXACT field; one Ampere loop
  const pts=[]; let guard=0; while(pts.length<72&&guard<5000){guard++; const x=rng()*1.2-0.6,y=rng()*1.2-0.6,z=rng()*1.2-0.6; let ok=true; for(const c of circuits){for(let s=0;s<c.pts.length;s+=6){const dx=c.pts[s]-x,dy=c.pts[s+1]-y,dz=c.pts[s+2]-z; if(dx*dx+dy*dy+dz*dz<0.0025){ok=false;break;}} if(!ok)break;} if(ok)pts.push(x,y,z);}
  const n=pts.length/3; const h=2e-4; const P=new Float64Array(n*6*3);
  for(let q=0;q<n;q++)for(let d=0;d<3;d++)for(let sgn=0;sgn<2;sgn++){const id=3*(q*6+2*d+sgn); P[id]=pts[3*q];P[id+1]=pts[3*q+1];P[id+2]=pts[3*q+2]; P[id+d]+=(sgn?-h:h);}
  const B=biotSavartPoints(circuits,P,a); const B0=biotSavartPoints(circuits,Float64Array.from(pts),a); let divMax=0,curlMax=0;
  for(let q=0;q<n;q++){const J=[[0,0,0],[0,0,0],[0,0,0]]; for(let d=0;d<3;d++){const ip=3*(q*6+2*d),im=3*(q*6+2*d+1); for(let c=0;c<3;c++)J[c][d]=(B[ip+c]-B[im+c])/(2*h);}
    const b=vlen(B0[3*q],B0[3*q+1],B0[3*q+2]); const div=Math.abs(J[0][0]+J[1][1]+J[2][2])*0.3/b; const curl=vlen(J[2][1]-J[1][2],J[0][2]-J[2][0],J[1][0]-J[0][1])*0.3/b; if(div>divMax)divMax=div; if(curl>curlMax)curlMax=curl;}
  // Ampère: loop of radius 3a around a segment of the first circuit
  // pick the segment (among 24 candidates on the first circuit) with the largest clearance from every other conductor, so the loop encloses exactly one current
  const c0=circuits[0]; const m0=c0.pts.length/3; let k=0, bestClr=-1;
  for(let cand=0;cand<24;cand++){const kk=Math.floor(m0*(cand+0.5)/24); const x=c0.pts[3*kk],y=c0.pts[3*kk+1],z=c0.pts[3*kk+2]; let clr=Infinity;
    for(let ci=0;ci<circuits.length;ci++){const p=circuits[ci].pts; const mm=p.length/3; for(let s=0;s<mm;s+=2){if(ci===0&&Math.abs(s-kk)<12)continue; const dx=p[3*s]-x,dy=p[3*s+1]-y,dz=p[3*s+2]-z; const d2=dx*dx+dy*dy+dz*dz; if(d2<clr)clr=d2;}}
    if(clr>bestClr){bestClr=clr;k=3*kk;}}
  const localClr=Math.sqrt(bestClr); const p1=[c0.pts[k],c0.pts[k+1],c0.pts[k+2]],p2=[c0.pts[k+3],c0.pts[k+4],c0.pts[k+5]]; const t=unit3([p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2]]);
  let trial=Math.abs(t[2])<0.9?[0,0,1]:[1,0,0]; const dd=dot(trial,t); const n1=unit3([trial[0]-dd*t[0],trial[1]-dd*t[1],trial[2]-dd*t[2]]); const n2=cross(t,n1); const rad=Math.max(1.05*a,Math.min(3*a,0.45*localClr)); const NL=720; const LP=new Float64Array(NL*3);
  for(let i=0;i<NL;i++){const th=2*Math.PI*i/NL; for(let c=0;c<3;c++)LP[3*i+c]=0.5*(p1[c]+p2[c])+rad*(Math.cos(th)*n1[c]+Math.sin(th)*n2[c]);}
  const BL=biotSavartPoints(circuits,LP,a); let circ=0; for(let i=0;i<NL;i++){const th=2*Math.PI*i/NL; const tg=[-Math.sin(th),Math.cos(th)]; for(let c=0;c<3;c++)circ+=BL[3*i+c]*rad*(tg[0]*n1[c]+tg[1]*n2[c]);} circ*=2*Math.PI/NL;
  return {divRelMax:divMax,curlRelMax:curlMax,ampere:circ/(4*Math.PI*1e-7*c0.sign),ampereRadius:rad,nPoints:n};
}
function wireDistanceGrid(circuits, n, half, stride){
  const N=n*n*n; const D=new Float32Array(N).fill(Infinity); const axis=new Float64Array(n); for(let i=0;i<n;i++)axis[i]=-half+2*half*i/(n-1);
  const P=[]; for(const c of circuits){for(let s=0;s<c.pts.length;s+=3*stride)P.push(c.pts[s],c.pts[s+1],c.pts[s+2]);}
  let idx=0; for(let i=0;i<n;i++)for(let j=0;j<n;j++)for(let k=0;k<n;k++){const x=axis[i],y=axis[j],z=axis[k]; let best=Infinity;
    for(let s=0;s<P.length;s+=3){const dx=P[s]-x,dy=P[s+1]-y,dz=P[s+2]-z; const d2=dx*dx+dy*dy+dz*dz; if(d2<best)best=d2;} D[idx++]=Math.sqrt(best);}
  return D;
}
function makeInterp(grid,n,half){
  const dx=2*half/(n-1);
  return function(x,y,z,out){
    let fx=(x+half)/dx,fy=(y+half)/dx,fz=(z+half)/dx; let i=Math.floor(fx),j=Math.floor(fy),k=Math.floor(fz);
    i=Math.min(Math.max(i,0),n-2);j=Math.min(Math.max(j,0),n-2);k=Math.min(Math.max(k,0),n-2);
    const wx=Math.min(Math.max(fx-i,0),1),wy=Math.min(Math.max(fy-j,0),1),wz=Math.min(Math.max(fz-k,0),1);
    out[0]=out[1]=out[2]=0;
    for(let a=0;a<2;a++){const ax=a?wx:1-wx; for(let b=0;b<2;b++){const ay=b?wy:1-wy; for(let c=0;c<2;c++){const az=c?wz:1-wz; const w=ax*ay*az; const id=3*(((i+a)*n+(j+b))*n+(k+c)); out[0]+=w*grid[id];out[1]+=w*grid[id+1];out[2]+=w*grid[id+2];}}}
  };
}
function makeInterpScalar(grid,n,half){
  const dx=2*half/(n-1);
  return function(x,y,z){
    let fx=(x+half)/dx,fy=(y+half)/dx,fz=(z+half)/dx; let i=Math.floor(fx),j=Math.floor(fy),k=Math.floor(fz);
    i=Math.min(Math.max(i,0),n-2);j=Math.min(Math.max(j,0),n-2);k=Math.min(Math.max(k,0),n-2);
    const wx=Math.min(Math.max(fx-i,0),1),wy=Math.min(Math.max(fy-j,0),1),wz=Math.min(Math.max(fz-k,0),1); let o=0;
    for(let a=0;a<2;a++){const ax=a?wx:1-wx; for(let b=0;b<2;b++){const ay=b?wy:1-wy; for(let c=0;c<2;c++){const az=c?wz:1-wz; o+=ax*ay*az*grid[((i+a)*n+(j+b))*n+(k+c)];}}} return o;
  };
}
function fieldStats(B,D,n,half){
  const axis=[]; for(let i=0;i<n;i++)axis.push(-half+2*half*i/(n-1)); let sum=0,cnt=0,bmax=0,divmax=0,bcore=[],bshell=[];
  let idx=0; for(let i=0;i<n;i++)for(let j=0;j<n;j++)for(let k=0;k<n;k++){const r=vlen(axis[i],axis[j],axis[k]); const b=vlen(B[3*idx],B[3*idx+1],B[3*idx+2]);
    if(r>=0.2&&r<=0.65&&D[idx]>=0.08){sum+=b*b;cnt++;} if(b>bmax)bmax=b; if(r<0.25&&D[idx]>0.05)bcore.push(b); if(Math.abs(r-0.75)<0.05&&D[idx]>0.05)bshell.push(b); idx++;}
  // divergence check by central differences on interior points away from conductors
  const dx=2*half/(n-1); let checked=0;
  for(let i=1;i<n-1;i++)for(let j=1;j<n-1;j++)for(let k=1;k<n-1;k++){const id=((i*n+j)*n+k); if(D[id]<0.08)continue;
    const d=(B[3*(((i+1)*n+j)*n+k)]-B[3*(((i-1)*n+j)*n+k)]+B[3*((i*n+j+1)*n+k)+1]-B[3*((i*n+j-1)*n+k)+1]+B[3*((i*n+j)*n+k+1)+2]-B[3*((i*n+j)*n+k-1)+2])/(2*dx);
    const b=vlen(B[3*id],B[3*id+1],B[3*id+2]); if(b>0){const rel=Math.abs(d)*dx/b; if(rel>divmax)divmax=rel; checked++;}}
  const med=a=>{if(!a.length)return NaN;a.sort((p,q)=>p-q);return a[a.length>>1];};
  return {brms:Math.sqrt(sum/Math.max(cnt,1)),bmax,divRelMax:divmax,divChecked:checked,bcore:med(bcore),bshell:med(bshell)};
}
function traceLines(B,D,n,half,seeds,wall,step,maxLen,wireClear){
  const interp=makeInterp(B,n,half), dist=makeInterpScalar(D,n,half); const out=[]; const tmp=[0,0,0];
  const tangent=(p,dir)=>{interp(p[0],p[1],p[2],tmp); const m=vlen(tmp[0],tmp[1],tmp[2])||1e-300; return [dir*tmp[0]/m,dir*tmp[1]/m,dir*tmp[2]/m];};
  for(let s=0;s<seeds.length;s+=3){const seed=[seeds[s],seeds[s+1],seeds[s+2]]; const res={seed,dirs:[]};
    for(const dir of [1,-1]){let p=seed.slice(); let L=0; const path=[p[0],p[1],p[2]]; interp(p[0],p[1],p[2],tmp); const b0=vlen(tmp[0],tmp[1],tmp[2]); let bmax=b0; let end='length';
      while(L<maxLen){const k1=tangent(p,dir),k2=tangent([p[0]+.5*step*k1[0],p[1]+.5*step*k1[1],p[2]+.5*step*k1[2]],dir),k3=tangent([p[0]+.5*step*k2[0],p[1]+.5*step*k2[1],p[2]+.5*step*k2[2]],dir),k4=tangent([p[0]+step*k3[0],p[1]+step*k3[1],p[2]+step*k3[2]],dir);
        p=[p[0]+step/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]),p[1]+step/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1]),p[2]+step/6*(k1[2]+2*k2[2]+2*k3[2]+k4[2])]; L+=step;
        if(!isFinite(p[0])){end='numerical';break;} const r=vlen(p[0],p[1],p[2]); if(r>=wall){end='wall';break;} if(Math.abs(p[0])>=half||Math.abs(p[1])>=half||Math.abs(p[2])>=half){end='grid-exit';break;}
        if(dist(p[0],p[1],p[2])<=wireClear){end='wire';break;} interp(p[0],p[1],p[2],tmp); const b=vlen(tmp[0],tmp[1],tmp[2]); if(b>bmax)bmax=b; path.push(p[0],p[1],p[2]);
        if(L>step*10&&vlen(p[0]-seed[0],p[1]-seed[1],p[2]-seed[2])<0.5*step){end='closed';break;}}
      res.dirs.push({end,L,bmax,path,b0});}
    const closed=res.dirs.some(d=>d.end==='closed')||(res.dirs[0].end==='length'&&res.dirs[1].end==='length');
    const R=closed?Infinity:Math.min(res.dirs[0].bmax,res.dirs[1].bmax)/res.dirs[0].b0; res.closed=closed; res.mirror=R; res.pred=closed?1:(R>1?Math.sqrt(1-1/R):0); out.push(res);}
  return out;
}
// --- Poincaré section on the EXACT field (batched RK4): punctures of the half-plane φ=0, x>0 ---
function poincare(circuits, a, seeds, turns, step, wireClear, wall){
  const n=seeds.length/3; const X=Float64Array.from(seeds); const alive=new Array(n).fill(true); const punct=Array.from({length:n},()=>[]); const turnsDone=new Int32Array(n);
  const phiPrev=new Float64Array(n); for(let i=0;i<n;i++)phiPrev[i]=Math.atan2(X[3*i+1],X[3*i]);
  // conductor proximity: coarse point list
  const cond=[]; for(const c of circuits){for(let s=0;s<c.pts.length;s+=6)cond.push(c.pts[s],c.pts[s+1],c.pts[s+2]);}
  const near=(x,y,z)=>{let best=Infinity; for(let s=0;s<cond.length;s+=3){const dx=cond[s]-x,dy=cond[s+1]-y,dz=cond[s+2]-z; const d=dx*dx+dy*dy+dz*dz; if(d<best)best=d;} return Math.sqrt(best);};
  const tangent=(P)=>{const B=biotSavartPoints(circuits,P,a); const m=P.length/3; for(let i=0;i<m;i++){const nb=vlen(B[3*i],B[3*i+1],B[3*i+2])||1e-300; B[3*i]/=nb;B[3*i+1]/=nb;B[3*i+2]/=nb;} return B;};
  const maxSteps=Math.floor(turns*2*Math.PI*1.2/step)+200; const ended=new Array(n).fill('complete');
  for(let it=0;it<maxSteps;it++){
    const idx=[]; for(let i=0;i<n;i++)if(alive[i])idx.push(i); if(!idx.length)break;
    const P=new Float64Array(idx.length*3); for(let j=0;j<idx.length;j++){P[3*j]=X[3*idx[j]];P[3*j+1]=X[3*idx[j]+1];P[3*j+2]=X[3*idx[j]+2];}
    const k1=tangent(P); const P2=new Float64Array(P.length); for(let i=0;i<P.length;i++)P2[i]=P[i]+.5*step*k1[i]; const k2=tangent(P2);
    for(let i=0;i<P.length;i++)P2[i]=P[i]+.5*step*k2[i]; const k3=tangent(P2); for(let i=0;i<P.length;i++)P2[i]=P[i]+step*k3[i]; const k4=tangent(P2);
    for(let j=0;j<idx.length;j++){const i=idx[j]; const px=P[3*j]+step/6*(k1[3*j]+2*k2[3*j]+2*k3[3*j]+k4[3*j]), py=P[3*j+1]+step/6*(k1[3*j+1]+2*k2[3*j+1]+2*k3[3*j+1]+k4[3*j+1]), pz=P[3*j+2]+step/6*(k1[3*j+2]+2*k2[3*j+2]+2*k3[3*j+2]+k4[3*j+2]);
      const phi=Math.atan2(py,px), pp=phiPrev[i];
      if(pp<0&&phi>=0&&Math.abs(phi-pp)<Math.PI){const f=-pp/(phi-pp); const qx=P[3*j]+f*(px-P[3*j]), qy=P[3*j+1]+f*(py-P[3*j+1]), qz=P[3*j+2]+f*(pz-P[3*j+2]); punct[i].push(Math.hypot(qx,qy),qz); turnsDone[i]++;}
      phiPrev[i]=phi; X[3*i]=px;X[3*i+1]=py;X[3*i+2]=pz;
      if(vlen(px,py,pz)>=wall){alive[i]=false;ended[i]='wall';} else if(it%4===0&&near(px,py,pz)<=wireClear){alive[i]=false;ended[i]='wire';} else if(turnsDone[i]>=turns){alive[i]=false;}}
  }
  return {punct,ended,turns:Array.from(turnsDone)};
}
function poincareAnalyse(circuits, a, R0, rTube, nSeeds, turns, step, wireClear){
  // pass 1: axis from a coarse fan (8 turns); pass 2: ray from the axis outboard
  const fan=[]; for(let i=0;i<12;i++){fan.push(R0+(-0.5+i/11)*rTube,0,0);} for(let i=0;i<6;i++){fan.push(R0,0,(-0.4+0.8*i/5)*rTube);}
  const p1=poincare(circuits,a,Float64Array.from(fan),8,step,wireClear,1.05);
  let axis=[R0,0], bestSpread=Infinity;
  p1.punct.forEach((pc,i)=>{ if(p1.ended[i]!=='complete'||pc.length<14)return; let mR=0,mz=0,k=pc.length/2; for(let j=0;j<pc.length;j+=2){mR+=pc[j];mz+=pc[j+1];} mR/=k;mz/=k; let sp=0; for(let j=0;j<pc.length;j+=2)sp+=Math.abs(pc[j]-mR)+Math.abs(pc[j+1]-mz); sp/=k; if(sp<bestSpread){bestSpread=sp;axis=[mR,mz];} });
  const reach=Math.max(R0+rTube-axis[0]-0.03,0.05); const seeds=[]; for(let i=0;i<nSeeds;i++){const s=0.03+0.92*i/(nSeeds-1); seeds.push(axis[0]+s*reach,0,axis[1]);}
  const p2=poincare(circuits,a,Float64Array.from(seeds),turns,step,wireClear,1.05);
  const lines=p2.punct.map((pc,i)=>{ const K=pc.length/2; if(p2.ended[i]!=='complete'||K<10) return {kind:'open',ended:p2.ended[i],turns:p2.turns[i],punct:pc};
    let th=[],r=[]; for(let j=0;j<pc.length;j+=2){th.push(Math.atan2(pc[j+1]-axis[1],pc[j]-axis[0])); r.push(Math.hypot(pc[j]-axis[0],pc[j+1]-axis[1]));}
    for(let j=1;j<th.length;j++){while(th[j]-th[j-1]>Math.PI)th[j]-=2*Math.PI; while(th[j]-th[j-1]<-Math.PI)th[j]+=2*Math.PI;}
    const iota=(th[th.length-1]-th[0])/(2*Math.PI*(K-1)); const rm=r.reduce((x,y)=>x+y,0)/K;
    // Fourier residual m<=4
    const M=9; const A=[]; for(let j=0;j<K;j++){const row=[1]; for(let m=1;m<=4;m++){row.push(Math.cos(m*th[j]),Math.sin(m*th[j]));} A.push(row);}
    // normal equations
    const AtA=Array.from({length:M},()=>new Array(M).fill(0)), Atb=new Array(M).fill(0);
    for(let j=0;j<K;j++)for(let p=0;p<M;p++){Atb[p]+=A[j][p]*r[j]; for(let q=0;q<M;q++)AtA[p][q]+=A[j][p]*A[j][q];}
    for(let p=0;p<M;p++)AtA[p][p]+=1e-9;
    // solve (Gauss)
    const c=Atb.slice(); const G=AtA.map(row=>row.slice());
    for(let p=0;p<M;p++){let piv=p; for(let q=p+1;q<M;q++)if(Math.abs(G[q][p])>Math.abs(G[piv][p]))piv=q; [G[p],G[piv]]=[G[piv],G[p]]; [c[p],c[piv]]=[c[piv],c[p]];
      for(let q=p+1;q<M;q++){const f=G[q][p]/G[p][p]; for(let k=p;k<M;k++)G[q][k]-=f*G[p][k]; c[q]-=f*c[p];}}
    for(let p=M-1;p>=0;p--){let s=c[p]; for(let k=p+1;k<M;k++)s-=G[p][k]*c[k]; c[p]=s/G[p][p];}
    let res=0; for(let j=0;j<K;j++){let fit=0; for(let p=0;p<M;p++)fit+=A[j][p]*c[p]; res+=(fit-r[j])**2;} res=Math.sqrt(res/K)/rm;
    return {kind:res<0.06?'surface':'island/chaotic',ended:'complete',turns:p2.turns[i],punct:pc,iota,r:rm,residual:res}; });
  return {axis,reach,lines};
}
// --- Boris pusher, code units: q/m = 1 (times charge sign), B in the units of the grid ---
function runParticles(B,D,n,half,params){
  const {x0,v0,q,tMax,dt,wall,wireR,omdt,sampleEvery}=params; const N=x0.length/3; const interp=makeInterp(B,n,half),dist=makeInterpScalar(D,n,half);
  const x=Float64Array.from(x0),v=Float64Array.from(v0); const alive=new Uint8Array(N).fill(1); const lossT=new Float64Array(N).fill(Infinity); const lossK=new Uint8Array(N); // 1 wall 2 wire 3 exit
  const steps=Math.ceil(tMax/dt); const frames=[]; const times=[]; const tmp=[0,0,0]; let t=0; let subMax=1;
  // mu diagnostics
  const mu0=new Float64Array(N),muMin=new Float64Array(N),muMax=new Float64Array(N);
  for(let i=0;i<N;i++){interp(x[3*i],x[3*i+1],x[3*i+2],tmp);const b=vlen(tmp[0],tmp[1],tmp[2])||1e-300;const vp=(v[3*i]*tmp[0]+v[3*i+1]*tmp[1]+v[3*i+2]*tmp[2])/b;const vv=v[3*i]*v[3*i]+v[3*i+1]*v[3*i+1]+v[3*i+2]*v[3*i+2];mu0[i]=Math.max(vv-vp*vp,0)/(2*b);muMin[i]=muMax[i]=mu0[i];}
  frames.push(Float32Array.from(x)); times.push(0); const aliveFrames=[Uint8Array.from(alive)];
  for(let s=0;s<steps;s++){
    for(let i=0;i<N;i++){ if(!alive[i])continue;
      interp(x[3*i],x[3*i+1],x[3*i+2],tmp); const bm=vlen(tmp[0],tmp[1],tmp[2]); const nsub=Math.min(Math.max(Math.ceil(Math.abs(q[i])*bm*dt/omdt),1),4096); if(nsub>subMax)subMax=nsub; const h=dt/nsub;
      for(let k=0;k<nsub;k++){ if(k>0)interp(x[3*i],x[3*i+1],x[3*i+2],tmp);
        const tx=q[i]*tmp[0]*0.5*h,ty=q[i]*tmp[1]*0.5*h,tz=q[i]*tmp[2]*0.5*h; const t2=tx*tx+ty*ty+tz*tz; const sx=2*tx/(1+t2),sy=2*ty/(1+t2),sz=2*tz/(1+t2);
        const vx=v[3*i],vy=v[3*i+1],vz=v[3*i+2]; const px=vx+(vy*tz-vz*ty),py=vy+(vz*tx-vx*tz),pz=vz+(vx*ty-vy*tx);
        const nx=vx+(py*sz-pz*sy),ny=vy+(pz*sx-px*sz),nz=vz+(px*sy-py*sx); v[3*i]=nx;v[3*i+1]=ny;v[3*i+2]=nz; x[3*i]+=nx*h;x[3*i+1]+=ny*h;x[3*i+2]+=nz*h; }
      const X=x[3*i],Y=x[3*i+1],Z=x[3*i+2]; const r=vlen(X,Y,Z);
      if(r>=wall){alive[i]=0;lossT[i]=t+dt;lossK[i]=1;continue;} if(Math.abs(X)>=half||Math.abs(Y)>=half||Math.abs(Z)>=half){alive[i]=0;lossT[i]=t+dt;lossK[i]=3;continue;}
      if(dist(X,Y,Z)<=wireR){alive[i]=0;lossT[i]=t+dt;lossK[i]=2;continue;}
      interp(X,Y,Z,tmp); const b=vlen(tmp[0],tmp[1],tmp[2])||1e-300; const vp=(v[3*i]*tmp[0]+v[3*i+1]*tmp[1]+v[3*i+2]*tmp[2])/b; const vv=v[3*i]*v[3*i]+v[3*i+1]*v[3*i+1]+v[3*i+2]*v[3*i+2]; const mu=Math.max(vv-vp*vp,0)/(2*b); if(mu<muMin[i])muMin[i]=mu; if(mu>muMax[i])muMax[i]=mu;
    }
    t+=dt; if((s+1)%sampleEvery===0||s===steps-1){frames.push(Float32Array.from(x));times.push(t);aliveFrames.push(Uint8Array.from(alive));}
  }
  let retained=0,nwall=0,nwire=0,nexit=0; for(let i=0;i<N;i++){if(alive[i])retained++;else if(lossK[i]===1)nwall++;else if(lossK[i]===2)nwire++;else nexit++;}
  const muVar=[]; for(let i=0;i<N;i++)muVar.push(muMax[i]>0?(muMax[i]-muMin[i])/muMax[i]:0); muVar.sort((a,b)=>a-b);
  return {frames,times,aliveFrames,lossT:Array.from(lossT),lossK:Array.from(lossK),retained,nwall,nwire,nexit,N,subMax,tMax:t,muVarMedian:muVar[muVar.length>>1]};
}
`;
