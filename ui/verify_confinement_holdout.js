'use strict';
const fs=require('fs'),vm=require('vm'),assert=require('assert'),path=require('path');
function element(){return {children:[],style:{},textContent:'',appendChild(x){this.children.push(x);return x;},replaceChildren(){this.children=[];}};}
function group(label,delta){return {comparison:label,complete_eligible_panels:true,conditions_consistent_across_panels:true,panels:[0,1,2].map(replicate=>({replicate,eligible:true,paired_markers_available:true,right_rmst:20,left_rmst:20+delta,difference:delta})),pooled:{right_rmst:20,left_rmst:20+delta,difference:delta,panels:3,marker_bootstrap_interval:{level:.95,lower:delta-1,upper:delta+1}}};}
const complete={expected_replicates:[0,1,2],nominee_locks:[{arm:'direct',nominee:'primary',design_sha256:'a'.repeat(64)},{arm:'construction',nominee:'primary',design_sha256:'b'.repeat(64)},{arm:'construction',nominee:'alternate',design_sha256:'c'.repeat(64)}],versus_baseline:[group('direct:primary minus baseline',1),group('construction:primary minus baseline',2),group('construction:alternate minus baseline',90)],primary_between_arms:group('construction:primary minus direct:primary',1)};
async function render(fixture){const nodes=new Map();const document={getElementById(id){if(!nodes.has(id))nodes.set(id,element());return nodes.get(id);},createElement:element};vm.runInNewContext(fs.readFileSync(path.join(__dirname,'confinement_holdout.js'),'utf8'),{document,fetch:async()=>({ok:true,json:async()=>fixture}),setTimeout:()=>{},Date});await new Promise(resolve=>setImmediate(resolve));return nodes;}
(async()=>{
  const good=await render(complete);
  assert.strictEqual(good.get('holdout-primary').children.length,8);
  assert.strictEqual(good.get('holdout-alternates').children.length,4);
  assert(good.get('primary-comparison').textContent.includes('+1.00 transits'));
  assert(!good.get('primary-comparison').textContent.includes('+90.00'));
  const incomplete=JSON.parse(JSON.stringify(complete));
  incomplete.versus_baseline[1].panels[1]={replicate:1,eligible:false,reasons:['left:numerical_qualification_missing_or_failed']};
  // Even accidentally retained pooled values and true summary flags cannot bypass a failed panel.
  const blocked=await render(incomplete);
  assert(blocked.get('primary-comparison').textContent.includes('withheld'));
  const construction=blocked.get('holdout-primary').children;
  assert.strictEqual(construction[7].children[4].textContent,'—');
  assert(construction[5].children[6].textContent.includes('Numerical qualification'));
  const absent=JSON.parse(JSON.stringify(complete));absent.nominee_locks=[];
  const unlocked=await render(absent);assert(unlocked.get('primary-comparison').textContent.includes('withheld'));
  console.log('Passed: three panel and aggregate values render; alternates cannot replace primary; failed/missing eligibility and missing locks suppress aggregate comparison.');
})();
