// UI logic check without a browser, simulation, or model call.
'use strict';
const fs=require('fs'),vm=require('vm'),assert=require('assert'),path=require('path');
function element(){return {children:[],style:{},value:'all',textContent:'',appendChild(x){this.children.push(x);return x;},replaceChildren(){this.children=[];},setAttribute(k,v){this[k]=v;},addEventListener(k,f){this[k]=f;}};}
const nodes=new Map();
const document={getElementById(id){if(!nodes.has(id))nodes.set(id,element());return nodes.get(id);},createElement:element,createElementNS:element};
const fixture={status:'holdout',arms:{direct:{used:4,limit:12,participants:['Astra','Fable']},construction:{used:1,limit:12}},evaluations:[
  {id:'qualified',arm:'direct',phase:'screen',accepted_ordinal:7,valid:true,numerically_qualified:true,primary:21,design:{components:[{family:'torushelix'}]}},
  {id:'unqualified',arm:'direct',phase:'screen',valid:true,numerically_qualified:false,primary:119},
  {id:'invalid',arm:'direct',phase:'screen',valid:false,numerically_qualified:true,primary:120},
  {id:'unknown-qualification',arm:'direct',phase:'screen',valid:true,primary:118},
  {id:'heldout',arm:'direct',phase:'holdout',valid:true,numerically_qualified:true,primary:99},
  {id:'null',arm:'construction',phase:'screen',valid:true,numerically_qualified:true,primary:null},
  {id:'sensitivity',arm:'construction',phase:'sensitivity',valid:true,numerically_qualified:true,primary:119},
  {id:'earlier-qualified',arm:'direct',phase:'screen',accepted_ordinal:3,valid:true,numerically_qualified:true,primary:17}
]};
vm.runInNewContext(fs.readFileSync(path.join(__dirname,'confinement_study.js'),'utf8'),{
  document,fetch:async()=>({ok:true,json:async()=>fixture}),setTimeout:()=>{},Date,console
});
setImmediate(()=>{
  assert.strictEqual(nodes.get('direct-best').textContent,'21.00');
  assert.strictEqual(nodes.get('construction-best').textContent,'—');
  assert(nodes.get('notice').textContent.includes('no winner'));
  const failed=nodes.get('rows').children.find(tr=>tr.children[0].textContent==='unqualified');
  assert(failed.children[2].children[0].textContent.includes('numerical check failed'));
  assert.strictEqual(failed.children[3].textContent,'119.00');
  const qualified=nodes.get('rows').children.find(tr=>tr.children[0].textContent==='qualified');
  assert.strictEqual(qualified.children[6].children[0].textContent,'View in 3D');
  const curve=nodes.get('chart').children.find(x=>x.d);
  assert(curve.d.includes('V201')); // RM120=21, never 119/120/118/99.
  assert(curve.d.startsWith('M305.5,')); // Accepted ordinal 3, despite reverse completion order.
  assert(curve.d.includes('H635.5')); // Accepted ordinal 7, preserving the submission gap.
  nodes.get('phase').value='sensitivity';nodes.get('phase').change();
  assert.strictEqual(nodes.get('rows').children.length,1);
  assert.strictEqual(nodes.get('rows').children[0].children[0].textContent,'sensitivity');
  console.log('Passed: qualification gates leaders/chart; raw failed results remain visible; holdouts/sensitivity never become screening wins; accepted submission order and gaps are preserved; sensitivity filter and torushelix link work.');
});
