'use strict';
(() => {
  const $ = id => document.getElementById(id);
  const labels = {direct:'Direct search',construction:'Build and explore'};
  const colors = {direct:'#efa86a',construction:'#71d4c3'};
  const previewFamilies=new Set(['circle','solenoid','precess','phase','codexhopf','hopfdrift','hopfcont','hopftorus','recursive','baseball','sphere','hopfmirror','hopfmesh','yinyang','picket','link','cone','torushelix']);
  let data = null;
  const numeric = x => typeof x === 'number' && Number.isFinite(x);
  const format = x => numeric(x) ? x.toFixed(2) : '—';
  const scalar = x => numeric(x) ? x : x && numeric(x.primary) ? x.primary : null;
  const percent = x => numeric(x) ? (100*x).toFixed(1)+'%' : '—';
  function cell(tr,value,cls) {const td=document.createElement('td');td.textContent=value;if(cls)td.className=cls;tr.appendChild(td);return td;}
  function survivor(row) {
    if(numeric(row.survival))return row.survival;
    const s=row.survival||{};
    return numeric(s.S120)?s.S120:numeric(s['120'])?s['120']:null;
  }
  function person(p) {return typeof p==='string'?p:p&&(p.name||p.actor||p.id||p.model)||'Participant';}
  function screens(arm) {
    return (data.evaluations||[]).filter(r=>r.arm===arm&&r.phase==='screen')
      .map((row,index)=>({row,index,ordinal:numeric(row.accepted_ordinal)&&row.accepted_ordinal>0?row.accepted_ordinal:index+1}))
      .sort((a,b)=>a.ordinal-b.ordinal||a.index-b.index)
      .map(item=>Object.assign({},item.row,{study_submission:item.ordinal}));
  }
  function renderChart() {
    const svg=$('chart');svg.replaceChildren();
    const ns='http://www.w3.org/2000/svg';
    function add(tag,attrs,text){const e=document.createElementNS(ns,tag);Object.entries(attrs).forEach(([k,v])=>e.setAttribute(k,v));if(text!==undefined)e.textContent=text;svg.appendChild(e);return e;}
    const max=Math.max(1,...Object.keys(labels).map(a=>Math.max(0,...screens(a).map(r=>r.study_submission),Number(data.arms?.[a]?.limit)||0)));
    const X=n=>58+990*n/max,Y=n=>236-200*n/120;
    for(const y of [0,30,60,90,120]){add('line',{x1:58,x2:1048,y1:Y(y),y2:Y(y),stroke:'#293743'});add('text',{x:45,y:Y(y)+4,fill:'#93a5b1','font-size':12,'text-anchor':'end'},y);}
    for(const x of [...new Set([0,Math.ceil(max/4),Math.ceil(max/2),Math.ceil(max*3/4),max])])add('text',{x:X(x),y:257,fill:'#93a5b1','font-size':12,'text-anchor':'middle'},x);
    add('text',{x:58,y:19,fill:'#93a5b1','font-size':12},'RM120 · transits');add('text',{x:1048,y:277,fill:'#93a5b1','font-size':12,'text-anchor':'end'},'Accepted screening submission');
    let plotted=false;
    for(const arm of Object.keys(labels)){
      let best=null,path='',lastX=null,lastY=null;
      screens(arm).forEach(r=>{
        if(r.valid===true&&r.numerically_qualified===true&&numeric(r.primary))best=best===null?r.primary:Math.max(best,r.primary);
        if(best===null)return;
        const x=X(r.study_submission),y=Y(Math.max(0,Math.min(120,best)));
        path+=lastX===null?`M${x},${y}`:`H${x}V${y}`;lastX=x;lastY=y;plotted=true;
      });
      if(path){add('path',{d:path,fill:'none',stroke:colors[arm],'stroke-width':2.5});add('circle',{cx:lastX,cy:lastY,r:4,fill:colors[arm]});}
    }
    $('chart-empty').hidden=plotted;
    $('chart-empty').style.display=plotted?'none':'grid';
  }
  function renderRows(){
    const tbody=$('rows');tbody.replaceChildren();
    const all=data.evaluations||[],phase=$('phase').value;
    const rows=all.filter(r=>phase==='all'||r.phase===phase).slice().reverse();
    $('count').textContent=all.length?'('+all.length+')':'';
    if(!rows.length){const tr=document.createElement('tr');const td=cell(tr,'No recorded evaluations in this phase.','muted');td.colSpan=7;tbody.appendChild(tr);return;}
    for(const r of rows){
      const tr=document.createElement('tr');
      cell(tr,r.name||r.design?.name||r.id||'Unnamed candidate','name');
      cell(tr,(labels[r.arm]||'Shared baseline')+(r.actor?' · '+r.actor:''));
      const phaseCell=cell(tr,'');const badge=document.createElement('span');badge.className='badge'+(r.valid===false||r.numerically_qualified===false?' bad':'');badge.textContent=(r.phase==='holdout'?'Held-out':r.phase==='screen'?'Screening':r.phase==='sensitivity'?'Sensitivity':r.phase||'Recorded')+(r.valid===false?' · invalid':r.numerically_qualified===false?' · numerical check failed':r.numerically_qualified!==true?' · qualification pending':'');phaseCell.appendChild(badge);
      cell(tr,r.valid===false?'Invalid':format(r.primary),'num');cell(tr,percent(survivor(r)),'num');cell(tr,numeric(r.seconds)?r.seconds.toFixed(1)+' s':'—','num');
      const td=cell(tr,'');
      if(r.design&&Array.isArray(r.design.components)&&r.id!==undefined){const a=document.createElement('a');const supported=r.design.components.every(c=>previewFamilies.has(c.family));a.textContent=supported?'View in 3D':'Inspect design JSON';a.title=supported?'Open this recorded candidate in the Lab':'This family is not supported by the current 3D preview; its complete JSON opens in the Lab editor.';a.href='chad_core_lab_standalone.html?study='+encodeURIComponent(String(r.id));a.target='_blank';a.rel='noopener';td.appendChild(a);}else td.textContent='—';
      tbody.appendChild(tr);
    }
  }
  function render(){
    const status=typeof data.status==='string'?data.status:'setup';
    $('status').textContent=({setup:'Preparing the study',search:'Search in progress',holdout:'Held-out validation',sensitivity:'Checking numerical sensitivity',complete:'Study complete'})[status]||status;
    $('dot').style.background=status==='complete'?'#93a5b1':'#71d4c3';
    $('notice').textContent=status==='holdout'?'Candidates are locked for held-out validation. Screening leads are not confirmed improvements; no winner is declared here.':status==='complete'?'Search and validation are recorded. Read the paired held-out results and study report before drawing a conclusion.':'RM120 is average particle confinement time, capped at 120 transits. It is separate from the Lab’s legacy Score. Screening leads still need held-out validation.';
    for(const arm of Object.keys(labels)){
      const info=data.arms?.[arm]||{};
      const values=screens(arm).filter(r=>r.valid===true&&r.numerically_qualified===true&&numeric(r.primary)).map(r=>r.primary);
      $(arm+'-best').textContent=values.length?format(Math.max(...values)):'—';
      const used=numeric(info.used)?info.used:screens(arm).length,limit=numeric(info.limit)?info.limit:null;
      $(arm+'-used').textContent=used+' / '+(limit===null?'—':limit);
      $(arm+'-fill').style.width=(limit>0?Math.min(100,100*used/limit):0)+'%';
      const participants=Array.isArray(info.participants)?info.participants:[];
      $(arm+'-people').textContent=participants.length?participants.map(person).join(' + '):'Participants not yet recorded';
    }
    const base=scalar(data.baseline);
    $('baseline').textContent=base===null?'Shared baseline: awaiting recorded evaluation.':'Shared baseline RM120: '+format(base)+' transits. Panel and phase are recorded in the study data.';
    renderChart();renderRows();
  }
  async function poll(){
    try{
      const response=await fetch('../results/confinement_study.json',{cache:'no-store'});
      if(!response.ok)throw new Error(response.status===404?'Waiting for the study to start.':'Study feed unavailable ('+response.status+').');
      const next=await response.json();if(!next||!Array.isArray(next.evaluations))throw new Error('Waiting for a complete study update.');
      data=next;render();$('updated').className='';$('updated').textContent='Last checked '+new Date().toLocaleTimeString()+'. Updates every 3 seconds.';
    }catch(error){$('updated').className='error';$('updated').textContent=error.message+(data?' Last recorded results remain visible.':'');if(!data)$('status').textContent='Waiting for study data';}
    setTimeout(poll,3000);
  }
  $('phase').addEventListener('change',()=>{if(data)renderRows();});
  poll();
})();
