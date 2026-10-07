'use strict';
(() => {
  const $=id=>document.getElementById(id);
  const labels={direct:'Direct search',construction:'Build and explore'};
  let lastAnalysis=null;
  const numeric=x=>typeof x==='number'&&Number.isFinite(x);
  const number=x=>numeric(x)?x.toFixed(2):'—';
  const signed=x=>numeric(x)?(x>0?'+':'')+x.toFixed(2):'—';
  function cell(tr,text,cls){const td=document.createElement('td');td.textContent=text;if(cls)td.className=cls;tr.appendChild(td);return td;}
  function interval(value){
    if(!value||!numeric(value.lower)||!numeric(value.upper))return 'Not available';
    const level=numeric(value.level)?Math.round(100*value.level)+'% ':'';
    return level+'['+signed(value.lower)+', '+signed(value.upper)+']';
  }
  function reasons(panel){
    if(panel.eligible===true)return panel.paired_markers_available===true?'Eligible · paired markers':'Eligible · descriptive only';
    const messages=(panel.reasons||[]).map(raw=>{
      const reason=String(raw);
      if(reason.includes('missing_result_file')||reason.includes('missing_expected_panel'))return 'Result missing';
      if(reason.includes('numerical_qualification'))return 'Numerical qualification failed or missing';
      if(reason.includes('locked_nominee_hash'))return 'Candidate differs from its locked design';
      if(reason.includes('conditions_differ'))return 'Evaluation conditions differ';
      if(reason.includes('invalid_design'))return 'Invalid geometry';
      if(reason.includes('diagnostic_or_nonprimary_grid'))return 'Diagnostic run; excluded from primary comparison';
      if(reason.includes('nonprimary_endpoint'))return 'Different endpoint';
      return reason.replace(/_/g,' ');
    });
    return [...new Set(messages)].join('; ')||'Eligibility not confirmed';
  }
  function qualifiedPooled(comparison,expected){
    if(!comparison||comparison.complete_eligible_panels!==true||comparison.conditions_consistent_across_panels!==true)return null;
    if(!expected.every(r=>(comparison.panels||[]).some(p=>p.replicate===r&&p.eligible===true)))return null;
    const pooled=comparison.pooled;
    return pooled&&numeric(pooled.left_rmst)&&numeric(pooled.right_rmst)&&numeric(pooled.difference)?pooled:null;
  }
  function addGroup(tbody,analysis,arm,nominee){
    const expected=analysis.expected_replicates;
    const label=arm+':'+nominee+' minus baseline';
    const lock=(analysis.nominee_locks||[]).find(x=>x.arm===arm&&(x.nominee||x.name)===nominee);
    const hash=lock&&typeof lock.design_sha256==='string'?lock.design_sha256:null;
    const comparison=hash?analysis.versus_baseline.find(c=>c.comparison===label):null;
    for(const replicate of expected){
      const panel=(comparison?.panels||[]).find(p=>p.replicate===replicate)||{replicate,eligible:false,reasons:['missing_expected_panel']};
      const tr=document.createElement('tr');
      const name=cell(tr,labels[arm]+' · '+nominee,'holdout-candidate');
      if(hash)name.title='Locked design SHA-256: '+hash;
      cell(tr,'Panel '+(replicate+1),'panel-label');
      cell(tr,panel.eligible===true?number(panel.right_rmst):'—','num');
      cell(tr,panel.eligible===true?number(panel.left_rmst):'—','num');
      cell(tr,panel.eligible===true?signed(panel.difference):'—','num');
      cell(tr,panel.eligible===true?interval(panel.marker_bootstrap_interval):'—','holdout-interval');
      cell(tr,reasons(panel),'validation-status');tbody.appendChild(tr);
    }
    const pooled=qualifiedPooled(comparison,expected);
    const tr=document.createElement('tr');tr.className='pooled-row';
    cell(tr,labels[arm]+' · '+nominee,'holdout-candidate');
    cell(tr,'Combined '+expected.length+' panels','panel-label');
    cell(tr,pooled?number(pooled.right_rmst):'—','num');
    cell(tr,pooled?number(pooled.left_rmst):'—','num');
    cell(tr,pooled?signed(pooled.difference):'—','num');
    cell(tr,pooled?interval(pooled.marker_bootstrap_interval):'—','holdout-interval');
    cell(tr,pooled?'All required panels eligible':comparison?.conditions_consistent_across_panels===false?'Conditions changed; no combined comparison':'Missing or ineligible panels; no combined comparison','validation-status');
    tbody.appendChild(tr);
    return !!pooled;
  }
  function render(analysis){
    $('holdout-primary').replaceChildren();$('holdout-alternates').replaceChildren();
    let completed=0;
    for(const arm of Object.keys(labels))if(addGroup($('holdout-primary'),analysis,arm,'primary'))completed++;
    const alternates=(analysis.nominee_locks||[]).filter(x=>(x.nominee||x.name)==='alternate'&&labels[x.arm]);
    $('alternate-section').hidden=alternates.length===0;
    for(const lock of alternates)addGroup($('holdout-alternates'),analysis,lock.arm,'alternate');
    $('holdout-state').textContent=completed===2?'Both locked primary candidates have complete eligible panel comparisons. Values below describe these candidates; this single pilot does not establish a superior general strategy.':'Held-out comparison is incomplete or has an eligibility failure. Available eligible panels remain visible; no missing result is treated as zero.';
    const between=analysis.primary_between_arms;
    const pooled=qualifiedPooled(between,analysis.expected_replicates);
    $('primary-comparison').hidden=false;
    if(pooled&&completed===2){
      const direction=String(between.comparison||'').replace('construction:primary','Build and explore primary').replace('direct:primary','Direct search primary');
      $('primary-comparison').textContent=direction+': '+signed(pooled.difference)+' transits. Conditional marker interval: '+interval(pooled.marker_bootstrap_interval)+'. This is a comparison of the locked candidates, not a causal strategy winner.';
    }else $('primary-comparison').textContent='An aggregate comparison between primary candidates is withheld until their required paired-baseline panels are eligible.';
    $('holdout-links').hidden=false;
  }
  async function poll(){
    try{
      const response=await fetch('../results/confinement_holdout_analysis.json',{cache:'no-store'});
      if(!response.ok){if(response.status===404&&!lastAnalysis){$('holdout-updated').textContent='Held-out analysis has not been published yet.';return;}throw new Error('Held-out comparison feed is unavailable.');}
      const data=await response.json();
      if(!data||!Array.isArray(data.versus_baseline)||!Array.isArray(data.nominee_locks)||!Array.isArray(data.expected_replicates)||!data.expected_replicates.length)throw new Error('Waiting for a complete comparison update.');
      lastAnalysis=data;render(data);$('holdout-updated').textContent='Comparison last checked '+new Date().toLocaleTimeString()+'. All values are in transits; combined means are weighted by marker count.';
    }catch(error){$('holdout-updated').textContent=error.message+(lastAnalysis?' Last recorded comparison remains visible.':'');}
    finally{setTimeout(poll,3000);}
  }
  poll();
})();
