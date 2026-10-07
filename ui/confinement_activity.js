'use strict';
(() => {
  const rooms={direct:'http://127.0.0.1:8788',construction:'http://127.0.0.1:8789'};
  const numeric=x=>typeof x==='number'&&Number.isFinite(x);
  let seen=false;
  function node(tag,text,cls){const element=document.createElement(tag);if(text!==undefined)element.textContent=text;if(cls)element.className=cls;return element;}
  function render(arms){
    for(const [arm,url] of Object.entries(rooms)){
      const record=arms[arm];if(!record)continue;
      const host=document.getElementById(arm+'-activity');host.replaceChildren();
      const head=node('div',undefined,'activity-heading');head.appendChild(node('span','Public activity · '+(record.status||'status unavailable')));
      const link=node('a','Open team room');link.href=url;link.target='_blank';link.rel='noopener';head.appendChild(link);host.appendChild(head);
      const participants=Array.isArray(record.participants)?record.participants:[];
      if(participants.length)host.appendChild(node('div',participants.map(p=>(p.name||'Participant')+(p.state?' · '+p.state:'')).join(' / ')));
      const messages=Array.isArray(record.last_public_messages)?record.last_public_messages:[];
      const ordered=messages.slice().sort((a,b)=>(Number(b.seq)||0)-(Number(a.seq)||0));
      if(ordered.length){
        const last=ordered[0];const text=String(last.text||'');
        const heading=(last.actor||'Participant')+(last.at?' · '+String(last.at):'');
        host.appendChild(node('div',heading,'activity-note'));
        host.appendChild(node('p',text.length>550?text.slice(0,547)+'…':text,'public-message'));
      }else host.appendChild(node('p','No public message recorded yet.','caption'));
      const usage=[];
      if(numeric(record.completed_responses))usage.push(record.completed_responses+' completed responses');
      usage.push(numeric(record.known_tokens)?record.known_tokens.toLocaleString()+' known tokens':'Token usage unavailable');
      if(numeric(record.usage_missing_activations)&&record.usage_missing_activations>0)usage.push(record.usage_missing_activations+' usage reports missing');
      host.appendChild(node('div',usage.join(' · '),'activity-usage'));
      if(numeric(record.instrument_experiment_seconds))host.appendChild(node('div','Recorded instrument runs: '+record.instrument_experiment_seconds.toFixed(1)+' s','activity-note'));
      const artifacts=Array.isArray(record.artifacts)?record.artifacts:[];
      if(artifacts.length){
        const details=node('details');details.appendChild(node('summary','Recent shared files'));
        const list=node('ul');for(const artifact of artifacts.slice(0,5))list.appendChild(node('li',String(artifact.path||'Shared file')+(artifact.revision!==undefined?' · revision '+artifact.revision:'')));
        details.appendChild(list);host.appendChild(details);
      }
      host.appendChild(node('div','Public messages and recorded files only. Token and instrument totals may be incomplete.','activity-note'));
    }
  }
  async function poll(){
    try{
      const response=await fetch('../results/confinement_activity.json',{cache:'no-store'});
      if(!response.ok){if(!seen&&response.status===404)return;throw new Error('Activity feed unavailable');}
      const data=await response.json();if(!data||!data.arms)return;
      render(data.arms);seen=true;
    }catch(error){
      if(seen)for(const arm of Object.keys(rooms)){
        const host=document.getElementById(arm+'-activity');
        if(!host.querySelector('.activity-stale'))host.appendChild(node('div','Activity feed unavailable; last recorded activity remains visible.','activity-note activity-stale'));
      }
    }finally{setTimeout(poll,3000);}
  }
  poll();
})();
