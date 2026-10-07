/* Recorded numerical sensitivity checks; never used to choose a primary. */
(() => {
  'use strict';
  const byId = id => document.getElementById(id);
  const number = value => typeof value === 'number' && Number.isFinite(value);
  const format = value => number(value) ? value.toFixed(2) : '—';
  let recorded = false;
  function render(data) {
    if (!data || !Array.isArray(data.rows)) throw new Error('Missing diagnostic rows');
    const body = byId('diagnostics-rows');
    body.replaceChildren();
    for (const row of data.rows) {
      const complete = row.complete === true && ['baseline', 'direct', 'construction'].every(key => number(row[key]));
      const tr = document.createElement('tr');
      const values = [row.setting || 'Unnamed check', Number.isInteger(row.particles) && row.particles > 0 ? String(row.particles) : '—', format(row.baseline), format(row.direct), format(row.construction), complete ? 'Recorded' : 'Incomplete · provisional values'];
      values.forEach((value, index) => {
        const td = document.createElement('td');
        td.textContent = value;
        if (index > 0 && index < 5) td.className = 'num';
        tr.appendChild(td);
      });
      body.appendChild(tr);
    }
    if (!data.rows.length) {
      const tr = document.createElement('tr'), td = document.createElement('td');
      td.colSpan = 6; td.className = 'muted'; td.textContent = 'No checks recorded yet.';
      tr.appendChild(td); body.appendChild(tr);
    }
    const complete = data.status === 'complete' && data.rows.length > 0 && data.rows.every(row => row.complete === true && ['baseline', 'direct', 'construction'].every(key => number(row[key])));
    byId('diagnostics-status').textContent = complete ? 'All listed checks recorded' : 'Checks in progress';
    byId('diagnostics-scope').textContent = (typeof data.scope === 'string' ? data.scope + ' ' : 'These diagnostics do not replace the locked primary comparison. ') + 'RM120 is measured in transits. Incomplete rows are provisional.';
    recorded = true;
  }
  async function poll() {
    try {
      const response = await fetch('../results/confinement_diagnostics.json', {cache: 'no-store'});
      if (!response.ok) throw new Error('Checks unavailable');
      render(await response.json());
    } catch (_) {
      byId('diagnostics-status').textContent = recorded ? 'Showing last recorded checks · update unavailable' : 'Awaiting recorded checks';
    }
    setTimeout(poll, 3000);
  }
  function loadFigure() {
    const probe = new Image();
    probe.onload = () => {
      byId('field-image').src = probe.src;
      byId('field-figure').hidden = false;
    };
    probe.onerror = () => { setTimeout(loadFigure, 15000); };
    probe.src = '../results/confinement_fields.png?check=' + Date.now();
  }
  poll();
  loadFigure();
})();
