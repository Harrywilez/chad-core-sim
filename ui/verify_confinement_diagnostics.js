'use strict';
const fs = require('fs'), vm = require('vm'), assert = require('assert'), path = require('path');
function element() { return {children: [], hidden: true, textContent: '', appendChild(x) { this.children.push(x); }, replaceChildren() { this.children = []; }}; }
async function render(data) {
  const nodes = new Map(), images = [];
  const document = {getElementById(id) { if (!nodes.has(id)) nodes.set(id, element()); return nodes.get(id); }, createElement: element};
  function Image() { images.push(this); }
  vm.runInNewContext(fs.readFileSync(path.join(__dirname, 'confinement_diagnostics.js'), 'utf8'), {document, Image, fetch: async () => ({ok: true, json: async () => data}), Date, setTimeout: () => {}});
  await new Promise(resolve => setImmediate(resolve));
  return {nodes, images, document};
}
(async () => {
  const row = {setting: 'Finer field grid', particles: 480, baseline: 20, direct: 30, construction: 40, complete: true};
  const good = await render({status: 'complete', rows: [row], scope: 'Matched particles; diagnostics do not replace primary results.'});
  assert.strictEqual(good.nodes.get('diagnostics-rows').children[0].children[4].textContent, '40.00');
  assert.strictEqual(good.nodes.get('diagnostics-status').textContent, 'All listed checks recorded');
  assert.strictEqual(good.document.getElementById('field-figure').hidden, true);
  good.images[0].onerror();
  assert.strictEqual(good.nodes.get('field-figure').hidden, true);
  good.images[0].onload();
  assert.strictEqual(good.nodes.get('field-figure').hidden, false);
  const partial = await render({status: 'complete', rows: [{...row, construction: null}, {...row, construction: 119, complete: false}]});
  assert.strictEqual(partial.nodes.get('diagnostics-status').textContent, 'Checks in progress');
  const rows = partial.nodes.get('diagnostics-rows').children;
  assert.strictEqual(rows[0].children[4].textContent, '—');
  assert.strictEqual(rows[1].children[4].textContent, '119.00');
  assert(rows.every(r => r.children[5].textContent.startsWith('Incomplete')));
  assert(rows.every(r => r.children.every(c => !c.className || c.className === 'num')));
  console.log('Passed: recorded values, missing versus zero, incomplete rows remain provisional without winner styling, and figure appears only after successful image load.');
})();
