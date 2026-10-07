// Offline render test of the Designs tab: node shot_designs.js [rowFilter] [out.png]   (NODE_PATH=$(npm root -g))
const { chromium } = require('playwright');
(async () => {
  const filter = process.argv[2] || '', out = process.argv[3] || 'shot_designs.png';
  const browser = await chromium.launch({ args: ['--use-gl=swiftshader','--enable-webgl','--ignore-gpu-blocklist'] });
  const page = await browser.newPage({ viewport: { width: 1500, height: 900 } });
  const logs=[]; page.on('console', m => logs.push(m.type()+': '+m.text())); page.on('pageerror', e => logs.push('PAGEERROR: '+e.message));
  await page.goto('file://' + process.cwd() + '/chad_core_lab_preview.html');
  await page.waitForTimeout(6000);
  await page.evaluate(() => document.querySelector('.tabs button[data-tab=designs]').click());
  await page.waitForTimeout(500);
  const info = await page.evaluate(f => { const el=document.getElementById('dfilter'); el.value=f; el.dispatchEvent(new Event('input')); const rows=[...document.querySelectorAll('#designs tr[data-id]')]; return {n:rows.length, first:rows[0]?rows[0].textContent.slice(0,120):null, feed:document.getElementById('feed-status').textContent}; }, filter);
  console.log(JSON.stringify(info));
  await page.evaluate(() => { const r=document.querySelector('#designs tr[data-id]'); if(r) r.click(); });
  await page.waitForTimeout(3000);
  await page.evaluate(() => { document.getElementById('gridn').value='17'; document.getElementById('compute').click(); });
  await page.waitForTimeout(25000);
  await page.evaluate(() => document.getElementById('lines').click());
  await page.waitForTimeout(8000);
  await page.screenshot({ path: out });
  await page.evaluate(() => document.querySelector('.tabs button[data-tab=designs]').click());
  await page.waitForTimeout(400);
  await page.screenshot({ path: out.replace('.png','_tab.png') });
  console.log(logs.slice(0,30).join('\n'));
  console.log(await page.evaluate(() => document.getElementById('status').textContent + ' | ' + document.getElementById('geo-note').textContent + ' | ' + document.getElementById('ro-design').textContent.slice(0,400)));
  await browser.close();
})();
