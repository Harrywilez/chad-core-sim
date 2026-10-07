// Offline render test: node shot.js [family] [out.png]   (NODE_PATH=$(npm root -g))
const { chromium } = require('playwright');
(async () => {
  const family = process.argv[2] || null, out = process.argv[3] || 'shot.png';
  const browser = await chromium.launch({ args: ['--use-gl=swiftshader','--enable-webgl','--ignore-gpu-blocklist'] });
  const page = await browser.newPage({ viewport: { width: 1500, height: 900 } });
  const logs=[]; page.on('console', m => logs.push(m.type()+': '+m.text())); page.on('pageerror', e => logs.push('PAGEERROR: '+e.message));
  await page.goto('file://' + process.cwd() + '/chad_core_lab_preview.html');
  await page.waitForTimeout(family ? 4000 : 25000);
  if (family) {
    await page.evaluate(f => { const s=document.getElementById('family'); s.value=f; s.dispatchEvent(new Event('change')); document.getElementById('gridn').value='17'; document.getElementById('build').click(); document.getElementById('compute').click(); }, family);
    await page.waitForTimeout(20000);
    await page.evaluate(() => document.getElementById('lines').click());
    await page.waitForTimeout(8000);
    await page.evaluate(() => { document.getElementById('npart').value='48'; document.getElementById('scale').value='reactor'; document.getElementById('species').value='D'; document.getElementById('run').click(); });
    await page.waitForTimeout(25000);
  }
  await page.screenshot({ path: out });
  console.log(logs.slice(0,30).join('\n'));
  console.log(await page.evaluate(() => document.getElementById('status').textContent + ' | ' + document.getElementById('geo-note').textContent + ' | ' + document.getElementById('ro-lines').textContent.slice(0,300)));
  await browser.close();
})();
