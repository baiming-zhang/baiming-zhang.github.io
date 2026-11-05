(() => {
  const number = /^[+−-]?(?:\d[\d,]*\.?\d*|\.\d+)(?:e[+−-]?\d+)?\s*%?$/i;
  function value(cell, part) {
    const text = (cell?.textContent || '').trim();
    const piece = (part == null ? text : text.split('/')[part] || '').trim();
    if (!piece || /^[-—–]$/.test(piece)) return null;
    return number.test(piece) ? Number(piece.replace(/[,\s%]/g, '').replace('−', '-')) : piece;
  }
  function init(table) {
    if (table.dataset.titleSortable) return;
    const headers = [...(table.tHead?.rows[table.tHead.rows.length - 1]?.cells || [])];
    if (!headers.length || !table.tBodies.length) return;
    table.dataset.titleSortable = 'true';
    headers.forEach((header, col) => {
      const label = header.textContent.trim();
      const paired = label.includes('平均') && label.includes('最差');
      header.textContent = '';
      if (paired) header.append(document.createTextNode(label.split(/[：:]/)[0] + '：'));
      (paired ? ['平均', '最差'] : [label]).forEach((name, part) => {
        if (part) header.append(document.createTextNode(' / '));
        const button = document.createElement('button');
        button.type = 'button'; button.className = 'm02-sort-title';
        button.textContent = name + ' ↕';
        button.title = '按' + label + (paired ? '（' + name + '）' : '') + '排序，再次点击切换升降序';
        button.setAttribute('aria-label', button.title);
        button.addEventListener('click', () => {
          const dir = button.dataset.direction === 'asc' ? 'desc' : 'asc';
          table.querySelectorAll('.m02-sort-title').forEach(b => {
            delete b.dataset.direction;
            b.textContent = b.textContent.replace(/ [↕↑↓]$/, ' ↕');
            b.removeAttribute('aria-pressed');
          });
          headers.forEach(h => h.removeAttribute('aria-sort'));
          button.dataset.direction = dir; button.setAttribute('aria-pressed', 'true');
          button.textContent = name + (dir === 'asc' ? ' ↑' : ' ↓');
          header.setAttribute('aria-sort', dir === 'asc' ? 'ascending' : 'descending');
          for (const body of table.tBodies) {
            const rows = [...body.rows];
            rows.map((row, index) => ({row, index, key: value(row.cells[col], paired ? part : null)}))
              .sort((a,b) => {
                if (a.key == null || b.key == null) return a.key == null ? (b.key == null ? a.index-b.index : 1) : -1;
                const cmp = typeof a.key === 'number' && typeof b.key === 'number' ? a.key-b.key : String(a.key).localeCompare(String(b.key), 'zh-CN', {numeric:true});
                return (dir === 'asc' ? cmp : -cmp) || a.index-b.index;
              }).forEach(item => body.append(item.row));
          }
        });
        header.append(button);
      });
    });
  }
  function scan() { document.querySelectorAll('table').forEach(init); }
  scan(); new MutationObserver(scan).observe(document.body, {childList:true, subtree:true});
})();
