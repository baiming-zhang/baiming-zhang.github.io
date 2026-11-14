"""Build a bounded Pages site while retaining large immutable resources in Git.

No source resources are deleted. Large binary downloads/field arrays are served
from their original public repository paths pinned to the deployment commit.
"""
from pathlib import Path
import argparse, json, os, re, shutil, urllib.parse

THRESHOLD = 1_000_000
MAX_SITE_BYTES = 900_000_000
LARGE_EXTENSIONS = {'.gz', '.pt', '.npz', '.zip', '.bin', '.ply'}

ROUTER = r'''(() => {
  'use strict';
  const external = __RESOURCE_MAP__;
  function route(input) {
    try {
      const url = new URL(String(input), document.baseURI);
      if (url.origin !== location.origin) return String(input);
      const path = decodeURIComponent(url.pathname).replace(/^\//, '');
      return external[path] || String(input);
    } catch (_) { return String(input); }
  }
  const nativeFetch = window.fetch;
  window.fetch = function(input, init) {
    const request = input instanceof Request;
    const oldUrl = request ? input.url : String(input);
    const target = route(oldUrl);
    if (target === oldUrl) return nativeFetch.call(this, input, init);
    const options = Object.assign({}, init, {mode:'cors', credentials:'omit'});
    const headers = new Headers(options.headers || (request ? input.headers : undefined));
    headers.delete('Authorization'); headers.delete('Cookie'); options.headers = headers;
    return nativeFetch.call(this, request ? new Request(target, input) : target, options);
  };
  const nativeOpen = XMLHttpRequest.prototype.open;
  XMLHttpRequest.prototype.open = function(method, url, ...rest) {
    return nativeOpen.call(this, method, route(url), ...rest);
  };
  function rewriteLinks(root) {
    const links = root.querySelectorAll ? root.querySelectorAll('a[href]') : [];
    for (const anchor of links) {
      const target = route(anchor.getAttribute('href'));
      if (target !== anchor.getAttribute('href')) anchor.setAttribute('href', target);
    }
    if (root.matches && root.matches('a[href]')) {
      const target = route(root.getAttribute('href'));
      if (target !== root.getAttribute('href')) root.setAttribute('href', target);
    }
  }
  function start() {
    rewriteLinks(document);
    new MutationObserver(records => {
      for (const record of records) {
        if (record.type === 'attributes') rewriteLinks(record.target);
        else for (const node of record.addedNodes) if (node.nodeType === 1) rewriteLinks(node);
      }
    }).observe(document.documentElement, {subtree:true, childList:true, attributes:true, attributeFilter:['href']});
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', start, {once:true});
  else start();
})();
'''

def build(source, output, repository, commit):
    source, output = Path(source).resolve(), Path(output).resolve()
    if source == output or source in output.parents:
        raise ValueError('Output must be outside source to prevent recursive packaging')
    if output.exists():
        raise ValueError('Use a new output directory; do not delete or overwrite unrelated files')
    if not re.fullmatch(r'[0-9a-f]{40}', commit):
        raise ValueError('A full immutable source commit is required')
    output.mkdir(parents=True)
    resources, copied, total_source = {}, [], 0
    for path in sorted(source.rglob('*')):
        relative = path.relative_to(source)
        if any(part in {'.git', '.github', '__pycache__'} for part in relative.parts):
            continue
        if path.is_symlink():
            raise ValueError('Pages artifact must not contain symlinks: ' + str(relative))
        if not path.is_file():
            continue
        size = path.stat().st_size
        total_source += size
        name = relative.as_posix()
        # Archive chunks must remain bounded individually for scientific Git
        # provenance; many sub-MB chunks can still exhaust Pages together.
        # Route ALL ZIPs by the same immutable, byte-preserving mechanism.
        # This report retains many small scientific assets. Keep its HTML in
        # Pages and serve its original data/figures immutably, without relaxing
        # the shared site's size guard or changing unrelated project resources.
        physical_report_asset = name.startswith('VariDrainage/edge-physics/') and path.suffix.lower() != '.html'
        # Scientific JSON ledgers are original downloadable datasets, too.
        # Reuse the existing anchor/fetch/XHR router and retain byte-exact Git
        # originals rather than letting adjacent reports exhaust the site.
        scientific_ledger = name.startswith('VariDrainage/') and path.suffix.lower() == '.json' and size >= 100_000
        # Keep the latest reference report and its inline PNGs in Pages. Its
        # downloadable evidence (sources, logs, audit ledgers and vector PDFs)
        # uses the same commit-pinned router, preserving all original bytes.
        reference_download = name.startswith('VariDrainage/retained-reference/') and path.suffix.lower() not in {'.html', '.png'}
        if physical_report_asset or scientific_ledger or reference_download or path.suffix.lower() == '.zip' or (size >= THRESHOLD and path.suffix.lower() in LARGE_EXTENSIONS):
            resources[name] = ('https://raw.githubusercontent.com/' + repository + '/' + commit + '/' + urllib.parse.quote(name, safe='/'))
            continue
        target = output / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        copied.append(target)
    tag = '<script src="/pages-resource-router.js?v=' + commit[:12] + '"></script>'
    html_count = 0
    for path in copied:
        if path.suffix.lower() != '.html':
            continue
        raw = path.read_bytes()
        text = raw.decode('utf-8')
        if path.relative_to(output).as_posix().startswith('VariDrainage/edge-physics/'):
            # img src is resolved before the anchor/fetch router runs. Rewrite
            # these report-local image URLs at build time as well.
            import posixpath
            parent = path.relative_to(output).parent.as_posix()
            def image_resource(match):
                name = posixpath.normpath(posixpath.join(parent, match.group(2)))
                return match.group(1) + resources.get(name, match.group(2)) + match.group(3)
            text = re.sub(r'(<img\b[^>]*\bsrc=["\x27])([^"\x27]+)(["\x27])', image_resource, text, flags=re.I)
        # Preserve exact source bytes apart from this routing tag.
        if re.search(r'<head\b[^>]*>', text, flags=re.I):
            text = re.sub(r'(<head\b[^>]*>)', lambda m: m.group(1) + tag, text, count=1, flags=re.I)
        else:
            text = tag + text
        path.write_bytes(text.encode('utf-8'))
        html_count += 1
    router = ROUTER.replace('__RESOURCE_MAP__', json.dumps(resources, ensure_ascii=True, separators=(',', ':')))
    (output / 'pages-resource-router.js').write_bytes(router.encode('utf-8'))
    (output / '.nojekyll').write_bytes(b'')
    report = dict(source_commit=commit, repository=repository, source_bytes=total_source,
                  external_resource_count=len(resources), html_pages=html_count,
                  resource_policy='Public original files remain in Git; immutable raw URLs; no numerical changes',
                  max_site_bytes=MAX_SITE_BYTES)
    report['site_bytes_before_report'] = sum(p.stat().st_size for p in output.rglob('*') if p.is_file())
    (output / 'pages-build-report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    total = sum(p.stat().st_size for p in output.rglob('*') if p.is_file())
    if total >= MAX_SITE_BYTES:
        raise ValueError(f'Pages payload {total:,} bytes exceeds {MAX_SITE_BYTES:,}-byte safety budget')
    print(json.dumps({**report, 'site_bytes':total}, indent=2))
    return report

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--repository', required=True)
    parser.add_argument('--commit', required=True)
    args = parser.parse_args()
    build(args.source, args.output, args.repository, args.commit)
