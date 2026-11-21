(() => {
  if (location.origin !== 'https://baiming-zhang.github.io' ||
      !['/', '/index.html'].includes(location.pathname) ||
      navigator.globalPrivacyControl || navigator.doNotTrack === '1' ||
      window.doNotTrack === '1') return;

  const recordView = () => {
    if (document.visibilityState !== 'visible') return;
    document.removeEventListener('visibilitychange', recordView);

    const collector = 'https://baiming-orbit-workspace-counters.zhangbaiming3918.chatgpt.site/api/';
    const post = (path, payload) => fetch(collector + path, {
      method: 'POST',
      mode: 'cors',
      credentials: 'omit',
      referrerPolicy: 'no-referrer',
      keepalive: true,
      body: JSON.stringify(payload),
    });

    const id = crypto.randomUUID();
    const started = performance.now();
    const maxElapsedMs = 12 * 60 * 60 * 1000;
    const device = {
      referrer: document.referrer ? new URL(document.referrer).origin : '',
      locale: navigator.language,
      timezone: Intl.DateTimeFormat().resolvedOptions().timeZone,
      screen: `${screen.width} × ${screen.height}`,
      viewport: `${window.innerWidth} × ${window.innerHeight}`,
    };
    const sections = [...document.querySelectorAll('.section-anchor[id]')].filter(el => el.getClientRects().length);
    let locationEstimate = {};
    let previousSection = '';
    let visibleSince = started;
    let activeMs = 0;
    let scrollPercent = 0;
    let lastSentMs = -1;
    let lastSentScroll = -1;
    let locationChanged = false;
    let sequence = 0;
    let framePending = false;
    const queue = [];

    function event(kind, target) {
      if (!target || sequence >= 120) return;
      queue.push({ sequence: sequence++, elapsedMs: Math.min(maxElapsedMs, Math.round(performance.now() - started)),
        kind, target: target.slice(0, 250) });
    }

    function measureScroll() {
      const height = document.documentElement.scrollHeight - window.innerHeight;
      const fraction = height > 0 ? window.scrollY / height : 0;
      scrollPercent = Math.max(scrollPercent, Math.min(100, Math.round(fraction * 100)));
      const section = sections.filter(el => el.getBoundingClientRect().top <= window.innerHeight / 3).at(-1);
      if (section && section.id !== previousSection) {
        previousSection = section.id;
        event('section', '#' + section.id);
      }
    }

    function flush() {
      const visibleMs = visibleSince === null ? 0 : performance.now() - visibleSince;
      const duration = Math.min(maxElapsedMs, Math.round(activeMs + visibleMs));
      if (!queue.length && !locationChanged && duration === lastSentMs && scrollPercent === lastSentScroll) return;
      lastSentMs = duration;
      lastSentScroll = scrollPercent;
      locationChanged = false;
      do {
        const events = queue.splice(0, 12);
        post('homepage-activity', { id, device, location: locationEstimate, activeMs: duration, scrollPercent, events })
          .then(response => {
            if (!response.ok) throw new Error('Activity delivery failed');
          }).catch(() => {
            queue.unshift(...events);
            lastSentMs = -1;
          });
      } while (queue.length);
    }

    function stopVisibleTime() {
      if (visibleSince !== null) activeMs += performance.now() - visibleSince;
      visibleSince = null;
      flush();
    }

    document.addEventListener('visibilitychange', () => {
      if (document.visibilityState === 'visible') {
        visibleSince ??= performance.now();
        measureScroll();
      } else {
        stopVisibleTime();
      }
    });
    window.addEventListener('pagehide', stopVisibleTime);
    window.addEventListener('pageshow', () => {
      if (document.visibilityState === 'visible') visibleSince ??= performance.now();
    });
    window.addEventListener('scroll', () => {
      if (framePending || visibleSince === null) return;
      framePending = true;
      requestAnimationFrame(() => {
        framePending = false;
        if (visibleSince !== null) measureScroll();
      });
    }, { passive: true });

    document.addEventListener('click', captureClick, true);
    document.addEventListener('auxclick', captureClick, true);
    function captureClick(click) {
      if (!click.isTrusted || (click.type === 'auxclick' && click.button !== 1) || !(click.target instanceof Element)) return;
      const control = click.target.closest('a[href], button');
      if (!control) return;
      if (control.tagName === 'BUTTON') {
        const label = control.getAttribute('aria-label') || control.textContent.trim();
        event('button', label);
      } else {
        const url = new URL(control.href);
        if (['mailto:', 'tel:'].includes(url.protocol)) event('contact', url.protocol);
        else if (['http:', 'https:'].includes(url.protocol)) {
          if (url.origin === location.origin && url.pathname === location.pathname && url.hash) {
            event('anchor', url.hash);
          } else {
            const destination = url.origin + url.pathname;
            const download = control.hasAttribute('download') || /\.(pdf|zip)$/i.test(url.pathname);
            event(download ? 'download' : url.origin === location.origin ? 'site' : 'outbound', destination);
          }
        }
      }
      flush();
    }

    measureScroll();
    flush();
    setInterval(() => { if (visibleSince !== null || queue.length) flush(); }, 15_000);
    fetch('https://get.geojs.io/v1/ip/geo.json', {
      credentials: 'omit',
      referrerPolicy: 'no-referrer',
      signal: AbortSignal.timeout(2500),
    }).then(async response => {
      if (!response.ok) return;
      const { ip, city, region, country_code: country } = await response.json();
      locationEstimate = { ip, city, region, country };
      locationChanged = true;
      flush();
    }).catch(() => {
      // The visit still appears when location is unavailable.
    }).finally(() => {
      const { city, country } = locationEstimate;
      post('homepage-view', { city, country }).catch(() => {});
    });
  };

  document.addEventListener('visibilitychange', recordView);
  recordView();
})();
