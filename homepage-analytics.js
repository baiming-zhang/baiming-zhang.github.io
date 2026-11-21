(() => {
  if (location.origin !== 'https://baiming-zhang.github.io' ||
      !['/', '/index.html'].includes(location.pathname) ||
      navigator.globalPrivacyControl || navigator.doNotTrack === '1' ||
      window.doNotTrack === '1') return;

  const recordView = async () => {
    if (document.visibilityState !== 'visible') return;
    document.removeEventListener('visibilitychange', recordView);

    let locationEstimate = {};
    try {
      const response = await fetch('https://get.geojs.io/v1/ip/geo.json', {
        credentials: 'omit',
        referrerPolicy: 'no-referrer',
        signal: AbortSignal.timeout(2500),
      });
      if (response.ok) {
        const { city, country_code: country } = await response.json();
        locationEstimate = { city: typeof city === 'string' ? city.trim().slice(0, 100) : '', country };
      }
    } catch {
      // A blocked or unavailable geolocation service leaves the city unknown.
    }

    fetch('https://baiming-orbit-workspace-counters.zhangbaiming3918.chatgpt.site/api/homepage-view', {
      method: 'POST',
      mode: 'cors',
      credentials: 'omit',
      referrerPolicy: 'no-referrer',
      keepalive: true,
      body: JSON.stringify(locationEstimate),
    }).catch(() => {});
  };

  document.addEventListener('visibilitychange', recordView);
  recordView();
})();
