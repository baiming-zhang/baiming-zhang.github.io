(() => {
  if (location.origin !== 'https://baiming-zhang.github.io' ||
      !['/', '/index.html'].includes(location.pathname) ||
      navigator.globalPrivacyControl || navigator.doNotTrack === '1' ||
      window.doNotTrack === '1') return;

  const recordView = () => {
    if (document.visibilityState !== 'visible') return;
    document.removeEventListener('visibilitychange', recordView);
    fetch('https://baiming-orbit-workspace-counters.zhangbaiming3918.chatgpt.site/api/homepage-view', {
      method: 'POST',
      mode: 'cors',
      credentials: 'omit',
      referrerPolicy: 'no-referrer',
      keepalive: true,
    }).catch(() => {});
  };

  document.addEventListener('visibilitychange', recordView);
  recordView();
})();
