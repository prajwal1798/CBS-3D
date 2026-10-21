window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex"
  }
};

document$.subscribe(() => {
  if (!window.MathJax || !MathJax.typesetPromise) {
    return;
  }

  if (MathJax.startup && MathJax.startup.output) {
    MathJax.startup.output.clearCache();
  }

  MathJax.typesetClear();
  MathJax.texReset();
  MathJax.typesetPromise();
});
