#!/bin/sh
# Build comment.pdf from comment.md via pandoc (HTML+MathML) + headless Chromium.
set -e
pandoc comment.md -s --mathml --embed-resources -c style.css \
  --metadata pagetitle="Comment on Tamo Tchomgui et al. (2026)" -o comment.html
# Chromium: prefer one on PATH, else the preinstalled Playwright build
CHROME=$(command -v chromium || command -v chromium-browser || command -v google-chrome || true)
[ -z "$CHROME" ] && CHROME=$(ls /opt/pw-browsers/chromium-*/chrome-linux/chrome 2>/dev/null | head -1)
"$CHROME" --headless --no-sandbox --disable-gpu --no-pdf-header-footer \
  --print-to-pdf=comment.pdf "file://$(pwd)/comment.html" 2>/dev/null
echo "built: $(ls -la comment.pdf)"
