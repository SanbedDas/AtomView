#!/bin/sh
# Adds AtomView to the application menu and opens .atomview sessions on double-click (current user only).
# Run from the unpacked AtomView folder:  ./install.sh      Remove again:  ./install.sh --uninstall
set -e
here=$(cd "$(dirname "$0")" && pwd)
data=${XDG_DATA_HOME:-$HOME/.local/share}
if [ "$1" = "--uninstall" ]; then
  rm -f "$data/applications/atomview.desktop" "$data/mime/packages/atomview.xml" "$data/icons/hicolor/256x256/apps/atomview.png"
else
  mkdir -p "$data/applications" "$data/mime/packages" "$data/icons/hicolor/256x256/apps"
  icon=$(find "$here" -path "*atomview/resources/icon.png" | head -n 1)
  [ -n "$icon" ] && cp "$icon" "$data/icons/hicolor/256x256/apps/atomview.png"
  cat > "$data/mime/packages/atomview.xml" <<XML
<?xml version="1.0" encoding="UTF-8"?>
<mime-info xmlns="http://www.freedesktop.org/standards/shared-mime-info">
  <mime-type type="application/x-atomview-session">
    <comment>AtomView session</comment>
    <sub-class-of type="application/json"/>
    <glob pattern="*.atomview"/>
    <icon name="atomview"/>
  </mime-type>
</mime-info>
XML
  cat > "$data/applications/atomview.desktop" <<DESKTOP
[Desktop Entry]
Type=Application
Name=AtomView
Comment=Visualize and analyze particle simulations
Exec="$here/AtomView" %f
Icon=atomview
Terminal=false
Categories=Science;Education;Graphics;
MimeType=application/x-atomview-session;chemical/x-xyz;chemical/x-pdb;
DESKTOP
fi
command -v update-mime-database >/dev/null && update-mime-database "$data/mime" || true
command -v update-desktop-database >/dev/null && update-desktop-database "$data/applications" || true
[ "$1" = "--uninstall" ] || { command -v xdg-mime >/dev/null && xdg-mime default atomview.desktop application/x-atomview-session || true; }
echo "AtomView ${1:+un}installed for $(whoami)"
