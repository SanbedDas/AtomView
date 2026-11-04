; Windows installer: Start menu entry, optional desktop icon, .atomview sessions open on double-click,
; AtomView under "Open with" for particle files. Installs per user, no administrator rights needed.
;   iscc /DVersion=0.2.0 packaging\atomview.iss    (after the PyInstaller build in dist\AtomView)
#ifndef Version
  #define Version "0.2.0"
#endif
[Setup]
AppId={{6F1C4B7E-2A0B-4C1D-9E3A-41A7D0B7C0DE}
AppName=AtomView
AppVersion={#Version}
AppPublisher=AtomView
DefaultDirName={localappdata}\Programs\AtomView
DefaultGroupName=AtomView
PrivilegesRequired=lowest
ChangesAssociations=yes
OutputDir=..
OutputBaseFilename=AtomView-{#Version}-windows-setup
SetupIconFile=icon.ico
UninstallDisplayIcon={app}\AtomView.exe
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
DisableProgramGroupPage=yes

[Tasks]
Name: desktopicon; Description: "Create a desktop icon"; Flags: unchecked

[Files]
Source: "..\dist\AtomView\*"; DestDir: "{app}"; Flags: recursesubdirs ignoreversion

[Icons]
Name: "{userprograms}\AtomView"; Filename: "{app}\AtomView.exe"
Name: "{userdesktop}\AtomView"; Filename: "{app}\AtomView.exe"; Tasks: desktopicon

[Registry]
Root: HKCU; Subkey: "Software\Classes\.atomview"; ValueType: string; ValueData: "AtomView.Session"; Flags: uninsdeletevalue
Root: HKCU; Subkey: "Software\Classes\AtomView.Session"; ValueType: string; ValueData: "AtomView session"; Flags: uninsdeletekey
Root: HKCU; Subkey: "Software\Classes\AtomView.Session\DefaultIcon"; ValueType: string; ValueData: "{app}\AtomView.exe,0"
Root: HKCU; Subkey: "Software\Classes\AtomView.Session\shell\open\command"; ValueType: string; ValueData: """{app}\AtomView.exe"" ""%1"""
Root: HKCU; Subkey: "Software\Classes\AtomView.Data"; ValueType: string; ValueData: "Particle data"; Flags: uninsdeletekey
Root: HKCU; Subkey: "Software\Classes\AtomView.Data\DefaultIcon"; ValueType: string; ValueData: "{app}\AtomView.exe,0"
Root: HKCU; Subkey: "Software\Classes\AtomView.Data\shell\open\command"; ValueType: string; ValueData: """{app}\AtomView.exe"" ""%1"""
Root: HKCU; Subkey: "Software\Classes\.xyz\OpenWithProgids"; ValueType: string; ValueName: "AtomView.Data"; ValueData: ""; Flags: uninsdeletevalue
Root: HKCU; Subkey: "Software\Classes\.extxyz\OpenWithProgids"; ValueType: string; ValueName: "AtomView.Data"; ValueData: ""; Flags: uninsdeletevalue
Root: HKCU; Subkey: "Software\Classes\.dump\OpenWithProgids"; ValueType: string; ValueName: "AtomView.Data"; ValueData: ""; Flags: uninsdeletevalue
Root: HKCU; Subkey: "Software\Classes\.lammpstrj\OpenWithProgids"; ValueType: string; ValueName: "AtomView.Data"; ValueData: ""; Flags: uninsdeletevalue
Root: HKCU; Subkey: "Software\Classes\.lmp\OpenWithProgids"; ValueType: string; ValueName: "AtomView.Data"; ValueData: ""; Flags: uninsdeletevalue
Root: HKCU; Subkey: "Software\Classes\.mat\OpenWithProgids"; ValueType: string; ValueName: "AtomView.Data"; ValueData: ""; Flags: uninsdeletevalue
Root: HKCU; Subkey: "Software\Classes\.pdb\OpenWithProgids"; ValueType: string; ValueName: "AtomView.Data"; ValueData: ""; Flags: uninsdeletevalue
Root: HKCU; Subkey: "Software\Classes\.gro\OpenWithProgids"; ValueType: string; ValueName: "AtomView.Data"; ValueData: ""; Flags: uninsdeletevalue
Root: HKCU; Subkey: "Software\Classes\.xtc\OpenWithProgids"; ValueType: string; ValueName: "AtomView.Data"; ValueData: ""; Flags: uninsdeletevalue
Root: HKCU; Subkey: "Software\Classes\.dcd\OpenWithProgids"; ValueType: string; ValueName: "AtomView.Data"; ValueData: ""; Flags: uninsdeletevalue
Root: HKCU; Subkey: "Software\Classes\.cif\OpenWithProgids"; ValueType: string; ValueName: "AtomView.Data"; ValueData: ""; Flags: uninsdeletevalue
Root: HKCU; Subkey: "Software\Classes\.vasp\OpenWithProgids"; ValueType: string; ValueName: "AtomView.Data"; ValueData: ""; Flags: uninsdeletevalue

[Run]
Filename: "{app}\AtomView.exe"; Description: "Start AtomView"; Flags: nowait postinstall skipifsilent
