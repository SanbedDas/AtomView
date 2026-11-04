#!/bin/zsh
# Double-click in Finder to start AtomView.
cd "$(dirname "$0")"
exec .venv/bin/python -m atomview "$@"
