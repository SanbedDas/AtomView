#!/usr/bin/env bash
# Start AtomView on Linux:  ./atomview.sh [file]   or   ./atomview.sh batch FOLDER
cd "$(dirname "$0")"
exec .venv/bin/python -m atomview "$@"
