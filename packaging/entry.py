import multiprocessing, os, sys
from atomview.cli import main

if __name__ == "__main__":
    multiprocessing.freeze_support()  # parallel batch workers in the frozen app
    os.environ.setdefault("QT_API", "pyside6")
    sys.exit(main())
