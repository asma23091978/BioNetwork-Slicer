from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]

REQUIRED_FILES = [
    ROOT / "slicing_script.py",
    ROOT / "requirements.txt",
    ROOT / "prolog" / "hgd0_04_02_reactions.pl",
]

REQUIRED_TO_RUN = [
    ROOT / "BioReSolvePositive.pl",
    ROOT / "output" / "hgd0_04-02.dot",
]

# Built this way so the checker does not flag its own source code.
FORBIDDEN_PATTERNS = [
    "/" + "Users/",
    "C:" + "\\\\Users\\\\",
]


def check_file(path: Path) -> bool:
    ok = path.exists()
    print(f"{'OK' if ok else 'MISSING'}  {path.relative_to(ROOT)}")
    return ok


def scan_for_private_paths() -> bool:
    bad = []
    for path in ROOT.rglob("*"):
        if path.is_file() and path.suffix in {".py", ".pl", ".md", ".txt"}:
            try:
                text = path.read_text(encoding="utf-8", errors="ignore")
            except Exception:
                continue
            for pattern in FORBIDDEN_PATTERNS:
                if pattern in text:
                    bad.append((path.relative_to(ROOT), pattern))
    if bad:
        print("\nPrivate/local paths found:")
        for path, pattern in bad:
            print(f"  {path}: {pattern}")
        return False
    print("OK  no private/local paths found")
    return True


def check_prolog_engine_file() -> bool:
    engine = ROOT / "BioReSolvePositive.pl"
    if not engine.exists():
        print("MISSING  BioReSolvePositive.pl")
        print("         This file must contain the Prolog slicing engine, including main_do(...).")
        return False
    text = engine.read_text(encoding="utf-8", errors="ignore")
    ok = "main_do" in text
    print(f"{'OK' if ok else 'MISSING'}  main_do predicate text found in BioReSolvePositive.pl")
    return ok


def main() -> None:
    print("Repository check\n================")
    ok_required = all(check_file(p) for p in REQUIRED_FILES)

    print("\nFiles needed to run the current script")
    ok_runtime_files = all(check_file(p) for p in REQUIRED_TO_RUN)

    print("\nProlog engine check")
    ok_engine = check_prolog_engine_file()

    print("\nPath hygiene check")
    ok_paths = scan_for_private_paths()

    print("\nExternal executable check")
    swipl = shutil.which("swipl")
    print(f"{'OK' if swipl else 'MISSING'}  SWI-Prolog executable 'swipl'")

    print("\nSummary")
    if ok_required and ok_runtime_files and ok_engine and ok_paths:
        print("OK: the repository structure looks runnable.")
    else:
        print("NOT READY TO RUN: see the missing items above.")
        print("You can still publish the repository as a clean draft, but add the missing runtime files before claiming it is reproducible.")


if __name__ == "__main__":
    main()
