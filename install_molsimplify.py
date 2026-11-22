"""Install only molSimplify's explicit-donor geometry backend."""
import os
from pathlib import Path
import subprocess
import sys


def main():
    base = Path(__file__).resolve().parent
    interpreter = base / ".molsimplify" / "Scripts" / "python.exe"
    if not interpreter.is_file():
        candidates = [[sys.executable]]
        candidates.extend([["py", f"-{version}"] for version in ("3.12", "3.11", "3.10")])
        for candidate in candidates:
            try:
                valid = subprocess.run(candidate + ["-c", "import sys; raise SystemExit(not ((3,10) <= sys.version_info[:2] <= (3,12)))"],
                                       capture_output=True).returncode == 0
                if valid:
                    subprocess.run(candidate + ["-m", "venv", str(base / ".molsimplify")], check=True)
                    break
            except OSError:
                continue
        else:
            raise RuntimeError("Please install Python 3.10-3.12 x64 first. The main application can keep Python 3.13.")
    pip = [str(interpreter), "-m", "pip", "install"]
    # PIP_INDEX_URL may be configured by the user; no mirror is hard-coded.
    subprocess.run(pip + ["-r", str(base / "requirements-molsimplify.txt")], check=True)
    # ML training, donor prediction and the discontinued upstream GUI are not used.
    subprocess.run(pip + ["--no-deps", "molSimplify==2.0.0"], check=True)
    env = os.environ.copy()
    env["PYTHONUTF8"] = "1"
    subprocess.run([str(interpreter), "-c", "from molSimplify.Scripts.generator import startgen_pythonic; from openbabel import openbabel; print('molSimplify geometry backend ready')"],
                   check=True, env=env)
    print("Done. Open start.bat and choose the molSimplify complex editor.")


if __name__ == "__main__":
    try:
        main()
    except (OSError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"Installation failed: {exc}", file=sys.stderr)
        raise SystemExit(1)
