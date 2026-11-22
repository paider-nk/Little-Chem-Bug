"""Launch external tools with their own DLLs on Windows."""
import os
import subprocess
import sys
import threading
from subprocess import DEVNULL, STDOUT, TimeoutExpired

CREATE_NO_WINDOW = getattr(subprocess, "CREATE_NO_WINDOW", 0)

_launch_lock = threading.Lock()


class Popen(subprocess.Popen):
    def __init__(self, *args, **kwargs):
        if os.name != "nt" or not getattr(sys, "frozen", False):
            super().__init__(*args, **kwargs)
            return
        import ctypes

        with _launch_lock:
            dll_dir = ctypes.windll.kernel32.SetDllDirectoryW
            dll_dir(None)
            try:
                super().__init__(*args, **kwargs)
            finally:
                dll_dir(sys._MEIPASS)


def run(args, *, input=None, capture_output=False, timeout=None, check=False, **kwargs):
    if not getattr(sys, "frozen", False):
        return subprocess.run(args, input=input, capture_output=capture_output,
                              timeout=timeout, check=check, **kwargs)
    if input is not None:
        kwargs["stdin"] = subprocess.PIPE
    if capture_output:
        kwargs.update(stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    with Popen(args, **kwargs) as process:
        try:
            stdout, stderr = process.communicate(input, timeout=timeout)
        except BaseException:
            process.kill()
            process.communicate()
            raise
        code = process.returncode
        if check and code:
            raise subprocess.CalledProcessError(code, args, stdout, stderr)
        return subprocess.CompletedProcess(args, code, stdout, stderr)
