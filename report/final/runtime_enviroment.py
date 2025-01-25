import os
import subprocess
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from postprocessing.routines import (
    write_settings,
    default_settings,
)
from preprocessing.generate_case import generate_case

class CFDWorker:
    def __init__(self, worker_id, base_folder, cfd_executable, av, shared_file, shared_lock):
        self.worker_id = worker_id
        self.base_folder = Path(base_folder)
        self.cfd_executable = cfd_executable
        self.shared_file = Path(shared_file)
        self.shared_lock = shared_lock

        self.av = av

    def setup_environment(self):
        self.worker_folder = self.base_folder / str(self.worker_id)
        self.worker_folder.mkdir(exist_ok=True)

        self.target_folder = self.worker_folder / self.av['casename']
        self.target_folder.mkdir(exist_ok=True)

        self.target_file = self.target_folder / f"input_{self.av['casename']}.txt"
        
        #write_settings(self.av, self.target_file)
        generate_case(self.av['casename'], self.worker_folder)

        # sleep
        time.sleep(1)

    def run(self):

        self.setup_environment()
        try:
            subprocess.run([self.cfd_executable, "--path", str(self.target_file).replace("\\", "/")], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Worker {self.worker_id} failed with error: {e}")
            return

        results = f"Worker {self.worker_id} results: Success\n"

        # Write results to the shared file in a thread-safe way
        with self.shared_lock:
            with open(self.shared_file, "a") as f:
                f.write(results)

        print(f"Worker {self.worker_id} finished.")

class CFDManager:
    def __init__(self, n_workers, base_folder, cfd_executable, avs):
        self.n_workers = n_workers
        self.base_folder = Path(base_folder)
        self.cfd_executable = cfd_executable
        self.shared_file = self.base_folder / "results.txt"
        self.shared_lock = threading.Lock()

        self.avs = avs

    def clear_shared_file(self):

        if self.shared_file.exists():
            self.shared_file.unlink()

    def start_workers(self):

        self.clear_shared_file()
        workers = [CFDWorker(i, self.base_folder, self.cfd_executable, self.avs[i-1], self.shared_file, self.shared_lock) for i in range(1, len(self.avs) + 1)]

        with ThreadPoolExecutor(max_workers=self.n_workers) as executor:
            futures = [executor.submit(worker.run) for worker in workers]
            for future in futures:
                future.result()

        print("All workers finished.")

if __name__ == "__main__":

    n_workers = 4
    base_folder = Path(os.getcwd()) / "case_env"
    cfd_executable = 'build/solverApp.exe'

    avs = [default_settings('bump')]

    manager = CFDManager(n_workers, base_folder, cfd_executable, avs)
    manager.start_workers()

