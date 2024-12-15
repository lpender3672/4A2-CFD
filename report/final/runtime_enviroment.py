import os
import subprocess
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from postprocessing.routines import (
    write_settings,
    default_settings,
    read_settings,
    read_conv,
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

        self.log_file = self.worker_folder / f"worker_{self.worker_id}.log"
        
        generate_case(self.av['casename'], self.worker_folder)
        print(self.av)
        write_settings(self.av, self.worker_folder) # overwrite default settings

        # sleep
        time.sleep(1)

    def parse_results(self):

        convpath = self.target_folder / str("convergence_history.txt")

        conv_hist = read_conv(convpath)
        
        if conv_hist['nstep'].shape[0] > 0:
            iterations = conv_hist['nstep'][-1]
            dro_max = conv_hist['dro_max'][-1]
            dro_avg = conv_hist['dro_avg'][-1]
        else:
            iterations = 1
            dro_max = np.inf
            dro_avg = np.inf

        converged = check_run_converged(lines)

        newrow = np.array([
            cfl, sfac, time, converged, iterations, dro_max, dro_avg, av['ni'], av['nj']
        ])

    def run(self):

        self.setup_environment()
        with open(self.log_file, "w") as log:
            try:
                parsed_file = str(self.target_file).replace("\\", "/")
                subprocess.run([self.cfd_executable, "--path", parsed_file],
                                check=True,
                                stdout=log)
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
        self.base_folder.mkdir(exist_ok=True)
        self.cfd_executable = cfd_executable
        self.shared_file = self.base_folder / "results.txt"
        self.shared_lock = threading.Lock()

        self.avs = avs

    def clear_shared_file(self):
        if self.shared_file.exists():
            self.shared_file.unlink()

    def clear_worker_folders(self):

        for folder in self.base_folder.iterdir():
            if folder.is_dir():
                for subfolder in folder.iterdir():
                    for file in subfolder.iterdir():
                        file.unlink()
                    subfolder.rmdir()
                folder.rmdir()

    def start_workers(self):

        self.clear_shared_file()
        workers = [CFDWorker(i, self.base_folder, self.cfd_executable, self.avs[i-1], self.shared_file, self.shared_lock) for i in range(1, len(self.avs) + 1)]

        with ThreadPoolExecutor(max_workers=self.n_workers) as executor:
            futures = [executor.submit(worker.run) for worker in workers]
            for future in futures:
                future.result()

        print("All workers finished.")

if __name__ == "__main__":

    n_workers = 2
    base_folder = "case_env"
    cfd_executable = 'build/solverApp.exe'

    av = read_settings('cases/bump/input_bump.txt')
    av['nsteps'] = 500
    avs = [av for _ in range(10)]

    manager = CFDManager(n_workers, base_folder, cfd_executable, avs)
    
    manager.clear_worker_folders()
    #manager.start_workers()

