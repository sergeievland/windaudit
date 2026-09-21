#!/usr/bin/env python3
"""Pinned, non-destructive fit launcher. Uses an already installed villa environment."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

PIN = '2dcfaf6a08c3bc796fde726c4fa32050c8fc90e7'


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def clean_env():
    return {k: v for k, v in os.environ.items()
            if not k.startswith(('FIT_SPIRAL_', 'WANDB_'))}


def prepare(repo, dataset, config, output):
    repo, dataset, config, output = [Path(p).resolve() for p in (repo, dataset, config, output)]
    head = subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip()
    if head != PIN:
        raise ValueError(f'villa revision must be {PIN}; found {head}')
    if subprocess.check_output(['git','status','--porcelain','--untracked-files=no'],cwd=repo,text=True).strip():
        raise ValueError('baseline requires a clean tracked villa checkout')
    if output.exists():
        raise ValueError('output already exists; choose a new directory (never overwritten)')
    if dataset == output or dataset in output.parents or output in dataset.parents:
        raise ValueError('output and dataset must be separate trees')
    settings = json.loads(config.read_text())
    if not isinstance(settings, dict):
        raise ValueError('config must be a JSON object')
    for key in ('z_begin','z_end','optimizer_num_training_steps','optimizer_random_seed'):
        if type(settings.get(key)) is not int:
            raise ValueError(f'explicit integer {key} is required')
    if not 0 <= settings['z_begin'] < settings['z_end'] or settings['optimizer_num_training_steps'] <= 0 or settings['optimizer_random_seed'] < 0:
        raise ValueError('invalid z range, steps or seed')
    required = ('spiral-scroll.json','umbilicus.json','abs_winding.json',
                'relative_windings.json','same_windings.json')
    hashes = {name:digest(dataset/name) for name in required}
    spec = json.loads((dataset/'spiral-scroll.json').read_text())
    if spec.get('schema_version') != 1 or spec.get('spiral_outward_sense','').upper() not in ('CW','ACW'):
        raise ValueError('invalid spiral-scroll.json schema or sense')
    if not isinstance(spec.get('voxel_size_um'), (int,float)) or spec['voxel_size_um'] <= 0:
        raise ValueError('positive voxel_size_um required')
    if not (dataset/'verified_patches').is_dir():
        raise ValueError('verified_patches directory missing')
    spiral = repo/'spiral-fitting'
    python = spiral/'.venv/bin/python'
    if not python.is_file():
        raise ValueError('run setup.sh first: villa Python environment is missing')
    manifest = dict(villa_commit=head, config=settings, annotation_sha256=hashes,
                    dataset=str(dataset), dataset_fully_hashed=False,
                    interpretation='Exploratory fit; training satisfaction is not held-out accuracy.',
                    lock_sha256=digest(spiral/'uv.lock'))
    return spiral, python, manifest, output


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--villa', required=True)
    p.add_argument('--dataset', required=True)
    p.add_argument('--config', default=str(Path(__file__).with_name('smoke.json')))
    p.add_argument('--out', required=True)
    p.add_argument('--dry-run', action='store_true')
    a = p.parse_args(argv)
    spiral, python, manifest, output = prepare(a.villa,a.dataset,a.config,a.out)
    cmd = [str(python),str(spiral/'fit_spiral.py'),'--dataset',str(Path(a.dataset).resolve())]
    manifest['command'] = cmd
    if a.dry_run:
        print(json.dumps(manifest,indent=2)); return 0
    # Validate upstream option names and allocate CUDA before creating a run.
    probe = ('import json,sys,torch; from config import Config; '
             'Config(json.loads(sys.argv[1])); '
             'assert torch.cuda.is_available(), "CUDA unavailable"; '
             'x=torch.empty(1,device="cuda"); torch.cuda.synchronize(); '
             'print(json.dumps({"torch":torch.__version__,"cuda":torch.version.cuda,'
             '"gpu":torch.cuda.get_device_name(0),"vram_bytes":torch.cuda.get_device_properties(0).total_memory}))')
    env=clean_env()
    probe_result=subprocess.check_output([str(python),'-c',probe,json.dumps(manifest['config'])],cwd=spiral,env=env,text=True)
    manifest['gpu_probe']=probe_result.strip()
    output.mkdir(parents=True,exist_ok=False)
    manifest['status']='running'
    manifest_path=output/'launch.json'
    manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
    env.update(FIT_SPIRAL_CONFIG_OVERRIDES=json.dumps(manifest['config']),
               FIT_SPIRAL_RUN_DIR=str(output/'fit'), WANDB_MODE='disabled',
               PYTHONUNBUFFERED='1')
    started=time.monotonic()
    code=1
    try:
        with open(output/'console.log','w') as log:
            child=subprocess.Popen(cmd,cwd=spiral,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
            try:
                for line in child.stdout:
                    print(line,end='',flush=True);log.write(line);log.flush()
                code=child.wait()
            finally:
                if child.poll() is None:
                    child.terminate()
                    try: child.wait(timeout=30)
                    except subprocess.TimeoutExpired: child.kill();child.wait()
                child.stdout.close()
    finally:
        manifest.update(status='completed' if code == 0 else 'failed_or_interrupted',
                        returncode=code,elapsed_seconds=round(time.monotonic()-started,2))
        manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
    return code

if __name__ == '__main__':
    try: raise SystemExit(main())
    except (ValueError,OSError,subprocess.CalledProcessError) as e:
        raise SystemExit(str(e))
