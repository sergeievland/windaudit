import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import pytest

spec=importlib.util.spec_from_file_location('launch',Path(__file__).parents[1]/'spiral_lab/launch.py')
launch=importlib.util.module_from_spec(spec)
spec.loader.exec_module(launch)


@pytest.fixture
def lab(tmp_path,monkeypatch):
    repo=tmp_path/'villa';repo.mkdir()
    spiral=repo/'spiral-fitting';spiral.mkdir()
    (spiral/'uv.lock').write_text('test lock')
    subprocess.run(['git','init','-q',str(repo)],check=True)
    subprocess.run(['git','add','.'],cwd=repo,check=True)
    subprocess.run(['git','-c','user.name=Test','-c','user.email=test@example.invalid','commit','-qm','fixture'],cwd=repo,check=True)
    monkeypatch.setattr(launch,'PIN',subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip())
    py=spiral/'.venv/bin/python';py.parent.mkdir(parents=True)
    py.write_text(f'#!{sys.executable}\nimport sys,os\nprint("fake GPU probe" if sys.argv[1]=="-c" else "fake fit")\nsys.exit(int(os.environ.get("TEST_FIT_EXIT","0")) if sys.argv[1]!="-c" else 0)\n')
    py.chmod(0o755)
    data=tmp_path/'data';data.mkdir();(data/'verified_patches').mkdir()
    for name in ['umbilicus.json','abs_winding.json','relative_windings.json','same_windings.json']:
        (data/name).write_text('{}')
    (data/'spiral-scroll.json').write_text(json.dumps(dict(schema_version=1,voxel_size_um=9.6,spiral_outward_sense='CW')))
    config=tmp_path/'config.json';config.write_text(json.dumps(dict(z_begin=0,z_end=10,optimizer_num_training_steps=1,optimizer_random_seed=1)))
    return repo,data,config,tmp_path/'out'


def args(lab):
    return [v for k,p in zip(['--villa','--dataset','--config','--out'],lab) for v in [k,str(p)]]


def test_dry_run_does_not_create_output(lab):
    assert launch.main(args(lab)+['--dry-run'])==0
    assert not lab[-1].exists()


@pytest.mark.parametrize('code',[0,7])
def test_records_success_and_failure(lab,monkeypatch,code):
    monkeypatch.setenv('TEST_FIT_EXIT',str(code))
    assert launch.main(args(lab))==code
    result=json.loads((lab[-1]/'launch.json').read_text())
    assert result['returncode']==code
    assert result['status']==('completed' if code==0 else 'failed_or_interrupted')
    assert 'fake fit' in (lab[-1]/'console.log').read_text()


def test_never_overwrites_output(lab):
    lab[-1].mkdir();(lab[-1]/'precious').write_text('keep')
    with pytest.raises(ValueError,match='already exists'):launch.main(args(lab))
    assert (lab[-1]/'precious').read_text()=='keep'


def test_removes_inherited_fit_controls(monkeypatch):
    monkeypatch.setenv('FIT_SPIRAL_RESUME_PATH','wrong.ckpt')
    monkeypatch.setenv('WANDB_MODE','online')
    assert 'FIT_SPIRAL_RESUME_PATH' not in launch.clean_env()
    assert 'WANDB_MODE' not in launch.clean_env()
