import os
import subprocess
import sys
from pathlib import Path


def make_stub_modules(base: Path, modules: dict[str, str]) -> None:
    for mod, content in modules.items():
        parts = mod.split('.')
        pkg = base
        for part in parts[:-1]:
            pkg /= part
            pkg.mkdir(exist_ok=True)
            init_file = pkg / '__init__.py'
            if not init_file.exists():
                init_file.write_text('')
        (pkg / (parts[-1] + '.py')).write_text(content)


def run_help(
    script: Path, modules: dict[str, str], tmp_path: Path
) -> subprocess.CompletedProcess:
    stub_dir = tmp_path / 'stubs'
    stub_dir.mkdir()
    make_stub_modules(stub_dir, modules)
    env = os.environ.copy()
    env['PYTHONPATH'] = os.pathsep.join([
        str(stub_dir),
        env.get('PYTHONPATH', ''),
    ])
    return subprocess.run(
        [sys.executable, str(script), '--help'],
        capture_output=True,
        text=True,
        env=env,
    )


def test_poasta_graphviz_region_help(tmp_path: Path) -> None:
    script = Path(__file__).resolve().parents[2] / 'contrib/poasta_tools/poasta_graphviz_region.py'
    result = run_help(script, {'pygraphviz': ''}, tmp_path)
    assert result.returncode == 0
    combined_output = result.stdout + result.stderr
    assert 'usage' in combined_output.lower()


def test_poasta_plot_help(tmp_path: Path) -> None:
    script = Path(__file__).resolve().parents[2] / 'contrib/poasta_tools/poasta_plot.py'
    modules = {
        'numpy': 'class ndarray: ...\n',
        'pandas': '',
        'matplotlib': '',
        'matplotlib.pyplot': '',
        'seaborn': '',
        'networkx': '',
        'matplotlib.animation': (
            'class FuncAnimation: ...\n'
            'class FFMpegWriter: ...\n'
        ),
    }
    result = run_help(script, modules, tmp_path)
    assert result.returncode == 0
    combined_output = result.stdout + result.stderr
    assert 'usage' in combined_output.lower()
