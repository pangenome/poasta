import sys
import types
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))


class DummySeries(list):
    def astype(self, typ):
        return DummySeries([typ(x) for x in self])


class DummyDataFrame(dict):
    def set_index(self, col):
        self['_index'] = list(self[col])
        return self

    def sort_index(self):
        return self

    def __getitem__(self, key):
        val = dict.__getitem__(self, key)
        if isinstance(val, list):
            return DummySeries(val)
        return val

    def __setitem__(self, key, val):
        if isinstance(val, DummySeries):
            val = list(val)
        dict.__setitem__(self, key, val)


def dummy_read_csv(fname, sep="\t", comment="#"):
    headers = None
    data = {}
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(comment):
                continue
            parts = line.split(sep)
            if headers is None:
                headers = parts
                for h in headers:
                    data[h] = []
            else:
                for h, v in zip(headers, parts):
                    data[h].append(v)
    return DummyDataFrame(data)


class DummyPandas(types.ModuleType):
    def __init__(self):
        super().__init__("pandas")
        self.read_csv = dummy_read_csv


@pytest.fixture(autouse=True)
def cleanup_modules():
    original = sys.modules.copy()
    yield
    for key in list(sys.modules.keys()):
        if key not in original:
            sys.modules.pop(key)


def test_load_astar_data(tmp_path, monkeypatch):
    def patch_deps():
        dummy_np = types.ModuleType("numpy")
        dummy_np.ndarray = object
        monkeypatch.setitem(sys.modules, "numpy", dummy_np)
        monkeypatch.setitem(sys.modules, "matplotlib", types.ModuleType("matplotlib"))
        mpl_anim = types.ModuleType("animation")
        mpl_anim.FuncAnimation = object
        mpl_anim.FFMpegWriter = object
        mpl_pyplot = types.ModuleType("pyplot")
        monkeypatch.setitem(sys.modules, "matplotlib.animation", mpl_anim)
        monkeypatch.setitem(sys.modules, "matplotlib.pyplot", mpl_pyplot)
        monkeypatch.setitem(sys.modules, "seaborn", types.ModuleType("seaborn"))
        monkeypatch.setitem(sys.modules, "networkx", types.ModuleType("networkx"))
        monkeypatch.setitem(sys.modules, "pandas", DummyPandas())

    patch_deps()

    tsv = tmp_path / "astar.tsv"
    tsv.write_text(
        "# seq_name: test - seq: AC - max_rank: 2\n"
        "node_id\toffset\tmatrix\tscore\n"
        "0\t0\tmatch\t1\n"
        "1\t1\tmatch\t2147482624\n"
    )

    from contrib.poasta_tools import poasta_plot as pp

    df, meta = pp.load_astar_data(tsv)

    assert meta == {"seq_name": "test", "seq": "AC", "max_rank": 2}

    mapping = {0: 0, 1: 1}
    df["rank"] = [mapping[int(n)] for n in df["node_id"]]
    df["offset"] = [int(o) for o in df["offset"]]

    assert df["rank"] == [0, 1]
    assert df["offset"] == [0, 1]
    assert df["score"] == [1, 2147482624]

