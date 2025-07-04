import importlib
import io
import sys
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import pytest


class DummyAGraph:
    def __init__(self, dot_str: str):
        self.dot = dot_str

    def to_string(self) -> str:
        return self.dot


class DummyPandas(types.ModuleType):
    def __init__(self):
        super().__init__('pandas')
        self.last_data = None

    def DataFrame(self, data):
        self.last_data = data
        return data


def patch_poasta_plot_deps(monkeypatch, nx_mod=None):
    dummy_np = types.ModuleType('numpy')
    dummy_np.ndarray = object
    monkeypatch.setitem(sys.modules, 'numpy', dummy_np)
    monkeypatch.setitem(sys.modules, 'matplotlib', types.ModuleType('matplotlib'))
    mpl_anim = types.ModuleType('animation')
    mpl_anim.FuncAnimation = object
    mpl_anim.FFMpegWriter = object
    mpl_pyplot = types.ModuleType('pyplot')
    monkeypatch.setitem(sys.modules, 'matplotlib.animation', mpl_anim)
    monkeypatch.setitem(sys.modules, 'matplotlib.pyplot', mpl_pyplot)
    monkeypatch.setitem(sys.modules, 'seaborn', types.ModuleType('seaborn'))
    monkeypatch.setitem(sys.modules, 'pandas', DummyPandas())
    if nx_mod is None:
        nx_mod = types.ModuleType('networkx')
    monkeypatch.setitem(sys.modules, 'networkx', nx_mod)


def make_dummy_networkx():
    class Graph:
        def __init__(self):
            self._nodes = {}
            self.adj = {}
            self.rev_adj = {}
            self.graph = {'graph': {}}

            class NodeView:
                def __init__(self, parent):
                    self.parent = parent

                def __getitem__(self, n):
                    return self.parent._nodes[n]

                def __call__(self, data=False):
                    return list(self.parent._nodes.items()) if data else list(self.parent._nodes)

            self.nodes = NodeView(self)

        def add_node(self, n, **attrs):
            self._nodes[n] = attrs
            self.adj.setdefault(n, {})
            self.rev_adj.setdefault(n, {})

        def nodes_data(self, data=False):
            return list(self._nodes.items()) if data else list(self._nodes)

        def nodes_iter(self, data=False):
            return self.nodes_data(data)

        def add_edge(self, u, v, **data):
            self.adj.setdefault(u, {})[v] = data
            self.rev_adj.setdefault(v, {})[u] = data

        def edges(self, data=False):
            res = []
            for u, targets in self.adj.items():
                for v, d in targets.items():
                    res.append((u, v, d) if data else (u, v))
            return res

        def remove_edges_from(self, edges):
            for u, v in edges:
                self.adj.get(u, {}).pop(v, None)
                self.rev_adj.get(v, {}).pop(u, None)

        def successors(self, n):
            return list(self.adj.get(n, {}))

        def out_degree(self, n):
            return len(self.adj.get(n, {}))

        def in_degree(self, n):
            return len(self.rev_adj.get(n, {}))

    def dfs_postorder_nodes(g):
        visited = set()
        order = []

        def dfs(n):
            visited.add(n)
            for nb in g.successors(n):
                if nb not in visited:
                    dfs(nb)
            order.append(n)

        for node in g.nodes_iter():
            if node not in visited:
                dfs(node)
        return order

    def topological_sort(g):
        indeg = {n: g.in_degree(n) for n in g.nodes_iter()}
        queue = [n for n, d in indeg.items() if d == 0]
        order = []
        while queue:
            n = queue.pop(0)
            order.append(n)
            for nb in g.successors(n):
                indeg[nb] -= 1
                if indeg[nb] == 0:
                    queue.append(nb)
        return order

    def graphviz_layout(g, prog='dot'):
        return {n: (i * 10, 0) for i, n in enumerate(g.nodes_iter())}

    nx_agraph = types.SimpleNamespace(read_dot=None, graphviz_layout=graphviz_layout)

    module = types.ModuleType('networkx')
    module.Graph = Graph
    module.nx_agraph = nx_agraph
    module.dfs_postorder_nodes = dfs_postorder_nodes
    module.topological_sort = topological_sort
    return module


@pytest.fixture(autouse=True)
def cleanup_modules():
    original = sys.modules.copy()
    yield
    for key in list(sys.modules.keys()):
        if key not in original:
            sys.modules.pop(key)


def test_parse_seq_meta(monkeypatch):
    dummy_pg = types.ModuleType('pygraphviz')
    dummy_pg.AGraph = DummyAGraph
    monkeypatch.setitem(sys.modules, 'pygraphviz', dummy_pg)
    from contrib.poasta_tools import poasta_graphviz_region as pgr

    name, node = pgr.parse_seq_meta(b'ref:Node0')
    assert name == 'ref'
    assert node == 'Node0'


def test_parse_poasta_graphviz(monkeypatch):
    dummy_pg = types.ModuleType('pygraphviz')
    dummy_pg.AGraph = DummyAGraph
    monkeypatch.setitem(sys.modules, 'pygraphviz', dummy_pg)

    import contrib.poasta_tools.poasta_graphviz_region as pgr

    dot = b"# seq:\tref:Node0\n" b"digraph {}"
    g, meta = pgr.parse_poasta_graphviz(io.BytesIO(dot))
    assert isinstance(g, DummyAGraph)
    assert meta == {'ref': (0, 'Node0\n')}


def test_contains():
    from contrib.poasta_tools import poasta_graphviz_region as pgr

    assert pgr.contains([1, 3, 5], 3) is True
    assert pgr.contains([1, 3, 5], 2) is False


def test_load_spoa_matrix(tmp_path, monkeypatch):
    pandas_mod = DummyPandas()
    monkeypatch.setitem(sys.modules, 'pandas', pandas_mod)
    patch_poasta_plot_deps(monkeypatch)

    tsv = tmp_path / 'mat.tsv'
    tsv.write_text('A\tB\nC\t0\t1\t2\nD\t1\t2147482624\t3\n')

    import contrib.poasta_tools.poasta_plot as pp

    xlabels, ylabels, df = pp.load_spoa_matrix(tsv)
    assert xlabels == ['0\nA', '1\nB']
    assert ylabels == ['C (0)', 'D (1)']
    assert df['rank'] == [0, 0, 1, 1]
    assert df['offset'] == [0, 1, 0, 1]
    assert df['score'][0] == 1
    assert df['score'][1] == 2
    assert df['score'][3] == 3 and str(df['score'][2]) == 'nan'


def test_load_graph_and_helpers(tmp_path, monkeypatch):
    nx_mod = make_dummy_networkx()
    patch_poasta_plot_deps(monkeypatch, nx_mod)

    def read_dot(path):
        g = nx_mod.Graph()
        with open(path) as f:
            for line in f:
                line = line.strip().strip(';')
                if not line or line.startswith('digraph') or line == '}':
                    continue
                if '->' in line:
                    u, rest = line.split('->')
                    u = u.strip()
                    if '[' in rest:
                        v, attr = rest.split('[', 1)
                        v = v.strip()
                        attr = attr.strip('[]')
                        attrs = {}
                        for item in attr.split(','):
                            if not item:
                                continue
                            k, v2 = item.split('=')
                            attrs[k.strip()] = v2.strip('"')
                    else:
                        v = rest.strip()
                        attrs = {}
                    g.add_edge(u, v, **attrs)
                else:
                    if '[' in line:
                        n, attr = line.split('[', 1)
                        n = n.strip()
                        attr = attr.strip('[]')
                        attrs = {}
                        for item in attr.split(','):
                            if not item:
                                continue
                            k, v2 = item.split('=')
                            attrs[k.strip()] = v2.strip('"')
                        g.add_node(n, **attrs)
        return g

    nx_mod.nx_agraph.read_dot = read_dot

    dot = """digraph {\n N0 [label="'A' (0)"];\n N1 [label="'B' (1)"];\n N2 [label="'C' (2)"];\n N0 -> N1;\n N1 -> N2 [style=dotted];\n}\n"""
    path = tmp_path / 'g.dot'
    path.write_text(dot)

    import contrib.poasta_tools.poasta_plot as pp

    g, mapping = pp.load_graph(path)
    assert mapping == {0: 0, 1: 1, 2: 2}
    assert list(g.edges()) == [('N0', 'N1')]
    layout = pp.poa_graph_layout(g)
    for n, data in g.nodes_iter(data=True):
        assert layout[n][1] == data['rank'] + 0.5

    lines = pp.poa_matrix_discontinuieties(g)
    assert lines == []

