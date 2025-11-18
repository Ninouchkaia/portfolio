from pathlib import Path
from typing import Dict, List, Tuple
import csv

def list_virus_dirs(rootdir: Path) -> List[Path]:
    """Retourne la liste des sous-dossiers (un par virus)."""
    return [
        d for d in rootdir.iterdir()
        if d.is_dir()
    ]

def load_nodes_edges(
    virus_dir: Path,
    nodes_filename: str = "nodes.csv",
    edges_filename: str = "edges.csv",
    delimiter_nodes: str = ",",
    delimiter_edges: str = ",",
    has_header: bool = False,
) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str], list[tuple[str, str]]]:
    """
    Lit nodes.csv et edges.csv et renvoie :
      - descr_dict_nodes: id -> nom
      - host_nodes: id -> nom (type 1)
      - virus_nodes: id -> nom (type 0)
      - edges: liste de tuples (u, v)
    """
    nodes_path = virus_dir / nodes_filename
    edges_path = virus_dir / edges_filename

    descr_dict_nodes: Dict[str, str] = {}
    host_nodes: Dict[str, str] = {}
    virus_nodes: Dict[str, str] = {}

    with nodes_path.open() as f:
        reader = csv.reader(f, delimiter=delimiter_nodes)
        if has_header:
            next(reader)
        for row in reader:
            if not row:
                continue
            node_id = row[0]
            node_name = row[1]
            node_type = row[2].strip()
            descr_dict_nodes[node_id] = node_name
            if node_type in {"1", "1.0"}:
                host_nodes[node_id] = node_name
            elif node_type in {"0", "0.0"}:
                virus_nodes[node_id] = node_name
            else:
                # on ignore les noeuds inconnus
                continue

    edges: list[tuple[str, str]] = []
    with edges_path.open() as f:
        reader = csv.reader(f, delimiter=delimiter_edges)
        if has_header:
            next(reader)
        for row in reader:
            if len(row) < 2:
                continue
            u, v = row[0], row[1]
            edges.append((u, v))

    return descr_dict_nodes, host_nodes, virus_nodes, edges
