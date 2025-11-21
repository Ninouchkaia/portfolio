from pathlib import Path
from typing import Dict, List
import logging

from .config import DATA_RAW, INTERACTORS_DIR
from .io_utils import list_virus_dirs, load_nodes_edges
from .network import build_graph, compute_order_interactors

logger = logging.getLogger(__name__)


def compute_all_orders_for_all_viruses(
    max_order: int = 2,
    nodes_filename: str = "nodes.csv",
    edges_filename: str = "edges.csv",
) -> None:
    """
    Pour chaque virus:
      - construit le graphe
      - calcule les interacteurs d'ordres 1..max_order
      - écrit des fichiers texte dans INTERACTORS_DIR/<virus>/ :
           direct_interactors.txt
           only_second_range_interactors.txt
           only_Third_range_interactors.txt
           ...
    """
    virus_dirs = list_virus_dirs(DATA_RAW)
    logger.info("Found %d virus folders", len(virus_dirs))

    for virus_dir in virus_dirs:
        virus_name = virus_dir.name
        logger.info("Processing virus %s", virus_name)

        descr, host_nodes, virus_nodes, edges = load_nodes_edges(
            virus_dir,
            nodes_filename=nodes_filename,
            edges_filename=edges_filename,
            delimiter_nodes=",",
            delimiter_edges=",",
            has_header=False,
        )

        G = build_graph(descr.keys(), edges)
        logger.info(
            "Graph for %s: n_nodes=%d, n_edges=%d",
            virus_name, G.number_of_nodes(), G.number_of_edges()
        )

        order_to_nodes = compute_order_interactors(
            G,
            virus_nodes=virus_nodes.keys(),
            host_nodes=set(host_nodes.keys()),
            max_order=max_order,
        )

        out_dir = INTERACTORS_DIR / virus_name
        out_dir.mkdir(parents=True, exist_ok=True)

        # ordre 1 : direct_interactors
        direct_path = out_dir / "direct_interactors.txt"
        _write_interactors_file(
            direct_path,
            virus_name,
            order_to_nodes.get(1, set()),
            descr,
            header_label="Direct interactors (order 1)",
        )

        # ordres suivants : only_second_range_interactors, etc.
        for order in range(2, max_order + 1):
            filename = {
                2: "only_second_range_interactors.txt",
                3: "only_Third_range_interactors.txt",
                4: "only_Fourth_range_interactors.txt",
                5: "only_Fifth_range_interactors.txt",
            }.get(order, f"only_{order}th_range_interactors.txt")
            out_path = out_dir / filename
            _write_interactors_file(
                out_path,
                virus_name,
                order_to_nodes.get(order, set()),
                descr,
                header_label=f"Only order {order} interactors",
            )


def _write_interactors_file(
    path: Path,
    virus_name: str,
    protein_ids: set[str],
    descr_dict: Dict[str, str],
    header_label: str,
) -> None:
    """
    Format inspiré des scripts legacy :
     - première ligne = virus_name
     - lignes suivantes = "protein_id<TAB>gene_symbol"
    """
    with path.open("w", encoding="utf-8") as f:
        f.write(f"{virus_name}\n")
        f.write(f"# {header_label}\n")
        for pid in sorted(protein_ids):
            gene = descr_dict.get(pid, "NA")
            f.write(f"{pid}\t{gene}\n")
    logger.info("Wrote %d interactors to %s", len(protein_ids), path)
