from collections import deque
from typing import Dict, Iterable, Set
import networkx as nx

def build_graph(
    node_ids: Iterable[str],
    edges: Iterable[tuple[str, str]],
) -> nx.Graph:
    G = nx.Graph()
    G.add_nodes_from(node_ids)
    G.add_edges_from(edges)
    return G

def compute_order_interactors(
    G: nx.Graph,
    virus_nodes: Iterable[str],
    host_nodes: set[str],
    max_order: int,
) -> Dict[int, Set[str]]:
    """
    Calcule les interacteurs d'ordre 1..max_order autour des virus_nodes.
    Retourne un dict: order -> set(node_ids).
    """
    virus_nodes = set(virus_nodes)
    visited = set(virus_nodes)
    frontier = set(virus_nodes)

    results: Dict[int, Set[str]] = {}

    for order in range(1, max_order + 1):
        next_frontier: set[str] = set()
        for u in frontier:
            for v in G.neighbors(u):
                if v not in visited:
                    visited.add(v)
                    next_frontier.add(v)
        # on ne garde que les noeuds humains pour le reporting
        results[order] = {n for n in next_frontier if n in host_nodes}
        frontier = next_frontier

    return results
