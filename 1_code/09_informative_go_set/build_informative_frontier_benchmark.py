#!/usr/bin/env python3
"""Build the GO informative-frontier benchmark from local OBO and human GAF.

This implementation follows the staged specification in which multi-parent
member candidates are retained until the final anchor set is known.  Anchor
confirmation (rule A), unique membership within that anchor set (rule B),
within-module ancestry pruning, and module-size filtering are deterministic.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import heapq
import json
import math
import re
import statistics
import sys
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping, Sequence


ROOT_ID = "GO:0008150"
BP_NAMESPACE = "biological_process"
DO_NOT_ANNOTATE = "gocheck_do_not_annotate"
DEFAULT_K_VALUES = (10, 20, 30, 50, 75, 100, 150, 200)
DEFAULT_MAX_GENE_SET = 500
MIN_MODULE_EXCLUSIVE = 3
SCRIPT_VERSION = "2.0.0"
OBO_URL = "https://purl.obolibrary.org/obo/go/go-basic.obo"
GAF_URL = "https://current.geneontology.org/annotations/gaf/HUMAN-uniprot.gaf.gz"

ENGLISH_STOPWORDS = frozenset(
    """a an and are as at be been being by can could did do does doing for from
    had has have having if in into is it its may might of on or such than that
    the their then there these this those through to under upon via was were
    when where which while with within without""".split()
)
GO_GENERIC_WORDS = frozenset(
    "regulation process positive negative response activity pathway involved cellular".split()
)
TOKEN_FILTERS = {
    "stopword_only": ENGLISH_STOPWORDS,
    "go_generic_removed": ENGLISH_STOPWORDS | GO_GENERIC_WORDS,
}


@dataclass(frozen=True)
class GOTerm:
    go_id: str
    name: str | None
    namespace: str | None
    definition: str | None
    obsolete: bool
    subsets: tuple[str, ...]
    is_a_parents: tuple[str, ...]
    part_of_parents: tuple[str, ...]
    alt_ids: tuple[str, ...]


@dataclass
class Graph:
    terms: dict[str, GOTerm]
    parents: dict[str, tuple[str, ...]]
    children: dict[str, tuple[str, ...]]
    edge_relations: dict[tuple[str, str], tuple[str, ...]]
    topo: tuple[str, ...]
    ancestors: dict[str, set[str]]
    min_depth: dict[str, int | None]


@dataclass
class GAFData:
    direct_genes: dict[str, set[str]]
    header: dict[str, str]
    stats: Counter[str]
    evidence_counts: Counter[str]
    raw_unique_pairs: int
    normalized_unique_pairs: int
    raw_unique_genes: int
    normalized_unique_genes: int


@dataclass
class SolverResult:
    anchors: tuple[str, ...]
    assignments: dict[str, tuple[str, ...]]
    skipped_rule_a: tuple[str, ...]
    ambiguous_final: tuple[str, ...]
    nested_removed_final: tuple[str, ...]
    dropped_modules: tuple[str, ...]
    anchor_memberships_removed: int
    iterations: list[dict[str, object]]
    converged: bool


def _field_values(stanza: Sequence[str], tag: str) -> list[str]:
    prefix = f"{tag}:"
    return [line[len(prefix):].strip() for line in stanza if line.startswith(prefix)]


def _first_field(stanza: Sequence[str], tag: str) -> str | None:
    values = _field_values(stanza, tag)
    return values[0] if values else None


def _quoted_text(value: str | None) -> str | None:
    if value is None:
        return None
    match = re.match(r'^"((?:\\.|[^"\\])*)"', value)
    if not match:
        return None
    text = match.group(1)
    replacements = {
        r"\n": "\n", r"\t": "\t", r'\"': '"', r"\\": "\\",
        r"\[": "[", r"\]": "]", r"\{": "{", r"\}": "}",
        r"\,": ",", r"\:": ":",
    }
    return re.sub(r"\\.", lambda m: replacements.get(m.group(0), m.group(0)[1:]), text).strip()


def _parent_id(value: str) -> str:
    return value.split("!", 1)[0].strip().split()[0]


def parse_obo(path: Path) -> tuple[dict[str, GOTerm], dict[str, str]]:
    terms: dict[str, GOTerm] = {}
    header: dict[str, str] = {}
    stanza_type: str | None = None
    stanza: list[str] = []

    def consume(kind: str | None, lines: Sequence[str]) -> None:
        if kind != "Term":
            return
        go_id = _first_field(lines, "id")
        if not go_id:
            raise ValueError("OBO [Term] stanza lacks id")
        if go_id in terms:
            raise ValueError(f"Duplicate GO ID: {go_id}")
        part_of: list[str] = []
        for value in _field_values(lines, "relationship"):
            fields = value.split("!", 1)[0].strip().split()
            if len(fields) >= 2 and fields[0] == "part_of":
                part_of.append(fields[1])
        terms[go_id] = GOTerm(
            go_id=go_id,
            name=_first_field(lines, "name"),
            namespace=_first_field(lines, "namespace"),
            definition=_quoted_text(_first_field(lines, "def")),
            obsolete=_first_field(lines, "is_obsolete") == "true",
            subsets=tuple(_field_values(lines, "subset")),
            is_a_parents=tuple(_parent_id(v) for v in _field_values(lines, "is_a")),
            part_of_parents=tuple(part_of),
            alt_ids=tuple(_field_values(lines, "alt_id")),
        )

    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith("[") and line.endswith("]"):
                consume(stanza_type, stanza)
                stanza_type = line[1:-1]
                stanza = []
            elif stanza_type is None:
                if ":" in line:
                    key, value = line.split(":", 1)
                    header.setdefault(key.strip(), value.strip())
            else:
                stanza.append(line)
    consume(stanza_type, stanza)
    return terms, header


def build_alt_map(terms: Mapping[str, GOTerm]) -> dict[str, str]:
    result: dict[str, str] = {}
    for primary, term in terms.items():
        for alt in term.alt_ids:
            previous = result.setdefault(alt, primary)
            if previous != primary:
                raise ValueError(f"alt_id {alt} maps to both {previous} and {primary}")
    return result


def build_full_graph(all_terms: Mapping[str, GOTerm]) -> Graph:
    terms = {
        go_id: term for go_id, term in all_terms.items()
        if term.namespace == BP_NAMESPACE and not term.obsolete
    }
    parent_sets = {go_id: set() for go_id in terms}
    child_sets = {go_id: set() for go_id in terms}
    relation_sets: dict[tuple[str, str], set[str]] = defaultdict(set)
    for child, term in terms.items():
        for relation, raw_parents in (
            ("is_a", term.is_a_parents), ("part_of", term.part_of_parents)
        ):
            for parent in raw_parents:
                if parent not in terms:
                    continue
                parent_sets[child].add(parent)
                child_sets[parent].add(child)
                relation_sets[(child, parent)].add(relation)

    indegree = {node: len(parent_sets[node]) for node in terms}
    ready = [node for node, degree in indegree.items() if degree == 0]
    heapq.heapify(ready)
    topo: list[str] = []
    while ready:
        node = heapq.heappop(ready)
        topo.append(node)
        for child in sorted(child_sets[node]):
            indegree[child] -= 1
            if indegree[child] == 0:
                heapq.heappush(ready, child)
    if len(topo) != len(terms):
        examples = sorted(node for node, degree in indegree.items() if degree > 0)[:10]
        raise ValueError(f"Full is_a/part_of graph is cyclic: {examples}")
    if ROOT_ID not in terms:
        raise ValueError(f"Missing BP root {ROOT_ID}")

    ancestors = {node: set() for node in terms}
    for node in topo:
        for parent in parent_sets[node]:
            ancestors[node].add(parent)
            ancestors[node].update(ancestors[parent])

    depths: dict[str, int | None] = {node: None for node in terms}
    depths[ROOT_ID] = 0
    queue = deque([ROOT_ID])
    while queue:
        parent = queue.popleft()
        assert depths[parent] is not None
        for child in child_sets[parent]:
            candidate = int(depths[parent]) + 1
            if depths[child] is None or candidate < int(depths[child]):
                depths[child] = candidate
                queue.append(child)

    return Graph(
        terms=terms,
        parents={node: tuple(sorted(parent_sets[node])) for node in terms},
        children={node: tuple(sorted(child_sets[node])) for node in terms},
        edge_relations={edge: tuple(sorted(v)) for edge, v in relation_sets.items()},
        topo=tuple(topo),
        ancestors=ancestors,
        min_depth=depths,
    )


def _parse_gaf_header(line: str, header: dict[str, str]) -> None:
    value = line[1:].strip()
    if ":" in value:
        key, content = value.split(":", 1)
        header.setdefault(key.strip(), content.strip())


def normalize_accession(accession: str) -> str:
    return re.sub(r"-\d+$", "", accession.strip())


def parse_gaf(
    path: Path,
    graph_ids: set[str],
    alt_map: Mapping[str, str],
    exclude_iea: bool,
) -> GAFData:
    direct: dict[str, set[str]] = defaultdict(set)
    header: dict[str, str] = {}
    stats: Counter[str] = Counter()
    evidence_counts: Counter[str] = Counter()
    raw_pairs: set[tuple[str, str]] = set()
    normalized_pairs: set[tuple[str, str]] = set()
    raw_genes: set[str] = set()
    normalized_genes: set[str] = set()

    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if raw.startswith("!"):
                _parse_gaf_header(raw.rstrip("\r\n"), header)
                stats["header_rows"] += 1
                continue
            stats["raw_data_rows"] += 1
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 17:
                raise ValueError(f"Malformed GAF row {line_number}: {len(fields)} columns")
            db, gene, qualifier, raw_go_id, evidence, aspect = (
                fields[0], fields[1], fields[3], fields[4], fields[6], fields[8]
            )
            if db != "UniProtKB":
                stats["F1_excluded_db_not_uniprotkb"] += 1
                continue
            stats["F1_pass_db_uniprotkb"] += 1
            if aspect != "P":
                stats["F2_excluded_aspect_not_p"] += 1
                continue
            stats["F2_pass_aspect_p"] += 1
            if "NOT" in qualifier.split("|"):
                stats["F3_excluded_qualifier_not"] += 1
                continue
            stats["F3_pass_qualifier"] += 1
            go_id = alt_map.get(raw_go_id, raw_go_id)
            if go_id != raw_go_id:
                stats["F4_alt_id_mapped_rows"] += 1
            if go_id not in graph_ids:
                stats["F4_excluded_go_id_not_in_full_graph"] += 1
                continue
            stats["F4_pass_go_id_in_full_graph"] += 1
            if exclude_iea and evidence == "IEA":
                stats["F5_excluded_iea_sensitivity"] += 1
                continue
            stats["F5_pass_evidence_setting"] += 1
            if not gene.strip():
                stats["F6_excluded_empty_gene_id"] += 1
                continue
            stats["retained_rows_before_pair_deduplication"] += 1
            normalized = normalize_accession(gene)
            if normalized != gene:
                stats["retained_rows_isoform_suffix_removed"] += 1
            raw_pairs.add((go_id, gene))
            normalized_pairs.add((go_id, normalized))
            raw_genes.add(gene)
            normalized_genes.add(normalized)
            direct[go_id].add(normalized)
            evidence_counts[evidence or "(missing)"] += 1

    stats["unique_term_gene_pairs_before_isoform_normalization"] = len(raw_pairs)
    stats["unique_term_gene_pairs_after_isoform_normalization"] = len(normalized_pairs)
    stats["duplicate_term_gene_rows_removed"] = (
        stats["retained_rows_before_pair_deduplication"] - len(normalized_pairs)
    )
    return GAFData(
        direct_genes=dict(direct), header=header, stats=stats,
        evidence_counts=evidence_counts, raw_unique_pairs=len(raw_pairs),
        normalized_unique_pairs=len(normalized_pairs), raw_unique_genes=len(raw_genes),
        normalized_unique_genes=len(normalized_genes),
    )


def propagate_genes(graph: Graph, direct: Mapping[str, set[str]]) -> dict[str, set[str]]:
    genes = {node: set(direct.get(node, set())) for node in graph.terms}
    for child in reversed(graph.topo):
        for parent in graph.parents[child]:
            genes[parent].update(genes[child])
    return genes


def informative_set(graph: Graph, n: Mapping[str, int], k: int) -> set[str]:
    return {
        term for term in graph.terms
        if n[term] >= k and all(n[child] < k for child in graph.children[term])
    }


def member_candidates_and_losses(
    all_terms: Mapping[str, GOTerm],
    graph: Graph,
    n: Mapping[str, int],
    informative: set[str],
    k: int,
    max_gene_set: int,
) -> tuple[set[str], dict[str, int]]:
    counts: dict[str, int] = {"starting_obo_term_stanzas": len(all_terms)}
    current = set(graph.terms)
    counts["M1_excluded_not_current_bp_full_graph"] = len(all_terms) - len(current)
    counts["M1_pass"] = len(current)

    next_set = {
        term for term in current
        if graph.terms[term].name and graph.terms[term].definition
    }
    counts["M2_excluded_missing_name_or_definition"] = len(current) - len(next_set)
    counts["M2_pass"] = len(next_set)
    current = next_set

    next_set = {term for term in current if DO_NOT_ANNOTATE not in graph.terms[term].subsets}
    counts["M3_excluded_do_not_annotate"] = len(current) - len(next_set)
    counts["M3_pass"] = len(next_set)
    current = next_set

    next_set = current - {ROOT_ID}
    counts["M4_excluded_bp_root"] = len(current) - len(next_set)
    counts["M4_pass"] = len(next_set)
    current = next_set

    next_set = {term for term in current if n[term] > 0}
    counts["M5_excluded_zero_propagated_genes"] = len(current) - len(next_set)
    counts["M5_pass"] = len(next_set)
    current = next_set

    next_set = current & informative
    counts["M6_excluded_not_informative"] = len(current) - len(next_set)
    counts["M6_pass"] = len(next_set)
    current = next_set

    next_set = {term for term in current if k <= n[term] <= max_gene_set}
    counts["M7_excluded_below_k"] = sum(n[term] < k for term in current)
    counts["M7_excluded_above_max_gene_set"] = sum(n[term] > max_gene_set for term in current)
    counts["M7_pass_member_candidates"] = len(next_set)
    return next_set, counts


def initial_anchor_candidates(
    graph: Graph, member_candidates: set[str], informative: set[str]
) -> dict[str, tuple[str, ...]]:
    result: dict[str, tuple[str, ...]] = {}
    for anchor in sorted(graph.terms):
        members = tuple(child for child in graph.children[anchor] if child in member_candidates)
        if len(members) > MIN_MODULE_EXCLUSIVE:
            if anchor in informative:
                raise AssertionError(f"P3 failed: anchor {anchor} is informative")
            result[anchor] = members
    return result


def _assignment_for_anchors(
    anchors: set[str], member_candidates: set[str], graph: Graph
) -> tuple[dict[str, set[str]], set[str], int]:
    assignments: dict[str, set[str]] = {anchor: set() for anchor in anchors}
    ambiguous: set[str] = set()
    anchor_memberships_removed = 0
    for member in member_candidates:
        parent_anchors = set(graph.parents[member]) & anchors
        if member in anchors:
            anchor_memberships_removed += len(parent_anchors)
            continue
        if len(parent_anchors) == 1:
            assignments[next(iter(parent_anchors))].add(member)
        elif len(parent_anchors) > 1:
            ambiguous.add(member)
    return assignments, ambiguous, anchor_memberships_removed


def _prune_nested(
    assignments: Mapping[str, set[str]], graph: Graph
) -> tuple[dict[str, set[str]], set[str]]:
    result: dict[str, set[str]] = {}
    removed: set[str] = set()
    for anchor, members in assignments.items():
        descendants = {
            member for member in members
            if any(other != member and other in graph.ancestors[member] for other in members)
        }
        result[anchor] = set(members) - descendants
        removed.update(descendants)
    return result, removed


def solve_modules(
    anchor_candidates: Mapping[str, Sequence[str]],
    member_candidates: set[str],
    graph: Graph,
) -> SolverResult:
    ordered = sorted(anchor_candidates, key=lambda p: (-len(anchor_candidates[p]), p))
    confirmed: list[str] = []
    skipped: list[str] = []
    for anchor in ordered:
        if any(anchor in anchor_candidates[prior] for prior in confirmed):
            skipped.append(anchor)
        else:
            confirmed.append(anchor)

    active = set(confirmed)
    dropped: set[str] = set()
    iteration_rows: list[dict[str, object]] = []
    previous_state: tuple[frozenset[str], frozenset[tuple[str, str]]] | None = None
    final_assignments: dict[str, set[str]] = {}
    final_ambiguous: set[str] = set()
    final_nested: set[str] = set()
    final_anchor_memberships_removed = 0
    converged = False

    for iteration in range(1, 100):
        assignments, ambiguous, anchor_memberships_removed = _assignment_for_anchors(
            active, member_candidates, graph
        )
        assignments, nested_removed = _prune_nested(assignments, graph)
        too_small = {
            anchor for anchor in active
            if len(assignments.get(anchor, set())) <= MIN_MODULE_EXCLUSIVE
        }
        active_after = active - too_small
        dropped.update(too_small)
        assignments_after = {
            anchor: assignments[anchor] for anchor in active_after
        }
        state = (
            frozenset(active_after),
            frozenset((anchor, member) for anchor, members in assignments_after.items() for member in members),
        )
        iteration_rows.append({
            "iteration": iteration,
            "anchors_before": len(active),
            "ambiguous_members": len(ambiguous),
            "nested_members_removed": len(nested_removed),
            "modules_dropped": len(too_small),
            "dropped_anchor_ids": ";".join(sorted(too_small)),
            "anchors_after": len(active_after),
            "assigned_members_after": sum(len(v) for v in assignments_after.values()),
            "state_unchanged": state == previous_state,
        })
        if not too_small and state == previous_state:
            converged = True
            final_assignments = assignments_after
            final_ambiguous = ambiguous
            final_nested = nested_removed
            final_anchor_memberships_removed = anchor_memberships_removed
            active = active_after
            break
        previous_state = state
        active = active_after
        final_assignments = assignments_after
        final_ambiguous = ambiguous
        final_nested = nested_removed
        final_anchor_memberships_removed = anchor_memberships_removed
    if not converged:
        raise AssertionError("Iterative anchor/member solver did not converge in <100 iterations")

    return SolverResult(
        anchors=tuple(sorted(active)),
        assignments={a: tuple(sorted(final_assignments[a])) for a in sorted(active)},
        skipped_rule_a=tuple(sorted(skipped)),
        ambiguous_final=tuple(sorted(final_ambiguous)),
        nested_removed_final=tuple(sorted(final_nested)),
        dropped_modules=tuple(sorted(dropped)),
        anchor_memberships_removed=final_anchor_memberships_removed,
        iterations=iteration_rows,
        converged=converged,
    )


def _tokenize(name: str, stopwords: frozenset[str]) -> frozenset[str]:
    return frozenset(
        token for token in re.findall(r"[a-z0-9]+", name.lower()) if token not in stopwords
    )


def _jaccard(left: frozenset[str], right: frozenset[str]) -> float:
    if not left and not right:
        return 0.0
    return len(left & right) / len(left | right)


def _quantile(values: Sequence[float | int], probability: float) -> float | None:
    if not values:
        return None
    ordered = sorted(float(v) for v in values)
    position = (len(ordered) - 1) * probability
    lower, upper = math.floor(position), math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def distribution(values: Sequence[float | int]) -> dict[str, float | int | None]:
    numbers = [float(v) for v in values]
    return {
        "count": len(numbers), "min": min(numbers) if numbers else None,
        "q25": _quantile(numbers, 0.25), "median": _quantile(numbers, 0.5),
        "mean": statistics.fmean(numbers) if numbers else None,
        "q75": _quantile(numbers, 0.75), "max": max(numbers) if numbers else None,
    }


def _auc_cliff(within: Sequence[float], between: Sequence[float]) -> tuple[float | None, float | None]:
    if not within or not between:
        return None, None
    negatives = sorted(between)
    score = 0.0
    for value in within:
        lower = bisect.bisect_left(negatives, value)
        upper = bisect.bisect_right(negatives, value)
        score += lower + 0.5 * (upper - lower)
    auc = score / (len(within) * len(between))
    return auc, 2 * auc - 1


def lexical_qc(result: SolverResult, graph: Graph) -> dict[str, object]:
    term_module = {
        member: anchor for anchor, members in result.assignments.items() for member in members
    }
    term_ids = sorted(term_module)
    output: dict[str, object] = {}
    for variant, stopwords in TOKEN_FILTERS.items():
        tokens = {term: _tokenize(graph.terms[term].name or "", stopwords) for term in term_ids}
        within: list[float] = []
        between: list[float] = []
        for index, left in enumerate(term_ids):
            for right in term_ids[index + 1:]:
                value = _jaccard(tokens[left], tokens[right])
                (within if term_module[left] == term_module[right] else between).append(value)
        auc, cliff = _auc_cliff(within, between)
        output[variant] = {
            "within_module_jaccard": distribution(within),
            "between_module_jaccard": distribution(between),
            "auc_within_greater_than_between": auc,
            "cliffs_delta_within_vs_between": cliff,
            "stopwords": sorted(stopwords),
        }
    return output


def membership_bias_qc(result: SolverResult, graph: Graph) -> dict[str, object]:
    retained = sorted(member for members in result.assignments.values() for member in members)
    ambiguous = list(result.ambiguous_final)
    universe = retained + ambiguous
    tokens = {
        term: _tokenize(graph.terms[term].name or "", TOKEN_FILTERS["go_generic_removed"])
        for term in universe
    }
    document_frequency = Counter(token for term_tokens in tokens.values() for token in term_tokens)
    rare_cutoff = max(1, math.floor(0.05 * len(universe))) if universe else 1

    def summarize(group: Sequence[str]) -> dict[str, object]:
        token_counts = [len(tokens[term]) for term in group]
        rare_fractions = [
            (sum(document_frequency[token] <= rare_cutoff for token in tokens[term]) / len(tokens[term]))
            if tokens[term] else 0.0
            for term in group
        ]
        return {
            "term_count": len(group),
            "mean_unique_name_tokens": statistics.fmean(token_counts) if token_counts else None,
            "median_unique_name_tokens": statistics.median(token_counts) if token_counts else None,
            "mean_rare_token_fraction": statistics.fmean(rare_fractions) if rare_fractions else None,
        }

    return {
        "token_filter": "ordinary stopwords plus the nine specified GO-generic words",
        "rare_token_definition": f"document frequency <= {rare_cutoff} terms (5% floor, minimum 1)",
        "retained_members": summarize(retained),
        "rule_b_excluded_ambiguous_members": summarize(ambiguous),
    }


def _upward_distances(node: str, parents: Mapping[str, Sequence[str]]) -> dict[str, int]:
    distances = {node: 0}
    queue = deque([node])
    while queue:
        child = queue.popleft()
        for parent in parents[child]:
            candidate = distances[child] + 1
            if parent not in distances or candidate < distances[parent]:
                distances[parent] = candidate
                queue.append(parent)
    return distances


def module_pair_qc(result: SolverResult, graph: Graph) -> tuple[dict[str, object], list[dict[str, object]]]:
    anchors = list(result.anchors)
    upward = {anchor: _upward_distances(anchor, graph.parents) for anchor in anchors}
    rows: list[dict[str, object]] = []
    ancestor_pairs = 0
    shared_grandparent_pairs = 0
    hard_negative_union = 0
    distances: list[int] = []
    for index, left in enumerate(anchors):
        for right in anchors[index + 1:]:
            ancestor_related = left in graph.ancestors[right] or right in graph.ancestors[left]
            shared_parents = sorted(set(graph.parents[left]) & set(graph.parents[right]))
            common = set(upward[left]) & set(upward[right])
            lca = ""
            left_distance: int | str = ""
            right_distance: int | str = ""
            total_distance: int | str = ""
            if common:
                lca = min(
                    common,
                    key=lambda term: (
                        upward[left][term] + upward[right][term],
                        -(graph.min_depth[term] if graph.min_depth[term] is not None else -1),
                        term,
                    ),
                )
                left_distance = upward[left][lca]
                right_distance = upward[right][lca]
                total_distance = int(left_distance) + int(right_distance)
                distances.append(total_distance)
            ancestor_pairs += int(ancestor_related)
            shared_grandparent_pairs += int(bool(shared_parents))
            hard_negative_union += int(ancestor_related or bool(shared_parents))
            rows.append({
                "anchor_a": left, "anchor_b": right,
                "ancestor_descendant_pair": ancestor_related,
                "shared_direct_parent": bool(shared_parents),
                "shared_direct_parent_ids": ";".join(shared_parents),
                "nearest_common_ancestor_id": lca,
                "nearest_common_ancestor_name": graph.terms[lca].name if lca else "",
                "anchor_a_to_lca": left_distance, "anchor_b_to_lca": right_distance,
                "total_lca_distance": total_distance,
            })
    return {
        "all_module_pairs": len(rows),
        "anchor_ancestor_descendant_pairs": ancestor_pairs,
        "anchors_sharing_direct_parent_pairs": shared_grandparent_pairs,
        "hard_negative_pair_union": hard_negative_union,
        "nearest_common_ancestor_total_distance_distribution": distribution(distances),
        "nearest_common_ancestor_total_distance_histogram": dict(sorted(Counter(distances).items())),
    }, rows


def validate(
    result: SolverResult,
    graph: Graph,
    informative: set[str],
    n: Mapping[str, int],
    k: int,
    max_gene_set: int,
) -> dict[str, object]:
    members = [member for values in result.assignments.values() for member in values]
    anchors = set(result.anchors)
    v1 = len(members) == len(set(members))
    v2 = not anchors.intersection(members)
    v3 = all(len(values) > MIN_MODULE_EXCLUSIVE for values in result.assignments.values())
    v4 = all(k <= n[member] <= max_gene_set for member in members)
    v5 = all(member in informative for member in members)
    v6 = all(anchor not in informative for anchor in anchors)
    v7 = all(
        left not in graph.ancestors[right] and right not in graph.ancestors[left]
        for values in result.assignments.values()
        for index, left in enumerate(values)
        for right in values[index + 1:]
    )
    v8 = result.converged
    v9 = all(len(set(graph.parents[member]) & anchors) == 1 for member in members)
    checks = {
        "V1_each_member_exactly_one_module": v1,
        "V2_no_anchor_in_mapa_input": v2,
        "V3_each_module_has_more_than_3_members": v3,
        "V4_member_gene_set_size_in_range": v4,
        "V5_all_members_informative": v5,
        "V6_all_anchors_non_informative": v6,
        "V7_no_within_module_ancestor_descendant_pair": v7,
        "V8_solver_reached_fixed_point": v8,
        "V9_each_member_has_exactly_one_parent_in_anchor_set": v9,
    }
    if not all(checks.values()):
        raise AssertionError(f"Validation failed: {checks}")
    return {"status": "passed", **checks}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_csv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]], delimiter: str = ",") -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


def write_k_dataset(
    output_dir: Path,
    k: int,
    max_gene_set: int,
    exclude_iea: bool,
    result: SolverResult,
    anchor_candidates: Mapping[str, Sequence[str]],
    member_candidates: set[str],
    term_losses: Mapping[str, int],
    graph: Graph,
    all_terms: Mapping[str, GOTerm],
    propagated: Mapping[str, set[str]],
    n: Mapping[str, int],
    informative: set[str],
    gaf: GAFData,
    obo_header: Mapping[str, str],
    obo_path: Path,
    gaf_path: Path,
) -> tuple[dict[str, object], dict[str, object]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    validation = validate(result, graph, informative, n, k, max_gene_set)
    lexical = lexical_qc(result, graph)
    bias = membership_bias_qc(result, graph)
    hard_negative, pair_rows = module_pair_qc(result, graph)
    members = sorted(member for values in result.assignments.values() for member in values)
    member_to_anchor = {
        member: anchor for anchor, values in result.assignments.items() for member in values
    }

    _write_csv(
        output_dir / "mapa_input.csv",
        ["go_id", "name", "definition", "gene_set_size"],
        ({
            "go_id": member, "name": graph.terms[member].name or "",
            "definition": graph.terms[member].definition or "", "gene_set_size": n[member],
        } for member in members),
    )
    _write_csv(
        output_dir / "ground_truth.csv",
        ["go_id", "term_name", "module_id", "anchor_go_id", "anchor_name"],
        ({
            "go_id": member, "term_name": graph.terms[member].name or "",
            "module_id": member_to_anchor[member],
            "anchor_go_id": member_to_anchor[member],
            "anchor_name": graph.terms[member_to_anchor[member]].name or "",
        } for member in members),
    )
    _write_csv(
        output_dir / "module_gene_sets.tsv",
        ["go_id", "gene_set_size", "uniprot_accessions"],
        ({
            "go_id": member, "gene_set_size": n[member],
            "uniprot_accessions": ";".join(sorted(propagated[member])),
        } for member in members), delimiter="\t",
    )
    module_fields = [
        "module_id", "anchor_go_id", "anchor_name", "anchor_definition",
        "module_size", "member_go_ids",
    ]
    module_rows = [
        {
            "module_id": anchor, "anchor_go_id": anchor,
            "anchor_name": graph.terms[anchor].name or "",
            "anchor_definition": graph.terms[anchor].definition or "",
            "module_size": len(result.assignments[anchor]),
            "member_go_ids": ";".join(result.assignments[anchor]),
        }
        for anchor in result.anchors
    ]
    _write_csv(output_dir / "modules.csv", module_fields, module_rows)
    _write_csv(output_dir / "module.csv", module_fields, module_rows)
    _write_csv(
        output_dir / "module_pair_lca.csv",
        ["anchor_a", "anchor_b", "ancestor_descendant_pair", "shared_direct_parent", "shared_direct_parent_ids", "nearest_common_ancestor_id", "nearest_common_ancestor_name", "anchor_a_to_lca", "anchor_b_to_lca", "total_lca_distance"],
        pair_rows,
    )
    _write_csv(
        output_dir / "excluded_ambiguous_members.csv",
        ["go_id", "name", "anchor_parent_ids", "anchor_parent_count", "gene_set_size"],
        ({
            "go_id": member, "name": graph.terms[member].name or "",
            "anchor_parent_ids": ";".join(sorted(set(graph.parents[member]) & set(result.anchors))),
            "anchor_parent_count": len(set(graph.parents[member]) & set(result.anchors)),
            "gene_set_size": n[member],
        } for member in result.ambiguous_final),
    )
    _write_csv(
        output_dir / "solver_iterations.csv",
        ["iteration", "anchors_before", "ambiguous_members", "nested_members_removed", "modules_dropped", "dropped_anchor_ids", "anchors_after", "assigned_members_after", "state_unchanged"],
        result.iterations,
    )
    _write_csv(
        output_dir / "anchor_candidates.csv",
        ["anchor_go_id", "anchor_name", "initial_member_count", "initial_member_go_ids", "rule_a_status", "final_status", "final_member_count"],
        ({
            "anchor_go_id": anchor, "anchor_name": graph.terms[anchor].name or "",
            "initial_member_count": len(anchor_candidates[anchor]),
            "initial_member_go_ids": ";".join(anchor_candidates[anchor]),
            "rule_a_status": "skipped_anchor_is_member" if anchor in result.skipped_rule_a else "confirmed",
            "final_status": "selected" if anchor in result.anchors else ("dropped_too_small" if anchor in result.dropped_modules else "not_confirmed_rule_a"),
            "final_member_count": len(result.assignments.get(anchor, ())),
        } for anchor in sorted(anchor_candidates)),
    )

    module_sizes = [len(result.assignments[a]) for a in result.anchors]
    gene_sizes = [n[member] for member in members]
    depths = [graph.min_depth[member] for member in members if graph.min_depth[member] is not None]
    relation_counts = Counter(rel for values in graph.edge_relations.values() for rel in values)
    qc = {
        "scale": {
            "modules": len(result.anchors), "member_terms": len(members),
            "module_size_histogram": dict(sorted(Counter(module_sizes).items())),
            "module_size_distribution": distribution(module_sizes),
            "member_gene_set_size_distribution": distribution(gene_sizes),
            "member_min_depth_from_GO_0008150_distribution": distribution(depths),
            "member_min_depth_histogram": dict(sorted(Counter(depths).items())),
        },
        "filter_losses": {
            "term_filters_M1_to_M7": dict(term_losses),
            "initial_anchor_candidates": len(anchor_candidates),
            "rule_A_skipped_anchor_count": len(result.skipped_rule_a),
            "rule_A_skipped_anchor_ids": list(result.skipped_rule_a),
            "rule_A_anchor_memberships_removed_final": result.anchor_memberships_removed,
            "rule_B_excluded_ambiguous_membership": len(result.ambiguous_final),
            "within_module_nested_descendants_removed": len(result.nested_removed_final),
            "module_size_recheck_dropped_modules": len(result.dropped_modules),
            "module_size_recheck_dropped_anchor_ids": list(result.dropped_modules),
            "solver_iterations_to_fixed_point": len(result.iterations),
            "gaf_filters": dict(sorted(gaf.stats.items())),
        },
        "hard_negative": hard_negative,
        "lexical_leakage": lexical,
        "rule_B_membership_bias": bias,
    }
    _write_json(output_dir / "qc_report.json", qc)

    manifest = {
        "sources": {
            "go_obo": {
                "url": OBO_URL, "path": str(obo_path), "sha256": _sha256(obo_path),
                "data_version": obo_header.get("data-version"),
                "format_version": obo_header.get("format-version"),
            },
            "human_uniprot_gaf": {
                "url": GAF_URL, "path": str(gaf_path), "sha256": _sha256(gaf_path),
                "gaf_version": gaf.header.get("gaf-version"),
                "date_generated": gaf.header.get("date-generated"),
                "go_version": gaf.header.get("go-version"),
            },
            "builder": {
                "version": SCRIPT_VERSION,
                "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            },
        },
        "rules": {
            "namespace": BP_NAMESPACE,
            "relationship": "is_a and part_of",
            "edge_direction": "child to parent",
            "gene_set_annotation_scope": "true-path propagated over is_a and part_of",
            "informative_threshold_k": k,
            "informative_definition": "n(t)>=k and every direct child in the full graph has n(child)<k",
            "member_definition": "informative direct children of an anchor",
            "module_size": "> 3",
            "anchor_rule_A": "anchor must not be a member of another module",
            "member_rule_B": "at most one direct parent within the confirmed anchor set",
            "within_module_nesting": "remove descendant member and retain ancestor member",
            "member_gene_set_size_inclusive": [k, max_gene_set],
            "evidence_codes": "all except IEA" if exclude_iea else "all (including IEA)",
            "exclude_iea": exclude_iea,
            "gene_unit": "deduplicated UniProtKB DB_Object_ID with terminal -<digits> isoform suffix removed",
            "reference_parent_excluded_from_mapa_input": True,
            "regulation_terms_excluded": False,
            "k_selection": "blind: dataset scale and lexical redundancy only; no method similarity or clustering result inspected",
        },
        "counts": {
            "all_obo_term_stanzas": len(all_terms),
            "full_graph_terms": len(graph.terms),
            "full_graph_edges": len(graph.edge_relations),
            "edge_relation_counts": dict(sorted(relation_counts.items())),
            "full_graph_terms_with_propagated_genes": sum(value > 0 for value in n.values()),
            "bp_root_propagated_gene_count": n[ROOT_ID],
            "gaf": dict(sorted(gaf.stats.items())),
            "retained_evidence_code_rows": dict(sorted(gaf.evidence_counts.items())),
            "unique_accessions_before_isoform_normalization": gaf.raw_unique_genes,
            "unique_accessions_after_isoform_normalization": gaf.normalized_unique_genes,
            "term_filters_M1_to_M7": dict(term_losses),
            "informative_terms_in_full_graph": len(informative),
            "member_candidates": len(member_candidates),
            "initial_anchor_candidates": len(anchor_candidates),
            "final_modules": len(result.anchors),
            "final_member_terms": len(members),
            "excluded_ambiguous_membership": len(result.ambiguous_final),
        },
        "validation": validation,
    }
    _write_json(output_dir / "manifest.json", manifest)
    return manifest, qc


def plot_scan(rows: Sequence[Mapping[str, object]], path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ks = [int(row["k"]) for row in rows]
    modules = [int(row["modules"]) for row in rows]
    members = [int(row["member_terms"]) for row in rows]
    aucs = [float(row["auc_go_generic_removed"]) if row["auc_go_generic_removed"] != "" else math.nan for row in rows]
    fig, left = plt.subplots(figsize=(8.2, 5.0), dpi=180)
    right = left.twinx()
    left.plot(ks, modules, marker="o", linewidth=2, color="#1f77b4", label="Modules")
    left.plot(ks, members, marker="s", linewidth=2, color="#2ca02c", label="Member terms")
    right.plot(ks, aucs, marker="^", linewidth=2, color="#d62728", label="Lexical AUC (filter B)")
    left.set_xlabel("Informative threshold k")
    left.set_ylabel("Dataset size (count)")
    right.set_ylabel("Within-vs-between name Jaccard AUC")
    right.set_ylim(0.45, 1.02)
    left.set_xticks(ks)
    left.grid(axis="y", color="#d9d9d9", linewidth=0.7, alpha=0.7)
    lines = left.get_lines() + right.get_lines()
    left.legend(lines, [line.get_label() for line in lines], loc="upper right", frameon=False)
    left.set_title("GO informative-frontier k scan")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def choose_k(rows: Sequence[Mapping[str, object]]) -> tuple[int, str]:
    adequate = [row for row in rows if row["adequate"] is True]
    if adequate:
        selected = min(
            adequate,
            key=lambda row: (
                float(row["auc_go_generic_removed"]),
                -int(row["member_terms"]),
                int(row["k"]),
            ),
        )
        reason = (
            "Among adequate datasets (>=8 modules and >=40 members), minimize lexical AUC "
            "after removing the nine specified GO-generic words; ties favor more members, "
            "then smaller k. Selection is blind to all method outputs."
        )
    else:
        selected = max(rows, key=lambda row: (int(row["modules"]), int(row["member_terms"]), -int(row["k"])))
        reason = "No dataset met adequacy; maximize modules, then members, then prefer smaller k."
    return int(selected["k"]), reason


def write_root_readme(output_dir: Path, rows: Sequence[Mapping[str, object]], selected_k: int, reason: str) -> None:
    lines = [
        "# GO informative-frontier benchmark scan", "",
        "The benchmark follows the attached staged specification. The full non-obsolete BP graph is built before eligibility filtering and contains only child-to-parent `is_a` and `part_of` edges.", "",
        "GO true-path propagation uses only these two relations. `regulates` relations are excluded because annotations do not propagate over them and because including regulation edges would create structural/gene-set inconsistency and severe name-template redundancy.", "",
        "`hard-negative pairs` is the union of anchor ancestor-descendant pairs and anchor pairs sharing at least one direct parent; the two component counts are retained separately in `k_scan_summary.csv` and each `qc_report.json`. AUC B removes exactly: `regulation, process, positive, negative, response, activity, pathway, involved, cellular`, in addition to ordinary English stopwords.", "",
        "| k | modules | members | hard-negative pairs | AUC A | AUC B | ambiguous removed | adequate |",
        "|---:|---:|---:|---:|---:|---:|---:|:---:|",
    ]
    for row in rows:
        def fmt(value: object) -> str:
            return "NA" if value == "" or value is None else f"{float(value):.3f}"
        lines.append(
            f"| {row['k']} | {row['modules']} | {row['member_terms']} | {row['hard_negative_pairs']} | "
            f"{fmt(row['auc_stopword_only'])} | {fmt(row['auc_go_generic_removed'])} | "
            f"{row['excluded_ambiguous_membership']} | {str(row['adequate']).lower()} |"
        )
    lines.extend([
        "", f"Recommended dataset: **`k_{selected_k}`**.", "", reason, "",
        "The k choice is blind: it uses only dataset size and name-token redundancy. No similarity matrix, clustering assignment, or performance result from MAPA or any comparator was inspected.", "",
        "Each `k_<value>/` directory contains the four required outputs plus module, LCA-pair, ambiguity, anchor-candidate, and fixed-point solver audit tables.", "",
    ])
    (output_dir / "README.md").write_text("\n".join(lines), encoding="utf-8")


def find_project_root(start: Path) -> Path:
    for candidate in (start, *start.parents):
        if (candidate / "mapa").is_dir() and (candidate / "mapa_manuscript").is_dir():
            return candidate
    raise RuntimeError("Could not find project root with mapa/ and mapa_manuscript/")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    root = find_project_root(Path(__file__).resolve().parent)
    data = root / "mapa/demo_data/GO_term_human_annotation"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--obo", type=Path, default=data / "go-basic.obo")
    parser.add_argument("--gaf", type=Path, default=data / "HUMAN-uniprot.gaf.gz")
    parser.add_argument("--output-dir", type=Path, default=root / "mapa_manuscript/3_data_analysis/informative_go_set")
    parser.add_argument("--k-values", type=int, nargs="+", default=list(DEFAULT_K_VALUES))
    parser.add_argument("--max-gene-set", type=int, default=DEFAULT_MAX_GENE_SET)
    parser.add_argument("--exclude-iea", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if any(k <= 0 for k in args.k_values) or len(args.k_values) != len(set(args.k_values)):
        raise ValueError("--k-values must contain unique positive integers")
    if args.max_gene_set <= 0:
        raise ValueError("--max-gene-set must be positive")
    obo_path, gaf_path, output_dir = args.obo.resolve(), args.gaf.resolve(), args.output_dir.resolve()
    for source in (obo_path, gaf_path):
        if not source.is_file():
            raise FileNotFoundError(source)
    output_dir.mkdir(parents=True, exist_ok=True)

    all_terms, obo_header = parse_obo(obo_path)
    alt_map = build_alt_map(all_terms)
    graph = build_full_graph(all_terms)
    gaf = parse_gaf(gaf_path, set(graph.terms), alt_map, args.exclude_iea)
    propagated = propagate_genes(graph, gaf.direct_genes)
    n = {term: len(genes) for term, genes in propagated.items()}

    scan_rows: list[dict[str, object]] = []
    for k in sorted(args.k_values):
        informative = informative_set(graph, n, k)
        member_candidates, term_losses = member_candidates_and_losses(
            all_terms, graph, n, informative, k, args.max_gene_set
        )
        anchor_candidates = initial_anchor_candidates(graph, member_candidates, informative)
        result = solve_modules(anchor_candidates, member_candidates, graph)
        _, qc = write_k_dataset(
            output_dir / f"k_{k}", k, args.max_gene_set, args.exclude_iea,
            result, anchor_candidates, member_candidates, term_losses, graph,
            all_terms, propagated, n, informative, gaf, obo_header, obo_path, gaf_path,
        )
        lexical = qc["lexical_leakage"]
        scale = qc["scale"]
        hard = qc["hard_negative"]
        losses = qc["filter_losses"]
        modules = int(scale["modules"])
        member_terms = int(scale["member_terms"])
        row = {
            "k": k, "modules": modules, "member_terms": member_terms,
            "hard_negative_pairs": hard["hard_negative_pair_union"],
            "anchor_ancestor_descendant_pairs": hard["anchor_ancestor_descendant_pairs"],
            "anchors_sharing_grandparent_pairs": hard["anchors_sharing_direct_parent_pairs"],
            "auc_stopword_only": lexical["stopword_only"]["auc_within_greater_than_between"],
            "auc_go_generic_removed": lexical["go_generic_removed"]["auc_within_greater_than_between"],
            "cliffs_delta_go_generic_removed": lexical["go_generic_removed"]["cliffs_delta_within_vs_between"],
            "median_gene_set_size": scale["member_gene_set_size_distribution"]["median"],
            "excluded_ambiguous_membership": losses["rule_B_excluded_ambiguous_membership"],
            "adequate": modules >= 8 and member_terms >= 40,
        }
        scan_rows.append(row)
        print(
            f"k={k}: modules={modules}, members={member_terms}, "
            f"hard-negative pairs={row['hard_negative_pairs']}, "
            f"AUC_B={row['auc_go_generic_removed']}"
        )

    summary_fields = [
        "k", "modules", "member_terms", "hard_negative_pairs",
        "anchor_ancestor_descendant_pairs", "anchors_sharing_grandparent_pairs",
        "auc_stopword_only", "auc_go_generic_removed", "cliffs_delta_go_generic_removed",
        "median_gene_set_size", "excluded_ambiguous_membership", "adequate",
    ]
    _write_csv(output_dir / "k_scan_summary.csv", summary_fields, scan_rows)
    plot_scan(scan_rows, output_dir / "k_scan.png")
    selected_k, reason = choose_k(scan_rows)
    _write_json(output_dir / "recommended_k.json", {
        "recommended_k": selected_k, "selection_reason": reason,
        "adequacy_definition": "modules >= 8 and member_terms >= 40",
        "primary_lexical_metric": "AUC after ordinary stopwords and the nine specified GO-generic words are removed",
        "blind_selection_statement": "No method similarity matrix or clustering result was inspected.",
        "selected_dataset_directory": f"k_{selected_k}",
    })
    write_root_readme(output_dir, scan_rows, selected_k, reason)
    print(f"Recommended k={selected_k}")
    print(f"Output: {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
