"""Deduplicação de registros bibliográficos por DOI e similaridade de título."""

import logging
import re
from typing import Dict, List, Tuple

from rapidfuzz import fuzz
from unidecode import unidecode

from src.models import BibRecord

logger = logging.getLogger(__name__)

# Prioridade de retenção: base com metadados mais completos
SOURCE_PRIORITY = {"pubmed": 0, "scopus": 1, "bvs": 2}

# Artigos/preposições a remover na normalização de título
STOPWORDS = {
    "the", "a", "an", "of", "in", "on", "at", "to", "for", "and", "or",
    "with", "by", "from", "as", "is", "was", "are", "were", "been",
    "o", "a", "os", "as", "um", "uma", "uns", "umas", "de", "do", "da",
    "dos", "das", "em", "no", "na", "nos", "nas", "por", "para", "com",
    "el", "la", "los", "las", "un", "una", "unos", "unas", "del", "en",
    "con", "por", "para", "que",
}


def normalize_doi(doi: str) -> str:
    """Normaliza DOI: lowercase, remove prefixos comuns."""
    doi = doi.strip().lower()
    for prefix in ["https://doi.org/", "http://doi.org/", "doi:", "doi.org/"]:
        if doi.startswith(prefix):
            doi = doi[len(prefix):]
    return doi


def normalize_title(title: str) -> str:
    """Normaliza título para comparação fuzzy."""
    text = unidecode(title).lower()
    text = re.sub(r"[^\w\s]", " ", text)
    words = [w for w in text.split() if w not in STOPWORDS]
    return " ".join(words)


def _choose_keeper(records: List[BibRecord]) -> BibRecord:
    """Escolhe o registro a manter: maior completude, depois prioridade da base."""
    return max(
        records,
        key=lambda r: (
            r.completeness_score(),
            -SOURCE_PRIORITY.get(r.source_db, 99),
        ),
    )


def deduplicate(
    records: List[BibRecord],
    fuzzy_threshold: int = 90,
) -> Tuple[List[BibRecord], List[BibRecord]]:
    """
    Deduplica registros em duas fases:
    1. Match exato por DOI
    2. Match fuzzy por título (+ verificação de ano e autor/periódico)

    Retorna (registros_únicos, registros_duplicados).
    """
    if not records:
        return [], []

    # === Fase 1: DOI exato ===
    doi_groups: Dict[str, List[int]] = {}
    for i, rec in enumerate(records):
        if rec.doi:
            ndoi = normalize_doi(rec.doi)
            if ndoi:
                doi_groups.setdefault(ndoi, []).append(i)

    duplicates_doi = set()
    for ndoi, indices in doi_groups.items():
        if len(indices) > 1:
            group = [records[i] for i in indices]
            keeper = _choose_keeper(group)
            for i in indices:
                if records[i] is not keeper:
                    records[i].is_duplicate = True
                    records[i].duplicate_of = keeper.source_id
                    duplicates_doi.add(i)

    logger.info("Fase 1 (DOI): %d duplicatas identificadas", len(duplicates_doi))

    # === Fase 2: Fuzzy por título ===
    remaining_indices = [
        i for i in range(len(records)) if i not in duplicates_doi
    ]

    # Pré-computar títulos normalizados e metadados para comparação
    normalized = {}
    for i in remaining_indices:
        rec = records[i]
        normalized[i] = {
            "title": normalize_title(rec.title),
            "year": rec.year,
            "first_author": rec.authors[0].lower() if rec.authors else "",
            "journal": rec.journal.lower() if rec.journal else "",
        }

    duplicates_fuzzy = set()
    checked = set()

    for idx, i in enumerate(remaining_indices):
        if i in duplicates_fuzzy:
            continue

        ni = normalized[i]
        if not ni["title"]:
            continue

        cluster = [i]

        for j in remaining_indices[idx + 1:]:
            if j in duplicates_fuzzy:
                continue

            nj = normalized[j]
            if not nj["title"]:
                continue

            # Verificação rápida: mesmo ano (se ambos tiverem ano)
            if ni["year"] and nj["year"] and ni["year"] != nj["year"]:
                continue

            # Similaridade de título
            score = fuzz.token_sort_ratio(ni["title"], nj["title"])
            if score < fuzzy_threshold:
                continue

            # Verificação adicional: mesmo primeiro autor OU mesmo periódico
            author_match = (
                ni["first_author"]
                and nj["first_author"]
                and fuzz.ratio(ni["first_author"], nj["first_author"]) > 80
            )
            journal_match = (
                ni["journal"]
                and nj["journal"]
                and fuzz.ratio(ni["journal"], nj["journal"]) > 80
            )

            if not author_match and not journal_match:
                continue

            cluster.append(j)

        if len(cluster) > 1:
            group = [records[i] for i in cluster]
            keeper = _choose_keeper(group)
            for ci in cluster:
                if records[ci] is not keeper:
                    records[ci].is_duplicate = True
                    records[ci].duplicate_of = keeper.source_id
                    duplicates_fuzzy.add(ci)

    logger.info("Fase 2 (fuzzy): %d duplicatas identificadas", len(duplicates_fuzzy))

    all_duplicate_indices = duplicates_doi | duplicates_fuzzy
    unique = [records[i] for i in range(len(records)) if i not in all_duplicate_indices]
    duplicates = [records[i] for i in all_duplicate_indices]

    logger.info(
        "Deduplicação concluída: %d únicos, %d duplicatas (DOI=%d, fuzzy=%d)",
        len(unique),
        len(duplicates),
        len(duplicates_doi),
        len(duplicates_fuzzy),
    )

    return unique, duplicates
