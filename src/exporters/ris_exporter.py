"""Exportador de registros bibliográficos para formato RIS."""

import logging
from pathlib import Path
from typing import List

import rispy

from src.models import BibRecord

logger = logging.getLogger(__name__)


def _record_to_ris_entry(rec: BibRecord) -> dict:
    """Converte BibRecord para dicionário no formato rispy."""
    entry = {
        "type_of_reference": "JOUR",
        "title": rec.title,
    }

    if rec.authors:
        entry["authors"] = list(rec.authors)

    if rec.journal:
        entry["secondary_title"] = rec.journal

    if rec.year:
        entry["year"] = str(rec.year)

    if rec.volume:
        entry["volume"] = rec.volume

    if rec.issue:
        entry["number"] = rec.issue

    if rec.pages:
        parts = rec.pages.split("-", 1)
        entry["start_page"] = parts[0].strip()
        if len(parts) > 1:
            entry["end_page"] = parts[1].strip()

    if rec.abstract:
        entry["abstract"] = rec.abstract

    if rec.doi:
        entry["doi"] = rec.doi

    if rec.keywords:
        entry["keywords"] = list(rec.keywords)

    if rec.mesh_terms:
        entry["notes_abstract"] = "MeSH: " + "; ".join(rec.mesh_terms)

    if rec.language:
        entry["language"] = rec.language

    if rec.issn:
        entry["issn"] = rec.issn

    if rec.url:
        entry["urls"] = [rec.url]

    if rec.source_id:
        entry["accession_number"] = rec.source_id

    if rec.source_db:
        entry["name_of_database"] = rec.source_db

    return entry


def export_ris(records: List[BibRecord], output_path: str) -> str:
    """
    Exporta registros para arquivo RIS (importável em Zotero, Mendeley, EndNote).

    Retorna o caminho do arquivo gerado.
    """
    entries = [_record_to_ris_entry(rec) for rec in records]

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        rispy.dump(entries, f)

    logger.info("RIS exportado: %s (%d registros)", output_path, len(records))
    return output_path
