"""Exportador de registros bibliográficos para CSV."""

import logging
from pathlib import Path
from typing import List

import pandas as pd

from src.models import BibRecord

logger = logging.getLogger(__name__)


def export_csv(records: List[BibRecord], output_path: str) -> str:
    """
    Exporta registros para CSV com encoding UTF-8-BOM (compatível com Excel Windows).

    Retorna o caminho do arquivo gerado.
    """
    rows = []
    for rec in records:
        rows.append({
            "source_db": rec.source_db,
            "source_id": rec.source_id,
            "doi": rec.doi or "",
            "title": rec.title,
            "authors": "; ".join(rec.authors),
            "journal": rec.journal,
            "year": rec.year or "",
            "volume": rec.volume or "",
            "issue": rec.issue or "",
            "pages": rec.pages or "",
            "issn": rec.issn or "",
            "abstract": rec.abstract,
            "keywords": "; ".join(rec.keywords),
            "mesh_terms": "; ".join(rec.mesh_terms),
            "language": rec.language or "",
            "publication_type": rec.publication_type or "",
            "url": rec.url or "",
        })

    df = pd.DataFrame(rows)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False, encoding="utf-8-sig", sep=";")

    logger.info("CSV exportado: %s (%d registros)", output_path, len(records))
    return output_path
