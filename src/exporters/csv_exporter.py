"""Exportador de registros bibliográficos para CSV e Excel."""

import csv
import logging
import re
from collections import Counter
from pathlib import Path
from typing import List

import pandas as pd

from src.models import BibRecord

logger = logging.getLogger(__name__)

# Limite de caracteres por célula no Excel
_EXCEL_CELL_LIMIT = 32_767

# Caracteres ilegais em XML 1.0 (usados internamente pelo xlsx)
_ILLEGAL_XML_RE = re.compile(
    r"[\x00-\x08\x0B\x0C\x0E-\x1F\x7F-\x84\x86-\x9F]"
)


def _sanitize_for_excel(value):
    """Remove caracteres ilegais e trunca para o limite do Excel."""
    if not isinstance(value, str):
        return value
    clean = _ILLEGAL_XML_RE.sub("", value)
    if len(clean) > _EXCEL_CELL_LIMIT:
        clean = clean[: _EXCEL_CELL_LIMIT - 3] + "..."
    return clean


def _build_dataframe(records: List[BibRecord]) -> pd.DataFrame:
    """Converte lista de BibRecord em DataFrame padronizado."""
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
            "matched_disasters_mesh": "; ".join(rec.matched_disasters_mesh),
            "matched_gestational_mesh": "; ".join(rec.matched_gestational_mesh),
            "matched_neonatal_mesh": "; ".join(rec.matched_neonatal_mesh),
            "language": rec.language or "",
            "publication_type": rec.publication_type or "",
            "url": rec.url or "",
        })
    return pd.DataFrame(rows)


def _build_mesh_stats(records: List[BibRecord]) -> pd.DataFrame:
    """Gera estatísticas de matching MeSH por categoria."""
    total = len(records)
    if total == 0:
        return pd.DataFrame()

    has_dis = [r for r in records if r.matched_disasters_mesh]
    has_ges = [r for r in records if r.matched_gestational_mesh]
    has_neo = [r for r in records if r.matched_neonatal_mesh]
    has_any = [r for r in records if r.matched_disasters_mesh or r.matched_gestational_mesh or r.matched_neonatal_mesh]
    has_none = total - len(has_any)

    # Combinações
    only_dis = [r for r in records if r.matched_disasters_mesh and not r.matched_gestational_mesh and not r.matched_neonatal_mesh]
    only_ges = [r for r in records if r.matched_gestational_mesh and not r.matched_disasters_mesh and not r.matched_neonatal_mesh]
    only_neo = [r for r in records if r.matched_neonatal_mesh and not r.matched_disasters_mesh and not r.matched_gestational_mesh]
    dis_ges = [r for r in records if r.matched_disasters_mesh and r.matched_gestational_mesh and not r.matched_neonatal_mesh]
    dis_neo = [r for r in records if r.matched_disasters_mesh and r.matched_neonatal_mesh and not r.matched_gestational_mesh]
    ges_neo = [r for r in records if r.matched_gestational_mesh and r.matched_neonatal_mesh and not r.matched_disasters_mesh]
    all_three = [r for r in records if r.matched_disasters_mesh and r.matched_gestational_mesh and r.matched_neonatal_mesh]

    def pct(n: int) -> str:
        return f"{n / total * 100:.1f}%"

    # Tabela 1: Resumo geral
    summary_rows = [
        {"Métrica": "Total de registros", "Contagem": total, "%": "100.0%"},
        {"Métrica": "Com algum matching MeSH", "Contagem": len(has_any), "%": pct(len(has_any))},
        {"Métrica": "Sem nenhum matching MeSH", "Contagem": has_none, "%": pct(has_none)},
        {"Métrica": "", "Contagem": "", "%": ""},
        {"Métrica": "--- Por categoria ---", "Contagem": "", "%": ""},
        {"Métrica": "Disasters (qualquer)", "Contagem": len(has_dis), "%": pct(len(has_dis))},
        {"Métrica": "Gestational (qualquer)", "Contagem": len(has_ges), "%": pct(len(has_ges))},
        {"Métrica": "Neonatal (qualquer)", "Contagem": len(has_neo), "%": pct(len(has_neo))},
        {"Métrica": "", "Contagem": "", "%": ""},
        {"Métrica": "--- Combinações exclusivas ---", "Contagem": "", "%": ""},
        {"Métrica": "Somente Disasters", "Contagem": len(only_dis), "%": pct(len(only_dis))},
        {"Métrica": "Somente Gestational", "Contagem": len(only_ges), "%": pct(len(only_ges))},
        {"Métrica": "Somente Neonatal", "Contagem": len(only_neo), "%": pct(len(only_neo))},
        {"Métrica": "Disasters + Gestational", "Contagem": len(dis_ges), "%": pct(len(dis_ges))},
        {"Métrica": "Disasters + Neonatal", "Contagem": len(dis_neo), "%": pct(len(dis_neo))},
        {"Métrica": "Gestational + Neonatal", "Contagem": len(ges_neo), "%": pct(len(ges_neo))},
        {"Métrica": "Todas as 3 categorias", "Contagem": len(all_three), "%": pct(len(all_three))},
    ]

    # Tabela 2: Frequência de cada termo MeSH individual
    dis_counter: Counter = Counter()
    ges_counter: Counter = Counter()
    neo_counter: Counter = Counter()
    for r in records:
        for t in r.matched_disasters_mesh:
            dis_counter[t] += 1
        for t in r.matched_gestational_mesh:
            ges_counter[t] += 1
        for t in r.matched_neonatal_mesh:
            neo_counter[t] += 1

    summary_rows.append({"Métrica": "", "Contagem": "", "%": ""})
    summary_rows.append({"Métrica": "=== Frequência por termo MeSH ===", "Contagem": "", "%": ""})
    summary_rows.append({"Métrica": "", "Contagem": "", "%": ""})

    summary_rows.append({"Métrica": "--- Disasters ---", "Contagem": "", "%": ""})
    for term, count in dis_counter.most_common():
        summary_rows.append({"Métrica": f"  {term}", "Contagem": count, "%": pct(count)})

    summary_rows.append({"Métrica": "", "Contagem": "", "%": ""})
    summary_rows.append({"Métrica": "--- Gestational ---", "Contagem": "", "%": ""})
    for term, count in ges_counter.most_common():
        summary_rows.append({"Métrica": f"  {term}", "Contagem": count, "%": pct(count)})

    summary_rows.append({"Métrica": "", "Contagem": "", "%": ""})
    summary_rows.append({"Métrica": "--- Neonatal ---", "Contagem": "", "%": ""})
    for term, count in neo_counter.most_common():
        summary_rows.append({"Métrica": f"  {term}", "Contagem": count, "%": pct(count)})

    return pd.DataFrame(summary_rows)


def export_csv(records: List[BibRecord], output_path: str) -> str:
    """
    Exporta registros para CSV com encoding UTF-8-BOM (compatível com Excel Windows).

    Retorna o caminho do arquivo gerado.
    """
    df = _build_dataframe(records)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False, encoding="utf-8-sig", sep=";",
              quoting=csv.QUOTE_ALL)

    logger.info("CSV exportado: %s (%d registros)", output_path, len(records))
    return output_path


def export_excel(records: List[BibRecord], output_path: str) -> str:
    """
    Exporta registros para Excel (.xlsx) com duas abas:
      - Resultados: dados bibliográficos completos
      - Estatísticas MeSH: resumo de matching por categoria e termo

    Retorna o caminho do arquivo gerado.
    """
    df_results = _build_dataframe(records)
    df_stats = _build_mesh_stats(records)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    # Sanitizar strings para compatibilidade com Excel/XML
    df_clean = df_results.apply(
        lambda col: col.map(_sanitize_for_excel) if col.dtype == object else col
    )

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        df_clean.to_excel(writer, sheet_name="Resultados", index=False)
        if not df_stats.empty:
            df_stats.to_excel(writer, sheet_name="Estatísticas MeSH", index=False)

    logger.info("Excel exportado: %s (%d registros, 2 abas)", output_path, len(records))
    return output_path
