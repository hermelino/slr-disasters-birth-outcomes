"""Modelo de dados normalizado para registros bibliográficos."""

from dataclasses import dataclass, field
from typing import List, Optional


@dataclass
class BibRecord:
    """Registro bibliográfico normalizado, comum a todas as bases de dados."""

    # Identidade
    source_db: str = ""            # "pubmed", "scopus", "bvs"
    source_id: str = ""            # PMID, Scopus EID ou LILACS ID
    doi: Optional[str] = None

    # Dados bibliográficos
    title: str = ""
    authors: List[str] = field(default_factory=list)
    journal: str = ""
    year: Optional[int] = None
    volume: Optional[str] = None
    issue: Optional[str] = None
    pages: Optional[str] = None
    issn: Optional[str] = None

    # Conteúdo
    abstract: str = ""
    keywords: List[str] = field(default_factory=list)
    mesh_terms: List[str] = field(default_factory=list)

    # Metadados
    language: Optional[str] = None
    publication_type: Optional[str] = None
    url: Optional[str] = None

    # Rastreamento de deduplicação
    is_duplicate: bool = False
    duplicate_of: Optional[str] = None  # source_id do registro mantido

    def completeness_score(self) -> int:
        """Retorna uma pontuação de completude dos metadados (mais campos = maior)."""
        score = 0
        if self.doi:
            score += 2
        if self.title:
            score += 1
        if self.authors:
            score += 1
        if self.abstract:
            score += 2
        if self.mesh_terms:
            score += 1
        if self.journal:
            score += 1
        if self.year:
            score += 1
        if self.keywords:
            score += 1
        return score
