"""Classe abstrata base para extratores de bases bibliográficas."""

from abc import ABC, abstractmethod
from typing import List

from src.models import BibRecord


class BaseExtractor(ABC):
    """Interface comum para todos os extratores de bases de dados."""

    def __init__(self, config):
        self.config = config
        self.search_query: str = ""
        self.total_results: int = 0
        self.records: List[BibRecord] = []

    @abstractmethod
    def build_query(self) -> str:
        """Constrói a string de busca booleana específica da base."""

    @abstractmethod
    def search(self) -> int:
        """Executa a busca e retorna a contagem total de resultados."""

    @abstractmethod
    def fetch_records(self) -> List[BibRecord]:
        """Recupera todos os registros e normaliza para BibRecord."""

    def get_prisma_stats(self) -> dict:
        """Retorna estatísticas para o diagrama PRISMA."""
        return {
            "database": self.__class__.__name__.replace("Extractor", ""),
            "query": self.search_query,
            "records_identified": self.total_results,
            "records_retrieved": len(self.records),
        }
