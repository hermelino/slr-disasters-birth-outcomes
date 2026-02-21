"""Extrator de registros da BVS/LILACS — modo manual com importação de RIS/CSV.

O portal BVS (pesquisa.bvsalud.org) utiliza proteção CDN (Bunny Shield) que
bloqueia acesso programático. Este módulo gera a query de busca formatada para
copiar/colar no portal e importa os arquivos exportados manualmente.

Fluxo:
  1. build_query()  → gera a query com descritores DeCS e texto livre
  2. O usuário executa a busca no navegador e exporta RIS/CSV
  3. import_from_file() → lê o arquivo exportado e converte para BibRecord
"""

import csv
import logging
import re
from pathlib import Path
from typing import List, Optional
from urllib.parse import quote_plus

import rispy

from src.config import Config
from src.extractors.base import BaseExtractor
from src.models import BibRecord

logger = logging.getLogger(__name__)

# Blocos de busca com descritores DeCS (mh:) e texto livre (tw:)
DISASTERS_DECS = [
    'mh:"Desastres"',
    'mh:"Desastres Naturais"',
    'mh:"Inundações"',
    'mh:"Terremotos"',
    'mh:"Tempestades Ciclônicas"',
    'mh:"Tsunamis"',
    'mh:"Erupções Vulcânicas"',
    'mh:"Incêndios Florestais"',
    'mh:"Secas"',
    'mh:"Liberação de Riscos Químicos"',
    'mh:"Liberação de Riscos Radioativos"',
    'mh:"Incidentes com Feridos em Massa"',
]

DISASTERS_FREE = [
    "tw:hurricane",
    "tw:typhoon",
    "tw:tornado",
    "tw:flood",
    "tw:earthquake",
    "tw:wildfire",
    "tw:tsunami",
    "tw:drought",
    "tw:cyclone",
    "tw:landslide",
    'tw:"oil spill"',
    'tw:"chemical spill"',
    'tw:"desastre ambiental"',
    'tw:"rompimento de barragem"',
    'tw:"desastre minerário"',
    'tw:"dam collapse"',
    'tw:"mining disaster"',
]

GESTATIONAL_DECS = [
    'mh:"Resultado da Gravidez"',
    'mh:"Complicações na Gravidez"',
    'mh:"Nascimento Prematuro"',
    'mh:"Recém-Nascido de Baixo Peso"',
    'mh:"Aborto Espontâneo"',
    'mh:"Natimorto"',
    'mh:"Pré-Eclâmpsia"',
    'mh:"Eclampsia"',
    'mh:"Morte Fetal"',
    'mh:"Peso ao Nascer"',
    'mh:"Mortalidade Materna"',
]

GESTATIONAL_FREE = [
    'tw:"nascimento prematuro"',
    'tw:"baixo peso ao nascer"',
    'tw:"preterm birth"',
    'tw:"low birth weight"',
    'tw:"pregnancy outcome"',
    "tw:miscarriage",
    "tw:stillbirth",
    "tw:preeclampsia",
]

NEONATAL_DECS = [
    'mh:"Mortalidade Infantil"',
    'mh:"Mortalidade Neonatal"',
    'mh:"Mortalidade Perinatal"',
    'mh:"Índice de Apgar"',
    'mh:"Terapia Intensiva Neonatal"',
    'mh:"Anormalidades Congênitas"',
]

NEONATAL_FREE = [
    'tw:"mortalidade infantil"',
    'tw:"mortalidade neonatal"',
    'tw:"infant mortality"',
    'tw:"neonatal mortality"',
    'tw:"neonatal death"',
    "tw:NICU",
    "tw:APGAR",
    'tw:"birth defect"',
    'tw:"congenital anomaly"',
]

BVS_PORTAL_URL = "https://pesquisa.bvsalud.org/portal/"


class BvsExtractor(BaseExtractor):
    """Extrator BVS/LILACS com exportação manual e importação de arquivos."""

    def __init__(self, config: Config):
        super().__init__(config)

    def build_query(self) -> str:
        block1 = "(" + " OR ".join(DISASTERS_DECS + DISASTERS_FREE) + ")"
        block2 = "(" + " OR ".join(GESTATIONAL_DECS + GESTATIONAL_FREE) + ")"
        block3 = "(" + " OR ".join(NEONATAL_DECS + NEONATAL_FREE) + ")"

        outcomes = f"({block2} OR {block3})"
        query = f"{block1} AND {outcomes}"

        self.search_query = query
        logger.info("Query BVS construída (%d caracteres)", len(query))
        return query

    def get_search_url(self) -> str:
        """Retorna a URL completa para busca manual no portal BVS."""
        if not self.search_query:
            self.build_query()
        params = f"?q={quote_plus(self.search_query)}&lang=en&filter[db][]=LILACS"
        return BVS_PORTAL_URL + params

    def print_manual_instructions(self):
        """Imprime instruções para exportação manual do portal BVS."""
        url = self.get_search_url()
        print()
        print("=" * 70)
        print("  BVS/LILACS — EXPORTAÇÃO MANUAL NECESSÁRIA")
        print("=" * 70)
        print()
        print("  O portal BVS bloqueia acesso programático (CDN Bunny Shield).")
        print("  Siga os passos abaixo para exportar manualmente:")
        print()
        print("  1. Abra o link abaixo no navegador:")
        print(f"     {url[:120]}...")
        print()
        print("  2. Verifique o total de resultados na página")
        print()
        print("  3. Selecione todos os registros e exporte como RIS:")
        print("     - Clique em 'Exportar' (ícone de download)")
        print("     - Formato: 'RIS (para Reference Manager, Zotero, etc.)'")
        print("     - Salve o arquivo .ris")
        print()
        print("  4. Execute novamente com --import-bvs:")
        print("     python main.py --import-bvs caminho/arquivo.ris")
        print()
        print("  DICA: A query completa foi salva no diretório de saída como")
        print("  'bvs_query.txt' para referência e reprodutibilidade.")
        print("=" * 70)
        print()

    def search(self) -> int:
        """Não executa busca automatizada. Imprime instruções manuais."""
        if not self.search_query:
            self.build_query()

        logger.info(
            "BVS/LILACS requer exportação manual. "
            "Use --import-bvs <arquivo.ris> para importar registros."
        )
        self.print_manual_instructions()
        self.total_results = 0
        return 0

    def fetch_records(self) -> List[BibRecord]:
        """Não recupera registros automaticamente."""
        self.records = []
        return self.records

    def save_query_file(self, output_dir: str):
        """Salva a query e URL de busca em arquivo para referência."""
        if not self.search_query:
            self.build_query()
        query_file = Path(output_dir) / "bvs_query.txt"
        url = self.get_search_url()
        content = (
            "BVS/LILACS — Query de Busca\n"
            "=" * 50 + "\n\n"
            "Query:\n"
            f"{self.search_query}\n\n"
            "URL de busca (abrir no navegador):\n"
            f"{url}\n\n"
            "Filtro: db=LILACS\n"
            "Instruções: exporte em formato RIS e use --import-bvs\n"
        )
        query_file.write_text(content, encoding="utf-8")
        logger.info("Query BVS salva em: %s", query_file)

    @staticmethod
    def import_from_file(file_path: str) -> List[BibRecord]:
        """Importa registros de um arquivo RIS ou CSV exportado da BVS.

        Retorna lista de BibRecord normalizados.
        """
        path = Path(file_path)
        if not path.exists():
            raise FileNotFoundError(f"Arquivo não encontrado: {file_path}")

        suffix = path.suffix.lower()
        if suffix == ".ris":
            return BvsExtractor._import_ris(path)
        elif suffix == ".csv":
            return BvsExtractor._import_csv(path)
        else:
            raise ValueError(
                f"Formato não suportado: {suffix}. Use .ris ou .csv exportado da BVS."
            )

    @staticmethod
    def _import_ris(path: Path) -> List[BibRecord]:
        """Importa registros de arquivo RIS exportado da BVS."""
        logger.info("Importando RIS da BVS: %s", path)

        with open(path, "r", encoding="utf-8") as f:
            entries = rispy.load(f)

        records = []
        for entry in entries:
            # Extrair DOI
            doi = entry.get("doi")
            if not doi:
                # Procurar DOI em URLs
                urls = entry.get("urls", [])
                for u in urls:
                    doi_match = re.search(r'10\.\d{4,}/[^\s]+', u)
                    if doi_match:
                        doi = doi_match.group(0)
                        break

            # Ano
            year = None
            year_str = entry.get("year", "")
            if year_str:
                try:
                    year = int(str(year_str)[:4])
                except ValueError:
                    pass

            # Páginas
            pages = None
            sp = entry.get("start_page", "")
            ep = entry.get("end_page", "")
            if sp:
                pages = f"{sp}-{ep}" if ep else sp

            # URL
            url = None
            urls = entry.get("urls", [])
            if urls:
                url = urls[0]

            rec = BibRecord(
                source_db="bvs",
                source_id=entry.get("accession_number", ""),
                doi=doi,
                title=entry.get("title", entry.get("primary_title", "")),
                authors=entry.get("authors", entry.get("first_authors", [])),
                journal=entry.get("secondary_title", entry.get("journal_name", "")),
                year=year,
                volume=entry.get("volume", None),
                issue=entry.get("number", None),
                pages=pages,
                issn=entry.get("issn", None),
                abstract=entry.get("abstract", ""),
                keywords=entry.get("keywords", []),
                mesh_terms=[],
                language=entry.get("language", None),
                publication_type=entry.get("type_of_reference", None),
                url=url,
            )
            records.append(rec)

        logger.info("BVS RIS: %d registros importados de %s", len(records), path.name)
        return records

    @staticmethod
    def _import_csv(path: Path) -> List[BibRecord]:
        """Importa registros de arquivo CSV exportado da BVS."""
        logger.info("Importando CSV da BVS: %s", path)

        records = []
        # Tentar diferentes encodings comuns em exports BVS
        for encoding in ["utf-8-sig", "utf-8", "latin-1", "cp1252"]:
            try:
                with open(path, "r", encoding=encoding) as f:
                    reader = csv.DictReader(f)
                    rows = list(reader)
                if rows:
                    break
            except (UnicodeDecodeError, UnicodeError):
                continue
        else:
            raise ValueError(f"Não foi possível ler {path} com nenhum encoding suportado")

        for row in rows:
            # Mapear campos (BVS CSV usa nomes variáveis)
            title = (
                row.get("Title", "")
                or row.get("Título", "")
                or row.get("ti", "")
                or ""
            )
            authors_str = (
                row.get("Authors", "")
                or row.get("Autores", "")
                or row.get("au", "")
                or ""
            )
            authors = [a.strip() for a in authors_str.split(";") if a.strip()]

            year = None
            year_str = (
                row.get("Year", "")
                or row.get("Ano", "")
                or row.get("da", "")
                or ""
            )
            if year_str:
                try:
                    year = int(str(year_str)[:4])
                except ValueError:
                    pass

            doi = row.get("DOI", "") or row.get("doi", "") or None

            rec = BibRecord(
                source_db="bvs",
                source_id=row.get("ID", "") or row.get("id", ""),
                doi=doi if doi else None,
                title=title,
                authors=authors,
                journal=row.get("Journal", "") or row.get("Revista", "") or row.get("ta", "") or "",
                year=year,
                volume=row.get("Volume", "") or row.get("volume", "") or None,
                issue=row.get("Issue", "") or row.get("Número", "") or None,
                pages=row.get("Pages", "") or row.get("Páginas", "") or None,
                issn=row.get("ISSN", "") or row.get("issn", "") or None,
                abstract=row.get("Abstract", "") or row.get("Resumo", "") or row.get("ab", "") or "",
                keywords=[],
                mesh_terms=[],
                language=row.get("Language", "") or row.get("Idioma", "") or None,
                publication_type=None,
                url=row.get("URL", "") or row.get("url", "") or None,
            )
            records.append(rec)

        logger.info("BVS CSV: %d registros importados de %s", len(records), path.name)
        return records
