"""Extrator de registros da BVS/LILACS via API iAHx (Solr)."""

import logging
import time
from typing import Any, Dict, List, Optional
from urllib.parse import quote

import requests

from src.config import Config
from src.extractors.base import BaseExtractor
from src.models import BibRecord

logger = logging.getLogger(__name__)

BATCH_SIZE = 100
MAX_RETRIES = 5
REQUEST_DELAY = 1.0  # segundos entre requisições

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

# URL base da API de pesquisa BVS
BVS_API_URL = "https://pesquisa.bvsalud.org/portal/"


class BvsExtractor(BaseExtractor):
    """Extrator de registros da BVS/LILACS via API iAHx."""

    def __init__(self, config: Config):
        super().__init__(config)
        self._session = requests.Session()
        self._session.headers.update({
            "User-Agent": "BibExtractor/1.0 (academic research)",
            "Accept": "application/json",
        })

    def build_query(self) -> str:
        block1 = "(" + " OR ".join(DISASTERS_DECS + DISASTERS_FREE) + ")"
        block2 = "(" + " OR ".join(GESTATIONAL_DECS + GESTATIONAL_FREE) + ")"
        block3 = "(" + " OR ".join(NEONATAL_DECS + NEONATAL_FREE) + ")"

        outcomes = f"({block2} OR {block3})"
        query = f"{block1} AND {outcomes}"

        self.search_query = query
        logger.info("Query BVS construída (%d caracteres)", len(query))
        return query

    def search(self) -> int:
        if not self.search_query:
            self.build_query()

        logger.info("Executando busca na BVS/LILACS...")
        params = {
            "q": self.search_query,
            "lang": "en",
            "filter[db][]": "LILACS",
            "output": "solr",
            "wt": "json",
            "rows": 0,
        }

        data = self._request_with_retry(params)
        if data:
            self.total_results = data.get("response", {}).get("numFound", 0)
        else:
            # Fallback: tentar formato alternativo
            self.total_results = self._search_fallback()

        logger.info("BVS/LILACS: %d registros encontrados", self.total_results)
        return self.total_results

    def _search_fallback(self) -> int:
        """Tenta formato alternativo de API se o padrão falhar."""
        params = {
            "q": self.search_query,
            "lang": "en",
            "filter[db][]": "LILACS",
            "format": "json",
            "count": 0,
        }
        data = self._request_with_retry(params)
        if data:
            return data.get("diapiResponse", {}).get("numFound", 0)
        return 0

    def fetch_records(self) -> List[BibRecord]:
        if not self.search_query:
            self.build_query()
        if self.total_results == 0:
            self.search()

        self.records = []
        total = min(self.total_results, self.config.search.max_results_per_db)
        logger.info("Recuperando %d registros da BVS em lotes de %d...", total, BATCH_SIZE)

        for start in range(0, total, BATCH_SIZE):
            batch = self._fetch_batch(start, min(BATCH_SIZE, total - start))
            self.records.extend(batch)
            logger.info(
                "BVS: lote %d-%d recuperado (%d registros acumulados)",
                start,
                start + len(batch),
                len(self.records),
            )
            time.sleep(REQUEST_DELAY)

        logger.info("BVS/LILACS: %d registros recuperados no total", len(self.records))
        return self.records

    def _fetch_batch(self, start: int, rows: int) -> List[BibRecord]:
        params = {
            "q": self.search_query,
            "lang": "en",
            "filter[db][]": "LILACS",
            "output": "solr",
            "wt": "json",
            "start": start,
            "rows": rows,
        }

        data = self._request_with_retry(params)
        if not data:
            return []

        docs = data.get("response", {}).get("docs", [])
        return [self._parse_doc(doc) for doc in docs]

    def _request_with_retry(self, params: dict) -> Optional[Dict[str, Any]]:
        for attempt in range(MAX_RETRIES):
            try:
                resp = self._session.get(BVS_API_URL, params=params, timeout=60)
                resp.raise_for_status()
                return resp.json()
            except requests.exceptions.JSONDecodeError:
                logger.warning(
                    "BVS: resposta não é JSON (tentativa %d/%d). Content-Type: %s",
                    attempt + 1,
                    MAX_RETRIES,
                    resp.headers.get("Content-Type", "unknown"),
                )
            except requests.RequestException as e:
                wait = 2 ** attempt
                logger.warning(
                    "BVS request falhou (tentativa %d/%d): %s. Aguardando %ds...",
                    attempt + 1,
                    MAX_RETRIES,
                    e,
                    wait,
                )
                time.sleep(wait)

        logger.error("BVS: falha após %d tentativas", MAX_RETRIES)
        return None

    def _parse_doc(self, doc: dict) -> BibRecord:
        # ID
        source_id = doc.get("id", "")

        # Título (priorizar português, depois inglês, depois espanhol)
        title = self._get_multilang_field(doc, "ti")

        # Autores
        authors = doc.get("au", [])
        if isinstance(authors, str):
            authors = [authors]

        # Periódico
        journal = doc.get("ta", "") or ""
        if isinstance(journal, list):
            journal = journal[0] if journal else ""

        # Ano
        year = None
        date_str = doc.get("da", "") or ""
        if isinstance(date_str, list):
            date_str = date_str[0] if date_str else ""
        if date_str and len(date_str) >= 4:
            try:
                year = int(date_str[:4])
            except ValueError:
                pass
        # Fallback: campo entry_date ou year
        if year is None:
            yr = doc.get("year_cluster", "") or doc.get("publication_year", "")
            if isinstance(yr, list):
                yr = yr[0] if yr else ""
            if yr:
                try:
                    year = int(str(yr)[:4])
                except ValueError:
                    pass

        # Abstract
        abstract = self._get_multilang_field(doc, "ab")

        # MeSH/DeCS terms
        mesh_terms = doc.get("mh", [])
        if isinstance(mesh_terms, str):
            mesh_terms = [mesh_terms]

        # DOI
        doi = None
        doi_field = doc.get("doi", None)
        if doi_field:
            if isinstance(doi_field, list):
                doi = doi_field[0] if doi_field else None
            else:
                doi = doi_field

        # Idioma
        language = None
        lang = doc.get("la", [])
        if isinstance(lang, list) and lang:
            language = lang[0]
        elif isinstance(lang, str):
            language = lang

        # Fonte completa (volume, issue, pages)
        volume = None
        issue = None
        pages = None
        full_source = doc.get("fo", "") or ""
        if isinstance(full_source, list):
            full_source = full_source[0] if full_source else ""

        # ISSN
        issn = None
        issn_field = doc.get("issn", None)
        if issn_field:
            if isinstance(issn_field, list):
                issn = issn_field[0] if issn_field else None
            else:
                issn = issn_field

        # URL
        url = None
        url_field = doc.get("ur", None) or doc.get("fulltext", None)
        if url_field:
            if isinstance(url_field, list):
                url = url_field[0] if url_field else None
            else:
                url = url_field

        return BibRecord(
            source_db="bvs",
            source_id=str(source_id),
            doi=doi,
            title=title,
            authors=authors,
            journal=journal,
            year=year,
            volume=volume,
            issue=issue,
            pages=pages,
            issn=issn,
            abstract=abstract,
            keywords=[],
            mesh_terms=mesh_terms,
            language=language,
            publication_type=None,
            url=url,
        )

    def _get_multilang_field(self, doc: dict, field_prefix: str) -> str:
        """Busca campo multilíngue priorizando pt > en > es > campo base."""
        for suffix in ["_pt", "_en", "_es", ""]:
            key = f"{field_prefix}{suffix}"
            value = doc.get(key)
            if value:
                if isinstance(value, list):
                    return value[0] if value else ""
                return str(value)
        return ""
