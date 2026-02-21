"""Extrator de registros do Scopus via pybliometrics."""

import logging
import math
from typing import List, Optional

from src.config import Config
from src.extractors.base import BaseExtractor
from src.models import BibRecord

logger = logging.getLogger(__name__)

SCOPUS_MAX_RESULTS = 5000

# Blocos de busca em sintaxe Scopus (TITLE-ABS-KEY)
DISASTERS_TERMS = [
    "disaster",
    '"natural disaster"',
    "flood",
    "earthquake",
    "hurricane",
    "typhoon",
    "tornado",
    "wildfire",
    "tsunami",
    "drought",
    "cyclone",
    '"chemical spill"',
    '"oil spill"',
    '"volcanic eruption"',
    "landslide",
    "mudslide",
    '"nuclear accident"',
    '"dam failure"',
    '"dam collapse"',
    '"mining disaster"',
    '"environmental disaster"',
    '"climate disaster"',
    '"mass casualty"',
]

OUTCOMES_TERMS = [
    '"pregnancy outcome"',
    '"pregnancy complication"',
    '"premature birth"',
    '"preterm birth"',
    '"preterm delivery"',
    '"low birth weight"',
    '"small for gestational age"',
    '"spontaneous abortion"',
    "stillbirth",
    "preeclampsia",
    "pre-eclampsia",
    "eclampsia",
    '"fetal death"',
    '"birth weight"',
    '"maternal mortality"',
    "miscarriage",
    '"gestational diabetes"',
    '"gestational hypertension"',
    '"intrauterine growth restriction"',
    '"infant mortality"',
    '"neonatal mortality"',
    '"perinatal mortality"',
    '"Apgar score"',
    '"congenital abnormalities"',
    '"congenital anomalies"',
    '"birth defect"',
    "NICU",
    '"neonatal intensive care"',
    '"neonatal death"',
    '"infant death"',
    '"neonatal outcome"',
]


def _setup_pybliometrics(config: Config):
    """Configura o pybliometrics escrevendo o config file completo e chamando init()."""
    import configparser
    from pathlib import Path

    cfg_dir = Path.home() / ".config"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    cfg_path = cfg_dir / "pybliometrics.cfg"

    # Recriar do zero para evitar entradas duplicadas
    if cfg_path.exists():
        cfg_path.unlink()
    cp = configparser.ConfigParser()

    # Authentication
    if "Authentication" not in cp:
        cp["Authentication"] = {}
    cp["Authentication"]["APIKey"] = config.api.scopus.api_key
    if config.api.scopus.institutional_token:
        cp["Authentication"]["InstToken"] = config.api.scopus.institutional_token

    # Directories
    if "Directories" not in cp:
        cp["Directories"] = {}
    cache_dir = Path.home() / ".cache" / "pybliometrics"
    cache_dir.mkdir(parents=True, exist_ok=True)
    for subdir in ["AbstractDir", "AffiliationDir", "AuthorDir", "CitationDir", "SearchDir", "SerialDir"]:
        d = cache_dir / subdir.replace("Dir", "").lower()
        d.mkdir(parents=True, exist_ok=True)
        cp["Directories"][subdir] = str(d)

    # Requests (exigido pelo pybliometrics v4+)
    if "Requests" not in cp:
        cp["Requests"] = {}
    cp["Requests"].setdefault("Timeout", "20")
    cp["Requests"].setdefault("Retries", "5")

    with open(cfg_path, "w") as f:
        cp.write(f)

    # Carregar o config no pybliometrics passando keys explicitamente
    import pybliometrics

    init_kwargs = {"config_path": str(cfg_path), "keys": [config.api.scopus.api_key]}
    if config.api.scopus.institutional_token:
        init_kwargs["inst_tokens"] = [config.api.scopus.institutional_token]
    pybliometrics.init(**init_kwargs)

    logger.info("pybliometrics configurado e inicializado em %s", cfg_path)


class ScopusExtractor(BaseExtractor):
    """Extrator de registros do Scopus usando pybliometrics."""

    def __init__(self, config: Config):
        super().__init__(config)
        if not config.api.scopus.api_key:
            raise ValueError(
                "Scopus requer API key. Configure api.scopus.api_key no config.yaml "
                "ou defina a variável de ambiente SCOPUS_API_KEY."
            )
        _setup_pybliometrics(config)

    def build_query(self) -> str:
        disasters_block = "TITLE-ABS-KEY(" + " OR ".join(DISASTERS_TERMS) + ")"
        outcomes_block = "TITLE-ABS-KEY(" + " OR ".join(OUTCOMES_TERMS) + ")"

        dr = self.config.search.date_range
        start_year = dr.start.split("/")[0]
        end_year = dr.end.split("/")[0]

        query = f"{disasters_block} AND {outcomes_block} AND PUBYEAR > {int(start_year) - 1} AND PUBYEAR < {int(end_year) + 1}"

        self.search_query = query
        logger.info("Query Scopus construída (%d caracteres)", len(query))
        return query

    def search(self) -> int:
        if not self.search_query:
            self.build_query()

        from pybliometrics.scopus import ScopusSearch

        logger.info("Executando busca no Scopus...")
        try:
            s = ScopusSearch(self.search_query, download=False, subscriber=False)
            self.total_results = s.get_results_size()
        except Exception as e:
            logger.warning("Scopus subscriber=False falhou: %s. Tentando com subscriber=True...", e)
            s = ScopusSearch(self.search_query, download=False)
            self.total_results = s.get_results_size()
        logger.info("Scopus: %d registros encontrados", self.total_results)
        return self.total_results

    def fetch_records(self) -> List[BibRecord]:
        if not self.search_query:
            self.build_query()

        from pybliometrics.scopus import ScopusSearch

        total = min(self.total_results or 999999, self.config.search.max_results_per_db)

        if total > SCOPUS_MAX_RESULTS:
            logger.info(
                "Scopus: %d resultados excedem o limite de %d. Dividindo por faixas de ano...",
                total,
                SCOPUS_MAX_RESULTS,
            )
            self.records = self._fetch_by_year_ranges()
        else:
            logger.info("Scopus: recuperando %d registros...", total)
            s = ScopusSearch(self.search_query, download=True, subscriber=False)
            raw_results = s.results or []
            self.records = [self._parse_result(r) for r in raw_results]

        logger.info("Scopus: %d registros recuperados no total", len(self.records))
        return self.records

    def _fetch_by_year_ranges(self) -> List[BibRecord]:
        """Divide a busca em faixas de 5 anos quando o total excede 5000."""
        from pybliometrics.scopus import ScopusSearch

        dr = self.config.search.date_range
        start_year = int(dr.start.split("/")[0])
        end_year = int(dr.end.split("/")[0])

        all_records = []
        step = 5
        for y_start in range(start_year, end_year + 1, step):
            y_end = min(y_start + step - 1, end_year)

            disasters_block = "TITLE-ABS-KEY(" + " OR ".join(DISASTERS_TERMS) + ")"
            outcomes_block = "TITLE-ABS-KEY(" + " OR ".join(OUTCOMES_TERMS) + ")"
            range_query = (
                f"{disasters_block} AND {outcomes_block} "
                f"AND PUBYEAR > {y_start - 1} AND PUBYEAR < {y_end + 1}"
            )

            logger.info("Scopus: buscando %d-%d...", y_start, y_end)
            try:
                s = ScopusSearch(range_query, download=True, subscriber=False)
                raw = s.results or []
                batch = [self._parse_result(r) for r in raw]
                all_records.extend(batch)
                logger.info("Scopus %d-%d: %d registros", y_start, y_end, len(batch))
            except Exception as e:
                logger.error("Scopus %d-%d falhou: %s", y_start, y_end, e)

        return all_records

    def _parse_result(self, result) -> BibRecord:
        # pybliometrics retorna namedtuples com campos acessíveis por atributo
        doi = getattr(result, "doi", None)
        title = getattr(result, "title", "") or ""

        author_names = getattr(result, "author_names", "") or ""
        authors = [a.strip() for a in author_names.split(";")] if author_names else []

        journal = getattr(result, "publicationName", "") or ""

        year = None
        cover_date = getattr(result, "coverDate", "") or ""
        if cover_date and len(cover_date) >= 4:
            try:
                year = int(cover_date[:4])
            except ValueError:
                pass

        abstract = getattr(result, "description", "") or ""

        auth_keywords = getattr(result, "authkeywords", "") or ""
        keywords = [k.strip() for k in auth_keywords.split("|")] if auth_keywords else []

        volume = getattr(result, "volume", None)
        issue = getattr(result, "issueIdentifier", None)
        pages = getattr(result, "pageRange", None)
        issn = getattr(result, "issn", None)
        eid = getattr(result, "eid", "") or ""

        source_id = eid

        return BibRecord(
            source_db="scopus",
            source_id=source_id,
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
            keywords=keywords,
            mesh_terms=[],
            language=None,
            publication_type=getattr(result, "subtypeDescription", None),
            url=f"https://doi.org/{doi}" if doi else None,
        )
