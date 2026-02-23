"""Extrator de registros do Scopus via pybliometrics ou exportação manual.

Suporta dois modos:
  1. API automática (pybliometrics) — requer API key Elsevier
  2. Exportação manual — o usuário busca pelo proxy CAPES e exporta CSV/RIS

Fluxo manual (proxy CAPES):
  1. build_query()  → gera a query
  2. O usuário acessa Scopus via proxy e executa a busca
  3. Exporta como CSV ou RIS
  4. import_from_file() → lê o arquivo e converte para BibRecord
"""

import csv
import logging
import math
import re
import time
from pathlib import Path
from typing import List, Optional

import rispy

from src.config import Config
from src.extractors.base import BaseExtractor
from src.models import BibRecord

logger = logging.getLogger(__name__)

SCOPUS_MAX_RESULTS = 5000
SCOPUS_PROXY_URL = "https://www-scopus-com.ez11.periodicos.capes.gov.br/search/form.uri?display=advanced"

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
    """Extrator de registros do Scopus (API ou manual)."""

    def __init__(self, config: Config):
        super().__init__(config)
        self._has_api_key = bool(config.api.scopus.api_key)
        if self._has_api_key:
            _setup_pybliometrics(config)

    def _disasters_terms_plain(self) -> str:
        """Termos de desastres sem prefixo TITLE-ABS-KEY."""
        return " OR ".join(DISASTERS_TERMS)

    def _outcomes_terms_plain(self) -> str:
        """Termos de desfechos sem prefixo TITLE-ABS-KEY."""
        return " OR ".join(OUTCOMES_TERMS)

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

    # ------------------------------------------------------------------
    # Modo manual (proxy institucional)
    # ------------------------------------------------------------------

    def print_manual_instructions(self):
        """Imprime instruções para exportação manual via proxy Scopus."""
        if not self.search_query:
            self.build_query()

        dr = self.config.search.date_range
        start_year = dr.start.split("/")[0]
        end_year = dr.end.split("/")[0]

        print()
        print("=" * 70)
        print("  SCOPUS — EXPORTAÇÃO MANUAL (PROXY CAPES)")
        print("=" * 70)
        print()
        print("  1. Abra o Scopus via proxy institucional:")
        print(f"     {SCOPUS_PROXY_URL}")
        print()
        print("  2. Na busca, preencha DUAS linhas:")
        print()
        print("     Linha 1:  Article title, Abstract, Keywords")
        print("               [colar conteúdo da LINHA 1 do scopus_query.txt]")
        print("     Operador: AND")
        print("     Linha 2:  Article title, Abstract, Keywords")
        print("               [colar conteúdo da LINHA 2 do scopus_query.txt]")
        print()
        print("  3. Clique em '+ Add date range' e defina:")
        print(f"     De: {start_year}  Até: {end_year}")
        print()
        print("  DICA: As linhas estão salvas em 'scopus_query.txt'")
        print("  no diretório de saída.")
        print()
        print("  4. Clique em 'Search' e verifique o total")
        print()
        print("  5. Exporte todos os registros:")
        print("     - Clique em 'Export' → 'CSV' ou 'RIS'")
        print("     - Selecione todos os campos (Citation + Abstract)")
        print("     - Salve o arquivo")
        print()
        print("  6. Execute novamente com --import-scopus:")
        print("     python main.py --import-scopus caminho/arquivo.csv")
        print("=" * 70)
        print()

    def save_query_file(self, output_dir: str):
        """Salva a query Scopus em arquivo para referência."""
        if not self.search_query:
            self.build_query()

        dr = self.config.search.date_range
        start_year = dr.start.split("/")[0]
        end_year = dr.end.split("/")[0]

        disasters = self._disasters_terms_plain()
        outcomes = self._outcomes_terms_plain()

        query_file = Path(output_dir) / "scopus_query.txt"
        content = (
            "Scopus — Query de Busca (Fielded Search)\n"
            + "=" * 60 + "\n\n"
            "LINHA 1 — Article title, Abstract, Keywords (copiar e colar):\n"
            + "-" * 60 + "\n"
            f"{disasters}\n\n"
            "LINHA 2 — Article title, Abstract, Keywords (operador AND):\n"
            + "-" * 60 + "\n"
            f"{outcomes}\n\n"
            f"DATE RANGE: {start_year} - {end_year}\n"
            "(usar botão '+ Add date range')\n\n"
            + "=" * 60 + "\n"
            "Instruções:\n"
            "  1. Abra: " + SCOPUS_PROXY_URL + "\n"
            "  2. Linha 1: 'Article title, Abstract, Keywords', cole LINHA 1\n"
            "  3. Operador AND, adicionar linha\n"
            "  4. Linha 2: 'Article title, Abstract, Keywords', cole LINHA 2\n"
            "  5. '+ Add date range': " + start_year + " a " + end_year + "\n"
            "  6. Clique 'Search'\n"
            "  7. Exporte como CSV (com Abstract) ou RIS\n"
            "  8. Use: python main.py --import-scopus arquivo.csv\n\n"
            + "=" * 60 + "\n"
            "Query completa (formato API, para referência):\n"
            f"{self.search_query}\n"
        )
        query_file.write_text(content, encoding="utf-8")
        logger.info("Query Scopus salva em: %s", query_file)

    @staticmethod
    def import_from_file(file_path: str) -> List[BibRecord]:
        """Importa registros de arquivo CSV ou RIS exportado do Scopus."""
        path = Path(file_path)
        if not path.exists():
            raise FileNotFoundError(f"Arquivo não encontrado: {file_path}")

        suffix = path.suffix.lower()
        if suffix == ".ris":
            return ScopusExtractor._import_ris(path)
        elif suffix == ".csv":
            return ScopusExtractor._import_csv(path)
        else:
            raise ValueError(
                f"Formato não suportado: {suffix}. Use .csv ou .ris exportado do Scopus."
            )

    @staticmethod
    def _import_ris(path: Path) -> List[BibRecord]:
        """Importa registros de arquivo RIS exportado do Scopus."""
        logger.info("Importando RIS do Scopus: %s", path)

        entries = None
        for encoding in ["utf-8-sig", "utf-8", "latin-1", "cp1252"]:
            try:
                with open(path, "r", encoding=encoding) as f:
                    entries = rispy.load(f)
                break
            except (UnicodeDecodeError, UnicodeError):
                continue

        if entries is None:
            raise ValueError(f"Não foi possível ler {path} com nenhum encoding suportado")

        records = []
        for entry in entries:
            doi = entry.get("doi")
            if not doi:
                urls = entry.get("urls", [])
                for u in urls:
                    doi_match = re.search(r'10\.\d{4,}/[^\s]+', u)
                    if doi_match:
                        doi = doi_match.group(0)
                        break

            year = None
            year_str = entry.get("year", "")
            if year_str:
                try:
                    year = int(str(year_str)[:4])
                except ValueError:
                    pass

            pages = None
            sp = entry.get("start_page", "")
            ep = entry.get("end_page", "")
            if sp:
                pages = f"{sp}-{ep}" if ep else sp

            url = None
            urls = entry.get("urls", [])
            if urls:
                url = urls[0]

            rec = BibRecord(
                source_db="scopus",
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

        logger.info("Scopus RIS: %d registros importados de %s", len(records), path.name)
        return records

    @staticmethod
    def _import_csv(path: Path) -> List[BibRecord]:
        """Importa registros de arquivo CSV exportado do Scopus."""
        logger.info("Importando CSV do Scopus: %s", path)

        records = []
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
            title = row.get("Title", "") or ""

            authors_str = row.get("Authors", "") or ""
            authors = [a.strip() for a in authors_str.split(",") if a.strip()]

            year = None
            year_str = row.get("Year", "") or ""
            if year_str:
                try:
                    year = int(str(year_str).strip()[:4])
                except ValueError:
                    pass

            doi = row.get("DOI", "") or None
            if doi:
                doi = doi.strip()

            abstract = row.get("Abstract", "") or ""

            keywords_str = row.get("Author Keywords", "") or ""
            keywords = [k.strip() for k in keywords_str.split(";") if k.strip()]

            # Index Keywords (Scopus-specific)
            idx_kw = row.get("Index Keywords", "") or ""
            if idx_kw:
                keywords.extend([k.strip() for k in idx_kw.split(";") if k.strip()])

            eid = row.get("EID", "") or ""

            rec = BibRecord(
                source_db="scopus",
                source_id=eid,
                doi=doi if doi else None,
                title=title,
                authors=authors,
                journal=row.get("Source title", "") or "",
                year=year,
                volume=row.get("Volume", "") or None,
                issue=row.get("Issue", "") or None,
                pages=row.get("Page start", "") or None,
                issn=row.get("ISSN", "") or None,
                abstract=abstract,
                keywords=keywords,
                mesh_terms=[],
                language=row.get("Language of Original Document", "") or None,
                publication_type=row.get("Document Type", "") or None,
                url=row.get("Link", "") or (f"https://doi.org/{doi.strip()}" if doi and doi.strip() else None),
            )
            records.append(rec)

        logger.info("Scopus CSV: %d registros importados de %s", len(records), path.name)
        return records

    # ------------------------------------------------------------------
    # API (modo automático)
    # ------------------------------------------------------------------

    def search(self) -> int:
        if not self.search_query:
            self.build_query()

        if not self._has_api_key:
            logger.info(
                "Scopus: sem API key configurada. "
                "Use --import-scopus <arquivo> para importar registros exportados manualmente."
            )
            self.print_manual_instructions()
            self.total_results = 0
            return 0

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

    def _search_with_fallback(self, query: str) -> list:
        """Executa ScopusSearch tentando subscriber=True (COMPLETE view) primeiro."""
        from pybliometrics.scopus import ScopusSearch

        try:
            s = ScopusSearch(query, download=True, subscriber=True)
            return s.results or []
        except Exception as e:
            logger.warning(
                "Scopus subscriber=True falhou: %s. Tentando subscriber=False...", e
            )
            s = ScopusSearch(query, download=True, subscriber=False)
            return s.results or []

    def fetch_records(self) -> List[BibRecord]:
        if not self._has_api_key:
            self.records = []
            return self.records

        if not self.search_query:
            self.build_query()

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
            raw_results = self._search_with_fallback(self.search_query)
            self.records = [self._parse_result(r) for r in raw_results]

        # Enriquecer registros sem abstract/autores via AbstractRetrieval
        self._enrich_records(self.records)

        logger.info("Scopus: %d registros recuperados no total", len(self.records))
        return self.records

    def _fetch_by_year_ranges(self) -> List[BibRecord]:
        """Divide a busca em faixas de 5 anos quando o total excede 5000."""
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
                raw = self._search_with_fallback(range_query)
                batch = [self._parse_result(r) for r in raw]
                all_records.extend(batch)
                logger.info("Scopus %d-%d: %d registros", y_start, y_end, len(batch))
            except Exception as e:
                logger.error("Scopus %d-%d falhou: %s", y_start, y_end, e)

        return all_records

    def _enrich_records(self, records: List[BibRecord]) -> None:
        """Enriquece registros sem abstract/autores via AbstractRetrieval (META_ABS)."""
        from pybliometrics.scopus import AbstractRetrieval

        to_enrich = [r for r in records if not r.abstract or not r.authors]
        if not to_enrich:
            return

        logger.info(
            "Scopus: %d/%d registros sem abstract/autores. "
            "Enriquecendo via AbstractRetrieval...",
            len(to_enrich), len(records),
        )

        enriched = 0
        for i, rec in enumerate(to_enrich):
            eid = rec.source_id
            if not eid:
                continue
            try:
                ab = AbstractRetrieval(eid, view="META_ABS")
                if not rec.abstract and ab.description:
                    rec.abstract = ab.description
                if not rec.authors and ab.authors:
                    rec.authors = [
                        f"{a.surname}, {a.initials}"
                        for a in ab.authors
                        if a.surname
                    ]
                if not rec.keywords and ab.authkeywords:
                    rec.keywords = list(ab.authkeywords)
                enriched += 1
            except Exception as e:
                logger.debug("AbstractRetrieval falhou para %s: %s", eid, e)
            time.sleep(0.1)
            if (i + 1) % 100 == 0:
                logger.info("Enriquecimento: %d/%d processados", i + 1, len(to_enrich))

        logger.info("Scopus: %d registros enriquecidos via AbstractRetrieval", enriched)

    def _parse_result(self, result) -> BibRecord:
        # pybliometrics retorna namedtuples com campos acessíveis por atributo
        doi = getattr(result, "doi", None)
        title = getattr(result, "title", "") or ""

        # COMPLETE view: author_names; STANDARD view: apenas creator (1º autor)
        author_names = getattr(result, "author_names", "") or ""
        if not author_names:
            creator = getattr(result, "creator", "") or ""
            author_names = creator
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
