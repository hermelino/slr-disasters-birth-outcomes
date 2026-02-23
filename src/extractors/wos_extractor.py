"""Extrator de registros do Web of Science.

Suporta dois modos de operação:
  1. API automática (Starter API) — requer api_key da Clarivate
  2. Exportação manual — o usuário busca pelo proxy institucional e exporta RIS/CSV

Fluxo manual (proxy CAPES):
  1. build_query()  → gera a query em sintaxe WoS
  2. O usuário acessa o portal WoS via proxy e executa a busca
  3. Exporta os resultados como arquivo (RIS, CSV ou Tab-delimited)
  4. import_from_file() → lê o arquivo e converte para BibRecord
"""

import csv
import logging
import re
import time
from pathlib import Path
from typing import List, Optional

import rispy

from src.config import Config
from src.extractors.base import BaseExtractor
from src.models import BibRecord

logger = logging.getLogger(__name__)

WOS_API_BASE = "https://api.clarivate.com/apis/wos-starter/v1"
WOS_PAGE_LIMIT = 50
MAX_RETRIES = 5

WOS_PROXY_URL = "https://www-webofscience-com.ez11.periodicos.capes.gov.br/wos/woscc/smart-search"

# Blocos de busca em sintaxe WoS (TS = Topic Search: título + abstract + keywords)
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


class WosExtractor(BaseExtractor):
    """Extrator de registros do Web of Science (API ou manual)."""

    def __init__(self, config: Config):
        super().__init__(config)
        self._has_api_key = bool(config.api.wos.api_key)
        if self._has_api_key:
            import requests
            self._api_key = config.api.wos.api_key
            self._session = requests.Session()
            self._session.headers.update({
                "X-ApiKey": self._api_key,
                "Accept": "application/json",
            })

    def _disasters_terms_plain(self) -> str:
        """Termos de desastres sem prefixo TS= (para Fielded Search)."""
        return " OR ".join(DISASTERS_TERMS)

    def _outcomes_terms_plain(self) -> str:
        """Termos de desfechos sem prefixo TS= (para Fielded Search)."""
        return " OR ".join(OUTCOMES_TERMS)

    def build_query(self) -> str:
        """Constrói query WoS em sintaxe de API (com TS= e PY=)."""
        disasters_block = "TS=(" + self._disasters_terms_plain() + ")"
        outcomes_block = "TS=(" + self._outcomes_terms_plain() + ")"

        dr = self.config.search.date_range
        start_year = dr.start.split("/")[0]
        end_year = dr.end.split("/")[0]

        query = f"{disasters_block} AND {outcomes_block} AND PY=({start_year}-{end_year})"

        self.search_query = query
        logger.info("Query WoS construída (%d caracteres)", len(query))
        return query

    # ------------------------------------------------------------------
    # Modo manual (proxy institucional)
    # ------------------------------------------------------------------

    def print_manual_instructions(self):
        """Imprime instruções para exportação manual via proxy WoS."""
        if not self.search_query:
            self.build_query()

        dr = self.config.search.date_range
        start_year = dr.start.split("/")[0]
        end_year = dr.end.split("/")[0]

        print()
        print("=" * 70)
        print("  WEB OF SCIENCE — EXPORTAÇÃO MANUAL (PROXY CAPES)")
        print("=" * 70)
        print()
        print("  1. Abra o portal WoS via proxy institucional:")
        print(f"     {WOS_PROXY_URL}")
        print()
        print("  2. Na aba 'FIELDED SEARCH', preencha DUAS linhas:")
        print()
        print("     Linha 1:  Topic  =  [colar conteúdo da LINHA 1 do wos_query.txt]")
        print("     Operador: AND")
        print("     Linha 2:  Topic  =  [colar conteúdo da LINHA 2 do wos_query.txt]")
        print()
        print("  3. Clique em '+ Add date range' e defina:")
        print(f"     De: {start_year}  Até: {end_year}")
        print()
        print("  DICA: As linhas estão salvas separadamente no arquivo")
        print("  'wos_query.txt' no diretório de saída.")
        print()
        print("  4. Clique em 'Search' e verifique o total de resultados")
        print()
        print("  5. Exporte todos os registros:")
        print("     - Clique em 'Export' → 'RIS' ou 'Tab-delimited'")
        print("     - Selecione 'Full Record' para incluir abstract")
        print("     - Salve o arquivo")
        print()
        print("  6. Execute novamente com --import-wos:")
        print("     python main.py --import-wos caminho/arquivo.ris")
        print("=" * 70)
        print()

    def save_query_file(self, output_dir: str):
        """Salva a query WoS em arquivo para referência."""
        if not self.search_query:
            self.build_query()

        dr = self.config.search.date_range
        start_year = dr.start.split("/")[0]
        end_year = dr.end.split("/")[0]

        disasters = self._disasters_terms_plain()
        outcomes = self._outcomes_terms_plain()

        query_file = Path(output_dir) / "wos_query.txt"
        content = (
            "Web of Science — Query de Busca (Fielded Search)\n"
            + "=" * 60 + "\n\n"
            "LINHA 1 — Topic (copiar e colar no 1º campo):\n"
            + "-" * 60 + "\n"
            f"{disasters}\n\n"
            "LINHA 2 — Topic (copiar e colar no 2º campo, operador AND):\n"
            + "-" * 60 + "\n"
            f"{outcomes}\n\n"
            f"DATE RANGE: {start_year} - {end_year}\n"
            "(usar botão '+ Add date range')\n\n"
            + "=" * 60 + "\n"
            "Instruções:\n"
            "  1. Abra: " + WOS_PROXY_URL + "\n"
            "  2. Aba 'FIELDED SEARCH'\n"
            "  3. Linha 1: dropdown 'Topic', cole a LINHA 1 acima\n"
            "  4. Clique '+ Add row', operador AND\n"
            "  5. Linha 2: dropdown 'Topic', cole a LINHA 2 acima\n"
            "  6. Clique '+ Add date range', defina " + start_year + " a " + end_year + "\n"
            "  7. Clique 'Search'\n"
            "  8. Exporte como RIS (Full Record)\n"
            "  9. Use: python main.py --import-wos arquivo.ris\n\n"
            + "=" * 60 + "\n"
            "Query completa (formato API, para referência):\n"
            f"{self.search_query}\n"
        )
        query_file.write_text(content, encoding="utf-8")
        logger.info("Query WoS salva em: %s", query_file)

    def search(self) -> int:
        if not self.search_query:
            self.build_query()

        if not self._has_api_key:
            logger.info(
                "WoS: sem API key configurada. "
                "Use --import-wos <arquivo> para importar registros exportados manualmente."
            )
            self.print_manual_instructions()
            self.total_results = 0
            return 0

        logger.info("Executando busca no Web of Science (API)...")
        resp = self._api_request("/documents", params={
            "db": "WOS",
            "q": self.search_query,
            "limit": 1,
        })

        self.total_results = resp.get("metadata", {}).get("total", 0)
        logger.info("Web of Science: %d registros encontrados", self.total_results)
        return self.total_results

    def fetch_records(self) -> List[BibRecord]:
        if not self._has_api_key:
            self.records = []
            return self.records

        if not self.search_query:
            self.build_query()

        total = min(self.total_results or 0, self.config.search.max_results_per_db)
        if total == 0:
            logger.warning("WoS: nenhum registro para recuperar")
            self.records = []
            return self.records

        logger.info("WoS: recuperando %d registros (páginas de %d)...", total, WOS_PAGE_LIMIT)
        self.records = []
        page = 1

        while len(self.records) < total:
            resp = self._api_request("/documents", params={
                "db": "WOS",
                "q": self.search_query,
                "limit": WOS_PAGE_LIMIT,
                "page": page,
            })

            hits = resp.get("hits", [])
            if not hits:
                break

            for doc in hits:
                self.records.append(self._parse_document(doc))

            logger.info(
                "WoS: página %d recuperada (%d registros acumulados)",
                page, len(self.records),
            )
            page += 1
            time.sleep(0.5)

        self.records = self.records[:total]
        logger.info("WoS: %d registros recuperados no total", len(self.records))
        return self.records

    # ------------------------------------------------------------------
    # Importação de arquivos exportados manualmente
    # ------------------------------------------------------------------

    @staticmethod
    def import_from_file(file_path: str) -> List[BibRecord]:
        """Importa registros de arquivo RIS ou Tab-delimited exportado do WoS."""
        path = Path(file_path)
        if not path.exists():
            raise FileNotFoundError(f"Arquivo não encontrado: {file_path}")

        suffix = path.suffix.lower()
        if suffix == ".ris":
            return WosExtractor._import_ris(path)
        elif suffix in (".txt", ".tsv", ".csv"):
            return WosExtractor._import_tabdelimited(path)
        else:
            raise ValueError(
                f"Formato não suportado: {suffix}. "
                "Use .ris ou .txt/.tsv (Tab-delimited) exportado do WoS."
            )

    @staticmethod
    def _import_ris(path: Path) -> List[BibRecord]:
        """Importa registros de arquivo RIS exportado do WoS."""
        logger.info("Importando RIS do WoS: %s", path)

        # Tentar diferentes encodings
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
                source_db="wos",
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

        logger.info("WoS RIS: %d registros importados de %s", len(records), path.name)
        return records

    @staticmethod
    def _import_tabdelimited(path: Path) -> List[BibRecord]:
        """Importa registros de arquivo Tab-delimited exportado do WoS."""
        logger.info("Importando Tab-delimited do WoS: %s", path)

        records = []
        for encoding in ["utf-8-sig", "utf-8", "latin-1", "cp1252"]:
            try:
                with open(path, "r", encoding=encoding) as f:
                    reader = csv.DictReader(f, delimiter="\t")
                    rows = list(reader)
                if rows:
                    break
            except (UnicodeDecodeError, UnicodeError):
                continue
        else:
            raise ValueError(f"Não foi possível ler {path} com nenhum encoding suportado")

        for row in rows:
            # Campos WoS Tab-delimited (Full Record)
            title = row.get("TI", "") or row.get("Article Title", "") or ""

            authors_str = row.get("AU", "") or row.get("Authors", "") or ""
            authors = [a.strip() for a in authors_str.split(";") if a.strip()]

            year = None
            year_str = row.get("PY", "") or row.get("Publication Year", "") or ""
            if year_str:
                try:
                    year = int(str(year_str).strip()[:4])
                except ValueError:
                    pass

            doi = row.get("DI", "") or row.get("DOI", "") or None
            if doi:
                doi = doi.strip()

            abstract = row.get("AB", "") or row.get("Abstract", "") or ""

            keywords_str = row.get("DE", "") or row.get("Author Keywords", "") or ""
            keywords = [k.strip() for k in keywords_str.split(";") if k.strip()]

            # Keywords Plus (WoS-specific)
            kw_plus = row.get("ID", "") or row.get("Keywords Plus", "") or ""
            if kw_plus:
                keywords.extend([k.strip() for k in kw_plus.split(";") if k.strip()])

            rec = BibRecord(
                source_db="wos",
                source_id=row.get("UT", "") or row.get("Accession Number", "") or "",
                doi=doi if doi else None,
                title=title,
                authors=authors,
                journal=row.get("SO", "") or row.get("Source Title", "") or "",
                year=year,
                volume=row.get("VL", "") or row.get("Volume", "") or None,
                issue=row.get("IS", "") or row.get("Issue", "") or None,
                pages=row.get("BP", "") or row.get("Beginning Page", "") or None,
                issn=row.get("SN", "") or row.get("ISSN", "") or None,
                abstract=abstract,
                keywords=keywords,
                mesh_terms=[],
                language=row.get("LA", "") or row.get("Language", "") or None,
                publication_type=row.get("DT", "") or row.get("Document Type", "") or None,
                url=f"https://doi.org/{doi.strip()}" if doi and doi.strip() else None,
            )
            records.append(rec)

        logger.info("WoS Tab-delimited: %d registros importados de %s", len(records), path.name)
        return records

    # ------------------------------------------------------------------
    # API (modo automático)
    # ------------------------------------------------------------------

    def _api_request(self, endpoint: str, params: dict) -> dict:
        """Executa requisição à API WoS com retry exponencial."""
        import requests as req

        url = WOS_API_BASE + endpoint

        for attempt in range(MAX_RETRIES):
            try:
                resp = self._session.get(url, params=params, timeout=30)
                if resp.status_code == 200:
                    return resp.json()
                elif resp.status_code == 429:
                    wait = 2 ** (attempt + 1)
                    logger.warning(
                        "WoS API rate limit (429). Aguardando %ds...", wait
                    )
                    time.sleep(wait)
                    continue
                else:
                    resp.raise_for_status()
            except req.exceptions.RequestException as e:
                wait = 2 ** attempt
                logger.warning(
                    "WoS API falhou (tentativa %d/%d): %s. Aguardando %ds...",
                    attempt + 1, MAX_RETRIES, e, wait,
                )
                time.sleep(wait)

        logger.error("WoS API: falha após %d tentativas", MAX_RETRIES)
        return {}

    def _parse_document(self, doc: dict) -> BibRecord:
        """Converte documento JSON da API WoS para BibRecord."""
        uid = doc.get("uid", "")
        title = doc.get("title", "")

        names = doc.get("names", {})
        authors_list = names.get("authors", [])
        authors = []
        for author in authors_list:
            display = author.get("displayName", "") or author.get("wosStandard", "")
            if display:
                authors.append(display)

        source = doc.get("source", {})
        journal = source.get("sourceTitle", "")
        year = source.get("publishYear", None)
        if year is not None:
            try:
                year = int(year)
            except (ValueError, TypeError):
                year = None
        volume = source.get("volume", None)
        issue = source.get("issue", None)
        pages_info = source.get("pages", {})
        pages = pages_info.get("range", None)

        identifiers = doc.get("identifiers", {})
        doi = identifiers.get("doi", None)
        issn = identifiers.get("issn", None)

        kw_data = doc.get("keywords", {})
        keywords = kw_data.get("authorKeywords", []) or []

        links = doc.get("links", {})
        url = links.get("record", None)

        return BibRecord(
            source_db="wos",
            source_id=uid,
            doi=doi,
            title=title,
            authors=authors,
            journal=journal,
            year=year,
            volume=volume,
            issue=issue,
            pages=pages,
            issn=issn,
            abstract="",  # API Starter não retorna abstract
            keywords=keywords,
            mesh_terms=[],
            language=None,
            publication_type=source.get("documentType", None),
            url=url or (f"https://doi.org/{doi}" if doi else None),
        )
