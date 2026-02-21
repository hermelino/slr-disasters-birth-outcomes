"""Extrator de registros do PubMed via E-utilities (Biopython)."""

import logging
import time
from typing import List

from Bio import Entrez

from src.config import Config
from src.extractors.base import BaseExtractor
from src.models import BibRecord

logger = logging.getLogger(__name__)

# Blocos da estratégia de busca
DISASTERS_MESH = [
    '"Disasters"[MeSH]',
    '"Natural Disasters"[MeSH]',
    '"Floods"[MeSH]',
    '"Earthquakes"[MeSH]',
    '"Cyclonic Storms"[MeSH]',
    '"Tsunamis"[MeSH]',
    '"Volcanic Eruptions"[MeSH]',
    '"Wildfires"[MeSH]',
    '"Droughts"[MeSH]',
    '"Chemical Hazard Release"[MeSH]',
    '"Radioactive Hazard Release"[MeSH]',
    '"Mass Casualty Incidents"[MeSH]',
]

DISASTERS_FREE = [
    "hurricane*[tiab]",
    "typhoon*[tiab]",
    "tornado*[tiab]",
    "flood*[tiab]",
    "earthquake*[tiab]",
    "wildfire*[tiab]",
    "tsunami*[tiab]",
    '"oil spill"[tiab]',
    '"chemical spill"[tiab]',
    '"nuclear accident"[tiab]',
    '"environmental disaster"[tiab]',
    '"climate disaster"[tiab]',
    "landslide*[tiab]",
    "mudslide*[tiab]",
    '"dam failure"[tiab]',
    '"dam collapse"[tiab]',
    '"mining disaster"[tiab]',
    '"volcanic eruption"[tiab]',
    "cyclone*[tiab]",
]

GESTATIONAL_MESH = [
    '"Pregnancy Outcome"[MeSH]',
    '"Pregnancy Complications"[MeSH]',
    '"Premature Birth"[MeSH]',
    '"Infant, Low Birth Weight"[MeSH]',
    '"Infant, Small for Gestational Age"[MeSH]',
    '"Gestational Age"[MeSH]',
    '"Abortion, Spontaneous"[MeSH]',
    '"Stillbirth"[MeSH]',
    '"Pre-Eclampsia"[MeSH]',
    '"Eclampsia"[MeSH]',
    '"Fetal Death"[MeSH]',
    '"Birth Weight"[MeSH]',
    '"Maternal Mortality"[MeSH]',
]

GESTATIONAL_FREE = [
    '"preterm birth"[tiab]',
    '"preterm delivery"[tiab]',
    '"low birth weight"[tiab]',
    "LBW[tiab]",
    "SGA[tiab]",
    '"gestational age"[tiab]',
    '"premature rupture"[tiab]',
    "miscarriage[tiab]",
    '"birth outcome"[tiab]',
    '"adverse birth"[tiab]',
    '"pregnancy loss"[tiab]',
    '"gestational diabetes"[tiab]',
    '"gestational hypertension"[tiab]',
    '"intrauterine growth restriction"[tiab]',
    "IUGR[tiab]",
]

NEONATAL_MESH = [
    '"Infant Mortality"[MeSH]',
    '"Perinatal Mortality"[MeSH]',
    '"Apgar Score"[MeSH]',
    '"Intensive Care, Neonatal"[MeSH]',
    '"Congenital Abnormalities"[MeSH]',
    '"Infant, Newborn, Diseases"[MeSH]',
    '"Fetal Growth Retardation"[MeSH]',
]

NEONATAL_FREE = [
    '"infant mortality"[tiab]',
    '"neonatal mortality"[tiab]',
    '"perinatal mortality"[tiab]',
    '"neonatal death"[tiab]',
    '"infant death"[tiab]',
    "APGAR[tiab]",
    "NICU[tiab]",
    '"congenital anomal*"[tiab]',
    '"congenital malformation*"[tiab]',
    '"birth defect*"[tiab]',
    '"neonatal outcome*"[tiab]',
]

BATCH_SIZE = 500
MAX_RETRIES = 5


class PubMedExtractor(BaseExtractor):
    """Extrator de registros do PubMed usando NCBI E-utilities."""

    def __init__(self, config: Config):
        super().__init__(config)
        Entrez.email = config.api.pubmed.email
        Entrez.tool = config.api.pubmed.tool_name
        if config.api.pubmed.api_key:
            Entrez.api_key = config.api.pubmed.api_key
        self._rate_delay = 0.11 if config.api.pubmed.api_key else 0.34
        self._webenv = None
        self._query_key = None

    def build_query(self) -> str:
        block1 = "(" + " OR ".join(DISASTERS_MESH + DISASTERS_FREE) + ")"
        block2 = "(" + " OR ".join(GESTATIONAL_MESH + GESTATIONAL_FREE) + ")"
        block3 = "(" + " OR ".join(NEONATAL_MESH + NEONATAL_FREE) + ")"

        outcomes = f"({block2} OR {block3})"
        query = f"{block1} AND {outcomes}"

        dr = self.config.search.date_range
        query += f' AND ("{dr.start}"[Date - Publication] : "{dr.end}"[Date - Publication])'

        self.search_query = query
        logger.info("Query PubMed construída (%d caracteres)", len(query))
        return query

    def search(self) -> int:
        if not self.search_query:
            self.build_query()

        logger.info("Executando busca no PubMed...")
        handle = Entrez.esearch(
            db="pubmed",
            term=self.search_query,
            usehistory="y",
            retmax=0,
        )
        results = Entrez.read(handle)
        handle.close()

        self.total_results = int(results["Count"])
        self._webenv = results["WebEnv"]
        self._query_key = results["QueryKey"]

        logger.info(
            "PubMed: %d registros encontrados (WebEnv=%s, QueryKey=%s)",
            self.total_results,
            self._webenv[:20] + "...",
            self._query_key,
        )
        return self.total_results

    def fetch_records(self) -> List[BibRecord]:
        if self._webenv is None:
            self.search()

        self.records = []
        total = min(self.total_results, self.config.search.max_results_per_db)
        logger.info("Recuperando %d registros do PubMed em lotes de %d...", total, BATCH_SIZE)

        for start in range(0, total, BATCH_SIZE):
            batch_records = self._fetch_batch(start, min(BATCH_SIZE, total - start))
            self.records.extend(batch_records)
            logger.info(
                "PubMed: lote %d-%d recuperado (%d registros acumulados)",
                start,
                start + len(batch_records),
                len(self.records),
            )
            time.sleep(self._rate_delay)

        logger.info("PubMed: %d registros recuperados no total", len(self.records))
        return self.records

    def _fetch_batch(self, retstart: int, retmax: int) -> List[BibRecord]:
        for attempt in range(MAX_RETRIES):
            try:
                handle = Entrez.efetch(
                    db="pubmed",
                    rettype="xml",
                    retmode="xml",
                    retstart=retstart,
                    retmax=retmax,
                    webenv=self._webenv,
                    query_key=self._query_key,
                )
                data = Entrez.read(handle)
                handle.close()
                return [self._parse_article(art) for art in data.get("PubmedArticle", [])]
            except Exception as e:
                wait = 2 ** attempt
                logger.warning(
                    "PubMed fetch falhou (tentativa %d/%d): %s. Aguardando %ds...",
                    attempt + 1,
                    MAX_RETRIES,
                    e,
                    wait,
                )
                time.sleep(wait)
        logger.error("PubMed: falha ao recuperar lote retstart=%d após %d tentativas", retstart, MAX_RETRIES)
        return []

    def _parse_article(self, article: dict) -> BibRecord:
        medline = article.get("MedlineCitation", {})
        art_data = medline.get("Article", {})
        journal_info = art_data.get("Journal", {})
        journal_issue = journal_info.get("JournalIssue", {})

        # PMID
        pmid = str(medline.get("PMID", ""))

        # Título
        title = str(art_data.get("ArticleTitle", ""))

        # Autores
        authors = []
        for author in art_data.get("AuthorList", []):
            last = author.get("LastName", "")
            initials = author.get("Initials", "")
            if last:
                authors.append(f"{last}, {initials}" if initials else last)

        # Periódico
        journal = str(journal_info.get("Title", ""))

        # Ano
        year = None
        pub_date = journal_issue.get("PubDate", {})
        year_str = pub_date.get("Year", "")
        if year_str:
            try:
                year = int(year_str)
            except ValueError:
                pass
        if year is None:
            medline_date = pub_date.get("MedlineDate", "")
            if medline_date and len(medline_date) >= 4:
                try:
                    year = int(medline_date[:4])
                except ValueError:
                    pass

        # Volume, Issue, Pages
        volume = str(journal_issue.get("Volume", "") or "")
        issue = str(journal_issue.get("Issue", "") or "")
        pagination = art_data.get("Pagination", {})
        pages = str(pagination.get("MedlinePgn", "") or "")

        # Abstract
        abstract_parts = art_data.get("Abstract", {}).get("AbstractText", [])
        abstract = " ".join(str(part) for part in abstract_parts)

        # MeSH terms
        mesh_terms = []
        for mesh in medline.get("MeshHeadingList", []):
            descriptor = mesh.get("DescriptorName", "")
            if descriptor:
                mesh_terms.append(str(descriptor))

        # Keywords
        keywords = []
        for kw_list in medline.get("KeywordList", []):
            for kw in kw_list:
                keywords.append(str(kw))

        # DOI
        doi = None
        article_ids = article.get("PubmedData", {}).get("ArticleIdList", [])
        for aid in article_ids:
            if hasattr(aid, "attributes") and aid.attributes.get("IdType") == "doi":
                doi = str(aid)
                break

        # Idioma
        languages = art_data.get("Language", [])
        language = str(languages[0]) if languages else None

        # ISSN
        issn_elem = journal_info.get("ISSN", "")
        issn = str(issn_elem) if issn_elem else None

        # Publication type
        pub_types = art_data.get("PublicationTypeList", [])
        pub_type = str(pub_types[0]) if pub_types else None

        return BibRecord(
            source_db="pubmed",
            source_id=pmid,
            doi=doi,
            title=title,
            authors=authors,
            journal=journal,
            year=year,
            volume=volume if volume else None,
            issue=issue if issue else None,
            pages=pages if pages else None,
            issn=issn,
            abstract=abstract,
            keywords=keywords,
            mesh_terms=mesh_terms,
            language=language,
            publication_type=pub_type,
            url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else None,
        )
