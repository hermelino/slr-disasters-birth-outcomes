"""Extrator de registros do PubMed via E-utilities (Biopython)."""

import logging
import re
import time
from typing import Dict, List, Set

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

# Environmental Pollution [N06.850.460] 


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


def _extract_mesh_name(term: str) -> str:
    """Extrai nome limpo do MeSH: '"Disasters"[MeSH]' → 'Disasters'."""
    m = re.search(r'"(.+?)"', term)
    return m.group(1) if m else term


DISASTERS_NAMES: Set[str] = {_extract_mesh_name(t) for t in DISASTERS_MESH}
GESTATIONAL_NAMES: Set[str] = {_extract_mesh_name(t) for t in GESTATIONAL_MESH}
NEONATAL_NAMES: Set[str] = {_extract_mesh_name(t) for t in NEONATAL_MESH}


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
        # PMIDs por categoria (preenchidos via sub-buscas com MeSH explosion)
        self._disaster_pmids: Set[str] = set()
        self._gestational_pmids: Set[str] = set()
        self._neonatal_pmids: Set[str] = set()

    def build_query(self) -> str:
        block1 = "(" + " OR ".join(DISASTERS_MESH) + ")"
        block2 = "(" + " OR ".join(GESTATIONAL_MESH) + ")"
        block3 = "(" + " OR ".join(NEONATAL_MESH) + ")"

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

    def _search_category_pmids(self) -> None:
        """Executa sub-buscas para classificar PMIDs por categoria MeSH (com explosion)."""
        dr = self.config.search.date_range
        date_filter = f'("{dr.start}"[Date - Publication] : "{dr.end}"[Date - Publication])'
        max_ret = self.config.search.max_results_per_db

        block1 = "(" + " OR ".join(DISASTERS_MESH) + ")"
        block2 = "(" + " OR ".join(GESTATIONAL_MESH) + ")"
        block3 = "(" + " OR ".join(NEONATAL_MESH) + ")"
        outcomes = f"({block2} OR {block3})"

        queries = {
            "disasters": f"{block1} AND {outcomes} AND {date_filter}",
            "gestational": f"{block1} AND {block2} AND {date_filter}",
            "neonatal": f"{block1} AND {block3} AND {date_filter}",
        }

        for label, query in queries.items():
            try:
                handle = Entrez.esearch(db="pubmed", term=query, retmax=max_ret)
                results = Entrez.read(handle)
                handle.close()
                pmids = set(results.get("IdList", []))
                if label == "disasters":
                    self._disaster_pmids = pmids
                elif label == "gestational":
                    self._gestational_pmids = pmids
                else:
                    self._neonatal_pmids = pmids
                logger.info("Sub-busca %s: %d PMIDs", label, len(pmids))
                time.sleep(self._rate_delay)
            except Exception as e:
                logger.warning("Sub-busca %s falhou: %s", label, e)

    def _fetch_mesh_trees_batch(self, term_names: Set[str]) -> Dict[str, List[str]]:
        """
        Obtém tree numbers para um conjunto de descritores MeSH via E-utilities.
        Retorna {term_name: [tree_numbers]}.
        """
        if not term_names:
            return {}

        result: Dict[str, List[str]] = {}
        names_list = sorted(term_names)

        # Batch esearch: grupos de 50 termos por query
        all_uids: List[str] = []
        es_batch = 50
        for i in range(0, len(names_list), es_batch):
            batch = names_list[i : i + es_batch]
            query = " OR ".join(f'"{n}"[MeSH Terms]' for n in batch)
            try:
                handle = Entrez.esearch(
                    db="mesh", term=query, retmax=len(batch) * 2
                )
                res = Entrez.read(handle)
                handle.close()
                all_uids.extend(res.get("IdList", []))
                time.sleep(self._rate_delay)
            except Exception as e:
                logger.warning("MeSH esearch batch falhou: %s", e)

        all_uids = list(set(all_uids))
        if not all_uids:
            return result

        # Batch efetch: grupos de 200 UIDs, formato texto legível
        # Formato: "1: TermName\n...\nTree Number(s): X.Y.Z, A.B.C\n..."
        heading_re = re.compile(r"^\d+:\s+(.+)$")
        tree_re = re.compile(r"^Tree Number\(s\):\s+(.+)$")

        ef_batch = 200
        for i in range(0, len(all_uids), ef_batch):
            batch = all_uids[i : i + ef_batch]
            try:
                handle = Entrez.efetch(
                    db="mesh", id=",".join(batch),
                    rettype="full", retmode="text",
                )
                text = handle.read()
                handle.close()
                if isinstance(text, bytes):
                    text = text.decode("utf-8")
                time.sleep(self._rate_delay)

                current_mh = None
                for line in text.split("\n"):
                    stripped = line.strip()
                    hm = heading_re.match(stripped)
                    if hm:
                        current_mh = hm.group(1).strip()
                        continue
                    tm = tree_re.match(stripped)
                    if tm and current_mh and current_mh in term_names:
                        trees = [t.strip() for t in tm.group(1).split(",")]
                        result[current_mh] = trees
            except Exception as e:
                logger.warning("MeSH efetch batch falhou: %s", e)

        logger.info(
            "Tree numbers obtidos para %d/%d termos MeSH",
            len(result), len(term_names),
        )
        return result

    def _resolve_mesh_explosions(self) -> None:
        """
        Substitui marcadores '(MeSH explosion)' pelos nomes reais dos
        descritores MeSH descendentes, usando a hierarquia de tree numbers.
        """
        # 1) Tree numbers dos nossos termos de busca
        search_names = DISASTERS_NAMES | GESTATIONAL_NAMES | NEONATAL_NAMES
        logger.info(
            "Resolvendo MeSH explosion: obtendo tree numbers para %d termos de busca...",
            len(search_names),
        )
        search_trees = self._fetch_mesh_trees_batch(search_names)

        # prefix → categoria
        prefix_map: Dict[str, str] = {}
        for name in DISASTERS_NAMES:
            for tn in search_trees.get(name, []):
                prefix_map[tn] = "disasters"
        for name in GESTATIONAL_NAMES:
            for tn in search_trees.get(name, []):
                prefix_map[tn] = "gestational"
        for name in NEONATAL_NAMES:
            for tn in search_trees.get(name, []):
                prefix_map[tn] = "neonatal"

        if not prefix_map:
            logger.warning("Nenhum tree number encontrado — explosion não resolvido")
            return

        # 2) Coletar termos únicos dos artigos com marcador explosion
        terms_to_resolve: Set[str] = set()
        for rec in self.records:
            has_expl = (
                "(MeSH explosion)" in rec.matched_disasters_mesh
                or "(MeSH explosion)" in rec.matched_gestational_mesh
                or "(MeSH explosion)" in rec.matched_neonatal_mesh
            )
            if has_expl:
                terms_to_resolve.update(rec.mesh_terms)
        terms_to_resolve -= search_names

        if not terms_to_resolve:
            return

        logger.info(
            "Obtendo tree numbers para %d termos de artigos...",
            len(terms_to_resolve),
        )
        article_trees = self._fetch_mesh_trees_batch(terms_to_resolve)

        # 3) Mapear cada termo de artigo → categorias que ele descende
        descendant_cats: Dict[str, Set[str]] = {}
        for term, trees in article_trees.items():
            for tn in trees:
                for prefix, category in prefix_map.items():
                    if tn == prefix or tn.startswith(prefix + "."):
                        descendant_cats.setdefault(term, set()).add(category)

        # 4) Atualizar registros
        resolved = 0
        for rec in self.records:
            for attr, cat in [
                ("matched_disasters_mesh", "disasters"),
                ("matched_gestational_mesh", "gestational"),
                ("matched_neonatal_mesh", "neonatal"),
            ]:
                matched = getattr(rec, attr)
                if "(MeSH explosion)" in matched:
                    names = [
                        t for t in rec.mesh_terms
                        if t in descendant_cats and cat in descendant_cats[t]
                    ]
                    if names:
                        setattr(rec, attr, names)
                        resolved += 1

        logger.info(
            "MeSH explosion resolvido: %d campos atualizados, %d termos descendentes identificados",
            resolved, len(descendant_cats),
        )

    def fetch_records(self) -> List[BibRecord]:
        if self._webenv is None:
            self.search()

        # Classificar PMIDs por categoria (usa MeSH explosion do PubMed)
        logger.info("Executando sub-buscas por categoria MeSH...")
        self._search_category_pmids()

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

        # Resolver "(MeSH explosion)" → nomes reais dos descendentes
        self._resolve_mesh_explosions()

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
        matched_disasters = []
        matched_gestational = []
        matched_neonatal = []
        for mesh in medline.get("MeshHeadingList", []):
            descriptor = mesh.get("DescriptorName", "")
            if descriptor:
                name = str(descriptor)
                mesh_terms.append(name)
                if name in DISASTERS_NAMES:
                    matched_disasters.append(name)
                if name in GESTATIONAL_NAMES:
                    matched_gestational.append(name)
                if name in NEONATAL_NAMES:
                    matched_neonatal.append(name)

        # Complementar com MeSH explosion (sub-busca por PMID)
        if not matched_disasters and pmid in self._disaster_pmids:
            matched_disasters = ["(MeSH explosion)"]
        if not matched_gestational and pmid in self._gestational_pmids:
            matched_gestational = ["(MeSH explosion)"]
        if not matched_neonatal and pmid in self._neonatal_pmids:
            matched_neonatal = ["(MeSH explosion)"]

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
            matched_disasters_mesh=matched_disasters,
            matched_gestational_mesh=matched_gestational,
            matched_neonatal_mesh=matched_neonatal,
            language=language,
            publication_type=pub_type,
            url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else None,
        )
