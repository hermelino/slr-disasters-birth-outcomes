"""
Extrator Bibliográfico — Desastres Ambientais e Indicadores Gestacionais/Neonatais

Pipeline de extração automatizada de registros científicos do PubMed, Scopus e BVS/LILACS,
com deduplicação e exportação em CSV e RIS para revisão sistemática.

Uso:
    python main.py                          # Executa com config.yaml padrão
    python main.py --dry-run                # Mostra queries sem executar
    python main.py --databases pubmed       # Apenas PubMed
    python main.py --databases pubmed,bvs   # PubMed e BVS
    python main.py --verbose                # Logs detalhados no console
    python main.py --skip-dedup             # Exporta sem deduplicação
    python main.py --import-bvs dados.ris   # Importa BVS exportado manualmente
"""

import argparse
import logging
import sys
from datetime import datetime
from pathlib import Path

from src.config import Config, load_config
from src.dedup.deduplicator import deduplicate
from src.exporters.csv_exporter import export_csv, export_excel
from src.exporters.ris_exporter import export_ris
from src.logging_prisma import (
    DatabaseStats,
    PrismaLog,
    print_summary,
    save_prisma_log,
    setup_logging,
)
from src.models import BibRecord

logger = logging.getLogger(__name__)

AVAILABLE_DATABASES = ["pubmed", "scopus", "bvs", "wos"]


def create_extractor(db_name: str, config: Config):
    """Instancia o extrator para a base de dados especificada."""
    if db_name == "pubmed":
        from src.extractors.pubmed_extractor import PubMedExtractor
        return PubMedExtractor(config)
    elif db_name == "scopus":
        from src.extractors.scopus_extractor import ScopusExtractor
        return ScopusExtractor(config)
    elif db_name == "bvs":
        from src.extractors.bvs_extractor import BvsExtractor
        return BvsExtractor(config)
    elif db_name == "wos":
        from src.extractors.wos_extractor import WosExtractor
        return WosExtractor(config)
    else:
        raise ValueError(f"Base de dados desconhecida: {db_name}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extrator bibliográfico para revisão sistemática sobre "
        "desastres ambientais e indicadores gestacionais/neonatais.",
    )
    parser.add_argument(
        "--config",
        default="config.yaml",
        help="Caminho para o arquivo de configuração YAML (padrão: config.yaml)",
    )
    parser.add_argument(
        "--databases",
        default=",".join(AVAILABLE_DATABASES),
        help="Bases de dados a consultar, separadas por vírgula (padrão: pubmed,scopus,bvs)",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Diretório de saída (sobrescreve config.yaml)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Apenas constrói e exibe as queries, sem executar buscas",
    )
    parser.add_argument(
        "--skip-dedup",
        action="store_true",
        help="Exporta sem deduplicação",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Ativa logs detalhados (DEBUG) no console",
    )
    parser.add_argument(
        "--import-bvs",
        default=None,
        metavar="ARQUIVO",
        help="Importa registros BVS de arquivo RIS ou CSV exportado manualmente",
    )
    parser.add_argument(
        "--import-scopus",
        default=None,
        metavar="ARQUIVO",
        help="Importa registros Scopus de arquivo CSV ou RIS exportado manualmente",
    )
    parser.add_argument(
        "--import-wos",
        default=None,
        metavar="ARQUIVO",
        help="Importa registros WoS de arquivo RIS ou Tab-delimited exportado manualmente",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Carregar configuração
    try:
        config = load_config(args.config)
    except (FileNotFoundError, ValueError) as e:
        print(f"Erro de configuração: {e}", file=sys.stderr)
        sys.exit(1)

    # Preparar diretório de saída
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = Path(args.output_dir or config.output.directory)
    if config.output.timestamp:
        output_dir = output_dir / timestamp
    output_dir.mkdir(parents=True, exist_ok=True)

    # Configurar logging
    log_file = str(output_dir / f"search_log_{timestamp}.log")
    setup_logging(log_file, verbose=args.verbose)

    logger.info("Iniciando extração bibliográfica")
    logger.info("Diretório de saída: %s", output_dir)

    # Inicializar PRISMA log
    prisma = PrismaLog(execution_start=datetime.now().isoformat())

    # Determinar bases de dados
    databases = [db.strip().lower() for db in args.databases.split(",")]
    for db in databases:
        if db not in AVAILABLE_DATABASES:
            logger.error("Base de dados desconhecida: %s. Opções: %s", db, AVAILABLE_DATABASES)
            sys.exit(1)

    # === Importação BVS manual ===
    bvs_imported = []
    if args.import_bvs:
        from src.extractors.bvs_extractor import BvsExtractor
        logger.info("=" * 50)
        logger.info("Importando BVS de: %s", args.import_bvs)
        logger.info("=" * 50)
        try:
            bvs_imported = BvsExtractor.import_from_file(args.import_bvs)
            prisma.databases.append(DatabaseStats(
                database_name="bvs",
                search_query="(importado de arquivo)",
                search_date=datetime.now().isoformat(),
                records_identified=len(bvs_imported),
                records_retrieved=len(bvs_imported),
            ))
        except Exception as e:
            logger.error("Falha ao importar BVS: %s", e)

    # === Importação Scopus manual ===
    scopus_imported = []
    if args.import_scopus:
        from src.extractors.scopus_extractor import ScopusExtractor
        logger.info("=" * 50)
        logger.info("Importando Scopus de: %s", args.import_scopus)
        logger.info("=" * 50)
        try:
            scopus_imported = ScopusExtractor.import_from_file(args.import_scopus)
            prisma.databases.append(DatabaseStats(
                database_name="scopus",
                search_query="(importado de arquivo)",
                search_date=datetime.now().isoformat(),
                records_identified=len(scopus_imported),
                records_retrieved=len(scopus_imported),
            ))
        except Exception as e:
            logger.error("Falha ao importar Scopus: %s", e)

    # === Importação WoS manual ===
    wos_imported = []
    if args.import_wos:
        from src.extractors.wos_extractor import WosExtractor
        logger.info("=" * 50)
        logger.info("Importando WoS de: %s", args.import_wos)
        logger.info("=" * 50)
        try:
            wos_imported = WosExtractor.import_from_file(args.import_wos)
            prisma.databases.append(DatabaseStats(
                database_name="wos",
                search_query="(importado de arquivo)",
                search_date=datetime.now().isoformat(),
                records_identified=len(wos_imported),
                records_retrieved=len(wos_imported),
            ))
        except Exception as e:
            logger.error("Falha ao importar WoS: %s", e)

    # === Extração ===
    all_records = []
    records_by_db = {}  # {db_name: [BibRecord, ...]}

    for db_name in databases:
        logger.info("=" * 50)
        logger.info("Processando: %s", db_name.upper())
        logger.info("=" * 50)

        # Se BVS já foi importado manualmente, apenas gerar query de referência
        if db_name == "bvs" and bvs_imported:
            logger.info("BVS já importado via --import-bvs (%d registros). Pulando.", len(bvs_imported))
            continue

        # Se Scopus já foi importado manualmente, pular
        if db_name == "scopus" and scopus_imported:
            logger.info("Scopus já importado via --import-scopus (%d registros). Pulando.", len(scopus_imported))
            continue

        # Se WoS já foi importado manualmente, pular
        if db_name == "wos" and wos_imported:
            logger.info("WoS já importado via --import-wos (%d registros). Pulando.", len(wos_imported))
            continue

        try:
            extractor = create_extractor(db_name, config)
        except ValueError as e:
            logger.warning("Ignorando %s: %s", db_name, e)
            continue

        # Construir query
        query = extractor.build_query()

        if args.dry_run:
            print(f"\n{'=' * 50}")
            print(f"  {db_name.upper()} — QUERY")
            print(f"{'=' * 50}")
            print(query)
            if db_name == "bvs":
                from src.extractors.bvs_extractor import BvsExtractor
                if isinstance(extractor, BvsExtractor):
                    extractor.save_query_file(str(output_dir))
                    print(f"\n  URL de busca:\n  {extractor.get_search_url()[:150]}...")
            if db_name == "scopus":
                from src.extractors.scopus_extractor import ScopusExtractor
                if isinstance(extractor, ScopusExtractor):
                    extractor.save_query_file(str(output_dir))
            if db_name == "wos":
                from src.extractors.wos_extractor import WosExtractor
                if isinstance(extractor, WosExtractor):
                    extractor.save_query_file(str(output_dir))
            print()
            continue

        # BVS: salvar query e instruções (sem tentar busca automática)
        if db_name == "bvs":
            from src.extractors.bvs_extractor import BvsExtractor
            if isinstance(extractor, BvsExtractor):
                extractor.save_query_file(str(output_dir))
                extractor.print_manual_instructions()
                prisma.databases.append(DatabaseStats(
                    database_name="bvs",
                    search_query=query,
                    search_date=datetime.now().isoformat(),
                    records_identified=0,
                    records_retrieved=0,
                ))
                continue

        # Scopus sem API key: salvar query e instruções para exportação manual
        if db_name == "scopus":
            from src.extractors.scopus_extractor import ScopusExtractor
            if isinstance(extractor, ScopusExtractor) and not extractor._has_api_key:
                extractor.save_query_file(str(output_dir))
                extractor.print_manual_instructions()
                prisma.databases.append(DatabaseStats(
                    database_name="scopus",
                    search_query=query,
                    search_date=datetime.now().isoformat(),
                    records_identified=0,
                    records_retrieved=0,
                ))
                continue

        # WoS sem API key: salvar query e instruções para exportação manual
        if db_name == "wos":
            from src.extractors.wos_extractor import WosExtractor
            if isinstance(extractor, WosExtractor) and not extractor._has_api_key:
                extractor.save_query_file(str(output_dir))
                extractor.print_manual_instructions()
                prisma.databases.append(DatabaseStats(
                    database_name="wos",
                    search_query=query,
                    search_date=datetime.now().isoformat(),
                    records_identified=0,
                    records_retrieved=0,
                ))
                continue

        # Executar busca (PubMed, Scopus, WoS)
        try:
            total = extractor.search()
        except Exception as e:
            logger.error("Falha na busca %s: %s", db_name, e)
            continue

        # Recuperar registros
        try:
            records = extractor.fetch_records()
        except Exception as e:
            logger.error("Falha ao recuperar registros %s: %s", db_name, e)
            records = extractor.records  # parcialmente recuperados

        all_records.extend(records)
        records_by_db[db_name] = records

        # Estatísticas PRISMA
        stats = extractor.get_prisma_stats()
        prisma.databases.append(DatabaseStats(
            database_name=db_name,
            search_query=query,
            search_date=datetime.now().isoformat(),
            records_identified=stats["records_identified"],
            records_retrieved=stats["records_retrieved"],
        ))

    if args.dry_run:
        logger.info("Modo dry-run: nenhuma busca executada.")
        return

    # Adicionar registros importados manualmente
    all_records.extend(bvs_imported)
    if bvs_imported:
        records_by_db["bvs"] = bvs_imported
    all_records.extend(scopus_imported)
    if scopus_imported:
        records_by_db["scopus"] = scopus_imported
    all_records.extend(wos_imported)
    if wos_imported:
        records_by_db["wos"] = wos_imported

    prisma.total_identified = sum(db.records_identified for db in prisma.databases)
    prisma.total_retrieved = len(all_records)

    logger.info("Total de registros recuperados: %d", len(all_records))

    # === Exportação individual por base ===
    for db_name, db_records in records_by_db.items():
        if not db_records:
            continue
        if "csv" in config.output.formats:
            db_csv = str(output_dir / f"{db_name}_{timestamp}.csv")
            export_csv(db_records, db_csv)
            db_xlsx = str(output_dir / f"{db_name}_{timestamp}.xlsx")
            export_excel(db_records, db_xlsx)
        if "ris" in config.output.formats:
            db_ris = str(output_dir / f"{db_name}_{timestamp}.ris")
            export_ris(db_records, db_ris)
        logger.info("%s: %d registros exportados individualmente", db_name.upper(), len(db_records))

    # === Deduplicação ===
    if args.skip_dedup:
        unique_records = all_records
        duplicate_records = []
        logger.info("Deduplicação ignorada (--skip-dedup)")
    else:
        logger.info("Iniciando deduplicação...")
        unique_records, duplicate_records = deduplicate(
            all_records,
            fuzzy_threshold=config.dedup.fuzzy_threshold,
        )

        # Contar duplicatas por tipo
        doi_dupes = sum(1 for r in duplicate_records if r.duplicate_of and r.doi)
        fuzzy_dupes = len(duplicate_records) - doi_dupes
        prisma.duplicates_removed_doi = doi_dupes
        prisma.duplicates_removed_fuzzy = fuzzy_dupes
        prisma.total_duplicates_removed = len(duplicate_records)

    prisma.records_after_dedup = len(unique_records)

    # === Exportação ===
    if "csv" in config.output.formats:
        csv_path = str(output_dir / f"results_{timestamp}.csv")
        export_csv(unique_records, csv_path)
        prisma.output_csv = csv_path

        xlsx_path = str(output_dir / f"results_{timestamp}.xlsx")
        export_excel(unique_records, xlsx_path)

    if "ris" in config.output.formats:
        ris_path = str(output_dir / f"results_{timestamp}.ris")
        export_ris(unique_records, ris_path)
        prisma.output_ris = ris_path

    # Exportar duplicatas separadamente para auditoria
    if duplicate_records:
        dupes_csv_path = str(output_dir / f"duplicates_{timestamp}.csv")
        export_csv(duplicate_records, dupes_csv_path)
        logger.info("Duplicatas exportadas: %s", dupes_csv_path)

    # === PRISMA log ===
    prisma.execution_end = datetime.now().isoformat()
    prisma_path = str(output_dir / f"prisma_log_{timestamp}.json")
    save_prisma_log(prisma, prisma_path)

    # === Resumo ===
    print_summary(prisma)

    logger.info("Extração concluída com sucesso.")


if __name__ == "__main__":
    main()
