# Estrategia de busca — Revisao sistematica

## Tema

Efeitos de desastres ambientais em indicadores gestacionais, neonatais e infantis.

A literatura esta dispersa entre epidemiologia, saude publica, economia da saude e ciencias ambientais.

## Bases de dados

### Primarias (utilizadas no pipeline)

| Base | Justificativa |
|---|---|
| PubMed/MEDLINE | Indispensavel para estudos epidemiologicos e clinicos sobre mortalidade materna, infantil e prematuridade associadas a desastres |
| Web of Science | Boa cobertura multidisciplinar, util para cruzar com literatura economica e ambiental |
| Scopus | Ampla cobertura em ciencias da saude e ambientais, com boas ferramentas de analise bibliometrica |
| LILACS (BVS) | Fundamental para literatura brasileira e latino-americana, incluindo o contexto do SUS |

### Complementares (busca manual) — *nao implementado*

| Base | Justificativa |
|---|---|
| EMBASE | Forte em farmacologia e saude publica, complementa o PubMed |
| SciELO | Acesso aberto com boa cobertura de periodicos brasileiros e regionais |
| REPIDISCA (PAHO/OPAS) | Repositorio especializado em saude ambiental nas Americas |
| Global Health (CAB Abstracts) | Foco em saude em paises de baixa e media renda |
| Disaster Lit (NLM) | Repositorio especializado em literatura sobre desastres e saude publica |
| ProQuest Dissertations | Teses e dissertacoes sobre o tema |

### Periodicos-chave para monitorar — *nao implementado*

- **Saude-desastres:** Disaster Medicine and Public Health Preparedness, Prehospital and Disaster Medicine, Environmental Health Perspectives
- **Saude materna e infantil:** The Lancet, Paediatric and Perinatal Epidemiology, BMC Pregnancy and Childbirth, Maternal and Child Health Journal
- **Epidemiologia ambiental:** Environmental Research, American Journal of Epidemiology, International Journal of Environmental Research and Public Health

### Literatura cinzenta — *nao implementado*

Repositorios da OMS, UNICEF, UNFPA, Ministerio da Saude (DATASUS) e Defesa Civil/CEMADEN. Muito do que existe sobre o Brasil esta em relatorios tecnicos, nao em periodicos indexados.

---

## Logica de busca

```
(Bloco 1: Desastres) AND (Bloco 2: Desfechos gestacionais e neonatais)
```

Filtro temporal: 2000-2025

## Bloco 1 — Desastres ambientais

disaster, natural disaster, flood, earthquake, hurricane, typhoon, tornado, wildfire, tsunami, drought, cyclone, chemical spill, oil spill, volcanic eruption, landslide, mudslide, nuclear accident, radioactive hazard release, dam failure, dam collapse, mining disaster, environmental disaster, climate disaster, mass casualty

## Bloco 2 — Desfechos gestacionais e neonatais

**Gestacionais:** pregnancy outcome, pregnancy complication, premature birth, preterm birth, preterm delivery, low birth weight (LBW), small for gestational age (SGA), gestational age, spontaneous abortion, stillbirth, preeclampsia, pre-eclampsia, eclampsia, fetal death, birth weight, maternal mortality, miscarriage, gestational diabetes, gestational hypertension, intrauterine growth restriction (IUGR), premature rupture, birth outcome, adverse birth, pregnancy loss

**Neonatais:** infant mortality, neonatal mortality, perinatal mortality, neonatal death, infant death, Apgar score, congenital abnormalities, congenital anomalies, congenital malformations, birth defect, NICU, neonatal intensive care, neonatal outcome, fetal growth retardation, infant newborn diseases

---

## Adaptacao por base

| Base | Campo de busca | Vocabulario controlado |
|---|---|---|
| PubMed | MeSH Terms + Title/Abstract [tiab] | MeSH (com explosion hierarquico) |
| Scopus | TITLE-ABS-KEY (titulo, abstract, keywords) | Nao |
| Web of Science | TS — Topic (titulo, abstract, keywords) | Nao |
| BVS/LILACS | DeCS/MeSH (pt/en/es) | DeCS |

---

## Resultados da busca (PRISMA)

Data da execucao: 2026-02-22/23

### Identificacao

| Base | Registros identificados |
|---|---|
| PubMed | 3.116 |
| Scopus | 1.439 |
| Web of Science | 518 |
| BVS/LILACS | 21 |
| **Total** | **5.094** |

### Deduplicacao

| Etapa | Registros |
|---|---|
| Total identificados | 5.094 |
| Duplicatas por DOI | 679 |
| Duplicatas por similaridade (fuzzy) | 68 |
| **Total duplicatas removidas** | **747** |
| **Registros unicos apos deduplicacao** | **4.347** |

### Fluxo PRISMA

```
Identificacao (n = 5.094)
    PubMed:  3.116
    Scopus:  1.439
    WoS:       518
    BVS:        21
        |
        v
Duplicatas removidas (n = 747)
    DOI exato:        679
    Titulo similar:    68
        |
        v
Registros para triagem (n = 4.347)
```

---

## Colunas do arquivo CSV de resultados

| Coluna | Descricao |
|---|---|
| source_db | Base de origem (pubmed, scopus, wos, bvs) |
| source_id | Identificador na base de origem (PMID, EID, UT, LILACS ID) |
| doi | Digital Object Identifier |
| title | Titulo do artigo |
| authors | Autores (separados por ponto e virgula) |
| journal | Nome do periodico |
| year | Ano de publicacao |
| volume | Volume do periodico |
| issue | Numero/fasciculo |
| pages | Paginas |
| issn | ISSN do periodico |
| abstract | Resumo do artigo |
| keywords | Palavras-chave dos autores (separadas por ponto e virgula) |
| mesh_terms | Descritores MeSH indexados (separados por ponto e virgula) |
| matched_disasters_mesh | Termos MeSH do artigo que correspondem ao Bloco 1 (desastres) |
| matched_gestational_mesh | Termos MeSH do artigo que correspondem aos desfechos gestacionais |
| matched_neonatal_mesh | Termos MeSH do artigo que correspondem aos desfechos neonatais |
| language | Idioma da publicacao |
| publication_type | Tipo de publicacao (journal article, review, etc.) |
| url | Link para o artigo |

---

## Estrutura do projeto

```
des_ambientais_e_ind_gestacionais/
|
|-- main.py                          # Ponto de entrada do pipeline
|-- config.yaml                      # Configuracao (API keys, filtros, output)
|-- config.example.yaml              # Exemplo de configuracao (sem credenciais)
|-- requirements.txt                 # Dependencias Python
|
|-- README.md                         # Este arquivo (estrategia de busca)
|-- ACESSO_BASES.md                  # Instrucoes de acesso a cada base de dados
|
|-- src/
|   |-- models.py                    # BibRecord — modelo normalizado de registro bibliografico
|   |-- config.py                    # Carregamento e validacao do config.yaml
|   |-- logging_prisma.py            # Logging, estatisticas PRISMA e resumo final
|   |
|   |-- extractors/
|   |   |-- base.py                  # BaseExtractor — classe abstrata (build_query, search, fetch_records)
|   |   |-- pubmed_extractor.py      # PubMed via E-utilities (API). MeSH matching e explosion
|   |   |-- scopus_extractor.py      # Scopus via pybliometrics (API) ou import manual (CSV/RIS)
|   |   |-- wos_extractor.py         # Web of Science via Starter API ou import manual (RIS/Tab-delimited)
|   |   |-- bvs_extractor.py         # BVS/LILACS — somente import manual (RIS)
|   |
|   |-- exporters/
|   |   |-- csv_exporter.py          # Exportacao CSV (UTF-8-BOM, separador ponto e virgula)
|   |   |-- ris_exporter.py          # Exportacao RIS (formato de intercambio bibliografico)
|   |
|   |-- dedup/
|       |-- deduplicator.py          # Deduplicacao em 2 fases: DOI exato + similaridade fuzzy de titulo
|
|-- output/                          # Resultados gerados pelo pipeline
|   |-- <base>_<timestamp>/          # Exportacao individual por execucao de base
|   |-- <timestamp>/                 # Execucao combinada (todas as bases + dedup)
|       |-- <base>_<ts>.csv/.ris     # Registros individuais por base
|       |-- results_<ts>.csv/.ris    # Registros combinados pos-deduplicacao
|       |-- duplicates_<ts>.csv      # Registros duplicados removidos
|       |-- search_log_<ts>.log      # Log completo da execucao
|
|-- tests/                           # Testes
```
