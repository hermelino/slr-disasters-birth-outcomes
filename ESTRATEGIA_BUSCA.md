# Estrategia de busca

## Logica

```
(Bloco 1: Desastres) AND (Bloco 2: Desfechos gestacionais e neonatais)
```

Filtro temporal: 2000-2025

## Bloco 1 — Desastres ambientais

disaster, natural disaster, flood, earthquake, hurricane, typhoon, tornado, wildfire, tsunami, drought, cyclone, chemical spill, oil spill, volcanic eruption, landslide, mudslide, nuclear accident, radioactive hazard release, dam failure, dam collapse, mining disaster, environmental disaster, climate disaster, mass casualty

## Bloco 2 — Desfechos gestacionais e neonatais

**Gestacionais:** pregnancy outcome, pregnancy complication, premature birth, preterm birth, preterm delivery, low birth weight (LBW), small for gestational age (SGA), gestational age, spontaneous abortion, stillbirth, preeclampsia, pre-eclampsia, eclampsia, fetal death, birth weight, maternal mortality, miscarriage, gestational diabetes, gestational hypertension, intrauterine growth restriction (IUGR), premature rupture, birth outcome, adverse birth, pregnancy loss

**Neonatais:** infant mortality, neonatal mortality, perinatal mortality, neonatal death, infant death, Apgar score, congenital abnormalities, congenital anomalies, congenital malformations, birth defect, NICU, neonatal intensive care, neonatal outcome, fetal growth retardation, infant newborn diseases

## Adaptacao por base

| Base | Campo de busca | Vocabulario controlado |
|---|---|---|
| PubMed | MeSH Terms + Title/Abstract [tiab] | MeSH (com explosion hierarquico) |
| Scopus | TITLE-ABS-KEY (titulo, abstract, keywords) | Nao |
| Web of Science | TS — Topic (titulo, abstract, keywords) | Nao |
| BVS/LILACS | DeCS/MeSH (pt/en/es) | DeCS |

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
| issue | Numero/fascículo |
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
