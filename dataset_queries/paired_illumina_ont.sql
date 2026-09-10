-- paired_illumina_ont.sql
WITH runs AS (
  SELECT
    biosample,
    organism,
    bioproject,
    librarysource,
    acc,
    platform,
    instrument,
    assay_type,
    librarylayout,
    mbases,
    avgspotlen,
    releasedate
  FROM `nih-sra-datastore.sra.metadata`
  WHERE consent = 'public'
    AND biosample IS NOT NULL
    AND platform IN ('ILLUMINA', 'OXFORD_NANOPORE')
    AND assay_type = 'WGS'
    AND librarysource IN ('GENOMIC', 'METAGENOMIC')
    AND mbases >= 100                      -- drop tiny/test runs
    -- AND organism = 'human gut metagenome'   -- taxid 408170; uncomment to match SRA-MetaGut
    -- AND releasedate <= TIMESTAMP('2025-12-31')  -- Illumina run must predate your Logan snapshot
),
per_sample AS (
  SELECT
    biosample,
    ANY_VALUE(organism)    AS organism,
    ANY_VALUE(bioproject)  AS bioproject,
    ANY_VALUE(librarysource) AS librarysource,
    -- largest Illumina run
    ARRAY_AGG(
      IF(platform = 'ILLUMINA',
         STRUCT(acc, instrument, librarylayout, mbases, avgspotlen, releasedate), NULL)
      IGNORE NULLS ORDER BY mbases DESC LIMIT 1
    )[SAFE_OFFSET(0)] AS ill,
    -- largest ONT run
    ARRAY_AGG(
      IF(platform = 'OXFORD_NANOPORE',
         STRUCT(acc, instrument, librarylayout, mbases, avgspotlen, releasedate), NULL)
      IGNORE NULLS ORDER BY mbases DESC LIMIT 1
    )[SAFE_OFFSET(0)] AS ont,
    COUNTIF(platform = 'ILLUMINA')        AS n_illumina_runs,
    COUNTIF(platform = 'OXFORD_NANOPORE') AS n_ont_runs
  FROM runs
  GROUP BY biosample
)
SELECT
  biosample,
  organism,
  bioproject,
  librarysource,
  ill.acc          AS illumina_acc,
  ill.instrument   AS illumina_instrument,
  ill.mbases       AS illumina_mbases,
  ill.avgspotlen   AS illumina_avgspotlen,
  ont.acc          AS ont_acc,
  ont.instrument   AS ont_instrument,
  ont.mbases       AS ont_mbases,
  ont.avgspotlen   AS ont_avgspotlen,
  ont.releasedate  AS ont_releasedate,
  n_illumina_runs,
  n_ont_runs
FROM per_sample
WHERE ill IS NOT NULL AND ont IS NOT NULL
ORDER BY organism, ont_mbases DESC;
