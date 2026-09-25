# NeRNA-Web

NeRNA-Web is an R Shiny application for generating and comparing synthetic negative RNA sequences. It supports miRNA, lncRNA, circRNA, and tRNA input in FASTA format.

Source repository: <https://github.com/DataRNA/NegativeRNAWeb>

Live application: <https://datarna.agu.edu.tr>

## Main functions

- Upload RNA sequences in FASTA format.
- Load a bundled example dataset with one click.
- Predict secondary structures and minimum free energy (MFE) with RNAfold.
- Generate negative sequences using NeRNA, dinucleotide-preserving shuffling, or random shuffling.
- Compare sequence identity, dinucleotide composition, GC content, RNAfold MFE, and paired-base statistics.
- Export sequence and analysis results as FASTA, CSV, or Excel files.

NeRNA-Web generates computational negative controls. The generated sequences should not automatically be interpreted as experimentally verified biological negatives.

## Load example data

On the landing page, choose an RNA type under **Try Example Data** and click **Load Example Data**. The same option is also available on **Upload FASTA**. The application loads the first 10 valid records from the corresponding bundled dataset:

- `cases/miRNA.fa`
- `cases/lncRNA_1000.fa`
- `cases/circRNA_1000.fa`
- `cases/tRNA.fasta`

The example loader is available for miRNA, lncRNA, circRNA, and tRNA. Selecting **Other** requires a user-supplied FASTA file.

## Negative-sequence methods

### NeRNA

NeRNA uses the RNAfold dot-bracket structure together with the RNA sequence, converts both into an octal representation, applies a circular bit shift, and decodes the shifted representation into a new sequence. Sequence length is preserved. RNAfold results are required before NeRNA generation.

### Dinucleotide-preserving shuffling

This method constructs an Eulerian path through the observed dinucleotide graph. It preserves sequence length and dinucleotide counts while changing the sequence order.

### Random shuffling

This method randomly permutes the nucleotides. It preserves sequence length and mononucleotide counts but does not preserve dinucleotide order.

## Reported metrics

Only directly calculated values are reported:

- **Sequence identity:** percentage of identical nucleotides at corresponding positions.
- **Dinucleotide similarity:** cosine similarity between the two named 16-dimensional dinucleotide-frequency vectors.
- **GC similarity:** `100 - absolute difference between GC percentages`.
- **MFE:** value returned by RNAfold.
- **Paired-base statistics:** calculated only from a valid RNAfold dot-bracket structure.

If RNAfold does not return a valid result, structure-dependent values are shown as **Not calculated** or `NA`. The application does not substitute fixed MFE values, estimated pairing ratios, or fabricated dot-bracket structures.

## Input limits

- Maximum sequence length: 1,000 nucleotides.
- Interactive analysis limit: 250 sequences per uploaded file.
- Files containing more than 250 valid sequences can be divided into downloadable batches.
- Standard RNA/DNA symbols are accepted; `T` is converted when required by the analysis.

## Repository structure

| Path | Purpose |
|---|---|
| `ui.R` | Shiny user interface |
| `server.R` | Shiny server logic and exports |
| `global.R` | Sequence processing, RNAfold integration, and helper functions |
| `cases/` | Bundled RNA datasets used by the example-data loader |
| `ML_models/` | Trained machine-learning models produced by the study workflows |
| `Results_data/` | Tables and analysis results produced for the study and revision |
| `utils/` | Local runtime utilities, including the Windows RNAfold executable when present |
| `scripts/` | Validation and maintenance scripts |

`ML_models/` and `Results_data/` remain part of the GitHub repository for reproducibility. They are listed in `.dockerignore` because the running web application does not need them.

## Run with Docker

```bash
docker build -t nerna-web .
docker run --rm -p 3838:3838 nerna-web
```

Open <http://localhost:3838>.

The Docker image installs ViennaRNA and verifies that `RNAfold` is available during the build.

## Run locally

Required software:

- R 4.4 or a compatible recent R version
- ViennaRNA/RNAfold
- R packages listed in the `Dockerfile`

From the repository directory:

```r
shiny::runApp(".", host = "127.0.0.1", port = 3838)
```

On Windows, the application first checks `RNAFOLD_PATH`, then the bundled executable under `utils/`, and finally the system `PATH`.

## Reproducibility notes

- Use an explicit random seed when generating shuffled sequences.
- Record the NeRNA shift value and software versions.
- RNAfold-derived values are reported only when a valid dot-bracket structure or finite MFE is available.
- Loading a new uploaded or example dataset clears results belonging to the previous dataset.

## License

See [LICENSE](LICENSE).
