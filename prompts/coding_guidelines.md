# Coding guidelines (read once at the start of every run)

These rules apply to every code cell you write in the notebook.

## Reproducibility
- Set `np.random.seed(42)` and `random.seed(42)` in the notebook setup cell.
- Pin all package versions in `manifest.json` via `conda run -n proteomics-agent pip freeze`.
- Use absolute paths everywhere (`os.path.abspath(...)`). The user may `cd` away mid-session.
- Hash every input file (sha256) and record in `manifest.json`.

## Caching (mandatory for all external API calls)
```python
import requests_cache
requests_cache.install_cache(
    f'{output_dir}/annotations/api_cache',
    backend='sqlite',
    expire_after=86400 * 30,  # 30 days
    allowable_methods=['GET', 'POST'],
)
```
This makes re-runs free and offline-capable. Verify cache works by re-running a query and checking `from_cache` attribute on the response.

## Retries (mandatory for all REST calls)
```python
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type
import requests

@retry(
    stop=stop_after_attempt(4),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    retry=retry_if_exception_type((requests.HTTPError, requests.ConnectionError, requests.Timeout)),
)
def safe_get(url, **kwargs):
    r = requests.get(url, timeout=30, **kwargs)
    r.raise_for_status()
    return r
```

## Rate limiting
- UniProt: ≤5 req/sec. Use batch endpoints (100 accessions per call).
- STRING: ≤1 req/sec. Always batch.
- KEGG REST: ≤3 req/sec. They block IPs that abuse it.
- AlphaFold: no published limit, but be polite — 2 req/sec.

## Notebook style
- One conceptual operation per cell. Don't put data loading + plotting + saving in one cell.
- Markdown cell **before** every code cell explaining the intent (this is what makes the notebook readable later).
- Markdown cell **after** every output cell interpreting what the result means biologically (1-3 sentences).
- Use `print()` for status, `display()` from IPython for DataFrames.

## Plotting
- Static: matplotlib, 300 DPI, save with `bbox_inches='tight'`. Always `plt.close()` after `savefig` to free memory.
- Interactive: plotly or pyvis, save as HTML.
- Color scheme: red (#D62728) for up, blue (#1F77B4) for down, grey (#AAAAAA) for non-significant. Use this consistently across all figures.
- Always label axes with units. Always include a colorbar legend if color encodes data.

## Error handling
- Wrap every external API call in try/except. On failure, log to `manifest.json:degraded_stages` and continue.
- For code execution errors: read traceback, attempt fix, retry up to 3 times. If still failing, pause and ask user.
- Never silently swallow exceptions. At minimum, print a warning.

## Memory
- For multi-file mode with large inputs: stream xlsx loading where possible, drop unused columns immediately, `del` large intermediates and `gc.collect()`.

## What NOT to do
- Don't compute new t-tests from `log2_mean` and `log2_std` — the user already did this and the input may have used a more appropriate test (Welch's, moderated, limma-like).
- Don't filter the input file on read — keep all rows visible in QC. Apply DE filters at Stage 4.
- Don't make up biology. If you're unsure what a gene does, query UniProt or say "function unclear" in the report.
- Don't paraphrase >15 words from any web source verbatim. Summarize in your own words and cite.
