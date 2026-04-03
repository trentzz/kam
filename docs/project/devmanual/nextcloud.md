# Nextcloud Storage

Large files that are not tracked in git are stored on Nextcloud. The `bigdata/` directory in the repository root mirrors the Nextcloud structure locally. Never commit files inside `bigdata/`.

## Shares

**Read-only share (public)**
URL: `https://nextcloudlocal.trentz.me/s/pTizAiSAJQsPcDo`
WebDAV root: `https://nextcloudlocal.trentz.me/public.php/dav/files/pTizAiSAJQsPcDo/`

**Edit share (developers only)**
Token: stored in `.env` as `NEXTCLOUD_EDIT_TOKEN`. Ask the project maintainer for access. See `.env.example`.

## Downloading

Use the public read-only WebDAV URL to download files:

```bash
curl -O "https://nextcloudlocal.trentz.me/public.php/dav/files/pTizAiSAJQsPcDo/<path>/<file>" \
     -u "pTizAiSAJQsPcDo:"
```

Replace `<path>/<file>` with the Nextcloud path from the table below.

## Uploading

Uploads require the edit share token. Copy `.env.example` to `.env` and fill in `NEXTCLOUD_EDIT_TOKEN`.

The upload scripts read the token automatically from `.env`. Run them from the repository root.

**Twist duplex ML dataset:**
```bash
./scripts/ml/upload_twist_duplex.sh [train|test|models|features|all]
```

The upload URL pattern is:
```
https://nextcloudlocal.trentz.me/public.php/dav/files/<NEXTCLOUD_EDIT_TOKEN>/<path>
```

## Structure

| Nextcloud path | Local path (`bigdata/`) | Contents | Reproduced by |
|---|---|---|---|
| `benchmarking/01-snvindel/` | `benchmarking/01-snvindel/` | varforge FASTQs and kam outputs for SNV/indel benchmark | `docs/benchmarking/01-snvindel/` scripts |
| `benchmarking/02-sv-core/` | `benchmarking/02-sv-core/` | varforge FASTQs for core SV types (large deletion, tandem duplication, inversion) | `docs/benchmarking/02-sv-core/` scripts |
| `benchmarking/03-sv-extended/` | `benchmarking/03-sv-extended/` | varforge FASTQs for extended SV types (insertion, large deletion, inversion-deletion) | `docs/benchmarking/03-sv-extended/` scripts |
| `benchmarking/04-comparison/` | `benchmarking/04-comparison/` | Alignment-baseline comparison data | `docs/benchmarking/04-comparison/` scripts |
| `benchmarking/05-public/` | `benchmarking/05-public/` | Public dataset downloads | `docs/benchmarking/05-public/` scripts |
| `benchmarking/06-runtime/` | `benchmarking/06-runtime/` | Runtime profiling outputs | `docs/benchmarking/06-runtime/` scripts |
| `experiments/01-ml-twist-duplex/` | `experiments/01-ml-twist-duplex/` | ML training simulations, extracted features, and trained models | `docs/project/experiments/01-ml-twist-duplex/` scripts |
| `paper/` | `paper/` | Compiled paper PDFs and supplementary materials | `docs/paper/` LaTeX source |
