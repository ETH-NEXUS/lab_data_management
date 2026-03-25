# UI (Nuxt)

For the recommended Docker-based workflow, follow the root `README.md`.

## Local development (without Docker)

From `ui/app`:

```bash
pnpm install
pnpm dev
```

## Environment

- `NUXT_PUBLIC_API_URL`:
  - if set, the UI will talk to `${NUXT_PUBLIC_API_URL}/api/v1`
  - if empty, it will default to `/api/v1` (intended for the Docker dev proxy)
- `NUXT_OPENAPI_SCHEMA` (optional):
  - if set, `nuxt-open-fetch` generates typed `$api` client from that schema URL/path
  - if unset, typed client generation is skipped and plain `$fetch` can be used
