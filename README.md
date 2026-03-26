# LDM User Guide

The LDM user guide id published [here](https://nexus-personalized-health-techno.gitbook.io/lab-data-management-app-user-g


# How to start the app in the development environment
## Getting started

1. Clone it from github

```bash
git clone git@github.com:ETH-NEXUS/lab_management.git
```

2. Create and edit the `.env` file

```bash
cd lab_management
cp .env.TEMPLATE .env
vi .env
```

&rarr; Replace all the `<>`.

3. Run it

```
docker compose  -f docker-compose.yml -f  docker-compose.dev.yml  up
```

# Production build assets

Before starting the production stack, prebuild `api`, `db`, and `ws` images, and generate UI/docs/static assets:

```bash
make build-assets
```

Then start production without building:

```bash
docker compose -f docker-compose.yml -f docker-compose.prod.yml up -d
```
