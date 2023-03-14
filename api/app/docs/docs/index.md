# Lab Management

## Getting started

1. Clone it from github

```bash
git clone git@github.com:ETH-NEXUS/lab_management.git
```

2. Create the `.env` file

```bash
cd lab_management
cp .env.TEMPLATE .env
vi .env
```

&rarr; Replace all the `<>`.

3. Run it

```
make up
```

## Import a compound library

In the api container run:

```bash
./manage.py import sdf -i /data/S230470.sdf -m /data/S230470_mapping.yml -r 16 -c 24
```
