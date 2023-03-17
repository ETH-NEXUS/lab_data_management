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

## Map echo output files

In the api container run:

```bash
./manage.py map echo -i /data/echo 
```

If the column names in your echo output files differ from the default ones, you should specify a mapping file:

```bash
./manage.py map echo -p /data/echo -m /data/echo_mapping.yml
```

## Map M1000 measurement files

In the api container run:

```bash
./manage.py map m1000 -p /data/m1000
```

If you need to apply a custom formula to the result values, e. g. to find the Acceptor/Donor ratio, you can specify the formula and provide the name of the result measurement value in the command as follows:
```bash
./manage.py map m1000 -p /data/ex1/m1000 -e Acceptor/Donor -n result_value
```

