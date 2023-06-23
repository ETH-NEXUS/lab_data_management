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

If you have not created barcode specifications for the plates in the ui, provide the experiment name so that they can be created during mapping:

```bash
./manage.py map echo -i /data/echo --experiment_name exp3
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
./manage.py map m1000 -p /data/m1000 -e Acceptor/Donor -n result_value
```

You can apply several formulas separating them with commas (in this case, you need to provide the same number of names for the result measurement values):
```bash
./manage.py map m1000 -p /data/m1000 -e Acceptor/Donor,Acceptor*Donor -n ratio,product
```

If you know that for some reason some destination plates were not created by the echo mapping, you can create them during the measurement mapping. 
For this, you need to add two parameters to the command: 

```bash
./manage.py map m1000 -p /data/ex3/M1000 --create_missing_plates --experiment_name exp3

```
Import control plate

```bash
./manage.py import library_plate -i /data/control_palte_example.csv --library_name controls_library --plate_barcode control12345

```
