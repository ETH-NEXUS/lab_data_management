-- We create the python langage
CREATE OR REPLACE LANGUAGE plpython3u;

-- A funciton to calculate the human readable position given a position and a number of columns
DROP FUNCTION IF EXISTS hr_position(integer, integer) CASCADE;

CREATE OR REPLACE FUNCTION hr_position(cols integer, pos integer) RETURNS text AS $$
import math
from string import ascii_uppercase
row = math.floor(pos / cols) + 1
col = pos - (row - 1) * cols + 1
try:
        letter = ""
        while row > 0:
                row, remainder = divmod(row - 1, 26)
                letter = ascii_uppercase[remainder] + letter
        return f"{letter}{col}"
except IndexError:
        return "?"
$$ LANGUAGE plpython3u;

-- A function to merge json objects with equal keys
DROP FUNCTION IF EXISTS json_merge(json) CASCADE;

CREATE OR REPLACE FUNCTION json_merge(jdata json) RETURNS json AS $$
import json

def merge_on_duplicate_keys(ordered_pairs):
        d = {}
        for k, v in ordered_pairs:
                if k in d:
                        d[k].update(v)
                else:
                        d[k] = v
        return d

d = json.loads(jdata, object_pairs_hook=merge_on_duplicate_keys)

ret = {}
for key, value in d.items():
        if key not in ret:
                ret[key] = {} 
        ret[key].update(value)
return json.dumps(ret)
$$ LANGUAGE plpython3u;

-- Caclulation of Median Absolute Deviation
DROP FUNCTION IF EXISTS _final_mad(double precision[]) CASCADE;

CREATE OR REPLACE FUNCTION _final_mad(vals double precision[]) RETURNS double precision AS $$
from scipy.stats import median_abs_deviation as mad
return mad(vals)
$$ LANGUAGE plpython3u;

CREATE OR REPLACE AGGREGATE MAD (
  SFUNC=array_append,
  BASETYPE=double precision,
  STYPE=double precision[],
  INITCOND='{}',
  FINALFUNC=_final_mad
);

-- Create the well details view
DROP VIEW IF EXISTS core_welldetail;
CREATE OR REPLACE VIEW core_welldetail AS
        SELECT 
                id,
                position,
                well_type,
                hr_position,
                plate_id,
                barcode,
                JSON_BUILD_OBJECT('id', sample_id, 'name', sample) AS sample,
                initial_amount,
                withdrawal,
                initial_amount - withdrawal AS amount,
                compounds,
                measurements
        FROM (
                SELECT 
                        w.id,
                        w.position,
                        wt.name AS well_type,
                        hr_position(d.cols, w.position) AS hr_position,
                        w.plate_id,
                        p.barcode,
                        w.sample_id,
                        s.name as sample,
                        COALESCE(( SELECT SUM(amount) FROM core_wellcompound WHERE well_id = w.id GROUP BY well_id ), 0) AS initial_amount,
                        COALESCE(( SELECT SUM(amount) FROM core_wellwithdrawal WHERE well_id = w.id GROUP BY well_id ), 0) AS withdrawal,
                        ( SELECT ARRAY(SELECT c.name from core_wellcompound AS wc INNER JOIN compoundlib_compound c ON wc.compound_id = c.id WHERE wc.well_id = w.id) ) AS compounds
                FROM core_well AS w
                INNER JOIN core_plate p ON w.plate_id = p.id
                INNER JOIN core_platedimension d ON p.dimension_id = d.id
                LEFT OUTER JOIN core_welltype wt ON w.type_id = wt.id
                LEFT OUTER JOIN core_sample s ON w.sample_id = s.id
                ORDER BY w.position
        ) s1
        LEFT JOIN (
                SELECT 
                        well_id,
                        JSON_OBJECT_AGG(label, value) AS measurements
                FROM core_measurement AS m 
                GROUP BY well_id
        ) s2 
        ON s1.id = s2.well_id;

-- Create the plate detail view
DROP VIEW IF EXISTS core_platedetail;
CREATE OR REPLACE VIEW core_platedetail AS
        SELECT
                id,
                barcode,
                JSON_BUILD_OBJECT('dimension', dimension, 'rows', MAX(rows), 'cols', MAX(cols), 'wells', MAX(rows) * MAX(cols)) AS dimension,
                ARRAY_AGG(DISTINCT label ORDER BY label) as measurement_labels,
                JSON_BUILD_OBJECT('min', MIN(min), 'max', MAX(max)) AS min_max,
                json_merge(COALESCE(JSON_OBJECT_AGG(label, JSON_BUILD_OBJECT(well_type, JSON_BUILD_OBJECT('min', min, 'max', max, 'mean', mean, 'median', median, 'stddev', stddev, 'mad', mad))) FILTER (WHERE label IS NOT NULL), '{}')) AS stats
        FROM (
                SELECT 
                        p.id,
                        p.barcode,
                        d.name AS dimension,
                        MIN(d.rows) AS rows,
                        MIN(d.cols) AS cols
                FROM core_plate AS p
                INNER JOIN core_well w ON w.plate_id = p.id
                INNER JOIN core_platedimension d ON p.dimension_id = d.id
                GROUP BY p.id, p.barcode, d.name
        ) s1
        LEFT JOIN (
                SELECT 
                        label,
                        wt.name AS well_type,
                        w.plate_id,
                        MIN(value),
                        MAX(value),
                        STDDEV_POP(value) AS stddev,
                        AVG(value) as mean,
                        PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY value) AS median,
                        MAD(value),
                        ARRAY_AGG(value)
                FROM core_measurement AS m 
                INNER JOIN core_well w ON m.well_id = w.id
                INNER JOIN core_welltype wt ON w.type_id = wt.id
                GROUP BY w.plate_id, label, wt.name
        ) s2 
        ON s1.id = s2.plate_id
        GROUP BY id, barcode, dimension
        
        
        
        
        
        
        
        
        