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