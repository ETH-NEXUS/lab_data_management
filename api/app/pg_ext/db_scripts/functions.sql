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
                        if isinstance(d[k], dict):
                                d[k].update(v)
                        elif isinstance(d[k], list):
                                d[k] = list(set(d[k] + v))
                        else:
                                d[k] += f", {v}"
                else:
                        d[k] = v
        return d

return json.dumps(json.loads(jdata, object_pairs_hook=merge_on_duplicate_keys))
$$ LANGUAGE plpython3u;

-- Caclulation of Median Absolute Deviation
DROP FUNCTION IF EXISTS _final_mad(double precision[]) CASCADE;

CREATE OR REPLACE FUNCTION _final_mad(vals double precision[]) RETURNS double precision AS $$
import statistics

def mad(values, constant=1.4826):
    median = statistics.median(values)
    deviations = [abs(x - median) for x in values]
    mad = statistics.median(deviations)
    return constant * mad

return mad(list(vals))

$$ LANGUAGE plpython3u;

CREATE OR REPLACE AGGREGATE MAD (
  SFUNC=array_append,
  BASETYPE=double precision,
  STYPE=double precision[],
  INITCOND='{}',
  FINALFUNC=_final_mad
);

