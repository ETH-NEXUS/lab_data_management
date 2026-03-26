-- A function to calculate the human readable position given a position and a number of columns.
CREATE OR REPLACE FUNCTION hr_position(cols integer, pos integer)
RETURNS text
LANGUAGE plpgsql
IMMUTABLE
AS $$
DECLARE
    row_num integer;
    col_num integer;
    letters text := '';
    remainder integer;
BEGIN
    IF cols IS NULL OR cols <= 0 OR pos IS NULL OR pos < 0 THEN
        RETURN '?';
    END IF;

    row_num := pos / cols + 1;
    col_num := pos - (row_num - 1) * cols + 1;

    WHILE row_num > 0 LOOP
        remainder := (row_num - 1) % 26;
        letters := chr(65 + remainder) || letters;
        row_num := (row_num - 1) / 26;
    END LOOP;

    RETURN letters || col_num::text;
EXCEPTION WHEN OTHERS THEN
    RETURN '?';
END;
$$;

-- A function to merge JSON objects with duplicate keys.
-- If values are objects: merge nested keys.
-- If values are arrays: set-like union.
-- Otherwise: concatenate scalar values as "a, b".
CREATE OR REPLACE FUNCTION json_merge(jdata json)
RETURNS json
LANGUAGE plpgsql
IMMUTABLE
AS $$
DECLARE
    pair record;
    result jsonb := '{}'::jsonb;
    existing jsonb;
    incoming jsonb;
BEGIN
    IF jdata IS NULL THEN
        RETURN NULL;
    END IF;

    FOR pair IN SELECT key, value FROM json_each(jdata)
    LOOP
        incoming := pair.value::jsonb;
        existing := result -> pair.key;

        IF existing IS NULL THEN
            result := result || jsonb_build_object(pair.key, incoming);
        ELSIF jsonb_typeof(existing) = 'object' AND jsonb_typeof(incoming) = 'object' THEN
            result := jsonb_set(result, ARRAY[pair.key], existing || incoming, true);
        ELSIF jsonb_typeof(existing) = 'array' AND jsonb_typeof(incoming) = 'array' THEN
            result := jsonb_set(
                result,
                ARRAY[pair.key],
                (
                    SELECT COALESCE(jsonb_agg(value), '[]'::jsonb)
                    FROM (
                        SELECT DISTINCT value
                        FROM jsonb_array_elements(existing || incoming) AS e(value)
                    ) dedup
                ),
                true
            );
        ELSE
            result := jsonb_set(
                result,
                ARRAY[pair.key],
                to_jsonb(trim(both '"' from existing::text) || ', ' || trim(both '"' from incoming::text)),
                true
            );
        END IF;
    END LOOP;

    RETURN result::json;
END;
$$;

-- Calculation of Median Absolute Deviation.
CREATE OR REPLACE FUNCTION _final_mad(vals double precision[])
RETURNS double precision
LANGUAGE sql
IMMUTABLE
AS $$
    WITH filtered AS (
        SELECT v
        FROM unnest(vals) AS t(v)
        WHERE v IS NOT NULL
    ),
    median_value AS (
        SELECT percentile_cont(0.5) WITHIN GROUP (ORDER BY v) AS med
        FROM filtered
    ),
    mad_value AS (
        SELECT percentile_cont(0.5) WITHIN GROUP (ORDER BY abs(f.v - m.med)) AS mad
        FROM filtered f
        CROSS JOIN median_value m
    )
    SELECT CASE
        WHEN EXISTS (SELECT 1 FROM filtered)
        THEN 1.4826 * (SELECT mad FROM mad_value)
        ELSE NULL
    END;
$$;

DO $$
BEGIN
    CREATE AGGREGATE MAD (
      SFUNC = array_append,
      BASETYPE = double precision,
      STYPE = double precision[],
      INITCOND = '{}',
      FINALFUNC = _final_mad
    );
EXCEPTION
    WHEN duplicate_function OR duplicate_object THEN
        NULL;
END $$;
