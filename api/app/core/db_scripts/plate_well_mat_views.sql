-- Create the well details view
DROP MATERIALIZED VIEW IF EXISTS core_welldetail;
CREATE MATERIALIZED VIEW core_welldetail AS
        SELECT 
                id,
                plate_id,
                type,
                status,
                position,
                hr_position,
                initial_amount,
                withdrawal,
                initial_amount - withdrawal AS amount,
                compounds,
                measurements
        FROM (
                SELECT 
                        w.id,
                        p.id as plate_id,
                        wt.name as type,
                        w.status,
                        w.position,
                        hr_position(d.cols, w.position) AS hr_position,
                        COALESCE(( SELECT SUM(amount) FROM core_wellcompound WHERE well_id = w.id GROUP BY well_id ), 0) AS initial_amount,
                        COALESCE(( SELECT SUM(amount) FROM core_wellwithdrawal WHERE well_id = w.id GROUP BY well_id ), 0) AS withdrawal,
                        ( SELECT ARRAY(SELECT c.name from core_wellcompound AS wc INNER JOIN compoundlib_compound c ON wc.compound_id = c.id WHERE wc.well_id = w.id) ) AS compounds
                FROM core_well AS w
                INNER JOIN core_plate p ON w.plate_id = p.id
                INNER JOIN core_platedimension d ON p.dimension_id = d.id
                INNER JOIN core_welltype wt ON w.type_id = wt.id
                ORDER BY w.position
        ) s1
        LEFT JOIN (
                SELECT 
                        well_id,
                        JSON_OBJECT_AGG(label, value) AS measurements
                FROM (
                        SELECT 
                                well_id, 
                                label,  
                                ARRAY_AGG(value ORDER BY measured_at) as value
                        FROM core_measurement
                        GROUP BY well_id, label     
                ) timeseries
                GROUP BY well_id
        ) s2 
        ON s1.id = s2.well_id;



-- Cerate indexes on core_welldetail
CREATE UNIQUE INDEX ON core_welldetail(id);
CREATE INDEX core_welldetail_plate_id_idx on core_welldetail(plate_id);

-- Create the plate detail view
DROP MATERIALIZED VIEW IF EXISTS core_platedetail;
CREATE MATERIALIZED VIEW core_platedetail AS
        SELECT
                s1.id,
                MAX(s1.rows) * MAX(s1.cols) AS num_wells,
                ARRAY_AGG(DISTINCT s2.label ORDER BY s2.label) as measurement_labels,
                -- MIN is a bit missleading it should just take the first
                MIN(s2.timestamp) as measurement_timestamps,
                json_merge(COALESCE(JSON_OBJECT_AGG(s2.label, JSON_BUILD_OBJECT(well_type, JSON_BUILD_OBJECT('min', s2.min, 'max', s2.max, 'mean', s2.mean, 'median', s2.median, 'std', s2.std, 'mad', s2.mad))) FILTER (WHERE s2.label IS NOT NULL), '{}')) AS stats,
                json_merge(COALESCE(JSON_OBJECT_AGG(s3.label, JSON_BUILD_OBJECT('min', s3.min, 'max', s3.max, 'mean', s3.mean, 'median', s3.median, 'std', s3.std, 'mad', s3.mad)) FILTER (WHERE s3.label IS NOT NULL), '{}')) AS overall_stats
        FROM (
                SELECT 
                        p.id,
                        d.rows,
                        d.cols
                FROM core_plate AS p
                INNER JOIN core_platedimension d ON p.dimension_id = d.id
        ) s1
        -- per well type stats
        LEFT JOIN (
                SELECT
                        label,
                        well_type,
                        plate_id,
                        ARRAY_AGG(measured_at ORDER BY measured_at) as timestamp,
                        ARRAY_AGG(min ORDER BY measured_at) as min,
                        ARRAY_AGG(max ORDER BY measured_at) as max,
                        ARRAY_AGG(std ORDER BY measured_at) as std,
                        ARRAY_AGG(mean ORDER BY measured_at) as mean,
                        ARRAY_AGG(median ORDER BY measured_at) as median,
                        ARRAY_AGG(mad ORDER BY measured_at) as mad
                FROM (
                        SELECT 
                                label,
                                wt.name AS well_type,
                                measured_at,
                                plate_id,
                                MIN(value),
                                MAX(value),
                                STDDEV_POP(value) AS std,
                                AVG(value) as mean,
                                PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY value) AS median,
                                MAD(value)
                        FROM core_measurement AS m 
                        INNER JOIN core_well w ON m.well_id = w.id
                        INNER JOIN core_welltype wt ON w.type_id = wt.id
                        GROUP BY plate_id, label, wt.name, measured_at
                ) s5
                GROUP BY plate_id, label, well_type
        ) s2
        ON s1.id = s2.plate_id
        -- overall stats
        LEFT JOIN (
                SELECT 
                        label,
                        plate_id,
                        --ARRAY_AGG(measured_at ORDER BY measured_at) as timestamp,
                        ARRAY_AGG(min ORDER BY measured_at) as min,
                        ARRAY_AGG(max ORDER BY measured_at) as max,
                        ARRAY_AGG(std ORDER BY measured_at) as std,
                        ARRAY_AGG(mean ORDER BY measured_at) as mean,
                        ARRAY_AGG(median ORDER BY measured_at) as median,
                        ARRAY_AGG(mad ORDER BY measured_at) as mad
                FROM (
                        SELECT 
                                label,
                                measured_at,
                                plate_id,
                                MIN(value),
                                MAX(value),
                                STDDEV_POP(value) AS std,
                                AVG(value) as mean,
                                PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY value) AS median,
                                MAD(value) as mad
                        FROM core_measurement AS m 
                        INNER JOIN core_well w ON m.well_id = w.id
                        GROUP BY plate_id, label, measured_at
                ) s4
                GROUP BY plate_id, label
        ) s3
        ON s1.id = s3.plate_id
        GROUP BY id;

-- Create indexes on core_platedetail
CREATE UNIQUE INDEX ON core_platedetail(id);
        
        
        
        
        
        
        
        
        