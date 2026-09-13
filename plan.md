# Part A: show all detail-view fields in the main inventory stock table

Each batch touches the same 3 files: `inventory-stock-table.columns.ts`,
`inventory-stock-table.values.ts`, `locales/en.json`.

- [x] Batch 1: favorite, storage_temperature (split out of itemType), brand, manufacturer, vendor
- [x] Batch 2: manufacturer_catalog_number, vendor_catalog_number, capacity, default_cost, is_active
- [x] Batch 3: description, serial_number, order_number, lifetime_days
- [x] Batch 4: quantity_in_base_units, minimum_quantity_in_base_units, safety_data_sheet
- [x] ~~Batch 5: created_at, updated_at~~ — not needed, skipped

# Part B: cover all material fields in the "Add new item" form

Feedback: all information needed to add a new item should be fillable in one
place (matching the stock table / Excel sheet), instead of having to go back
and edit the item from the table afterward. All of these fields stay optional
— only Quantity, Material, Room/Sector, and Stock unit remain required.

Same pattern each batch: extend the "Additional material details" section in
`InventoryAddItemModal.vue`, its draft fields in `inventoryAddItemForm.utils.ts`
and `useInventoryAddItemFormState.ts`, and its optional-PATCH payload builder
in `useInventoryAddItemForm.ts`.

- [x] Batch 1: brand, default_cost
- [x] Batch 2: manufacturer, vendor, manufacturer_catalog_number, vendor_catalog_number (currently read-only display -> made editable)
- [x] Batch 3: capacity_value, capacity_unit
- [x] Batch 4: description
- [x] Batch 5: serial_number, order_number, lifetime_days, is_active

# Part C: personalize inventory favorites

- [x] Backend and database: store favorites per user and return the current user's value.
- [x] UI: keep the existing favorite actions and remove the outdated creation checkbox.

# Part D: inventory activity and check-in/check-out cards

- [x] Backend: include the source order and its project in check-in history.
- [x] Backend: paginate inventory history for dashboard requests.
- [x] Backend: filter check-in and check-out history records.
- [x] UI data: request the first five history records.
- [x] UI data: request the first five check-in and check-out records.
- [x] Dashboard: show color-coded activities and check-in/check-out tiles.

# Part E: paginated recent activities table

- [x] UI data: add a server-paginated history query.
- [x] UI: add the history table, page navigation, and dashboard link.

# Part F: paginated check-in/check-out table

- [x] UI data: add a server-paginated check-in/check-out query.
- [x] UI: add the check-in/check-out table, page navigation, and dashboard link.

# Part G: inventory dashboard cleanup

- [x] Backend: fix favorite sorting, avoid history favorite N+1 queries, and add regression tests.
- [x] UI: open Material Usage from history usage records and show dashboard history request errors.
- [x] UI: merge the duplicate history workspaces into one variant-based workspace.
- [x] Backend: limit awaiting-check-in orders and recent project usages to five records server-side.
- [x] UI: use the new limited dashboard requests.
- [x] UI: format the remaining Add Item files so the inventory lint check passes.

# Part H: recently linked LDM experiments

- [x] Backend: return the five latest Harvest-project and LDM-experiment usages.
- [x] UI: add the LDM experiment dashboard card.

# Part I: personalized inventory dashboard tiles

- [x] Backend: store available tiles and each user's selected tiles in the database.
- [x] UI: let the user choose any dashboard tiles and render only the saved selection.

# Part J: dashboard tile reliability and ordering

- [x] Backend: make first-load preferences conflict-safe and validate duplicate or unknown keys.
- [x] UI: allow any number of tiles, handle an empty dashboard, and let the user choose card positions.
- [x] UI: render cards in the saved order for both visual and keyboard navigation.
- [x] UI: do not load data for hidden dashboard tiles.
- [x] UI: reset unsaved tile selection when the settings dialog is reopened.
- [x] UI: load each stock preview as one server-paginated request of five items.
- [x] UI: show preview loading failures explicitly and move stock preview rendering into its own component.
- [x] UI: move the device preview and its queries into a dedicated dashboard component.
- [x] UI: move the awaiting-check-in preview and its query into a dedicated dashboard component.
- [x] UI: render project and experiment usage previews with one shared component.

# Part K: organize inventory UI components

- [x] UI: group dashboard, history, stock-table, and add-item components into thematic directories and update imports.

# Problematic wells (messages page) — volume threshold not reported

- [x] Step 1 (backend): fix `find_problems mark_empty_wells` — decide the plate status
      after all wells are checked, treat 0 as "empty" instead of "no data", and read the
      newest withdrawal by `created_at`.
- [x] Step 2 (backend): store `current_amount` in µL everywhere (variant A) — the
      library-copy path now carries the last reported fill level over in µL, or None
      when it was never reported. No copy rows to migrate in the dev database.
- [x] Step 3 (backend): the threshold check in `Plate.map` now runs for every mapped
      well, not only when an existing `WellWithdrawal` is updated, and both it and
      `find_problems` share one check in `core/thresholds.py`.
- [x] Step 4 (UI): the well panel now says explicitly when the instrument reported
      nothing, distinguishes that from a reported zero, and labels every volume with
      its unit (µL for the reported fill level, nL for the stored amount).
- [x] Audit follow-up 1: `/api/compoundlib/redflag/` and `/api/compoundlib/recalculate_status/`
      now require a logged in user (they were reachable without any authentication).
- [x] Audit follow-up 2a: the test suite runs again — it now uses Postgres instead of
      the broken sqlite test settings, the pytest config moved into `api/pyproject.toml`
      (the only one the api image sees) and the stale `Compound.identifier` usages in
      `core/tests.py` are gone. The 10 tests that failed before are skipped (see 2a-1).
- [x] `StatisticsTest` deleted: it tested `plate.z_prime` / `z_factor` / `z_scores`, which no
      longer exist.
- [ ] Audit follow-up 2a-1: rewrite the 7 skipped tests that cover existing code —
      `MapperTests` (missing test data files) and `InventoryMaterialReagentTests`
      (never authenticates, gets 403).
- [x] Audit follow-up 2b: unit tests for `core/thresholds.py` (pure function, no database).
- [x] Audit follow-up 2c: tests for `find_problems` and for the flagging in `Plate.map`.
      Writing them surfaced a regression from step 3: `Mapping` defaulted its reported
      values to 0 ("empty") instead of None ("not reported"), so a plate copy or a csv
      mapping of a library plate would have marked every source well as empty.
- [x] Audit follow-up 3a (API): `/api/compoundlib/redflag/` returns each marked well with
      its last reported values and the thresholds it is below.
- [x] Audit follow-up 3b (UI): the problematic plates card shows those values, the value
      below the threshold in red and the reason in words.
- [x] Audit follow-up 4: the recalculation runs on POST instead of GET, so DRF checks the
      CSRF token and a link or another site cannot start it.
- [x] Audit follow-up 5a (model): `Threshold` has units in its help text, rejects negative
      values and a DMSO above 100 %, and is read everywhere through `Threshold.current()`,
      which creates the defaults when the row is missing instead of silently marking nothing.
- [x] Refactor: `core/views.py` (1042 lines) is now the `core/views/` package, one file per
      topic, re-exported from `__init__.py`. Pure move: every definition is unchanged and the
      467 URL routes resolve to the same views.
- [x] Audit follow-up 5b (API only): through the API the threshold can be read and changed by
      every logged in user, but no longer created or deleted. The admin stays as it is.
- [x] Audit follow-up 5c (UI): the threshold form refuses a DMSO above 100 % and shows the
      error when the API refuses a value.
- [x] `/api/refresh/` (refreshing the materialized views) requires a logged in user and a POST
      with a CSRF token; it was reachable without any authentication and ran on GET.

- [x] Before release, step 1: recalculation and Echo import only mark plates without a status
      (None or ""). A status set by hand, like "disposed" (6 plates on production), is kept.
- [x] Tests centralized in `api/app/tests/<app>/`; the large `tests.py` files of core, importer
      and inventory are split by class, the 7 empty `tests.py` stubs are gone, and pytest only
      looks in `tests/` (`testpaths`). Same 130 tests, code unchanged.
- [ ] Before release, step 2: data migration that turns the 0 uL / 100 % fill levels written by
      the old library copy (2484 wells on Drug01_C … Drug08_C in production) into "unknown".

## Deferred to a separate audit (found, not changed — the code has worked for years)

- Echo import crashes with `KeyError: 'DMSO'` on files from the newer Echo software, which
  write `Fluid Composition` / `Fluid Units` / `Fluid Type` instead of `% DMSO`
  (example: `data_temp/test/_data_examples_echo_Testrun1_*_Transfer_*.csv`). The column
  schema in `core/config.py` is strict, so a fix touches the schema, `ldm.yaml` and the mapper.
- `MapperTests` (4 skipped) fail with `FileNotFoundError: ./temp/M1000/20210902-131750_BAF210901_1.asc`.
- `InventoryMaterialReagentTests` (skipped) never authenticates, so the API answers 403.
- Dev database: 1 library plate has no dimension. If it ever gets flagged, `RedFlagView`
  would crash on `well.hr_position`.
- Duplicate withdrawals with a target well: 0 in dev. `Plate.map` does not catch
  `MultipleObjectsReturned`, and there is no unique constraint on (well, target_well).
- Other audit items not started: re-importing a file doubles withdrawals, archived plates are
  listed on the messages page, `find_problems` ignores an unknown argument silently,
  recalculation runs synchronously in the request, the `Problem` model is unused.

- [x] `%_COPY%` check on production (2026-09-13): plate copies never wrote a withdrawal there,
      so no data migration is needed and Recalculate cannot raise false alarms from them.
