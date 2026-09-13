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
- [ ] Audit follow-up 2a-1: unskip and fix the 10 tests that were failing before the
      suite could run again — `StatisticsTest` (calls `plate.z_prime` / `z_factor` /
      `z_scores`, which the model no longer has), `MapperTests` (missing test data
      files) and `InventoryMaterialReagentTests` (never authenticates, gets 403).
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
- [ ] Open: run the `%_COPY%` check on production and add a data migration if it
      returns rows.
