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
