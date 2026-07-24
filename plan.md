# Part A: show all detail-view fields in the main inventory stock table

Each batch touches the same 3 files: `inventory-stock-table.columns.ts`,
`inventory-stock-table.values.ts`, `locales/en.json`.

- [x] Batch 1: favorite, storage_temperature (split out of itemType), brand, manufacturer, vendor
- [x] Batch 2: manufacturer_catalog_number, vendor_catalog_number, capacity, default_cost, is_active
- [x] Batch 3: description, serial_number, order_number, lifetime_days
- [x] Batch 4: quantity_in_base_units, minimum_quantity_in_base_units, safety_data_sheet
- [ ] Batch 5: created_at, updated_at

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
- [ ] Batch 5: serial_number, order_number, lifetime_days, is_active
