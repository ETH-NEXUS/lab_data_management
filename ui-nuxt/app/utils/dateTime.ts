/**
 * Formats a date-like input into a localized human-readable date/time string.
 *
 * Accepted data examples:
 * - `dateValue = '2026-03-14T10:45:00Z'`
 * - `dateValue = new Date('2026-03-14T10:45:00Z')`
 * - `dateValue = null` (returns fallback)
 *
 * Returned data examples:
 * - `'Mar 14, 2026, 11:45 AM'` (depends on user locale/timezone)
 * - `'-'` when no date is provided
 * - original input string when parsing fails (for example `'invalid-date'`)
 */
export const formatDateTime = (
  dateValue?: string | Date | null,
  options: Intl.DateTimeFormatOptions = { dateStyle: 'medium', timeStyle: 'short' },
  fallback = '-',
): string => {
  if (!dateValue) return fallback

  const date = dateValue instanceof Date ? dateValue : new Date(dateValue)
  if (Number.isNaN(date.getTime())) return String(dateValue)

  return new Intl.DateTimeFormat(undefined, options).format(date)
}
