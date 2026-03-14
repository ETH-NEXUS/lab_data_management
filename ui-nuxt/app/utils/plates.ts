/**
 * Returns a safe plate identifier string from a dynamic route param.
 *
 * Accepted input examples:
 * - `'42'`
 * - `42`
 * - `['42']`
 * - `undefined`
 *
 * Returned output examples:
 * - `'42'`
 * - `'-'` (fallback when value is missing)
 */
export const getPlateRouteIdLabel = (routeParam: unknown): string => {
  if (Array.isArray(routeParam)) {
    const first = routeParam[0]
    return typeof first === 'string' || typeof first === 'number' ? String(first) : '-'
  }

  if (typeof routeParam === 'string' || typeof routeParam === 'number') {
    return String(routeParam)
  }

  return '-'
}

